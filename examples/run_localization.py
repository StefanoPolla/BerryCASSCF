#!/usr/bin/env python
"""Locate a degeneracy using loop transport alone, by bisection and triangulation.

`docs/todo.md` §2. The gap scan finds an intersection by evaluating the gap everywhere and
looking for the smallest. Loop transport never evaluates a gap -- it returns one bit per
loop. This driver turns those bits into a position:

  1. **bisect** the loop radius about a fixed centre until the phase turns over. The
     transition radius is the elliptical distance from that centre to the degeneracy, so one
     centre confines it to an ellipse;
  2. **triangulate** by repeating from other centres and intersecting the ellipses.

Formaldimine is the validation case and is run first, because its answer is known
independently: an FCI gap scan puts the intersection at alpha = 132.61 deg
(`docs/results.md`). Anything the method claims there can be checked.

It is also the sharpest possible test of the method's value. At CAS(2,2) the SA-CASSCF gap
scan does not merely misplace this intersection -- it **finds no minimum in the region at
all**, reporting a spurious one 18 deg away. If bisection at the same CAS(2,2) recovers a
position near 132.6, then loop transport yields a *geometric* answer where the comparator
at equal cost yields a wrong one.

    python examples/run_localization.py formaldimine --cas 2 2
    python examples/run_localization.py butadiene --cas 2 2
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ContinuationConfig
from berrycasscf.butadiene import butadiene_geom
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.localize import (
    bisect_radius,
    elliptical_radius,
    merge_centre_records,
    resume_bisections,
    triangulate,
)
from berrycasscf.runlog import JobLog
from berrycasscf.store import load_json, save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SYSTEMS = {
    "formaldimine": {
        "geom_fn": formaldimine_geom,
        "basis": "sto-3g",
        "shape": (10.0, 10.0),
        # Three centres, deliberately NOT collinear: two on the phi = 90 line to fix alpha,
        # one displaced in phi to break the mirror ambiguity that two centres always leave.
        "centres": [(130.0, 89.9), (137.0, 90.0), (133.0, 97.0)],
        "scale_hi": 1.0,
        "scale_lo": 0.05,
        "tol": 0.03,
        "max_probes": 11,
        # FCI gap scan, docs/results.md. phi = 90 by the symmetry of the coordinate.
        "reference": (132.61, 90.0),
        "reference_note": "FCI gap scan (docs/results.md); phi = 90 by symmetry",
    },
    "butadiene": {
        "geom_fn": butadiene_geom,
        "basis": "6-31g*",
        "shape": (12.0, 18.0),
        # Every centre must ENCLOSE the target at scale 1 and exclude it at scale_lo, or the
        # bisection has nothing to bracket. Taking CAS(2,2)'s searched intersection
        # (89.97, 105.67) as the target, the elliptical radii are 0.21, 0.87 and 0.78 -- all
        # inside 1. The third is deliberately off the tw = 90 line, because two centres on a
        # line leave a mirror ambiguity that nothing else here can break.
        "centres": [(90.0, 101.85321091497578), (90.0, 90.0), (99.0, 101.85321091497578)],
        "scale_hi": 1.0,
        "scale_lo": 0.08,
        "tol": 0.05,
        "max_probes": 9,
        "reference": None,
        "reference_note": "no exact reference exists for butadiene",
    },
}


def report(results, centres, shape, ref, spec, args, ne, ncas, out, wall_time):
    """Print the brackets, triangulate, and write the record.

    Shared by a single-job run and by --merge so that a localization assembled from
    per-centre jobs is byte-for-byte the same kind of record as one run in series.
    """
    payload = {
        "system": args.system, "cas": [ne, ncas], "basis": spec["basis"],
        "shape": list(shape), "centres": [list(c) for c in centres],
        "reference": list(ref) if ref else None,
        "reference_note": spec["reference_note"],
        "bisections": [r.to_dict() for r in results],
        "complete": True,
        "wall_time": wall_time,
    }

    print("=" * 78)
    for r in results:
        print("  " + r.summary())
        if ref and r.bracketed:
            expected = elliptical_radius(ref, r.centre, shape)
            print(f"      reference rho {expected:.4f} -> "
                  f"{'INSIDE' if r.lo <= expected <= r.hi else 'OUTSIDE'} the bracket "
                  f"(off by {abs(expected - r.rho):.4f})")

    bracketed = [r for r in results if r.bracketed]
    if len(bracketed) >= 2:
        tri = triangulate(bracketed, prefer=ref)
        payload["triangulation"] = tri.to_dict()
        print(f"\n  triangulated position: "
              f"({tri.chosen[0]:.3f}, {tri.chosen[1]:.3f})" if tri.chosen else
              "\n  triangulation failed")
        if tri.candidates:
            print(f"    candidates: "
                  + ", ".join(f"({p[0]:.3f}, {p[1]:.3f})" for p in tri.candidates))
            print(f"    residual {tri.residual:.5f} (in units of the loop semi-axes)")
            print(f"    note: {tri.note}")
        if ref and tri.chosen:
            err = np.hypot(tri.chosen[0] - ref[0], tri.chosen[1] - ref[1])
            payload["error_vs_reference"] = float(err)
            print(f"\n  ERROR VS REFERENCE: {err:.3f} deg "
                  f"(loop transport {tri.chosen[0]:.2f} vs reference {ref[0]:.2f})")
    else:
        print("\n  fewer than two centres bracketed; no triangulation")

    payload["total_micro"] = sum(r.total_micro for r in results)
    print(f"\n  total cost: {payload['total_micro']} micro-iterations, {wall_time:.0f} s")
    save_json(payload, out)
    print(f"-> {os.path.relpath(out, ROOT)}")
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS))
    ap.add_argument("--cas", nargs=2, type=int, default=(2, 2), metavar=("NE", "NCAS"))
    ap.add_argument("--centres", type=int, default=3, help="how many centres to use")
    ap.add_argument("--centre-xy", nargs="*", default=None, metavar="X,Y",
                    help="override the centres, e.g. --centre-xy 90,101.85 90,90 99,110. "
                         "A centre is only useful if its full-size loop encloses the target "
                         "without grazing it: one that grazes is refused and costs a bisection "
                         "for nothing (docs/todo.md §9).")
    ap.add_argument("--only-centre", type=int, default=None, metavar="K",
                    help="run just centre K (0-based) and save it on its own. Bisections "
                         "about different centres share nothing, so at a large active space "
                         "they run as concurrent jobs and are combined afterwards with "
                         "--merge, which turns days in series into a day in parallel.")
    ap.add_argument("--merge", action="store_true",
                    help="combine the per-centre records written by --only-centre into the "
                         "full record, and triangulate")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    spec = SYSTEMS[args.system]
    ne, ncas = args.cas
    tag = f"{args.system}_cas{ne}-{ncas}"
    out = os.path.join(ROOT, "results", "localize", f"{tag}.json")
    if args.only_centre is not None:
        out = os.path.join(ROOT, "results", "localize", f"{tag}_centre{args.only_centre}.json")
    existing = load_json(out) if os.path.exists(out) else None
    # Records written before incremental saving carry no flag and are complete by
    # construction, so a missing key means complete.
    if existing is not None and existing.get("complete", True) and not args.force:
        print(f"SKIP: {os.path.relpath(out, ROOT)} exists; --force to redo")
        return 0

    if args.only_centre is not None and args.merge:
        raise SystemExit("--only-centre and --merge are alternatives, not a combination")

    cas = CasConfig(basis=spec["basis"], ncas=ncas, nelecas=ne)
    cont = ContinuationConfig()
    if args.centre_xy:
        centres = [tuple(float(v) for v in c.split(",")) for c in args.centre_xy]
    else:
        centres = spec["centres"][: args.centres]
    shape = spec["shape"]
    ref = spec["reference"]

    if args.merge:
        paths = [os.path.join(ROOT, "results", "localize", f"{tag}_centre{k}.json")
                 for k in range(len(centres))]
        missing = [os.path.relpath(q, ROOT) for q in paths if not os.path.exists(q)]
        if missing:
            raise SystemExit("cannot merge, these per-centre records are missing:\n  "
                             + "\n  ".join(missing))
        results = merge_centre_records([load_json(q) for q in paths], centres)
        print(f"Merged {len(results)} per-centre records for {args.system} "
              f"{cas.cas_label}/{spec['basis']}")
        report(results, centres, shape, ref, spec, args, ne, ncas, out,
               wall_time=sum(load_json(q).get("wall_time", 0.0) for q in paths))
        return 0

    print(f"Locating the degeneracy with loop transport alone: {args.system} "
          f"{cas.cas_label}/{spec['basis']}")
    print(f"  loop shape {shape}, {len(centres)} centres")
    if ref:
        print(f"  reference position {ref}  ({spec['reference_note']})")
        for c in centres:
            print(f"    centre {c}: reference rho would be "
                  f"{elliptical_radius(ref, c, shape):.4f}")
    print()

    # One log per *task*: with --only-centre the three centres of a rung run concurrently,
    # and a shared log interleaves three bisections into something no one can follow.
    log_tag = tag if args.only_centre is None else f"{tag}_centre{args.only_centre}"
    log = JobLog(f"localize_{log_tag}",
                 total=(len(centres) if args.only_centre is None else 1) * spec["max_probes"])
    t0 = time.time()
    # In --only-centre mode the record is written once, when that centre finishes, so an
    # incomplete file cannot exist and the completeness check above has already skipped a
    # finished one. Only the all-centres path has anything to resume.
    saved = (existing or {}).get("bisections", []) if args.only_centre is None else []
    results = [] if args.force else resume_bisections(saved, centres)
    if results:
        print(f"  resuming: {len(results)} of {len(centres)} centres already saved")
        for r in results:
            print("    " + r.summary())

    def payload_now(complete: bool) -> dict:
        return {
            "system": args.system, "cas": [ne, ncas], "basis": spec["basis"],
            "shape": list(shape), "centres": [list(c) for c in centres],
            "reference": list(ref) if ref else None,
            "reference_note": spec["reference_note"],
            "bisections": [r.to_dict() for r in results],
            "complete": complete,
            "wall_time": time.time() - t0,
        }

    wanted = range(len(centres)) if args.only_centre is None else [args.only_centre]
    for i, centre in enumerate(centres):
        if i not in wanted or i < len(results):
            continue
        print(f"--- centre {i+1}/{len(centres)}: {centre} ---")
        res = bisect_radius(
            centre, shape, cas=cas, cont=cont, geom_fn=spec["geom_fn"],
            scale_hi=spec["scale_hi"], scale_lo=spec["scale_lo"],
            tol=spec["tol"], max_probes=spec["max_probes"],
            name=f"{args.system[:3].upper()}{i}", progress=log,
        )
        results.append(res)
        # Save before starting the next centre: hours of bisection should not depend on the
        # job surviving to the end.
        done = len(results) if args.only_centre is None else 1
        want = len(centres) if args.only_centre is None else 1
        save_json(payload_now(complete=done == want), out)
        print(f"    saved {done}/{want} centres -> {os.path.relpath(out, ROOT)}")
        print()

    if args.only_centre is not None:
        print(f"centre {args.only_centre} done: {results[-1].summary()}")
        print(f"-> {os.path.relpath(out, ROOT)}")
        print(f"   merge with: python examples/run_localization.py {args.system} "
              f"--cas {ne} {ncas} --merge"
              + (f" --centre-xy {' '.join(f'{c[0]},{c[1]}' for c in centres)}"
                 if args.centre_xy else ""))
        return 0

    log.done()
    report(results, centres, shape, ref, spec, args, ne, ncas, out,
           wall_time=time.time() - t0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
