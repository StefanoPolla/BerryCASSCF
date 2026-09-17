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
from berrycasscf.localize import bisect_radius, elliptical_radius, triangulate
from berrycasscf.runlog import JobLog
from berrycasscf.store import save_json

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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS))
    ap.add_argument("--cas", nargs=2, type=int, default=(2, 2), metavar=("NE", "NCAS"))
    ap.add_argument("--centres", type=int, default=3, help="how many centres to use")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    spec = SYSTEMS[args.system]
    ne, ncas = args.cas
    tag = f"{args.system}_cas{ne}-{ncas}"
    out = os.path.join(ROOT, "results", "localize", f"{tag}.json")
    if os.path.exists(out) and not args.force:
        print(f"[skip] {os.path.relpath(out, ROOT)} exists; --force to redo")
        return 0

    cas = CasConfig(basis=spec["basis"], ncas=ncas, nelecas=ne)
    cont = ContinuationConfig()
    centres = spec["centres"][: args.centres]
    shape = spec["shape"]
    ref = spec["reference"]

    print(f"Locating the degeneracy with loop transport alone: {args.system} "
          f"{cas.cas_label}/{spec['basis']}")
    print(f"  loop shape {shape}, {len(centres)} centres")
    if ref:
        print(f"  reference position {ref}  ({spec['reference_note']})")
        for c in centres:
            print(f"    centre {c}: reference rho would be "
                  f"{elliptical_radius(ref, c, shape):.4f}")
    print()

    log = JobLog(f"localize_{tag}", total=len(centres) * spec["max_probes"])
    t0 = time.time()
    results = []
    for i, centre in enumerate(centres):
        print(f"--- centre {i+1}/{len(centres)}: {centre} ---")
        res = bisect_radius(
            centre, shape, cas=cas, cont=cont, geom_fn=spec["geom_fn"],
            scale_hi=spec["scale_hi"], scale_lo=spec["scale_lo"],
            tol=spec["tol"], max_probes=spec["max_probes"],
            name=f"{args.system[:3].upper()}{i}", progress=log,
        )
        results.append(res)
        print()

    payload = {
        "system": args.system, "cas": [ne, ncas], "basis": spec["basis"],
        "shape": list(shape), "reference": list(ref) if ref else None,
        "reference_note": spec["reference_note"],
        "bisections": [r.to_dict() for r in results],
        "wall_time": time.time() - t0,
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

    total_micro = sum(r.total_micro for r in results)
    payload["total_micro"] = total_micro
    print(f"\n  total cost: {total_micro} micro-iterations, {time.time()-t0:.0f} s")
    log.done()
    save_json(payload, out)
    print(f"-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
