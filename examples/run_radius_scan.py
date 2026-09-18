#!/usr/bin/env python
"""Does a loop stay walkable as it shrinks onto a conical intersection?

Everything about locating a degeneracy by loop transport turns on this, and it had never been
measured. Bisection works by making a loop *graze* the seam, which is the expensive, failure-prone
configuration; an iterative scheme that re-centres on the intersection and shrinks would instead
keep the loop **concentric**, where every point is the same distance from the seam.

Two competing expectations, and this settles which one holds:

* near a conical intersection the adiabatic state depends on the *azimuthal angle*, so the
  wavefunction turns by pi over any enclosing loop regardless of its size. The discretization
  requirement is then **angular** and a small concentric loop should cost the same as a large one;
* but as the radius shrinks every point on the loop becomes near-degenerate, and the tracked state
  is a *state-specific* CASSCF ground state with S1 right on top of it. Root flipping should
  eventually defeat it, putting a floor on the usable radius.

The scan reports, per radius: the verdict, how many points the adaptive walk needed, its worst
adjacent overlap, and its cost. Flat cost with radius favours the first reading; a knee where the
verdicts turn to refusals measures the floor, which is the number an iterative scheme has to live
within.

    python examples/run_radius_scan.py ethylene --cas 2 2 --centre 90 111.0
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
from berrycasscf.ethylene import ethylene_geom
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.localize import evaluate_radius
from berrycasscf.runlog import JobLog
from berrycasscf.store import load_json, save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GEOM = {"formaldimine": (formaldimine_geom, "sto-3g"),
        "ethylene": (ethylene_geom, "6-31g*"),
        "butadiene": (butadiene_geom, "6-31g*")}

# Halving, then finer than any loop the study has used. The largest is comparable to the loops
# the ladder walks; the smallest is well below the resolution any bisection here has reached.
DEFAULT_RADII = [4.0, 2.0, 1.0, 0.5, 0.25, 0.1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(GEOM))
    ap.add_argument("--cas", nargs=2, type=int, default=(2, 2), metavar=("NE", "NCAS"))
    ap.add_argument("--centre", nargs=2, type=float, required=True, metavar=("X", "Y"),
                    help="the point to centre the loops on -- the best available estimate of "
                         "the intersection, since the question is what happens when a loop is "
                         "concentric with it")
    ap.add_argument("--radii", nargs="*", type=float, default=DEFAULT_RADII,
                    help="loop radii in degrees, circular (default: 4 2 1 0.5 0.25 0.1)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    geom_fn, basis = GEOM[args.system]
    ne, ncas = args.cas
    tag = f"{args.system}_cas{ne}-{ncas}" + (f"_{args.label}" if args.label else "")
    out = os.path.join(ROOT, "results", "radius_scan", f"{tag}.json")
    done = {}
    if os.path.exists(out) and not args.force:
        record = load_json(out)
        if record.get("complete"):
            print(f"SKIP: {os.path.relpath(out, ROOT)} exists; --force to redo")
            return 0
        # Each radius is independent, so a killed scan resumes at the first one missing.
        done = {float(p["radius"]): p for p in record.get("probes", [])}
        if done:
            print(f"resuming: {len(done)} radii already measured")

    cas = CasConfig(basis=basis, ncas=ncas, nelecas=ne)
    cont = ContinuationConfig()
    centre = tuple(args.centre)
    log = JobLog(f"radius_scan_{tag}", total=len(args.radii))
    print(f"Concentric loops about {centre}: {args.system} {cas.cas_label}/{basis}")
    print(f"  radii (deg): {args.radii}\n")

    t0 = time.time()
    probes = []
    for r in args.radii:
        if float(r) in done:
            probes.append(done[float(r)])
            print(f"  [cached] r = {r}")
            continue
        probe = evaluate_radius(centre, (r, r), 1.0, cas=cas, cont=cont,
                                geom_fn=geom_fn, name=f"R{r:g}", progress=log)
        rec = probe.to_dict()
        rec["radius"] = float(r)
        probes.append(rec)
        save_json({"system": args.system, "cas": [ne, ncas], "basis": basis,
                   "centre": list(centre), "radii": list(args.radii),
                   "probes": probes, "complete": len(probes) == len(args.radii),
                   "wall_time": time.time() - t0}, out)

    print("\n" + "=" * 78)
    print(f"{'radius':>8} {'verdict':>14} {'points':>14} {'worst overlap':>14} "
          f"{'micro':>9} {'wall (s)':>9}")
    print("-" * 78)
    for rec in probes:
        runs = rec.get("runs", [])
        pts = "/".join(str(r["n_points"]) for r in runs) if runs else "-"
        worst = min((r["min_overlap"] for r in runs if np.isfinite(r["min_overlap"])),
                    default=float("nan"))
        print(f"{rec['radius']:>8.3g} {rec['verdict']:>14} {pts:>14} {worst:>14.3f} "
              f"{rec['cost_micro']:>9} {rec['wall_time']:>9.1f}")
    log.done()
    print(f"\n-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
