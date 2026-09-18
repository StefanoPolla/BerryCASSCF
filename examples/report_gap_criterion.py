#!/usr/bin/env python
"""When is a state-averaged minimum gap small enough to call a conical intersection?

`docs/todo.md` §6. The project has no such criterion, and its absence has already produced
a wrong statement: butadiene CAS(12,12) bottoms out at 1.5 mHa where every other rung reaches
0.003-0.006, and that was written up as "no intersection in this plane" -- which over-reads
it, because 1.5 mHa is 0.04 eV, inside the error of the model itself.

**This does not invent a threshold.** A tolerance picked after seeing which answer it gives is
worthless, and any fixed number is system- and basis-dependent. What it does is report the
*ingredients*, so that a reader can apply their own:

``floor``      the smallest gap a direct 2D search reached, in mHa. An upper bound on the
               true minimum, not the minimum.
``slope``      how fast the gap opens away from the apex along the scanned cut, mHa/deg,
               from the cone fit. This is what converts an energy into a distance.
``miss``       ``floor / slope``: how far the cut passes from the apex **if the surface really
               is a cone**, in degrees. The honest reading of a floor -- an intersection whose
               apex the cut misses by ``miss`` degrees would present exactly this floor.
``noise``      the convergence floor, from the SA-CASSCF threshold and the measured mirror
               asymmetries. A floor below this is not distinguishable from zero.
``radius``     the loop's semi-axis in the same coordinate, for scale: a miss distance small
               against the radius means the loop still encloses whatever is there.

The comparison that matters is ``miss`` against ``radius``, not ``floor`` against a tolerance:
loop transport asks whether a degeneracy is inside a loop, and a cut that misses an apex by
0.4 deg has found an intersection the loop encloses, whatever the gap at the cut says.

    python examples/report_gap_criterion.py butadiene
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf.refine import cone_apex
from berrycasscf.scan import ScanResult
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The SA-CASSCF energy threshold is 1e-9 Ha, and gaps are differences of two converged
# energies, so nothing below ~1e-6 mHa is meaningful. The measured mirror asymmetries of the
# ethylene scans (1e-11 to 3e-5 mHa) bracket the same figure from the other side.
NOISE_MHA = 1e-3

SYSTEMS = {
    "butadiene": {"cut": "butadiene_refinepyr_cas{ne}-{ncas}.npz",
                  "search": "gap_minimum_search.json",
                  "radius": 18.0, "coord": "pyr"},
    "ethylene": {"cut": "ethylene_scan_cas{ne}-{ncas}_*.npz",
                 "search": None, "radius": 12.0, "coord": "phi"},
}
LADDER = [(2, 2), (4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]


def searched_floor(system: str, ne: int, ncas: int) -> float | None:
    """The smallest gap a direct 2D search reached at this rung, in mHa."""
    spec = SYSTEMS[system]
    if not spec["search"]:
        return None
    path = os.path.join(ROOT, "results", system, spec["search"])
    if not os.path.exists(path):
        return None
    for run in json.load(open(path))["runs"]:
        if tuple(run["cas"]) == (ne, ncas):
            return float(run["gap"])
    return None


def cut_fit(system: str, ne: int, ncas: int):
    """Cone fit along this rung's finest cut: (position, closest approach, slope, residual)."""
    pattern = os.path.join(ROOT, "results", system,
                           SYSTEMS[system]["cut"].format(ne=ne, ncas=ncas))
    hits = sorted(glob.glob(pattern))
    if not hits:
        return None
    res = ScanResult.load(hits[0])
    row = int(np.nanargmin(np.nanmin(res.gap, axis=1)))
    return cone_apex(res.phis, res.gap[row] * 1e3, window=3)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS), nargs="?", default="butadiene")
    args = ap.parse_args()
    spec = SYSTEMS[args.system]

    rows = []
    for ne, ncas in LADDER:
        fit = cut_fit(args.system, ne, ncas)
        if fit is None:
            continue
        floor = searched_floor(args.system, ne, ncas)
        best = floor if floor is not None else fit.closest_approach
        slope = abs(fit.slope)
        miss = best / slope if slope > 0 else float("nan")
        rows.append({
            "cas": [ne, ncas],
            "floor_mHa": floor,
            "fitted_closest_approach_mHa": float(fit.closest_approach),
            "slope_mHa_per_deg": float(slope),
            "fit_residual": float(fit.residual),
            "miss_deg": float(miss),
            "miss_over_radius": float(miss / spec["radius"]),
            "above_noise": bool(best > NOISE_MHA),
        })

    if not rows:
        print(f"no refined cuts for {args.system}; run the refine stage first")
        return 1

    print(f"{args.system}: is the gap-scan minimum an intersection?   "
          f"(loop semi-axis in {spec['coord']} = {spec['radius']:.0f} deg, "
          f"convergence noise {NOISE_MHA:g} mHa)\n")
    print(f"{'CAS':>8} {'floor (mHa)':>12} {'slope (mHa/deg)':>16} {'miss (deg)':>11} "
          f"{'miss/radius':>12} {'above noise':>12}")
    print("-" * 76)
    for r in rows:
        floor = r["floor_mHa"]
        shown = f"{floor:.4f}" if floor is not None else f"{r['fitted_closest_approach_mHa']:.4f}*"
        print(f"{('(%d,%d)' % tuple(r['cas'])):>8} {shown:>12} "
              f"{r['slope_mHa_per_deg']:>16.3f} {r['miss_deg']:>11.4f} "
              f"{r['miss_over_radius']:>12.4f} {str(r['above_noise']):>12}")
    print("\n* fitted closest approach, where no direct search exists for that rung.")
    print("\nRead the miss distance, not the floor: a cut that passes this far from the apex")
    print("of a true cone would show exactly the floor observed. Compare it against the loop")
    print("semi-axis -- that is the question loop transport actually answers.")

    out = os.path.join(ROOT, "results", args.system, "gap_criterion.json")
    save_json({"system": args.system, "noise_mHa": NOISE_MHA,
               "loop_semi_axis_deg": spec["radius"], "coordinate": spec["coord"],
               "rows": rows}, out)
    print(f"\n-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
