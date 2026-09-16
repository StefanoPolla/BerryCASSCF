#!/usr/bin/env python
"""Butadiene active-space ladder, on the plane located by the search phase.

Same design as the ethylene study: everything except the active space is held fixed, and the
loops are the same for every rung. Run ``search_butadiene_ci.py`` first.

Two differences from ethylene, both consequences of the molecule being bigger:

* **There is no exact in-basis reference.** The full valence space is CAS(22,22). The largest
  affordable rung is the best available, so "converged" here means "the answer has stopped
  moving", never "the answer is exact".
* **The grids are thin (5 x 13).** Ethylene's intersection was pinned to ``tau = 90`` by
  symmetry; butadiene has no such symmetry, so ``tw`` has to be sampled rather than assumed.
  Five rows resolve it while keeping the cost of the expensive rungs bounded; ``pyr``, along
  which the gap varies fastest, gets thirteen.

    python examples/run_butadiene_study.py scan
    python examples/run_butadiene_study.py berry
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ContinuationConfig, ScanConfig, run_loop, scan_gap
from berrycasscf.butadiene import (
    DEFAULT_BASIS,
    SEARCH_PLANES,
    mirror_partner,
    plane_geom_fn,
)
from berrycasscf.casscf import apply_singlet_constraint, build_mol, run_rhf
from berrycasscf.geometry import Loop
from berrycasscf.scan import ScanResult
from berrycasscf.store import berry_record_exists, save_berry_run

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(ROOT, "results", "butadiene")

# The plane the search selected, and the region around the intersection it found.
PLANE = "tw_pyr"
CI_REGION = Loop("butadiene_ci", centre=(90.0, 110.0), radius=(20.0, 30.0))

CAS_LADDER = [(4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]
REFERENCE_CAS = (12, 12)
GRID = (5, 13)                 # tw in {70, 80, 90, 100, 110}; pyr in 5 deg steps, 80..140

# CAS(12,12) measured out at >2.4 min per cold point (853k determinants on 68 AOs), so a full
# 5 x 13 grid would be ~2.6 h for one rung. It is run instead on the single tw = 90 row, which
# is where every other rung's two-dimensional minimum lies and which carries the only quantity
# compared across the ladder -- the refined pyr position. The cost is that its tw is assumed
# rather than resolved; see docs/butadiene.md.
GRID_OVERRIDE: dict[tuple[int, int], tuple[int, int]] = {(12, 12): (1, 13)}
REGION_OVERRIDE: dict[tuple[int, int], Loop] = {}

# The Berry stage stops below CAS(12,12): a state-specific solve there costs minutes, so three
# loops at two discretizations would run to several hours for no change in the conclusion, which
# is already established over four rungs.
BERRY_LADDER = [(4, 4), (6, 6), (8, 8), (10, 10)]
LOOP_RADIUS = (12.0, 18.0)     # degrees in (tw, pyr)
NPOINTS = [13, 21]


def grid_for(ne: int, ncas: int) -> tuple[int, int]:
    return GRID_OVERRIDE.get((ne, ncas), GRID)


def region_for(ne: int, ncas: int) -> Loop:
    """The scanned region; a one-row grid is centred on tw = 90 by giving it zero tw radius."""
    if grid_for(ne, ncas)[0] == 1:
        return Loop(CI_REGION.name, (90.0, CI_REGION.centre[1]), (0.0, CI_REGION.radius[1]))
    return CI_REGION


def scan_path(ne: int, ncas: int) -> str:
    g = grid_for(ne, ncas)
    return os.path.join(RESULT_DIR, f"butadiene_scan_cas{ne}-{ncas}_{g[0]}x{g[1]}.npz")


def berry_path(loop: str, ne: int, ncas: int, n: int) -> str:
    return os.path.join(RESULT_DIR, f"butadiene_{loop}_cas{ne}-{ncas}_N{n}.json")


def geom_fn():
    (cx, cy), fixed, _ = SEARCH_PLANES[PLANE]
    return plane_geom_fn(cx, cy, **fixed)


def symmetry_spot_check(ne: int, ncas: int, basis: str, points, geom) -> float:
    """Max |gap(tw, pyr) - gap(-tw, -pyr)| in mHa over ``points``.

    Butadiene has no ``tw -> 180 - tw`` symmetry (the two methylene hydrogens are inequivalent),
    so the ethylene-style check on the grid itself does not apply. The exact symmetry is
    reflection through the molecular plane, which maps outside the scanned region -- hence a
    handful of extra solves rather than a property of the grid.
    """
    from pyscf import mcscf

    worst = 0.0
    for tw, pyr in points:
        gaps = []
        for x, y in ((tw, pyr), mirror_partner(tw, pyr)):
            mol = build_mol(geom(x, y), basis)
            mc = mcscf.CASSCF(run_rhf(mol), ncas, ne)
            apply_singlet_constraint(mc)
            mc = mc.state_average_([0.5, 0.5])
            mc.verbose = 0
            mc.kernel()
            gaps.append(float(mc.e_states[1] - mc.e_states[0]) * 1000.0)
        worst = max(worst, abs(gaps[0] - gaps[1]))
    return worst


def refined_pyr(res: ScanResult) -> float:
    """Sub-grid intersection position along pyr, on the row holding the 2D minimum.

    Ethylene's analogue used the symmetry line; here the row is chosen from the data, since
    nothing forces the intersection to sit at tw = 90.
    """
    i = int(np.unravel_index(np.nanargmin(res.gap), res.gap.shape)[0])
    row = res.gap[i]
    j = int(np.nanargmin(row))
    if j in (0, len(row) - 1):
        return float(res.phis[j])
    y0, y1, y2 = row[j - 1], row[j], row[j + 1]
    d = y0 - 2 * y1 + y2
    if abs(d) < 1e-18:
        return float(res.phis[j])
    return float(res.phis[j] + 0.5 * (y0 - y2) / d * (res.phis[1] - res.phis[0]))


def do_scans(args) -> int:
    os.makedirs(RESULT_DIR, exist_ok=True)
    fn = geom_fn()
    for ne, ncas in CAS_LADDER:
        path = scan_path(ne, ncas)
        if os.path.exists(path) and not args.force:
            res = ScanResult.load(path)
            if np.isfinite(res.e_states).all():
                print(f"[skip] CAS({ne},{ncas}): pyr = {refined_pyr(res):7.2f}   "
                      f"min {res.min_gap_point()[2]*1e3:8.3f} mHa")
                continue
        g = grid_for(ne, ncas)
        print(f"\n=== scan CAS({ne},{ncas})/{args.basis}  {g[0]}x{g[1]} ===")
        t0 = time.time()
        res = scan_gap(
            region_for(ne, ncas),
            cas=CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
            scan=ScanConfig(n_alpha=g[0], n_phi=g[1], margin=0.0),
            geom_fn=fn,
            progress=None if args.quiet else print,
            checkpoint=path,
        )
        res.save(path)
        sym = symmetry_spot_check(ne, ncas, args.basis,
                                  [(90.0, 100.0), (90.0, 120.0)], fn)
        print(f"  pyr = {refined_pyr(res):.2f}   min {res.min_gap_point()[2]*1e3:.3f} mHa   "
              f"plane-reflection check {sym:.1e} mHa   [{time.time()-t0:.0f} s]")
    return 0


def study_loops(centre: tuple[float, float]) -> dict[str, Loop]:
    t0, p0 = centre
    return {
        "B_x": Loop("B_x", (t0, p0), LOOP_RADIUS),
        "B_1": Loop("B_1", (t0 - 2.5 * LOOP_RADIUS[0], p0), LOOP_RADIUS),
        "B_2": Loop("B_2", (t0 + 2.5 * LOOP_RADIUS[0], p0), LOOP_RADIUS),
    }


def do_berry(args) -> int:
    os.makedirs(RESULT_DIR, exist_ok=True)
    if args.centre:
        centre = tuple(args.centre)
    else:
        path = scan_path(*REFERENCE_CAS)
        if not os.path.exists(path):
            raise SystemExit(f"Reference scan {os.path.relpath(path, ROOT)} not found; "
                             "run the scan stage first or pass --centre TW PYR.")
        centre = (90.0, refined_pyr(ScanResult.load(path)))
    fn = geom_fn()
    loops = study_loops(centre)
    print(f"Loops centred on (tw, pyr) = ({centre[0]:.2f}, {centre[1]:.2f}), "
          f"radius {LOOP_RADIUS}")
    for name, lp in loops.items():
        print(f"  {name}: centre ({lp.centre[0]:7.2f}, {lp.centre[1]:7.2f})  "
              f"encloses: {lp.encloses(*centre)}")

    for ne, ncas in BERRY_LADDER:
        for name, loop in loops.items():
            for n in NPOINTS:
                path = berry_path(name, ne, ncas, n)
                if berry_record_exists(path) and not args.force:
                    print(f"[skip] {os.path.basename(path)}")
                    continue
                print(f"\n=== {name}  CAS({ne},{ncas})/{args.basis}  N={n} ===")
                try:
                    res, trav = run_loop(
                        loop.with_n_points(n),
                        cas=CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
                        cont=ContinuationConfig(),
                        geom_fn=fn,
                        progress=None if args.quiet else print,
                    )
                except Exception as exc:                        # noqa: BLE001
                    print(f"  RUN FAILED: {type(exc).__name__}: {exc}")
                    continue
                print(res.summary())
                save_berry_run(res, trav, path)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["scan", "berry"])
    ap.add_argument("--basis", default=DEFAULT_BASIS)
    ap.add_argument("--centre", nargs=2, type=float, default=None, metavar=("TW", "PYR"))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    return do_scans(args) if args.stage == "scan" else do_berry(args)


if __name__ == "__main__":
    raise SystemExit(main())
