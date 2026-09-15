#!/usr/bin/env python
"""Fulvene, stage 2: Berry phase on loops built around the intersection located in stage 1.

Reads the stage-1 scan, places a CI-enclosing loop on the minimum-gap point and two control
loops displaced along the bond-length coordinate, then runs the continuation on all three.

Intended for the cluster; see slurm/fulvene_berry.sbatch and docs/followup.md.

    python examples/run_fulvene_berry.py --scan results/fulvene/<stage1>.npz --basis 6-31g*
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf import CasConfig, ContinuationConfig, run_loop
from berrycasscf.fulvene import fulvene_geom, loops_around
from berrycasscf.scan import ScanResult
from berrycasscf.store import save_berry_run, berry_record_exists

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "fulvene")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", required=True, help="stage-1 .npz produced by run_fulvene_scan.py")
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--cas", default="6,6")
    ap.add_argument("--npoints", nargs="*", type=int, default=[17, 25])
    ap.add_argument("--radius", nargs=2, type=float, default=[0.06, 15.0],
                    metavar=("DR", "DTHETA"),
                    help="loop radii in (Angstrom, degrees)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    scan = ScanResult.load(args.scan)
    r_ci, t_ci, gap = scan.min_gap_point()
    print(f"Stage-1 CI candidate from {os.path.basename(args.scan)}:")
    print(f"  r = {r_ci:.4f} Angstrom, theta = {t_ci:.2f} deg, gap = {gap:.6f} Ha")
    if gap > 5e-3:
        print(f"  WARNING: the stage-1 minimum gap is {gap:.6f} Ha, which is large for a "
              "conical intersection. The loops below may not enclose anything; widen the "
              "search region or enlarge the active space before trusting stage 2.")

    ne, ncas = (int(x) for x in args.cas.split(","))
    loops = loops_around((r_ci, t_ci), radius=tuple(args.radius))
    os.makedirs(RESULT_DIR, exist_ok=True)
    cas = CasConfig(basis=args.basis, ncas=ncas, nelecas=ne)
    cont = ContinuationConfig()

    print(f"\nLoops (radius {args.radius[0]} A, {args.radius[1]} deg):")
    for name, lp in loops.items():
        print(f"  {name}: centre (r, theta) = ({lp.centre[0]:.4f}, {lp.centre[1]:.2f})")

    for name, loop in loops.items():
        for n in args.npoints:
            tag = f"fulvene_{name}_cas{ne}-{ncas}_{args.basis.replace('*','s')}_N{n}"
            path = os.path.join(RESULT_DIR, tag + ".json")
            if berry_record_exists(path) and not args.force:
                print(f"[skip] {os.path.basename(path)}")
                continue
            print(f"\n=== {name}  CAS({ne},{ncas})/{args.basis}  N={n} ===")
            try:
                res, trav = run_loop(loop.with_n_points(n), cas=cas, cont=cont,
                                     geom_fn=fulvene_geom,
                                     progress=None if args.quiet else print)
            except Exception as exc:                        # noqa: BLE001
                print(f"  RUN FAILED: {type(exc).__name__}: {exc}")
                continue
            print(res.summary())
            save_berry_run(res, trav, path)
            print(f"  -> {os.path.relpath(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
