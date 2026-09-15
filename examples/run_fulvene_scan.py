#!/usr/bin/env python
"""Fulvene, stage 1: locate the S1/S0 intersection in the (r, theta) plane.

Equal-weight two-state SA-CASSCF over the exocyclic-bond-length / methylene-torsion plane.
The conical intersection position is *not* assumed anywhere: this scan finds it, and stage 2
(``run_fulvene_berry.py``) builds the loops around whatever is found.

Intended for the cluster; see slurm/fulvene_scan.sbatch and docs/followup.md.

    python examples/run_fulvene_scan.py --basis 6-31g* --grid 21 21 --cas 6,6
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf import CasConfig, ScanConfig, scan_gap
from berrycasscf.fulvene import SEARCH_REGION, fulvene_geom, DEFAULT_NCAS, DEFAULT_NELECAS

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "fulvene")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--cas", default="6,6", help="'nelec,norb' (default: the full pi manifold)")
    ap.add_argument("--grid", nargs=2, type=int, default=[21, 21],
                    metavar=("N_R", "N_THETA"))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    ne, ncas = (int(x) for x in args.cas.split(","))
    nr, nt = args.grid
    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = f"fulvene_scan_cas{ne}-{ncas}_{args.basis.replace('*','s')}_{nr}x{nt}"
    path = os.path.join(RESULT_DIR, tag + ".npz")

    cas = CasConfig(basis=args.basis, ncas=ncas, nelecas=ne)
    # ScanConfig's fields are named for formaldimine's (alpha, phi); here they carry
    # (r, theta). The scan machinery is agnostic to what the two coordinates mean.
    scan = ScanConfig(n_alpha=nr, n_phi=nt, margin=0.0)

    print(f"=== fulvene stage 1: SA-CASSCF({ne},{ncas})/{args.basis}, {nr}x{nt} grid ===")
    lo_r, hi_r, lo_t, hi_t = SEARCH_REGION.bounding_box()
    print(f"    r     in [{lo_r:.3f}, {hi_r:.3f}] Angstrom")
    print(f"    theta in [{lo_t:.1f}, {hi_t:.1f}] degrees")

    res = scan_gap(SEARCH_REGION, cas=cas, scan=scan, geom_fn=fulvene_geom,
                   progress=None if args.quiet else print, checkpoint=path)
    res.save(path)
    print()
    print(res.summary().replace("alpha", "r    ").replace("phi", "theta"))
    r_ci, t_ci, gap = res.min_gap_point()
    print(f"\n  CI candidate: r = {r_ci:.4f} Angstrom, theta = {t_ci:.2f} deg, "
          f"gap = {gap:.6f} Ha")
    print(f"  -> {os.path.relpath(path)}")
    print("\nNext: python examples/run_fulvene_berry.py "
          f"--scan {os.path.relpath(path)} --basis {args.basis} --cas {args.cas}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
