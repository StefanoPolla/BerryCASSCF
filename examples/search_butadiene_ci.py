#!/usr/bin/env python
"""Locate an S1/S0 intersection in butadiene, before committing to a plane.

Nothing downstream assumes where the intersection is. This scans several candidate planes of
the four rigid coordinates at a cheap active space, reports where each one's gap is smallest,
and leaves the choice of plane to the evidence.

    python examples/search_butadiene_ci.py                    # all candidate planes, coarse
    python examples/search_butadiene_ci.py --planes tw_pyr --grid 15 15
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ScanConfig, scan_gap
from berrycasscf.butadiene import (
    DEFAULT_BASIS,
    DEFAULT_NCAS,
    DEFAULT_NELECAS,
    SEARCH_PLANES,
    plane_geom_fn,
)
from berrycasscf.scan import ScanResult
from berrycasscf.runlog import JobLog

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "butadiene")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--planes", nargs="*", default=list(SEARCH_PLANES))
    ap.add_argument("--grid", nargs=2, type=int, default=[13, 13], metavar=("NX", "NY"))
    ap.add_argument("--basis", default=DEFAULT_BASIS)
    ap.add_argument("--cas", default=f"{DEFAULT_NELECAS},{DEFAULT_NCAS}")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    ne, ncas = (int(x) for x in args.cas.split(","))
    nx, ny = args.grid
    os.makedirs(RESULT_DIR, exist_ok=True)

    summary = []
    for name in args.planes:
        (cx, cy), fixed, region = SEARCH_PLANES[name]
        path = os.path.join(RESULT_DIR, f"search_{name}_cas{ne}-{ncas}_{nx}x{ny}.npz")
        if os.path.exists(path) and not args.force:
            res = ScanResult.load(path)
            if np.isfinite(res.e_states).all():
                print(f"SKIP: {os.path.basename(path)}")
                summary.append((name, cx, cy, fixed, res))
                continue
        lo_x, hi_x, lo_y, hi_y = region.bounding_box()
        print(f"\n=== plane {name}: {cx} in [{lo_x:.0f}, {hi_x:.0f}], "
              f"{cy} in [{lo_y:.0f}, {hi_y:.0f}], fixed {fixed} ===")
        print(f"    SA-CASSCF({ne},{ncas})/{args.basis}, {nx}x{ny} points")
        res = scan_gap(
            region,
            cas=CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
            scan=ScanConfig(n_alpha=nx, n_phi=ny, margin=0.0),
            geom_fn=plane_geom_fn(cx, cy, **fixed),
            progress=None if args.quiet else print,
            checkpoint=path,
        )
        res.save(path)
        summary.append((name, cx, cy, fixed, res))

    # Reference S0 energy at the undistorted planar geometry, so that strain is visible.
    from berrycasscf.butadiene import butadiene_coords, _to_string
    from berrycasscf.casscf import apply_singlet_constraint, build_mol, run_rhf
    from pyscf import mcscf as _mcscf

    mol0 = build_mol(_to_string(butadiene_coords()), args.basis)
    mc0 = _mcscf.CASSCF(run_rhf(mol0), ncas, ne)
    apply_singlet_constraint(mc0)
    mc0 = mc0.state_average_([0.5, 0.5])
    mc0.verbose = 0
    mc0.kernel()
    e_ref = float(mc0.e_states[0])

    HARTREE_KCAL = 627.509

    print("\n" + "=" * 96)
    print(f"{'plane':>10} {'min gap (mHa)':>14} {'at':>26} {'strain (kcal/mol)':>18}  fixed")
    print("-" * 96)
    for name, cx, cy, fixed, res in summary:
        x, y, g = res.min_gap_point()
        i = int(np.argmin(np.abs(res.alphas - x)))
        j = int(np.argmin(np.abs(res.phis - y)))
        strain = (res.e_states[i, j, 0] - e_ref) * HARTREE_KCAL
        print(f"{name:>10} {g*1e3:14.3f} {f'{cx}={x:.1f}, {cy}={y:.1f}':>26} "
              f"{strain:18.1f}  {fixed}")

    print("\nA gap below ~5 mHa is worth refining; above ~20 mHa the plane probably does not")
    print("contain an intersection. **Check the strain column**: a grid minimum can sit at an")
    print("absurdly distorted geometry (a methylene folded back onto its own C-C bond, say)")
    print("where the two states are degenerate but the structure is chemically meaningless.")
    print("Anything much above ~200 kcal/mol deserves inspection before it is built on.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
