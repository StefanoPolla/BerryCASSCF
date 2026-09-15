#!/usr/bin/env python
"""Exact (FCI) S0/S1 gap reference in the same basis, to calibrate the truncated active spaces.

This is the same quantity the gap map of arXiv:2304.06070 Fig. 1a reports. Expensive:
~30 s per point for formaldimine/STO-3G. Checkpoints after every row.

    python examples/run_fci_reference.py --line          # alpha line at phi = 89.9 (default)
    python examples/run_fci_reference.py --grid 9 5      # coarse 2D map
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf import CasConfig, ScanConfig, SCAN_REGIONS
from berrycasscf.geometry import Loop
from berrycasscf.scan import scan_gap_fci, ScanResult

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "scan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loop", default="C_x", help="region whose bounding box is scanned")
    ap.add_argument("--grid", nargs=2, type=int, default=None, metavar=("N_ALPHA", "N_PHI"))
    ap.add_argument("--line", action="store_true",
                    help="scan alpha only, at the loop's centre phi (default)")
    ap.add_argument("--margin", type=float, default=0.0)
    ap.add_argument("--basis", default="sto-3g")
    args = ap.parse_args()

    region = SCAN_REGIONS[args.loop]
    if args.grid is None or args.line:
        na, npi = 17, 1
        # a degenerate 'region' one point wide in phi, centred on the loop's phi
        region = Loop(region.name, region.centre, (region.radius[0], 0.0),
                      region.phase, region.n_points)
    else:
        na, npi = args.grid

    os.makedirs(RESULT_DIR, exist_ok=True)
    path = os.path.join(RESULT_DIR, f"{region.name}_fci_{na}x{npi}.npz")
    cas = CasConfig(basis=args.basis)
    scan = ScanConfig(n_alpha=na, n_phi=npi, margin=args.margin)

    print(f"=== FCI reference: {region.name}, {na}x{npi} points, basis {args.basis} ===")
    print(f"    (~30 s/point; {na*npi} points)")
    res = scan_gap_fci(region, cas=cas, scan=scan, progress=print, checkpoint=path)
    res.save(path)
    print()
    print(res.summary())
    print(f"  -> {os.path.relpath(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
