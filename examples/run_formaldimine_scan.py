#!/usr/bin/env python
"""State-averaged CASSCF two-dimensional gap scan over the region enclosed by a loop.

Writes ``results/scan/<loop>_cas<ne>-<ncas>_<na>x<np>.npz``. The scan checkpoints after every
grid row, so an interrupted run resumes where it stopped.

    python examples/run_formaldimine_scan.py --loops C_x --cas 4,4 --grid 21 21
    python examples/run_formaldimine_scan.py                 # default sweep
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf import SCAN_REGIONS, CasConfig, ScanConfig, scan_gap
from berrycasscf.scan import ScanResult
from berrycasscf.runlog import JobLog

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "scan")


def record_path(loop_name: str, ne: int, ncas: int, na: int, npi: int) -> str:
    return os.path.join(RESULT_DIR, f"{loop_name}_cas{ne}-{ncas}_{na}x{npi}.npz")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loops", nargs="*", default=["C_x", "C_1", "C_2"],
                    help=f"scan regions: {list(SCAN_REGIONS)} "
                         "('overview' covers all three loops in one map)")
    ap.add_argument("--cas", nargs="*", default=["2,2", "4,4"],
                    help="active spaces as 'nelec,norb'")
    ap.add_argument("--grid", nargs=2, type=int, default=[21, 21],
                    metavar=("N_ALPHA", "N_PHI"))
    ap.add_argument("--margin", type=float, default=1.0,
                    help="degrees of padding around the loop bounding box")
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    os.makedirs(RESULT_DIR, exist_ok=True)
    na, npi = args.grid

    for loop_name in args.loops:
        loop = SCAN_REGIONS[loop_name]
        for spec in args.cas:
            ne, ncas = (int(x) for x in spec.split(","))
            path = record_path(loop_name, ne, ncas, na, npi)
            if os.path.exists(path) and not args.force:
                res = ScanResult.load(path)
                if res.converged.all():
                    print(f"[skip] {os.path.basename(path)}")
                    print(res.summary())
                    continue
            cas = CasConfig(basis=args.basis, ncas=ncas, nelecas=ne)
            scan = ScanConfig(n_alpha=na, n_phi=npi, margin=args.margin)
            print(f"\n=== scan {loop_name}  {cas.cas_label}/{args.basis}  {na}x{npi} ===")
            res = scan_gap(loop, cas=cas, scan=scan,
                           progress=None if args.quiet else print,
                           checkpoint=path)
            res.save(path)
            print(res.summary())
            print(f"  -> {os.path.relpath(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
