#!/usr/bin/env python
"""Berry-phase continuation for formaldimine: all loops x active spaces x discretizations.

Every completed (loop, CAS, N) combination is written to ``results/berry/`` as JSON and is
skipped on a re-run, so the sweep is restartable.

    python examples/run_formaldimine_berry.py                 # default sweep
    python examples/run_formaldimine_berry.py --loops C_x --cas 2,2 --npoints 9 25
    python examples/run_formaldimine_berry.py --force         # recompute everything
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf import BENCHMARK_LOOPS, CasConfig, ContinuationConfig, run_loop
from berrycasscf.store import save_berry_run, berry_record_exists
from berrycasscf.runlog import JobLog

RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "berry")

# The sweep: active spaces tested, smallest chemically sensible first.
DEFAULT_CAS = [(2, 2), (4, 4), (6, 6)]
# Discretizations, for the robustness-to-N half of the stability criterion.
DEFAULT_NPOINTS = [9, 17, 25, 41]


def record_path(loop_name: str, ne: int, ncas: int, n: int) -> str:
    return os.path.join(RESULT_DIR, f"{loop_name}_cas{ne}-{ncas}_N{n}.json")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loops", nargs="*", default=list(BENCHMARK_LOOPS),
                    help=f"loop names (default: all of {list(BENCHMARK_LOOPS)})")
    ap.add_argument("--cas", nargs="*", default=None,
                    help="active spaces as 'nelec,norb' (default: 2,2 4,4 6,6)")
    ap.add_argument("--npoints", nargs="*", type=int, default=DEFAULT_NPOINTS)
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--force", action="store_true", help="recompute existing records")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    cas_list = DEFAULT_CAS if args.cas is None else [
        tuple(int(x) for x in spec.split(",")) for spec in args.cas
    ]
    os.makedirs(RESULT_DIR, exist_ok=True)
    cont = ContinuationConfig()

    for loop_name in args.loops:
        loop = BENCHMARK_LOOPS[loop_name]
        for ne, ncas in cas_list:
            for n in args.npoints:
                path = record_path(loop_name, ne, ncas, n)
                if berry_record_exists(path) and not args.force:
                    print(f"[skip] {os.path.basename(path)}")
                    continue
                cas = CasConfig(basis=args.basis, ncas=ncas, nelecas=ne)
                print(f"\n=== {loop_name}  {cas.cas_label}/{args.basis}  N={n} ===")
                try:
                    res, trav = run_loop(
                        loop.with_n_points(n), cas=cas, cont=cont,
                        progress=None if args.quiet else print,
                    )
                except Exception as exc:                      # noqa: BLE001
                    print(f"  RUN FAILED: {type(exc).__name__}: {exc}")
                    continue
                print(res.summary())
                save_berry_run(res, trav, path)
                print(f"  -> {os.path.relpath(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
