#!/usr/bin/env python
"""Single-update versus fully-converged continuation: the cost/accuracy trade-off.

arXiv:2304.06070 takes **one** parameter update per loop point rather than optimizing to
convergence, on the argument that updates are the currency and, at fixed total budget, more points
with fewer updates each is the better spend. This package converges every point instead. The two
choices are compared here on the same loops.

**Cost** is reported three ways, because the right one depends on where the work runs:

* ``micro`` — CASSCF micro-iterations, the closest analogue of the parameter update the quantum
  algorithm counts, and the currency in which the paper's argument is framed;
* ``macro`` — macro-iterations, each additionally containing a CI diagonalization;
* ``wall`` — seconds, which also carries the per-point overhead (integrals, transformation, the
  mean-field solve) that does *not* shrink when fewer updates are taken per point. That overhead
  is why the classical and quantum cost models differ.

**Accuracy** is reported as two independent figures of merit:

* ``1 - |omega|`` — how far the transported state is from returning to itself after a closed loop,
  measured at a geometry identical to the start. Exactly zero for converged continuation by
  construction; for single-update it measures the accumulated lag, and is the quantity that
  decides whether the endpoint estimator is interpretable at all;
* ``1 - |Pi|`` — the discretization error of the cyclic product.

The Z2 answer itself is recorded separately: it is what both methods actually exist to produce,
and it can be right while both figures of merit are poor.

    python examples/run_stepping_comparison.py
    python examples/run_stepping_comparison.py --cas 6,6 --npoints 9 17 33
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import BENCHMARK_LOOPS, CasConfig, ContinuationConfig
from berrycasscf.berry import analyse
from berrycasscf.continuation import traverse_loop
from berrycasscf.runlog import JobLog
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepping", "stepping_comparison.json")

# (label, macro-iterations allowed per point, require convergence, use the fallback ladder)
MODES = [
    ("converged", 200, True, True),
    ("two-step", 2, False, False),
    ("single-step", 1, False, False),
]
DEFAULT_NPOINTS = [9, 13, 17, 25, 33, 49, 65, 97]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loops", nargs="*", default=["C_x", "C_1"])
    ap.add_argument("--cas", nargs="*", default=["2,2"])
    ap.add_argument("--npoints", nargs="*", type=int, default=DEFAULT_NPOINTS)
    ap.add_argument("--basis", default="sto-3g")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    log = JobLog("stepping_comparison", echo=not args.quiet)

    print(f"{'loop':>5} {'CAS':>9} {'mode':>12} {'N':>4} {'micro':>7} {'macro':>7} "
          f"{'wall(s)':>8} {'1-|om|':>9} {'1-|Pi|':>8} {'min|ovl|':>9}  phase")
    print("-" * 104)
    for spec in args.cas:
        ne, ncas = (int(v) for v in spec.split(","))
        for loop_name in args.loops:
            loop = BENCHMARK_LOOPS[loop_name]
            for label, maxmac, require_conv, fallback in MODES:
                for n in args.npoints:
                    cas = CasConfig(basis=args.basis, ncas=ncas, nelecas=ne,
                                    max_cycle_macro=maxmac)
                    cont = ContinuationConfig(require_converged=require_conv,
                                              use_fallback=fallback)
                    t0 = time.time()
                    try:
                        trav = traverse_loop(loop.with_n_points(n), cas=cas, cont=cont)
                    except Exception as exc:                    # noqa: BLE001
                        log.note(f"{loop_name} {label} N={n} failed: {exc}")
                        continue
                    res = analyse(trav, cont)
                    wall = time.time() - t0
                    om = abs(res.endpoint_estimator) if res.endpoint_estimator else np.nan
                    row = {
                        "loop": loop_name, "cas": [ne, ncas], "mode": label, "n_points": n,
                        "micro": trav.total_micro, "macro": trav.total_macro, "wall": wall,
                        "endpoint_abs": float(om), "product_abs": abs(res.product_estimator),
                        "product": res.product_estimator,
                        "min_abs_overlap": res.min_abs_adjacent_overlap,
                        "phase": res.berry_phase, "status": res.status,
                        "all_converged": bool(trav.all_converged),
                    }
                    rows.append(row)
                    print(f"{loop_name:>5} {'(%d,%d)' % (ne, ncas):>9} {label:>12} {n:>4} "
                          f"{trav.total_micro:>7} {trav.total_macro:>7} {wall:>8.2f} "
                          f"{1-om:>9.2e} {1-abs(res.product_estimator):>8.4f} "
                          f"{res.min_abs_adjacent_overlap:>9.4f}  "
                          f"{res.berry_phase.split()[0]}")
                    log(f"{loop_name} {label} N={n}: micro={trav.total_micro} "
                        f"1-|omega|={1-om:.2e} phase={res.berry_phase}")
                    save_json({"basis": args.basis, "modes": [m[0] for m in MODES],
                               "runs": rows}, OUT)
    log.done()

    # Did the Z2 answer ever come out wrong, in either mode?
    expected = {"C_x": "non-trivial", "C_1": "trivial", "C_2": "trivial"}
    wrong = [r for r in rows if r["phase"].split()[0] != expected.get(r["loop"], "")
             and r["phase"] != "undetermined"]
    print(f"\n{len(rows)} runs. Z2 answer wrong in {len(wrong)} of them.")
    for r in wrong:
        print(f"  {r['loop']} {r['mode']} N={r['n_points']}: {r['phase']}")
    print(f"-> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
