#!/usr/bin/env python
"""Confirm the fitted intersection positions by minimizing the gap directly.

The positions quoted elsewhere come from a *model* — a parabola in ``gap**2``, exact for an ideal
cone. This starts from each fitted position and minimizes the gap with a derivative-free search,
so the fit is checked rather than believed. Two outcomes are informative:

* the search barely moves and finds a similar gap — the fit was sound;
* the search moves, or finds a markedly lower gap — the fit was biased, and by how much.

Each evaluation is a full cold SA-CASSCF, so the budget is deliberately small: this refines a good
guess over a few tens of solves. Expensive rungs are skipped by default for that reason.

    python examples/run_gap_minimum_search.py                       # CAS(2,2)..(10,10)
    python examples/run_gap_minimum_search.py --cas 6,6 --budget 60
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig
from berrycasscf.butadiene import plane_geom_fn
from berrycasscf.minsearch import minimize_gap
from berrycasscf.refine import cone_apex
from berrycasscf.runlog import JobLog
from berrycasscf.scan import ScanResult
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(ROOT, "results", "butadiene")
OUT = os.path.join(RESULT_DIR, "gap_minimum_search.json")

# The cheap end of the ladder. CAS(12,12) is excluded by default: at >2 min per solve a 40-solve
# search is well over an hour for one rung.
DEFAULT_CAS = ["2,2", "4,4", "6,6", "8,8", "10,10"]


def fitted_start(ne: int, ncas: int) -> tuple[tuple[float, float], float] | None:
    """The cone-fit position from the fine 1-deg cut, with the tw it was taken at."""
    path = os.path.join(RESULT_DIR, f"butadiene_refinepyr_cas{ne}-{ncas}.npz")
    if not os.path.exists(path):
        return None
    r = ScanResult.load(path)
    if not np.isfinite(r.e_states).all():
        return None
    fit = cone_apex(r.phis, r.gap[0] * 1e3, window=None)
    return (float(r.alphas[0]), fit.position), fit.closest_approach


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cas", nargs="*", default=DEFAULT_CAS)
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--budget", type=int, default=40, help="max SA-CASSCF solves per rung")
    ap.add_argument("--step", nargs=2, type=float, default=[1.0, 1.0],
                    metavar=("D_TW", "D_PYR"), help="initial simplex size")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    fn = plane_geom_fn("tw", "pyr", tc=0.0, bend=0.0)
    rows = []
    print(f"{'CAS':>10} {'fitted (tw, pyr)':>22} {'fit gap':>9} "
          f"{'searched (tw, pyr)':>22} {'gap':>9} {'moved':>7} {'evals':>6}")
    print("-" * 92)
    for spec in args.cas:
        ne, ncas = (int(v) for v in spec.split(","))
        seed = fitted_start(ne, ncas)
        if seed is None:
            print(f"{'CAS(%d,%d)' % (ne, ncas):>10}   no fine cut; run the refine stage first")
            continue
        start, fit_gap = seed
        log = JobLog(f"butadiene_gapsearch_cas{ne}-{ncas}", echo=not args.quiet)
        log.note(f"start from cone fit at (tw={start[0]:.3f}, pyr={start[1]:.3f}), "
                 f"fitted closest approach {fit_gap:.4f} mHa")
        t0 = time.time()
        res = minimize_gap(fn, CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
                           start=start, step=tuple(args.step),
                           max_evaluations=args.budget, progress=log)
        log.done()
        rows.append({
            "cas": [ne, ncas], "start": list(res.start), "start_gap": res.start_gap,
            "fitted_closest_approach": fit_gap,
            "found": [res.x, res.y], "gap": res.gap, "moved": res.moved,
            "n_evaluations": res.n_evaluations, "converged": res.converged,
            "message": res.message, "wall_time": time.time() - t0,
            "trace": res.trace,
        })
        print(f"{'CAS(%d,%d)' % (ne, ncas):>10} "
              f"{'(%.2f, %.2f)' % res.start:>22} {res.start_gap:9.4f} "
              f"{'(%.2f, %.2f)' % (res.x, res.y):>22} {res.gap:9.4f} "
              f"{res.moved:7.3f} {res.n_evaluations:6d}")
        save_json({"basis": args.basis, "budget": args.budget, "runs": rows}, OUT)

    if rows:
        worst = max(rows, key=lambda r: r["moved"])
        print(f"\nLargest correction to a fitted position: {worst['moved']:.3f} deg "
              f"at CAS({worst['cas'][0]},{worst['cas'][1]}).")
        deeper = [r for r in rows if r["gap"] < 0.5 * r["start_gap"]]
        if deeper:
            print(f"{len(deeper)} rung(s) reached a gap below half the fitted value — "
                  "those fits were optimistic about how close the cut passes.")
        print(f"-> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
