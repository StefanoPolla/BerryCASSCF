#!/usr/bin/env python
"""Does loop transport fail loudly, and what does that conservatism cost?

A method that returns a wrong answer confidently is worse than one that refuses. This script
probes the refusal behaviour directly, by shrinking the loop toward a degeneracy and varying the
discretization independently.

It is written to be able to embarrass the method, not to flatter it. Three things are recorded
for every run, including the ones that are refused:

* the phase it **would** have reported had the checks not run — a refused run that would have been
  right is a cost of the conservatism, and is reported as such;
* the minimum adjacent overlap, against the 0.80 continuity threshold;
* whether refinement in N rescues a loop that failed at coarse N, which distinguishes "too coarse"
  from "genuinely too close to the seam".

A displaced loop is included as a control: it encloses nothing but passes near the seam, so it
tests whether the checks refuse indiscriminately or only when continuity is actually lost.

    python examples/run_failure_modes.py
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, run_loop
from berrycasscf.butadiene import plane_geom_fn
from berrycasscf.geometry import Loop
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "butadiene", "failure_modes.json")

CENTRE = (90.0, 101.85)        # the enclosing loop's centre, from the CAS(12,12) reference
# Loops that pass close to a degeneracy converge slowly (the solver's fallback ladder retries
# several strategies per point), so the sweep is kept deliberately small. N = 61 is spent only
# where it answers something: on the radii that were refused at coarser N, to separate
# "discretization too coarse" from "genuinely too close to the seam".
RADII = [6.0, 8.0, 12.0, 18.0]
NPOINTS = [15, 31]
RESCUE_N = 61
RESCUE_RADII = [6.0, 8.0]
CAS = (8, 8)                   # the rung whose gap-scan minimum leaves the loop


def classify(res) -> str:
    """What the run reports, and what it would have said if the checks had not run."""
    would = "pi" if res.product_estimator < 0 else "0"
    if res.status == "OK":
        return f"reported {would}"
    return f"REFUSED (would have said {would})"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--basis", default="6-31g*")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    fn = plane_geom_fn("tw", "pyr", tc=0.0, bend=0.0)
    cas = CasConfig(basis=args.basis, ncas=CAS[1], nelecas=CAS[0])
    rows = []

    def checkpoint():
        save_json({"centre": CENTRE, "cas": list(CAS), "basis": args.basis, "runs": rows}, OUT)

    print(f"Loop transport at CAS{CAS}/{args.basis}, centre (tw, pyr) = {CENTRE}")
    print(f"{'radius':>7} {'reach':>7} {'N':>4} {'min|ovl|':>9} {'product':>9} "
          f"{'|endpoint|':>11}  outcome")
    print("-" * 78)
    for rp in RADII:
        for n in NPOINTS:
            loop = Loop("probe", CENTRE, (12.0, rp), n_points=n)
            t0 = time.time()
            try:
                res, _ = run_loop(loop, cas=cas, geom_fn=fn,
                                  progress=None if args.quiet else print)
            except Exception as exc:                       # noqa: BLE001
                print(f"{rp:7.0f} {CENTRE[1]+rp:7.1f} {n:4d}   RUN FAILED: {exc}")
                continue
            rows.append({
                "kind": "enclosing", "radius_pyr": rp, "n_points": n,
                "min_abs_overlap": res.min_abs_adjacent_overlap,
                "product": res.product_estimator,
                "endpoint": res.endpoint_estimator,
                "status": res.status, "phase": res.berry_phase,
                "would_have_said": "pi" if res.product_estimator < 0 else "0",
                "messages": res.messages, "wall_time": time.time() - t0,
            })
            print(f"{rp:7.0f} {CENTRE[1]+rp:7.1f} {n:4d} {res.min_abs_adjacent_overlap:9.3f} "
                  f"{res.product_estimator:+9.4f} {abs(res.endpoint_estimator):11.4f}  "
                  f"{classify(res)}", flush=True)
            checkpoint()

    # Spend the expensive discretization only where a refusal has to be explained.
    print("\nRescue attempts: refine N on the radii that were refused")
    for rp in RESCUE_RADII:
        loop = Loop("probe", CENTRE, (12.0, rp), n_points=RESCUE_N)
        res, _ = run_loop(loop, cas=cas, geom_fn=fn,
                          progress=None if args.quiet else print)
        rows.append({
            "kind": "enclosing", "radius_pyr": rp, "n_points": RESCUE_N,
            "min_abs_overlap": res.min_abs_adjacent_overlap,
            "product": res.product_estimator, "endpoint": res.endpoint_estimator,
            "status": res.status, "phase": res.berry_phase,
            "would_have_said": "pi" if res.product_estimator < 0 else "0",
            "messages": res.messages,
        })
        print(f"{rp:7.0f} {CENTRE[1]+rp:7.1f} {RESCUE_N:4d} "
              f"{res.min_abs_adjacent_overlap:9.3f} {res.product_estimator:+9.4f} "
              f"{abs(res.endpoint_estimator):11.4f}  {classify(res)}", flush=True)
        checkpoint()

    # Control: a loop that encloses nothing but passes near the seam. If the checks refuse this
    # too, they are refusing proximity rather than lost continuity.
    print("\nControl: displaced loop, encloses nothing, still passes near the seam")
    for n in NPOINTS:
        loop = Loop("displaced", (CENTRE[0], CENTRE[1] - 26.0), (12.0, 10.0), n_points=n)  # noqa: E501
        res, _ = run_loop(loop, cas=cas, geom_fn=fn, progress=None if args.quiet else print)
        rows.append({
            "kind": "displaced_control", "radius_pyr": 10.0, "n_points": n,
            "min_abs_overlap": res.min_abs_adjacent_overlap,
            "product": res.product_estimator, "endpoint": res.endpoint_estimator,
            "status": res.status, "phase": res.berry_phase,
            "would_have_said": "pi" if res.product_estimator < 0 else "0",
            "messages": res.messages,
        })
        print(f"{'—':>7} {'—':>7} {n:4d} {res.min_abs_adjacent_overlap:9.3f} "
              f"{res.product_estimator:+9.4f} {abs(res.endpoint_estimator):11.4f}  "
              f"{classify(res)}", flush=True)
        checkpoint()

    checkpoint()
    print(f"\n-> {os.path.relpath(OUT, ROOT)}")

    refused = [r for r in rows if r["status"] != "OK"]
    rescued = [r for r in rows if r["status"] == "OK" and any(
        q["status"] != "OK" and q["radius_pyr"] == r["radius_pyr"] and q["n_points"] < r["n_points"]
        for q in rows)]
    print(f"\n{len(refused)} of {len(rows)} runs refused. "
          f"{len(rescued)} were rescued by refining N alone.")
    agree = [r for r in refused if r["would_have_said"] == "pi"]
    print(f"Of the refused runs, {len(agree)} would have given the same sign as the "
          "well-resolved ones — the cost of erring toward refusal.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
