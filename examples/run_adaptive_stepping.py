#!/usr/bin/env python
"""Adaptive versus uniform discretization of a loop: does step control pay?

`docs/todo.md` §1. A uniform walk has to be fine enough for the hardest arc of the loop and
so oversamples everywhere else. The adaptive walk steers its step size by the continuity it
measures (`berrycasscf.adaptive`). This driver asks whether that actually buys anything, on
loops whose answers are already known.

**The comparison is deliberately a family against a family.** Adaptive stepping has knobs
(`d_max`, `target_mismatch`), so comparing one adaptive setting against one N would prove
nothing -- either method can be made to look good by choosing its parameter well. Both are
swept and compared on the *frontier*: for each method, the cheapest run that passes every
check, and the quality reached at a given cost.

Cost is counted in CASSCF **micro-iterations** -- parameter updates, the currency of
arXiv:2304.06070 -- because wall time carries per-point overhead that differs between
machines. Rejected adaptive trials are charged in full; hiding them would flatter the method.

    python examples/run_adaptive_stepping.py formaldimine
    python examples/run_adaptive_stepping.py butadiene
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ContinuationConfig, AdaptiveConfig
from berrycasscf.adaptive import traverse_loop_adaptive
from berrycasscf.berry import analyse
from berrycasscf.continuation import traverse_loop
from berrycasscf.geometry import LOOP_CI, LOOP_CONTROL_UPPER, formaldimine_geom
from berrycasscf.butadiene import butadiene_geom
from berrycasscf.geometry import Loop
from berrycasscf.runlog import JobLog
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

UNIFORM_N = [9, 13, 17, 25, 33, 49]
# Sweep the two knobs that matter. d_max caps the coarsest step (so it sets the cost on an
# easy loop); target_mismatch sets how much margin the walk carries on a hard one.
ADAPTIVE_SETTINGS = [
    {"d_max": 0.20, "target_mismatch": 0.02},
    {"d_max": 0.10, "target_mismatch": 0.02},
    {"d_max": 0.10, "target_mismatch": 0.005},
    {"d_max": 0.05, "target_mismatch": 0.02},
]

SYSTEMS = {
    "formaldimine": {
        "geom_fn": formaldimine_geom,
        "basis": "sto-3g",
        "cas": [(2, 2), (6, 6)],
        "loops": {"C_x": LOOP_CI, "C_2": LOOP_CONTROL_UPPER},
        "expected": {"C_x": "non-trivial (pi)", "C_2": "trivial (0)"},
    },
    "butadiene": {
        "geom_fn": butadiene_geom,
        "basis": "6-31g*",
        "cas": [(2, 2)],
        "loops": {
            "B_x": Loop("B_x", (90.0, 101.85321091497578), (12.0, 18.0)),
            "B_2": Loop("B_2", (90.0, 101.85321091497578 + 45.0), (12.0, 18.0)),
        },
        "expected": {"B_x": "non-trivial (pi)", "B_2": None},
    },
}


def row(res, trav, method: str, setting: str) -> dict:
    ad = trav.adaptive or {}
    return {
        "method": method,
        "setting": setting,
        "loop": res.loop_name,
        "cas": res.cas_label,
        "n_points": trav.n_points,
        "n_rejected": ad.get("n_rejected", 0),
        "micro": trav.total_micro,
        "macro": trav.total_macro,
        "wall": trav.wall_time,
        "min_overlap": float(res.min_abs_adjacent_overlap),
        "one_minus_abs_pi": float(1.0 - abs(res.product_estimator)),
        "endpoint": res.endpoint_estimator,
        "status": res.status,
        "phase": res.berry_phase,
        "closed": ad.get("closed", True),
        "hit_step_floor": ad.get("hit_step_floor", False),
        "step_min": float(np.min([e["d"] for e in ad["events"] if e["accepted"]]))
                    if ad.get("events") else None,
        "step_max": float(np.max([e["d"] for e in ad["events"] if e["accepted"]]))
                    if ad.get("events") else None,
        "messages": res.messages,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    spec = SYSTEMS[args.system]
    out = os.path.join(ROOT, "results", "adaptive", f"{args.system}_stepping.json")
    if os.path.exists(out) and not args.force:
        print(f"[skip] {os.path.relpath(out, ROOT)} exists; --force to redo")
        return 0

    cont = ContinuationConfig()
    rows: list[dict] = []
    total = len(spec["cas"]) * len(spec["loops"]) * (len(UNIFORM_N) + len(ADAPTIVE_SETTINGS))
    log = JobLog(f"adaptive_{args.system}", total=total)

    for ne, ncas in spec["cas"]:
        cas = CasConfig(basis=spec["basis"], ncas=ncas, nelecas=ne)
        for name, loop in spec["loops"].items():
            for n in UNIFORM_N:
                trav = traverse_loop(loop.with_n_points(n), cas=cas, cont=cont,
                                     geom_fn=spec["geom_fn"])
                r = row(analyse(trav, cont), trav, "uniform", f"N={n}")
                rows.append(r)
                log(f"{name} CAS({ne},{ncas}) uniform N={n}: {r['status']} "
                    f"{r['phase']} micro={r['micro']}")
            for setting in ADAPTIVE_SETTINGS:
                adcfg = AdaptiveConfig(**setting)
                trav = traverse_loop_adaptive(loop, cas=cas, cont=cont, adaptive=adcfg,
                                              geom_fn=spec["geom_fn"])
                label = f"dmax={setting['d_max']},tgt={setting['target_mismatch']}"
                r = row(analyse(trav, cont), trav, "adaptive", label)
                rows.append(r)
                log(f"{name} CAS({ne},{ncas}) adaptive {label}: {r['status']} "
                    f"{r['phase']} micro={r['micro']} pts={r['n_points']} "
                    f"rej={r['n_rejected']}")
    log.done()

    print(f"\n{'loop':>5} {'CAS':>10} {'method':>9} {'setting':>26} {'pts':>4} {'rej':>4} "
          f"{'micro':>7} {'wall(s)':>8} {'minovl':>7} {'status':>7}  phase")
    print("-" * 110)
    for r in rows:
        print(f"{r['loop']:>5} {r['cas']:>10} {r['method']:>9} {r['setting']:>26} "
              f"{r['n_points']:>4} {r['n_rejected']:>4} {r['micro']:>7} {r['wall']:>8.1f} "
              f"{r['min_overlap']:>7.4f} {r['status']:>7}  {r['phase']}")

    # --- the frontier: cheapest run of each method that passes every check --------------
    print("\nCheapest TRUSTWORTHY run of each method (passes every check), by micro-iterations:")
    print(f"{'loop':>5} {'CAS':>10} {'uniform':>28} {'adaptive':>34} {'saving':>8}")
    print("-" * 92)
    summary = []
    for cas_label in sorted({r["cas"] for r in rows}):
        for loop in sorted({r["loop"] for r in rows}):
            sel = [r for r in rows if r["cas"] == cas_label and r["loop"] == loop
                   and r["status"] == "OK"]
            u = min((r for r in sel if r["method"] == "uniform"),
                    key=lambda r: r["micro"], default=None)
            a = min((r for r in sel if r["method"] == "adaptive"),
                    key=lambda r: r["micro"], default=None)
            us = f"{u['setting']} {u['micro']} micro" if u else "none passes"
            as_ = f"{a['setting']} {a['micro']} micro" if a else "none passes"
            saving = f"{u['micro']/a['micro']:.2f}x" if (u and a) else "--"
            print(f"{loop:>5} {cas_label:>10} {us:>28} {as_:>34} {saving:>8}")
            summary.append({"loop": loop, "cas": cas_label,
                            "uniform": u, "adaptive": a, "saving": saving})

    # --- the fairer comparison: cost to reach a given quality -------------------------
    # "Cheapest run that passes" rewards whichever method happens to skate closest to the
    # threshold, and §5's calibration shows runs near the threshold are the dangerous ones.
    # Matching on achieved continuity instead asks what each method costs for the SAME
    # margin of safety.
    print("\nCost (micro-iterations) to reach a given worst-case adjacent overlap:")
    print(f"{'loop':>5} {'CAS':>10} {'target':>8} {'uniform':>18} {'adaptive':>18} {'ratio':>7}")
    print("-" * 72)
    quality = []
    for cas_label in sorted({r["cas"] for r in rows}):
        for loop in sorted({r["loop"] for r in rows}):
            for target in (0.90, 0.95, 0.98):
                sel = [r for r in rows if r["cas"] == cas_label and r["loop"] == loop
                       and r["status"] == "OK" and r["min_overlap"] >= target]
                u = min((r for r in sel if r["method"] == "uniform"),
                        key=lambda r: r["micro"], default=None)
                a_ = min((r for r in sel if r["method"] == "adaptive"),
                         key=lambda r: r["micro"], default=None)
                ratio = f"{u['micro']/a_['micro']:.2f}x" if (u and a_) else "--"
                print(f"{loop:>5} {cas_label:>10} {target:>8.2f} "
                      f"{(str(u['micro']) + ' (' + u['setting'] + ')') if u else 'none':>18} "
                      f"{(str(a_['micro']) + ' (' + a_['setting'].split(',')[0] + ')') if a_ else 'none':>18} "
                      f"{ratio:>7}")
                quality.append({"loop": loop, "cas": cas_label, "target": target,
                                "uniform_micro": u["micro"] if u else None,
                                "adaptive_micro": a_["micro"] if a_ else None,
                                "uniform_setting": u["setting"] if u else None,
                                "adaptive_setting": a_["setting"] if a_ else None})

    save_json({"system": args.system, "rows": rows, "frontier": summary,
               "matched_quality": quality,
               "uniform_n": UNIFORM_N, "adaptive_settings": ADAPTIVE_SETTINGS}, out)
    print(f"\n-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
