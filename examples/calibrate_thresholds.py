#!/usr/bin/env python
"""Calibrate the loop-transport acceptance thresholds against every run on disk.

`docs/todo.md` §5 asks for two changes -- promote "continuation chain broken" from a
message to a check, and revisit the 0.80 adjacent-overlap threshold -- and insists they be
*calibrated* rather than guessed. A threshold picked after seeing which answer it gives is
worthless; a threshold picked by how well it separates self-consistent runs from
self-contradictory ones is not.

The calibration set is every Berry-phase record in ``results/``. Runs are grouped into
**questions**: same system, same loop (centre and radius), same active space, differing
only in the discretization N. Physics does not depend on N, so within a question every
trustworthy run must report the same phase. That gives a label without knowing any right
answer:

    contradiction  : a question where two or more runs pass the checks and *disagree*.
                     At least one accepted answer is wrong.
    lone dissenter : a question where exactly ONE run passes, while some *refused* run of
                     the same question reports the opposite sign. The multi-N stability
                     criterion cannot fire, so a user running that single discretization
                     gets a confident answer with nothing to warn them. This is the
                     radius-6 case of `docs/findings.md` §5, and it is the error mode the
                     corpus actually contains -- there are no contradictions at any
                     threshold.
    decidable      : a question where two or more runs pass and agree -- the stability
                     criterion can be met, so the question gets an answer.

A threshold that admits contradictions or lone dissenters is unsafe; one that leaves
nothing decidable is useless.

    python examples/calibrate_thresholds.py
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "diagnostics", "threshold_calibration.json")

NONTRIVIAL = "non-trivial (pi)"
TRIVIAL = "trivial (0)"


def load_sweep_runs() -> list[dict]:
    """The radius sweep from `run_failure_modes.py`, which is stored as one summary file.

    It is the most informative part of the corpus -- it is the only place where runs of the
    same question disagree -- so excluding it for a format reason would gut the calibration.
    Its records carry no per-point strategies, but the messages do name the broken chain.
    """
    path = os.path.join(ROOT, "results", "butadiene", "failure_modes.json")
    if not os.path.exists(path):
        return []
    rec = json.load(open(path))
    runs = []
    for r in rec["runs"]:
        msgs = " ".join(r.get("messages", []))
        n_broken = msgs.count("continuation chain broken")
        sign = -1 if float(r["product"]) < 0 else +1
        runs.append({
            "path": "results/butadiene/failure_modes.json",
            "system": "butadiene",
            "loop": f"sweep-{r['kind']}",
            "cas": f"CAS({rec['cas'][0]},{rec['cas'][1]})",
            "centre": tuple(round(float(c), 4) for c in rec["centre"]),
            "radius": (12.0, round(float(r["radius_pyr"]), 4)),
            "n_points": r["n_points"],
            "phase": r["phase"],
            "status": r["status"],
            "sign": sign,
            "endpoint": r.get("endpoint"),
            "min_adj": float(r["min_abs_overlap"]),
            # These runs were all solved with convergence required, so a FAILED status with
            # no convergence message means continuity or endpoint, not convergence.
            "all_converged": "did not converge" not in msgs,
            "n_broken": n_broken,
            "n_degraded": msgs.count("weaker warm start"),
            "broken": [],
        })
    return runs


def load_runs() -> list[dict]:
    """Every Berry record on disk, flattened to the fields the calibration needs."""
    runs = load_sweep_runs()
    for path in sorted(glob.glob(os.path.join(ROOT, "results", "**", "*.json"), recursive=True)):
        try:
            rec = json.load(open(path))
        except json.JSONDecodeError:
            continue
        if not (isinstance(rec, dict) and "result" in rec and "traversal" in rec):
            continue
        res, trav = rec["result"], rec["traversal"]
        pts = trav["points"]
        # Points reached by a cold start (or by nothing at all) were not reached by
        # continuation: the chain is broken there. Point 0 is cold by definition.
        broken = [p["index"] for p in pts
                  if p.get("strategy") in ("cold", "none-converged") and p["index"] != 0]
        degraded = [p["index"] for p in pts
                    if p.get("strategy") in ("warm-mo", "warm-mo-loose")]
        adj = [p["abs_overlap_with_prev"] for p in pts
               if p.get("abs_overlap_with_prev") is not None]
        runs.append({
            "path": os.path.relpath(path, ROOT),
            "system": os.path.basename(os.path.dirname(path)),
            "loop": res["loop_name"],
            "cas": res["cas_label"],
            "centre": tuple(round(float(c), 4) for c in res["loop_centre"]),
            "radius": tuple(round(float(r), 4) for r in res["loop_radius"]),
            "n_points": res["n_points"],
            "phase": res["berry_phase"],
            "status": res["status"],
            "sign": -1 if float(res["product_estimator"]) < 0 else +1,
            "endpoint": res.get("endpoint_estimator"),
            "min_adj": float(min(adj)) if adj else float("nan"),
            "all_converged": all(p["converged"] for p in pts),
            "n_broken": len(broken),
            "n_degraded": len(degraded),
            "broken": broken,
        })
    return runs


def verdict(run: dict, min_adj: float, require_chain: bool) -> tuple[bool, int]:
    """Re-decide one run under candidate thresholds, from its stored diagnostics.

    Returns (passes, sign). Mirrors berry.analyse: convergence, continuity, endpoint
    magnitude, estimator agreement -- plus the proposed chain check.
    """
    ok = run["all_converged"]
    ok = ok and np.isfinite(run["min_adj"]) and run["min_adj"] >= min_adj
    ep = run["endpoint"]
    if ep is not None:
        ok = ok and abs(ep) >= 0.90
        ok = ok and (np.sign(ep) == np.sign(run["sign"]))
    if require_chain:
        ok = ok and run["n_broken"] == 0
    return bool(ok), run["sign"]


def group_key(run: dict):
    return (run["system"], run["loop"], run["cas"], run["centre"], run["radius"])


def evaluate(runs: list[dict], min_adj: float, require_chain: bool) -> dict:
    groups: dict = defaultdict(list)
    for r in runs:
        groups[group_key(r)].append(r)

    n_pass = contradictions = decidable = singletons = dissenters = undecided = 0
    contradicting: list[str] = []
    dissenting: list[str] = []
    for key, members in groups.items():
        passing = [m for m in members if verdict(m, min_adj, require_chain)[0]]
        signs = {m["sign"] for m in passing}
        all_signs = {m["sign"] for m in members}
        label = f"{key[0]}/{key[1]}/{key[2]} r={key[4]}"
        n_pass += len(passing)
        if len(passing) >= 2 and len(signs) > 1:
            contradictions += 1
            contradicting.append(label)
        elif len(passing) >= 2:
            decidable += 1
        elif len(passing) == 1:
            singletons += 1
            # A refused run of the same question disagreeing with the single accepted one.
            if len(all_signs) > 1:
                dissenters += 1
                dissenting.append(label)
        else:
            undecided += 1
    return {
        "min_abs_overlap": min_adj,
        "require_chain": require_chain,
        "runs_passing": n_pass,
        "runs_total": len(runs),
        "questions": len(groups),
        "contradictions": contradictions,
        "lone_dissenters": dissenters,
        "decidable": decidable,
        "singletons": singletons,
        "undecided": undecided,
        "contradicting_questions": sorted(set(contradicting)),
        "dissenting_questions": sorted(set(dissenting)),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", type=float, nargs=3, default=(0.50, 0.98, 0.02),
                    metavar=("LO", "HI", "STEP"))
    args = ap.parse_args()

    runs = load_runs()
    print(f"{len(runs)} Berry records, "
          f"{len({group_key(r) for r in runs})} distinct questions "
          f"(same loop + CAS, differing N)\n")

    # How common are the two diagnostics, and do they track wrong answers?
    n_broken = sum(1 for r in runs if r["n_broken"])
    n_degraded = sum(1 for r in runs if r["n_degraded"])
    print(f"runs with a broken continuation chain : {n_broken:3d} / {len(runs)}")
    print(f"runs with a degraded warm start only  : {n_degraded:3d} / {len(runs)}")
    print(f"min adjacent overlap: worst {min(r['min_adj'] for r in runs):.4f}, "
          f"median {np.median([r['min_adj'] for r in runs]):.4f}\n")

    lo, hi, step = args.grid
    thresholds = [round(x, 4) for x in np.arange(lo, hi + 1e-9, step)]
    rows = []
    for chain in (False, True):
        for t in thresholds:
            rows.append(evaluate(runs, t, chain))

    print(f"{'chain?':>7} {'min_ovl':>8} {'pass':>6} {'decidable':>10} {'contra':>7} "
          f"{'dissent':>8} {'undecided':>10}")
    print("-" * 62)
    for row in rows:
        flag = "req" if row["require_chain"] else "-"
        star = "  <-- UNSAFE" if row["contradictions"] or row["lone_dissenters"] else ""
        print(f"{flag:>7} {row['min_abs_overlap']:8.2f} {row['runs_passing']:6d} "
              f"{row['decidable']:10d} {row['contradictions']:7d} "
              f"{row['lone_dissenters']:8d} {row['undecided']:10d}{star}")

    safe = [r for r in rows if r["contradictions"] == 0 and r["lone_dissenters"] == 0]
    best = max(safe, key=lambda r: (r["decidable"], -r["min_abs_overlap"])) if safe else None
    print()
    if best:
        print(f"Widest safe setting: min_abs_overlap={best['min_abs_overlap']:.2f}, "
              f"chain check {'ON' if best['require_chain'] else 'off'} "
              f"-> {best['decidable']} decidable, 0 contradictions, 0 lone dissenters.")
    else:
        print("Nothing in the swept range is free of both failure modes.")
    for row in rows:
        if row["dissenting_questions"] and row["min_abs_overlap"] in (0.80, 0.8):
            print(f"  at 0.80 / chain {'ON ' if row['require_chain'] else 'off'}: "
                  f"dissenting {row['dissenting_questions']}")

    save_json({"runs": runs, "sweep": rows, "best": best}, OUT)
    print(f"-> {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
