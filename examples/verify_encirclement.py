#!/usr/bin/env python
"""Does loop transport encircle the intersection the gap scan found?

The check the project should have run first. Every comparison so far has been indirect --
locate the state-specific object by bisection, locate the state-averaged one by scanning, and
compare two positions each with its own error bar. This asks the question in one loop:

    put a small circular loop ON the gap scan's intersection and walk it.

    pi    the state-specific transport encircles it too: one object, and the two methods agree
          about what they are looking at
    0     it does not: whatever carries the phase is somewhere else, and the two methods are
          not sensing the same thing
    refused  the loop could not be walked there, which is itself localizing -- the seam is
          about that far away

The radius has to be small enough to contain nothing else and large enough to be walkable; the
floor measured in notebooks/probing_by_small_loops.ipynb is about 1 deg, so 2 deg is the
default. A zero means an *even* count, so a radius enclosing two degeneracies would also read
0; that is why the radius is reported with every verdict.

    python examples/verify_encirclement.py butadiene
    python examples/verify_encirclement.py ethylene --cas 2,2 4,4
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ContinuationConfig
from berrycasscf.butadiene import butadiene_geom
from berrycasscf.ethylene import ethylene_geom
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.localize import evaluate_radius
from berrycasscf.refine import cone_apex
from berrycasscf.runlog import JobLog
from berrycasscf.scan import ScanResult
from berrycasscf.store import load_json, save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SYSTEMS = {
    "formaldimine": {"geom_fn": formaldimine_geom, "basis": "sto-3g",
                     "scan": "results/scan/C_x_cas{ne}-{ncas}_25x25.npz"},
    "ethylene": {"geom_fn": ethylene_geom, "basis": "6-31g*",
                 "scan": "results/ethylene/ethylene_scan_cas{ne}-{ncas}_*.npz"},
    "butadiene": {"geom_fn": butadiene_geom, "basis": "6-31g*",
                  "search": "results/butadiene/gap_minimum_search.json",
                  "scan": "results/butadiene/butadiene_scan_cas{ne}-{ncas}_*.npz"},
}
LADDER = [(2, 2), (4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]


def sa_intersection(system: str, ne: int, ncas: int):
    """Where the state-averaged methods put this rung's intersection, and how it was found."""
    spec = SYSTEMS[system]
    search = spec.get("search")
    if search and os.path.exists(os.path.join(ROOT, search)):
        for run in json.load(open(os.path.join(ROOT, search)))["runs"]:
            if tuple(run["cas"]) == (ne, ncas):
                return tuple(run["found"]), "direct 2D search"
    hits = sorted(glob.glob(os.path.join(ROOT, spec["scan"].format(ne=ne, ncas=ncas))))
    if not hits:
        return None, None
    res = ScanResult.load(hits[0])
    row = int(np.nanargmin(np.nanmin(res.gap, axis=1)))
    fit = cone_apex(res.phis, res.gap[row] * 1e3, window=3)
    return (float(res.alphas[row]), float(fit.position)), "grid scan, cone-refined"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS))
    ap.add_argument("--cas", nargs="*", default=None, metavar="NE,NCAS")
    ap.add_argument("--radius", type=float, default=2.0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    rungs = ([tuple(int(v) for v in c.split(",")) for c in args.cas] if args.cas else LADDER)
    spec = SYSTEMS[args.system]
    out = os.path.join(ROOT, "results", "encirclement",
                       f"{args.system}_r{args.radius:g}.json")
    done = {}
    if os.path.exists(out) and not args.force:
        record = load_json(out)
        done = {tuple(r["cas"]): r for r in record.get("rungs", [])}
        if len(done) >= len(rungs):
            print(f"SKIP: {os.path.relpath(out, ROOT)} covers every rung; --force to redo")
            return 0
        if done:
            print(f"resuming: {len(done)} rungs already probed")

    log = JobLog(f"encirclement_{args.system}", total=len(rungs))
    t0 = time.time()
    results = []
    for ne, ncas in rungs:
        if (ne, ncas) in done:
            results.append(done[(ne, ncas)])
            continue
        point, how = sa_intersection(args.system, ne, ncas)
        if point is None:
            print(f"  CAS({ne},{ncas}): no gap-scan record; skipped")
            continue
        print(f"\n=== CAS({ne},{ncas}): probing ({point[0]:.2f}, {point[1]:.2f}) "
              f"[{how}] with r = {args.radius} deg ===")
        probe = evaluate_radius(point, (args.radius, args.radius), 1.0,
                                cas=CasConfig(basis=spec["basis"], ncas=ncas, nelecas=ne),
                                cont=ContinuationConfig(), geom_fn=spec["geom_fn"],
                                name=f"V{ne}", progress=log)
        rec = probe.to_dict()
        rec.update({"cas": [ne, ncas], "sa_intersection": list(point), "how": how})
        results.append(rec)
        save_json({"system": args.system, "radius": args.radius, "rungs": results,
                   "wall_time": time.time() - t0}, out)

    print("\n" + "=" * 80)
    print(f"{args.system}: does loop transport encircle the gap scan's intersection?"
          f"   (r = {args.radius} deg)\n")
    print(f"{'CAS':>9} {'gap-scan intersection':>24} {'verdict':>14} {'encircles?':>12} "
          f"{'wall (s)':>9}")
    print("-" * 80)
    for rec in results:
        p = rec["sa_intersection"]
        answer = {"pi": "YES", "zero": "NO", "undetermined": "refused"}[rec["verdict"]]
        print(f"{('(%d,%d)' % tuple(rec['cas'])):>9} "
              f"{f'({p[0]:.2f}, {p[1]:.2f})':>24} {rec['verdict']:>14} {answer:>12} "
              f"{rec['wall_time']:>9.1f}")
    log.done()
    print(f"\n-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
