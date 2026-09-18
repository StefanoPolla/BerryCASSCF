#!/usr/bin/env python
"""Which of these candidate centres has a degeneracy inside it?

The primitive an iterative search needs. Bisection asks "how far is it?" by making a loop
graze the seam, which is the expensive and failure-prone configuration. This asks the cheap
question instead -- "is it inside *this* small loop?" -- at several centres, and lets the
pattern of yes/no say where to look next.

Each centre gets one circular loop of the given radius, walked with both step-control
settings. The verdicts mean:

  pi            an odd number of degeneracies inside -- normally one, and the search has found
                the neighbourhood
  zero          an even number, normally none. In a system with a mirror symmetry a loop
                straddling the symmetry line can enclose a pair and report zero, so a zero on
                the line and a zero off it do not mean the same thing
  undetermined  the loop could not be walked, which for a small loop usually means it passes
                very close to the seam. That is localizing information, not a failure

    python examples/run_centre_probe.py ethylene --cas 2 2 --radius 2 \\
           --centres 90,105 90,107 90,109 90,113 90,115
"""

from __future__ import annotations

import argparse
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
from berrycasscf.runlog import JobLog
from berrycasscf.store import load_json, save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GEOM = {"formaldimine": (formaldimine_geom, "sto-3g"),
        "ethylene": (ethylene_geom, "6-31g*"),
        "butadiene": (butadiene_geom, "6-31g*")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(GEOM))
    ap.add_argument("--cas", nargs=2, type=int, default=(2, 2), metavar=("NE", "NCAS"))
    ap.add_argument("--centres", nargs="+", required=True, metavar="X,Y")
    ap.add_argument("--radius", type=float, required=True,
                    help="loop radius in degrees, circular and the same at every centre, so "
                         "the verdicts are comparable")
    ap.add_argument("--label", required=True, help="names the record; use one per sweep")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    geom_fn, basis = GEOM[args.system]
    ne, ncas = args.cas
    centres = [tuple(float(v) for v in c.split(",")) for c in args.centres]
    tag = f"{args.system}_cas{ne}-{ncas}_{args.label}"
    out = os.path.join(ROOT, "results", "centre_probe", f"{tag}.json")

    done = {}
    if os.path.exists(out) and not args.force:
        record = load_json(out)
        if record.get("complete"):
            print(f"SKIP: {os.path.relpath(out, ROOT)} exists; --force to redo")
            return 0
        done = {tuple(p["centre"]): p for p in record.get("probes", [])}
        if done:
            print(f"resuming: {len(done)} centres already probed")

    cas = CasConfig(basis=basis, ncas=ncas, nelecas=ne)
    cont = ContinuationConfig()
    log = JobLog(f"centre_probe_{tag}", total=len(centres))
    print(f"{args.system} {cas.cas_label}/{basis}: radius {args.radius} deg at "
          f"{len(centres)} centres\n")

    t0 = time.time()
    probes = []
    for c in centres:
        if c in done:
            probes.append(done[c])
            print(f"  [cached] {c}")
            continue
        probe = evaluate_radius(c, (args.radius, args.radius), 1.0, cas=cas, cont=cont,
                               geom_fn=geom_fn, name=f"C{c[1]:g}", progress=log)
        rec = probe.to_dict()
        rec["centre"] = list(c)
        probes.append(rec)
        save_json({"system": args.system, "cas": [ne, ncas], "basis": basis,
                   "radius": args.radius, "centres": [list(x) for x in centres],
                   "probes": probes, "complete": len(probes) == len(centres),
                   "wall_time": time.time() - t0}, out)

    print("\n" + "=" * 70)
    print(f"{'centre':>18} {'verdict':>14} {'points':>10} {'worst overlap':>14} {'wall (s)':>9}")
    print("-" * 70)
    for rec in probes:
        runs = rec.get("runs", [])
        pts = "/".join(str(r["n_points"]) for r in runs) if runs else "-"
        worst = min((r["min_overlap"] for r in runs if np.isfinite(r["min_overlap"])),
                    default=float("nan"))
        c = rec["centre"]
        print(f"{f'({c[0]:g}, {c[1]:g})':>18} {rec['verdict']:>14} {pts:>10} "
              f"{worst:>14.3f} {rec['wall_time']:>9.1f}")
    enclosing = [p["centre"] for p in probes if p["verdict"] == "pi"]
    print(f"\n  encloses at: {enclosing if enclosing else 'none of these centres'}")
    log.done()
    print(f"-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
