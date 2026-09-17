#!/usr/bin/env python
"""Time a single CASSCF solve, so cluster jobs can be sized from measurement.

Every runtime estimate in `docs/compute.md` is a measured per-point cost multiplied by a
point count. That per-point cost is machine-dependent -- and, for the larger active spaces,
thread-dependent in a way that is not worth predicting. This prints it.

    python examples/bench_point.py butadiene --cas 12 12
    OMP_NUM_THREADS=8 python examples/bench_point.py butadiene --cas 12 12 --sa

`--sa` times a state-averaged solve (the gap-scan cost) instead of the state-specific one
(the loop-transport cost); they differ by roughly a factor of two and the ladders use both.
Nothing is saved: this is a measurement of the machine, not a result about the molecule.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf.butadiene import butadiene_geom
from berrycasscf.casscf import build_mol, run_casscf, run_rhf
from berrycasscf.ethylene import ethylene_geom
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.scan import run_sa_casscf

# One geometry per system, at or near the intersection each study is about -- the expensive
# region, so an estimate built on this does not flatter itself.
SYSTEMS = {
    "formaldimine": (formaldimine_geom, "sto-3g", (130.0, 89.9)),
    "ethylene": (ethylene_geom, "6-31g*", (90.0, 110.0)),
    "butadiene": (butadiene_geom, "6-31g*", (90.0, 101.85)),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", choices=sorted(SYSTEMS))
    ap.add_argument("--cas", nargs=2, type=int, default=(12, 12), metavar=("NE", "NCAS"))
    ap.add_argument("--sa", action="store_true", help="time a state-averaged solve instead")
    ap.add_argument("--repeat", type=int, default=1)
    args = ap.parse_args()

    geom_fn, basis, point = SYSTEMS[args.system]
    ne, ncas = args.cas
    mol = build_mol(geom_fn(*point), basis)

    print(f"{args.system} CAS({ne},{ncas})/{basis} at {point}")
    print(f"  {mol.nao} AOs, {mol.nelectron} electrons; "
          f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'unset')}")

    t0 = time.time()
    mf = run_rhf(mol)
    print(f"  RHF                     {time.time() - t0:8.2f} s   E = {mf.e_tot:.8f}")

    times = []
    for i in range(args.repeat):
        t0 = time.time()
        if args.sa:
            e_states, _, converged = run_sa_casscf(mol, ncas, ne, (0.5, 0.5), mf=mf)
            energy, extra = e_states[0], f"gap = {(e_states[1]-e_states[0])*1e3:.3f} mHa"
        else:
            wfn = run_casscf(mol, ncas, ne, mf=mf)
            energy, converged = wfn.energy, wfn.converged
            extra = f"{wfn.ci.size} determinants"
        dt = time.time() - t0
        times.append(dt)
        kind = "SA-CASSCF" if args.sa else "CASSCF"
        print(f"  {kind:<22}  {dt:8.2f} s   E = {energy:.8f}  "
              f"converged={converged}  {extra}")

    if args.repeat > 1:
        print(f"  mean {np.mean(times):.2f} s over {args.repeat} solves")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
