#!/usr/bin/env python
"""What is in the region whose loops report pi at ethylene CAS(2,2)?

Small-loop probing left an arithmetic impossibility: loops of radius 6 and above about
(90, 110.9) report pi, loops of radius 3 and below report 0, nothing on the mirror line
encloses anything, and the only candidates off the line come in exact mirror pairs -- which a
loop centred on the line must enclose two of, or none. An even count cannot give pi.

So the region is worth looking at directly, with two cheap maps that need no loop at all:

``incas_gap``  the S0-to-second-root gap *within the tracked active space*, at the
               state-specific orbitals. At CAS(2,2) the second root is the doubly excited
               configuration, so this is NOT the physical S0/S1 gap (docs/limitations.md) --
               it is a picture of the solution the continuation is tracking. Raggedness means
               neighbouring geometries converged to different CASSCF solutions.
``ring``       cold-start energies around circles of the radii that gave pi and 0. Each point
               is solved independently, so a jump is the cold solver landing on a different
               solution: it maps where competing solutions exist, not what a warm-started walk
               does.

Neither proves what the loops are doing. Together they say whether the region is one smooth
sheet with a degeneracy in it, or several sheets meeting -- and a sign picked up crossing from
one sheet to another is not a Berry phase.

    python examples/diagnose_ethylene_solutions.py
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf.casscf import build_mol, casci_roots, run_casscf
from berrycasscf.ethylene import ethylene_geom
from berrycasscf.runlog import JobLog

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "diagnostics")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cas", nargs=2, type=int, default=(2, 2), metavar=("NE", "NCAS"))
    ap.add_argument("--centre", nargs=2, type=float, default=(90.0, 110.9))
    ap.add_argument("--radii", nargs="*", type=float, default=[6.0, 3.0])
    ap.add_argument("--n-ring", type=int, default=72)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    ne, ncas = args.cas
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"ethylene_cas{ne}-{ncas}_solutions.npz")
    if os.path.exists(path) and not args.force:
        print(f"SKIP: {os.path.relpath(path, ROOT)} exists; --force to redo")
        return 0

    def energy(tw, pyr, want_gap=False):
        wfn = run_casscf(build_mol(ethylene_geom(tw, pyr), "6-31g*"), ncas, ne)
        if not want_gap:
            return wfn.energy, np.nan
        roots = casci_roots(wfn, nroots=2)
        return wfn.energy, float((roots[1] - roots[0]) * 1e3)

    tws = np.arange(84.0, 96.5, 1.0)
    pyrs = np.arange(104.0, 118.5, 1.0)
    log = JobLog(f"ethylene_solutions_cas{ne}-{ncas}",
                 total=tws.size * pyrs.size + len(args.radii) * args.n_ring)
    gap = np.full((tws.size, pyrs.size), np.nan)
    ener = np.full_like(gap, np.nan)
    t0 = time.time()
    for i, tw in enumerate(tws):
        for j, p in enumerate(pyrs):
            try:
                ener[i, j], gap[i, j] = energy(tw, p, want_gap=True)
            except Exception as exc:                            # noqa: BLE001
                log.note(f"   ({tw}, {p}) failed: {type(exc).__name__}")
            log(f"grid ({tw:.0f}, {p:.0f})")

    rings = {}
    th = np.linspace(0, 2 * np.pi, args.n_ring, endpoint=False)
    for r in args.radii:
        es = np.full(args.n_ring, np.nan)
        for k, t in enumerate(th):
            tw = args.centre[0] + r * np.cos(t)
            p = args.centre[1] + r * np.sin(t)
            try:
                es[k], _ = energy(tw, p)
            except Exception:                                   # noqa: BLE001
                pass
            log(f"ring r={r:g} theta={np.degrees(t):.0f}")
        rings[f"ring_{r:g}"] = es

    np.savez(path, tws=tws, pyrs=pyrs, incas_gap=gap, energy=ener, theta=th,
             centre=np.array(args.centre), radii=np.array(args.radii), **rings)
    log.done()
    print(f"\n  in-CAS gap: min {np.nanmin(gap):.1f} mHa, max {np.nanmax(gap):.1f} mHa")
    for r in args.radii:
        es = rings[f"ring_{r:g}"]
        d = np.abs(np.diff(np.concatenate([es, es[:1]])))
        print(f"  ring r={r:g}: median step {np.median(d):.2e} Ha, largest {np.max(d):.2e}, "
              f"steps above 10x median: {int((d > 10 * np.median(d)).sum())}")
    print(f"  [{time.time() - t0:.0f} s]  -> {os.path.relpath(path, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
