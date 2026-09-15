#!/usr/bin/env python
"""How large an active space does each method actually need?

Runs both workflows over a ladder of active spaces on one system, holding **everything else
fixed** -- same basis, same geometry family, and above all the *same loops*. Varying the loops
with the active space would mean the methods were answering different questions.

Two stages:

  scan   equal-weight SA-CASSCF gap scans over a fixed region, one per active space.
         Answers "does this CAS put a conical intersection inside the loop?"
  berry  Berry-phase continuation on three fixed loops, one set per active space.
         Answers "does this CAS give the right topology?"

The loops are centred on the intersection located by the *reference* (largest) active space,
so the ladder is judged against a single fixed target.

    python examples/run_active_space_study.py scan
    python examples/run_active_space_study.py berry
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from berrycasscf import CasConfig, ContinuationConfig, ScanConfig, run_loop, scan_gap
from berrycasscf.ethylene import CI_REGION, DEFAULT_BASIS, ethylene_geom
from berrycasscf.geometry import Loop
from berrycasscf.scan import ScanResult
from berrycasscf.store import save_berry_run, berry_record_exists

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(ROOT, "results", "ethylene")

# The ladder. CAS(12,12) spans the full valence space (only the two C 1s orbitals stay core)
# and is the in-basis reference.
CAS_LADDER = [(2, 2), (4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]
REFERENCE_CAS = (12, 12)

GRID = (11, 11)
# "anchor" costs two SA-CASSCF solves per grid point but is path-independent (see
# berrycasscf.config.ScanConfig). CAS(12,12) at ~30 s/solve would double a one-hour scan, so it
# uses cold-only; the mirror-symmetry check printed below verifies that this was safe.
STRATEGY = {(12, 12): "cold"}
DEFAULT_STRATEGY = "anchor"

LOOP_RADIUS = (12.0, 12.0)       # degrees in (tau, phi)
# Three discretizations: N=13 is deliberately coarse enough to fail the continuity test for
# the larger active spaces, which is itself part of the result.
NPOINTS = [13, 21, 31]


def scan_path(ne: int, ncas: int) -> str:
    return os.path.join(RESULT_DIR, f"ethylene_scan_cas{ne}-{ncas}_{GRID[0]}x{GRID[1]}.npz")


def berry_path(loop: str, ne: int, ncas: int, n: int) -> str:
    return os.path.join(RESULT_DIR, f"ethylene_{loop}_cas{ne}-{ncas}_N{n}.json")


def mirror_asymmetry(res: ScanResult) -> float:
    """Max |gap(tau) - gap(180-tau)| in mHa.

    Geometries at tau and 180-tau are exact mirror images, so a correct scan is symmetric.
    Any asymmetry is active-space drift along the scan path, not physics.
    """
    g = res.gap * 1000.0
    return float(np.nanmax(np.abs(g - g[::-1, :])))


def do_scans(args) -> int:
    os.makedirs(RESULT_DIR, exist_ok=True)
    for ne, ncas in CAS_LADDER:
        path = scan_path(ne, ncas)
        if os.path.exists(path) and not args.force:
            res = ScanResult.load(path)
            if np.isfinite(res.e_states).all():
                t, p, g = res.min_gap_point()
                print(f"[skip] CAS({ne},{ncas}): min {g*1e3:8.3f} mHa at "
                      f"(tau={t:6.2f}, phi={p:6.2f})   mirror asymmetry "
                      f"{mirror_asymmetry(res):.2e} mHa")
                continue
        strategy = STRATEGY.get((ne, ncas), DEFAULT_STRATEGY)
        print(f"\n=== scan CAS({ne},{ncas})/{args.basis}  {GRID[0]}x{GRID[1]}  "
              f"strategy={strategy} ===")
        t0 = time.time()
        res = scan_gap(
            CI_REGION,
            cas=CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
            scan=ScanConfig(n_alpha=GRID[0], n_phi=GRID[1], margin=0.0, strategy=strategy),
            geom_fn=ethylene_geom,
            progress=None if args.quiet else print,
            checkpoint=path,
        )
        res.save(path)
        t, p, g = res.min_gap_point()
        print(f"  min gap {g*1e3:.3f} mHa at (tau={t:.2f}, phi={p:.2f})   "
              f"mirror asymmetry {mirror_asymmetry(res):.2e} mHa   "
              f"[{time.time()-t0:.0f} s]")
    return 0


def reference_ci(override: tuple[float, float] | None = None) -> tuple[float, float]:
    """Intersection position defining the fixed loops, from the reference active space.

    ``override`` lets the Berry stage start before the (slow) reference scan has finished;
    the value used is always re-checked against the reference scan by ``verify_centre``
    once it exists, so a wrong guess cannot pass silently.
    """
    if override is not None:
        return override
    path = scan_path(*REFERENCE_CAS)
    if not os.path.exists(path):
        raise SystemExit(
            f"Reference scan {os.path.relpath(path, ROOT)} not found. "
            "Run `run_active_space_study.py scan` first, or pass --centre TAU PHI."
        )
    t, p, _ = ScanResult.load(path).min_gap_point()
    return t, p


def verify_centre(used: tuple[float, float]) -> None:
    """Check the loop centre against the reference scan, once that scan exists."""
    path = scan_path(*REFERENCE_CAS)
    if not os.path.exists(path):
        print("  NOTE: reference scan not yet available; loop centre unverified.")
        return
    res = ScanResult.load(path)
    if not np.isfinite(res.e_states).all():
        print("  NOTE: reference scan incomplete; loop centre unverified.")
        return
    t, p, _ = res.min_gap_point()
    if abs(t - used[0]) < 1e-6 and abs(p - used[1]) < 1e-6:
        print(f"  loop centre verified against CAS{REFERENCE_CAS}: ({t:.2f}, {p:.2f})")
    else:
        print(f"  WARNING: loops were centred on {used}, but CAS{REFERENCE_CAS} puts the "
              f"intersection at ({t:.2f}, {p:.2f}). Re-run this stage with --force.")


def study_loops(centre: tuple[float, float]) -> dict[str, Loop]:
    """One loop enclosing the reference CI and two controls displaced along tau."""
    t0, p0 = centre
    return {
        "E_x": Loop("E_x", (t0, p0), LOOP_RADIUS),
        "E_1": Loop("E_1", (t0 - 2.5 * LOOP_RADIUS[0], p0), LOOP_RADIUS),
        "E_2": Loop("E_2", (t0 + 2.5 * LOOP_RADIUS[0], p0), LOOP_RADIUS),
    }


def do_berry(args) -> int:
    os.makedirs(RESULT_DIR, exist_ok=True)
    centre = reference_ci(tuple(args.centre) if args.centre else None)
    loops = study_loops(centre)
    print(f"Reference CI from CAS{REFERENCE_CAS}: (tau, phi) = "
          f"({centre[0]:.2f}, {centre[1]:.2f})")
    for name, lp in loops.items():
        print(f"  {name}: centre ({lp.centre[0]:7.2f}, {lp.centre[1]:7.2f})  "
              f"radius {lp.radius}  encloses reference CI: {lp.encloses(*centre)}")
    verify_centre(centre)

    for ne, ncas in CAS_LADDER:
        for name, loop in loops.items():
            for n in NPOINTS:
                path = berry_path(name, ne, ncas, n)
                if berry_record_exists(path) and not args.force:
                    print(f"[skip] {os.path.basename(path)}")
                    continue
                print(f"\n=== {name}  CAS({ne},{ncas})/{args.basis}  N={n} ===")
                try:
                    res, trav = run_loop(
                        loop.with_n_points(n),
                        cas=CasConfig(basis=args.basis, ncas=ncas, nelecas=ne),
                        cont=ContinuationConfig(),
                        geom_fn=ethylene_geom,
                        progress=None if args.quiet else print,
                    )
                except Exception as exc:                      # noqa: BLE001
                    print(f"  RUN FAILED: {type(exc).__name__}: {exc}")
                    continue
                print(res.summary())
                save_berry_run(res, trav, path)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["scan", "berry"])
    ap.add_argument("--basis", default=DEFAULT_BASIS)
    ap.add_argument("--centre", nargs=2, type=float, default=None, metavar=("TAU", "PHI"),
                    help="loop centre, to start the berry stage before the reference scan "
                         "finishes; it is re-checked against that scan afterwards")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    return do_scans(args) if args.stage == "scan" else do_berry(args)


if __name__ == "__main__":
    raise SystemExit(main())
