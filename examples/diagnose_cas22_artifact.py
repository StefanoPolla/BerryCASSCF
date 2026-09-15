#!/usr/bin/env python
"""Why the CAS(2,2) state-averaged comparator disagrees with the Berry phase.

The SA-CASSCF(2,2) gap scan reports a near-degeneracy inside the *trivial* control loop C_2.
This script recomputes the gap at that point with progressively larger active spaces to
establish whether the degeneracy is real. Results are written to
``results/diagnostics/cas22_artifact.json`` and read back by the results notebook.

    python examples/diagnose_cas22_artifact.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pyscf import mcscf

from berrycasscf.casscf import apply_singlet_constraint, build_mol, run_rhf
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.scan import ScanResult
from berrycasscf.store import save_json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ACTIVE_SPACES = [(2, 2), (4, 4), (6, 6), (8, 8)]


def gap_at(alpha: float, phi: float, nelecas: int, ncas: int, basis: str = "sto-3g") -> dict:
    mol = build_mol(formaldimine_geom(alpha, phi), basis)
    mf = run_rhf(mol)
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    apply_singlet_constraint(mc)
    mc = mc.state_average_([0.5, 0.5])
    mc.verbose = 0
    mc.kernel()
    return {
        "cas": f"CAS({nelecas},{ncas})",
        "e0": float(mc.e_states[0]),
        "e1": float(mc.e_states[1]),
        "gap": float(mc.e_states[1] - mc.e_states[0]),
        "converged": bool(mc.converged),
    }


def main() -> int:
    # Locate the spurious minimum from the saved CAS(2,2) scan of the trivial loop C_2.
    scan_path = os.path.join(ROOT, "results", "scan", "C_2_cas2-2_25x25.npz")
    if os.path.exists(scan_path):
        alpha, phi, gap22 = ScanResult.load(scan_path).min_gap_point()
        source = os.path.relpath(scan_path, ROOT)
    else:
        alpha, phi, gap22 = 150.92, 89.9, float("nan")
        source = "hardcoded fallback (run the C_2 CAS(2,2) scan first)"

    print(f"Spurious CAS(2,2) minimum taken from {source}")
    print(f"  (alpha, phi) = ({alpha:.2f}, {phi:.2f}), CAS(2,2) gap = {gap22:.6f} Ha\n")
    print("Recomputing the gap there with larger active spaces:\n")

    rows = []
    for ne, ncas in ACTIVE_SPACES:
        row = gap_at(alpha, phi, ne, ncas)
        rows.append(row)
        note = "  <-- claims a conical intersection" if row["gap"] < 1e-3 else ""
        print(f"  {row['cas']:9s} gap = {row['gap']:.6f} Ha   "
              f"E0 = {row['e0']:.6f}  converged={row['converged']}{note}")

    verdict = (
        "ARTIFACT: the larger active spaces agree there is no degeneracy here"
        if min(r["gap"] for r in rows[1:]) > 1e-2
        else "INCONCLUSIVE: larger active spaces also find a small gap"
    )
    print(f"\n{verdict}")

    out = os.path.join(ROOT, "results", "diagnostics", "cas22_artifact.json")
    save_json(
        {
            "point": {"alpha": alpha, "phi": phi},
            "source_scan": source,
            "rows": rows,
            "verdict": verdict,
        },
        out,
    )
    print(f"-> {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
