"""Configuration dataclasses.

Every scientific choice in this package is a named field here, so that nothing is
buried in code. These serialize to/from JSON for the run records.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field, fields
from typing import Any


def _from_dict(cls, data: dict):
    known = {f.name for f in fields(cls)}
    return cls(**{k: v for k, v in data.items() if k in known})


@dataclass
class CasConfig:
    """Electronic-structure settings shared by the Berry-phase and scan workflows."""

    basis: str = "sto-3g"
    ncas: int = 2
    nelecas: int = 2
    charge: int = 0
    spin: int = 0
    conv_tol: float = 1e-10
    # PySCF's own default is sqrt(conv_tol). Tightening this past ~1e-5 makes a *warm*
    # start fail its convergence test while sitting on the right answer to 1e-10: the
    # predictor lands so close that the macro iterations rattle around below the energy
    # threshold without ever meeting an over-tight gradient threshold.
    conv_tol_grad: float | None = 1e-5
    max_cycle_macro: int = 200
    fix_spin: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "CasConfig":
        return _from_dict(cls, data)

    @property
    def cas_label(self) -> str:
        n = self.nelecas if isinstance(self.nelecas, int) else sum(self.nelecas)
        return f"CAS({n},{self.ncas})"


@dataclass
class ContinuationConfig:
    """Thresholds that decide whether a loop traversal is trustworthy.

    These are the 'stability criteria' referred to in the documentation. A run that
    violates any of them is reported as FAILED rather than yielding a Berry phase.
    """

    # Smallest acceptable |<Psi_{k-1}|Psi_k>| anywhere on the loop. Below this the
    # discretization is too coarse or the solution jumped to a different branch.
    min_abs_overlap: float = 0.80
    # |<Psi_0|Psi_N>| at the closing point, where the geometry is identical to the start.
    # Well below 1 means the endpoint is not the same physical state.
    min_endpoint_overlap_magnitude: float = 0.90
    # Every CASSCF point must report convergence.
    require_converged: bool = True
    # Record the spin-adapted CASCI S0/S1 gap *within the tracked active space* at the
    # state-specific orbitals. This is a root-flipping risk indicator, not the physical
    # S0/S1 gap -- see berrycasscf.casscf.casci_roots.
    compute_casci_gap: bool = True
    # Warn if the tracked solution comes within this of another root of its own CI problem.
    small_gap_warning: float = 0.02

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ContinuationConfig":
        return _from_dict(cls, data)


@dataclass
class ScanConfig:
    """State-averaged CASSCF two-dimensional gap scan."""

    n_alpha: int = 21
    n_phi: int = 21
    margin: float = 2.0          # degrees of padding around the loop's bounding box
    nroots: int = 2
    weights: tuple[float, ...] = (0.5, 0.5)
    conv_tol: float = 1e-9
    max_cycle_macro: int = 200
    warm_start: bool = True      # reuse the previous grid point's orbitals

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ScanConfig":
        return _from_dict(cls, data)
