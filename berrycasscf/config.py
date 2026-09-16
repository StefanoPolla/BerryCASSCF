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
    #
    # 0.80 was once suspected of being too permissive, because the most suspicious run in
    # the study passed at 0.844. Calibrating against all 116 saved runs
    # (`examples/calibrate_thresholds.py`) says otherwise, and identifies the real culprit:
    # with `require_unbroken_chain` ON, no threshold from 0.70 to 0.92 admits a single
    # contradiction or lone dissenter. With it OFF, 0.80 admits a lone dissenter and 0.70
    # admits an outright contradiction, and the threshold has to go to 0.86 to be safe --
    # which costs 5 decidable questions against the chain check's 1. The threshold is not
    # where the danger was.
    min_abs_overlap: float = 0.80
    # |<Psi_0|Psi_N>| at the closing point, where the geometry is identical to the start.
    # Well below 1 means the endpoint is not the same physical state.
    min_endpoint_overlap_magnitude: float = 0.90
    # Every CASSCF point must report convergence.
    require_converged: bool = True
    # Refuse a run in which the solver fallback ladder fell through to a *cold* start at any
    # point after the first. Such a point was never reached by continuation, so the chain is
    # broken there. Gauge fixing repairs a cold start's arbitrary overall sign but NOT a
    # change of branch, and a modest branch change still clears min_abs_overlap -- which is
    # exactly how butadiene radius 6, N=31 came to pass every check at 0.844 and report a
    # phase its neighbours contradict (docs/findings.md §5). This was reported as a warning
    # for most of the project and is now binding; see the calibration above for the cost.
    require_unbroken_chain: bool = True
    # Walk the solver fallback ladder when a warm start fails to converge. Turn this OFF to
    # study the single-update regime of arXiv:2304.06070, where a point is *deliberately* not
    # converged and retrying would defeat the point of the exercise.
    use_fallback: bool = True
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
class AdaptiveConfig:
    """Step-size control for walking a loop, in units of the loop parameter t in [0, 1].

    Uniform discretization sizes the whole loop for its worst arc. The per-step difficulty
    is very non-uniform -- measured adjacent overlaps on the harder loops span 0.85 to 0.99
    -- so a uniform schedule spends most of its solves where they are not needed. The error
    indicator is free: |<Psi_prev|Psi_new>| is already computed at every step, and costs
    ~0.1 ms against a solve of 30 ms to 3 min.

    Unequal spacing costs nothing here, because a topological readout needs the chain to be
    continuous and closed, not evenly sampled.
    """

    # Largest allowed step. 1/10 of the loop reproduces a coarse uniform walk if nothing
    # ever has to shrink, so a run that never rejects costs about what N=10 costs.
    d_max: float = 0.1
    # Smallest allowed step. Reaching it is a RESULT, not a safeguard: if shrinking restores
    # continuity the problem was discretization, and if it does not, the loop is passing
    # through something. 1e-3 of a loop of radius ~12 deg is a nuclear displacement of
    # ~0.08 deg, far below any resolution the chemistry justifies.
    d_min: float = 1e-3
    # Target for the mismatch m = 1 - |<Psi_prev|Psi_new>|. The controller steers toward
    # this rather than merely rejecting at a floor: a hard accept/reject at the failure
    # threshold oscillates and wastes solves, and a rejection costs a full solve with no
    # cheap predictor available. 0.02 is |overlap| = 0.98, comfortably inside the 0.80 that
    # ContinuationConfig refuses, so the walk carries margin rather than skating the limit.
    target_mismatch: float = 0.02
    # Reject a trial step whose mismatch exceeds this. Kept well above target_mismatch so
    # that ordinary controller overshoot does not trigger an expensive rejection.
    accept_mismatch: float = 0.10
    # Step-size response m ~ d^2 gives d_new = d * sqrt(target/m). These clip that factor.
    growth_cap: float = 2.0
    max_shrink: float = 0.1
    # Guard against a pathological walk. Hitting either is reported, never silently ignored.
    max_points: int = 500
    max_rejections: int = 200
    # Avoids division by zero when a step is so good the measured mismatch underflows.
    mismatch_floor: float = 1e-12

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "AdaptiveConfig":
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
    # How each grid point is started. Unlike the Berry-phase loop, a scan is NOT a
    # continuation problem, so nothing here needs a consistent gauge -- the only question is
    # which SA-CASSCF stationary point each solve lands on.
    #   "warm" : reuse the previous grid point's orbitals. Fastest, but for larger active
    #            spaces the active space can drift along the scan path, making the result
    #            path-dependent. Measured for ethylene CAS(8,8)/6-31G*: two geometries that
    #            are exact mirror images returned 19.1 and 29.4 mHa, against 26.9 mHa cold.
    #   "cold" : fresh RHF orbitals at every point. Path-independent and deterministic, but
    #            can miss a lower-energy SA-CASSCF solution that a warm start finds.
    #   "best" : warm and cold, keep the lower SA energy. Finds better solutions than cold,
    #            but inherits warm's path dependence -- measured 7.9 mHa mirror asymmetry for
    #            ethylene CAS(4,4).
    #   "anchor": cold, plus a warm start transferred from ONE fixed anchor geometry solved
    #            cold at the centre of the region; keep the lower SA energy. Path-independent,
    #            and reaches a lower SA energy than cold -- but the anchor generally does not sit
    #            on a symmetry element, so the transferred guess is NOT equivariant.
    #
    # Default is "cold". It is the only strategy that respects molecular symmetry for free:
    # mirror-image geometries give integrals related by a signed permutation and the minao guess
    # transforms the same way, so the two solves are one calculation in different coordinates.
    # Measured: cold reproduces ethylene's mirror symmetry to 1e-11 mHa against anchor's 4e-04.
    # It is also the cheapest (one solve per point), embarrassingly parallel (no anchor to solve
    # first), and defines the active space by the same prescription at every geometry rather than
    # by whatever was carried in from elsewhere.
    #
    # The cost is that cold is not always the lowest-energy solution -- at ethylene CAS(4,4) a
    # warm sweep reached 1.3e-03 Ha lower. That did not move the observable (the intersection
    # position and gap agree across all three strategies at every rung tested), and a lower SA
    # energy is not unambiguously better anyway: distinct minima can correspond to different
    # active spaces, so the variationally lowest solution need not be the intended one.
    strategy: str = "cold"
    warm_start: bool = True      # deprecated alias; ignored unless strategy is None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ScanConfig":
        return _from_dict(cls, data)
