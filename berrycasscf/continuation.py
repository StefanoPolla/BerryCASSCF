"""Continuous tracking of one CASSCF stationary solution around a closed nuclear loop.

This is deliberately *not* a series of independent CASSCF calculations. Each point is
warm-started from its predecessor in both the orbitals (carried through the OAO frame)
and the CI vector, and every step is checked with an exact nonorthogonal overlap.

Gauge convention
----------------
CASSCF returns a CI vector whose overall sign is arbitrary. We fix it sequentially: the
wavefunction at point k is flipped, if necessary, so that ``<Psi_{k-1}|Psi_k> > 0``. After
that every adjacent overlap on the open chain is positive by construction and all the
topological information sits in the single closing overlap.

Two Berry-phase estimators come out of one traversal (see :mod:`berrycasscf.berry`).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, asdict
from typing import Callable, Sequence

import numpy as np

from .casscf import (
    build_mol,
    run_casscf,
    run_rhf,
    transfer_mo,
    orthonormality_error,
    casci_roots,
)
from .config import CasConfig, ContinuationConfig
from .geometry import Loop, formaldimine_geom
from .overlap import CasWavefunction, cas_overlap


@dataclass
class PointRecord:
    """Everything worth knowing about one point of the traversal."""

    index: int
    t: float
    alpha: float
    phi: float
    energy: float
    converged: bool
    # <Psi_{k-1} | Psi_k~> with Psi_{k-1} gauge-fixed and Psi_k~ as the solver returned it.
    raw_overlap_with_prev: float | None
    gauge_sign: int
    # |<Psi_{k-1}|Psi_k>| after gauge fixing; the continuity diagnostic.
    abs_overlap_with_prev: float | None
    casci_gap: float | None
    guess_orthonormality_error: float | None
    max_mo_change: float | None
    strategy: str
    wall_time: float

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class LoopTraversal:
    """Result of walking one loop once."""

    loop_name: str
    loop_centre: tuple[float, float]
    loop_radius: tuple[float, float]
    n_points: int
    cas_label: str
    basis: str
    points: list[PointRecord]
    # <Psi_{N-1}|Psi_0>: closes the cycle without re-solving. Sign is the product estimator.
    closing_overlap: float
    # <Psi_0|Psi_N> where R_N == R_0 exactly, after an extra continuation step.
    endpoint_overlap: float | None
    endpoint_converged: bool | None
    wall_time: float
    warnings: list[str] = field(default_factory=list)

    @property
    def adjacent_abs_overlaps(self) -> np.ndarray:
        return np.array(
            [p.abs_overlap_with_prev for p in self.points if p.abs_overlap_with_prev is not None]
        )

    @property
    def energies(self) -> np.ndarray:
        return np.array([p.energy for p in self.points])

    @property
    def all_converged(self) -> bool:
        conv = all(p.converged for p in self.points)
        if self.endpoint_converged is not None:
            conv = conv and self.endpoint_converged
        return conv

    def to_dict(self) -> dict:
        d = asdict(self)
        d["points"] = [p.to_dict() for p in self.points]
        return d


# Fallback ladder for a point that will not converge. Each rung keeps as much of the
# continuation as it can; the final rung abandons the warm start entirely and is flagged,
# because it breaks the continuation chain. Every rung is still overlap-checked afterwards,
# so a fallback that lands on a different branch is caught by the continuity test rather
# than silently accepted.
SOLVE_STRATEGIES: tuple[tuple[str, bool, bool, float], ...] = (
    #  name                 warm MO  warm CI  conv_tol_grad multiplier
    ("warm-mo+warm-ci",     True,    True,    1.0),
    ("warm-mo",             True,    False,   1.0),
    ("warm-mo-loose",       True,    False,   20.0),
    ("cold",                False,   False,   1.0),
)


def _solve_point(
    geom: str,
    cas: CasConfig,
    previous: CasWavefunction | None,
    label: str,
) -> tuple[CasWavefunction, dict]:
    """Solve CASSCF at one geometry, warm-started from ``previous`` when given.

    Walks :data:`SOLVE_STRATEGIES` until one converges, and reports which rung was used.
    """
    mol = build_mol(geom, cas.basis, charge=cas.charge, spin=cas.spin)
    info: dict = {"guess_orthonormality_error": None, "max_mo_change": None,
                  "strategy": "initial"}

    if previous is None:
        mf = run_rhf(mol)
        wfn = run_casscf(
            mol, cas.ncas, cas.nelecas,
            conv_tol=cas.conv_tol, conv_tol_grad=cas.conv_tol_grad,
            max_cycle_macro=cas.max_cycle_macro, fix_spin=cas.fix_spin,
            label=label, mf=mf,
        )
        return wfn, info

    mo_guess = transfer_mo(previous.mol, previous.mo_coeff, mol)
    info["guess_orthonormality_error"] = orthonormality_error(mol, mo_guess)
    # RHF here only supplies integrals; warm-start it from the transferred orbitals so it
    # converges in a couple of cycles.
    mf = run_rhf(mol, dm0=_density_from_mo(mo_guess, previous.ncore))

    attempted: list[str] = []
    wfn = None
    for name, warm_mo, warm_ci, grad_mult in SOLVE_STRATEGIES:
        attempted.append(name)
        grad = None if cas.conv_tol_grad is None else cas.conv_tol_grad * grad_mult
        candidate = run_casscf(
            mol, cas.ncas, cas.nelecas,
            mo_guess=mo_guess if warm_mo else None,
            ci0=np.asarray(previous.ci) if warm_ci else None,
            conv_tol=cas.conv_tol, conv_tol_grad=grad,
            max_cycle_macro=cas.max_cycle_macro, fix_spin=cas.fix_spin,
            label=label, mf=mf,
        )
        wfn = candidate if wfn is None else wfn
        if candidate.converged:
            wfn = candidate
            info["strategy"] = name
            break
    else:
        # Nothing converged: keep the first attempt and let the verdict record the failure.
        info["strategy"] = "none-converged"

    info["strategies_tried"] = attempted
    nocc = wfn.ncore + wfn.ncas
    info["max_mo_change"] = float(
        np.abs(np.abs(mo_guess[:, :nocc]) - np.abs(wfn.mo_coeff[:, :nocc])).max()
    )
    return wfn, info


def _density_from_mo(mo: np.ndarray, ncore: int) -> np.ndarray:
    occ = mo[:, :ncore]
    return 2.0 * occ @ occ.T


def traverse_loop(
    loop: Loop,
    cas: CasConfig | None = None,
    cont: ContinuationConfig | None = None,
    geom_fn: Callable[[float, float], str] = formaldimine_geom,
    solve_endpoint: bool = True,
    progress: Callable[[str], None] | None = None,
) -> LoopTraversal:
    """Track the state-specific CASSCF ground state once around ``loop``.

    Returns the full diagnostic record; interpretation is left to
    :func:`berrycasscf.berry.analyse`.
    """
    cas = cas or CasConfig()
    cont = cont or ContinuationConfig()
    say = progress or (lambda _msg: None)
    t_start = time.time()

    pts = loop.points()
    t_vals = loop.t_values
    records: list[PointRecord] = []
    wfns: list[CasWavefunction] = []
    warnings: list[str] = []

    previous: CasWavefunction | None = None
    for k, ((alpha, phi), t) in enumerate(zip(pts, t_vals)):
        tic = time.time()
        wfn, info = _solve_point(geom_fn(alpha, phi), cas, previous, label=f"{loop.name}[{k}]")

        if previous is None:
            raw, gauge, absov = None, 1, None
        else:
            raw = float(cas_overlap(previous, wfn))
            gauge = 1 if raw >= 0.0 else -1
            if gauge == -1:
                wfn = wfn.scaled(-1.0)       # gauge fixing
            absov = abs(raw)

        gap = None
        if cont.compute_casci_gap:
            roots = casci_roots(wfn, nroots=2)
            if roots.size >= 2:
                gap = float(roots[1] - roots[0])

        records.append(
            PointRecord(
                index=k,
                t=float(t),
                alpha=float(alpha),
                phi=float(phi),
                energy=wfn.energy,
                converged=wfn.converged,
                raw_overlap_with_prev=raw,
                gauge_sign=gauge,
                abs_overlap_with_prev=absov,
                casci_gap=gap,
                guess_orthonormality_error=info["guess_orthonormality_error"],
                max_mo_change=info["max_mo_change"],
                strategy=info.get("strategy", "initial"),
                wall_time=time.time() - tic,
            )
        )
        wfns.append(wfn)
        previous = wfn
        say(
            f"  [{loop.name}] point {k:3d}/{loop.n_points}  "
            f"alpha={alpha:7.3f} phi={phi:7.3f}  E={wfn.energy:.10f}  "
            f"conv={wfn.converged}  "
            + ("overlap=  n/a" if absov is None else f"overlap={raw:+.6f}")
            + ("" if gap is None else f"  gap={gap:.4f}")
        )

    # --- close the cycle without re-solving: <Psi_{N-1} | Psi_0> --------------------
    closing_overlap = float(cas_overlap(wfns[-1], wfns[0]))

    # --- endpoint estimator: one more continuation step onto the identical geometry --
    endpoint_overlap = None
    endpoint_converged = None
    if solve_endpoint:
        alpha_c, phi_c = loop.closing_point()
        wfn_n, _ = _solve_point(geom_fn(alpha_c, phi_c), cas, wfns[-1], label=f"{loop.name}[close]")
        raw_n = float(cas_overlap(wfns[-1], wfn_n))
        if raw_n < 0.0:
            wfn_n = wfn_n.scaled(-1.0)
        endpoint_overlap = float(cas_overlap(wfns[0], wfn_n))
        endpoint_converged = bool(wfn_n.converged)
        say(
            f"  [{loop.name}] closing point: last->close overlap={raw_n:+.6f}, "
            f"endpoint <0|N>={endpoint_overlap:+.6f}"
        )
        if abs(raw_n) < cont.min_abs_overlap:
            warnings.append(
                f"closing continuation step has |overlap|={abs(raw_n):.3f} "
                f"< {cont.min_abs_overlap}"
            )

    fallbacks = [p.index for p in records if p.strategy in ("cold", "none-converged")
                 and p.index != 0]
    if fallbacks:
        warnings.append(
            f"continuation chain broken at point(s) {fallbacks}: the warm start was "
            "abandoned there, so those points were not reached by continuation. The "
            "overlap checks still apply, but treat the result with suspicion."
        )
    degraded = [p.index for p in records if p.strategy in ("warm-mo", "warm-mo-loose")]
    if degraded:
        warnings.append(
            f"fell back to a weaker warm start at point(s) {degraded} "
            "(orbitals still continued, CI guess dropped or gradient threshold relaxed)"
        )

    if cont.compute_casci_gap:
        gaps = [p.casci_gap for p in records if p.casci_gap is not None]
        if gaps and min(gaps) < cont.small_gap_warning:
            warnings.append(
                f"root-flipping risk: the tracked solution comes within "
                f"{min(gaps):.4f} Ha of the next root of its own active-space CI problem "
                f"(threshold {cont.small_gap_warning}). This is an in-CAS diagnostic, not "
                f"the physical S0/S1 gap."
            )

    return LoopTraversal(
        loop_name=loop.name,
        loop_centre=tuple(float(x) for x in loop.centre),
        loop_radius=tuple(float(x) for x in loop.radius),
        n_points=loop.n_points,
        cas_label=cas.cas_label,
        basis=cas.basis,
        points=records,
        closing_overlap=closing_overlap,
        endpoint_overlap=endpoint_overlap,
        endpoint_converged=endpoint_converged,
        wall_time=time.time() - t_start,
        warnings=warnings,
    )
