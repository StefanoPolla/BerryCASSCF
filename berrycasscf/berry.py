"""Berry-phase extraction and the pass/fail verdict.

Under the real-wavefunction assumption the Berry phase is a Z2 quantity read off from a
*sign*. Two independent estimators are formed from one traversal:

**Product estimator** -- the cyclic product of adjacent overlaps around the closed loop,

    Pi = ( prod_{k=1}^{N-1} <Psi_{k-1}|Psi_k> ) * <Psi_{N-1}|Psi_0>

After the sequential gauge fixing in :mod:`berrycasscf.continuation` every factor except
the last is positive, so ``sign(Pi)`` is carried entirely by the closing overlap. Each
wavefunction appears once as a bra and once as a ket around the cycle, so the arbitrary
per-point signs enter squared and cancel: the estimator is gauge-invariant by construction.
``|Pi|`` carries no phase information and is used only as a continuity diagnostic.

**Endpoint estimator** -- one further continuation step onto the geometry identical to the
start, then ``omega = <Psi_0|Psi_N>``. Because the two geometries coincide exactly, this
overlap involves no change of AO basis.

The two must agree in sign. They cannot disagree for physical reasons, so a disagreement is
reported as a failure of the calculation, not as a result.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np

from .config import ContinuationConfig
from .continuation import LoopTraversal

TRIVIAL = "trivial (0)"
NONTRIVIAL = "non-trivial (pi)"
UNDETERMINED = "undetermined"


@dataclass
class BerryResult:
    """Interpretation of one loop traversal."""

    loop_name: str
    cas_label: str
    basis: str
    n_points: int
    loop_centre: tuple[float, float]
    loop_radius: tuple[float, float]

    product_estimator: float
    endpoint_estimator: float | None
    closing_overlap: float

    berry_phase: str
    status: str                      # "OK" or "FAILED"
    checks: dict[str, bool]
    messages: list[str] = field(default_factory=list)

    min_abs_adjacent_overlap: float = np.nan
    mean_abs_adjacent_overlap: float = np.nan
    min_casci_gap: float | None = None
    energy_range: tuple[float, float] = (np.nan, np.nan)
    wall_time: float = np.nan

    @property
    def ok(self) -> bool:
        return self.status == "OK"

    @property
    def is_nontrivial(self) -> bool | None:
        if self.berry_phase == NONTRIVIAL:
            return True
        if self.berry_phase == TRIVIAL:
            return False
        return None

    def to_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        lines = [
            f"Loop {self.loop_name}  centre={self.loop_centre}  radius={self.loop_radius}  "
            f"N={self.n_points}  {self.cas_label}/{self.basis}",
            f"  product estimator   Pi = {self.product_estimator:+.6f}",
            (
                "  endpoint estimator  w  = n/a"
                if self.endpoint_estimator is None
                else f"  endpoint estimator  w  = {self.endpoint_estimator:+.6f}"
            ),
            f"  min |adjacent overlap| = {self.min_abs_adjacent_overlap:.4f}",
        ]
        if self.min_casci_gap is not None:
            lines.append(f"  min CASCI S0/S1 gap    = {self.min_casci_gap:.4f} Ha")
        lines.append(f"  BERRY PHASE: {self.berry_phase}   [{self.status}]")
        for m in self.messages:
            lines.append(f"    ! {m}")
        return "\n".join(lines)


def analyse(trav: LoopTraversal, cont: ContinuationConfig | None = None) -> BerryResult:
    """Turn a traversal into a Berry phase plus an explicit pass/fail verdict."""
    cont = cont or ContinuationConfig()

    adj = trav.adjacent_abs_overlaps
    product = float(np.prod(adj) * trav.closing_overlap) if adj.size else float(trav.closing_overlap)

    messages: list[str] = list(trav.warnings)
    checks: dict[str, bool] = {}

    # 1. every CASSCF point converged.
    #    Skipped when the configuration says not to require it, which is the single-update
    #    regime of arXiv:2304.06070: there a point is *deliberately* left unconverged and the
    #    verdict must rest on continuity and the endpoint alone.
    if cont.require_converged:
        checks["all_points_converged"] = bool(trav.all_converged)
        if not checks["all_points_converged"]:
            bad = [p.index for p in trav.points if not p.converged]
            messages.append(f"CASSCF did not converge at point(s) {bad}")
    else:
        checks["all_points_converged"] = True
        unconverged = sum(1 for p in trav.points if not p.converged)
        if unconverged:
            messages.append(
                f"{unconverged}/{len(trav.points)} points left unconverged by design "
                "(convergence not required); the verdict rests on continuity and the endpoint"
            )

    # 2. continuity: no adjacent overlap may collapse
    min_abs = float(adj.min()) if adj.size else float("nan")
    checks["continuity"] = bool(adj.size and min_abs >= cont.min_abs_overlap)
    if adj.size and not checks["continuity"]:
        worst = int(np.argmin(adj)) + 1
        messages.append(
            f"continuity lost: |overlap| = {min_abs:.3f} < {cont.min_abs_overlap} "
            f"at step {worst}; increase n_points or the solution jumped branch"
        )

    # 3. the endpoint must be the same physical state as the start
    if trav.endpoint_overlap is None:
        checks["endpoint_is_same_state"] = True
    else:
        mag = abs(trav.endpoint_overlap)
        checks["endpoint_is_same_state"] = bool(mag >= cont.min_endpoint_overlap_magnitude)
        if not checks["endpoint_is_same_state"]:
            messages.append(
                f"endpoint |<Psi_0|Psi_N>| = {mag:.3f} < "
                f"{cont.min_endpoint_overlap_magnitude}: the loop did not return to the "
                "same state, so its sign is not interpretable"
            )

    # 4. the two estimators must agree in sign
    if trav.endpoint_overlap is None:
        checks["estimators_agree"] = True
    else:
        agree = np.sign(product) == np.sign(trav.endpoint_overlap)
        checks["estimators_agree"] = bool(agree)
        if not agree:
            messages.append(
                f"estimators disagree: product {product:+.4f} vs endpoint "
                f"{trav.endpoint_overlap:+.4f}. Under the stated assumptions these must "
                "match; this indicates a gauge or continuity bug, not physics"
            )

    status = "OK" if all(checks.values()) else "FAILED"
    if status == "OK":
        phase = NONTRIVIAL if product < 0 else TRIVIAL
    else:
        phase = UNDETERMINED

    gaps = [p.casci_gap for p in trav.points if p.casci_gap is not None]
    energies = trav.energies

    return BerryResult(
        loop_name=trav.loop_name,
        cas_label=trav.cas_label,
        basis=trav.basis,
        n_points=trav.n_points,
        loop_centre=trav.loop_centre,
        loop_radius=trav.loop_radius,
        product_estimator=product,
        endpoint_estimator=trav.endpoint_overlap,
        closing_overlap=trav.closing_overlap,
        berry_phase=phase,
        status=status,
        checks=checks,
        messages=messages,
        min_abs_adjacent_overlap=min_abs,
        mean_abs_adjacent_overlap=float(adj.mean()) if adj.size else float("nan"),
        min_casci_gap=float(min(gaps)) if gaps else None,
        energy_range=(float(energies.min()), float(energies.max())),
        wall_time=trav.wall_time,
    )


def run_loop(loop, cas=None, cont=None, geom_fn=None, progress=None, solve_endpoint=True):
    """Convenience: traverse a loop and analyse it in one call.

    Returns ``(BerryResult, LoopTraversal)``.
    """
    from .continuation import traverse_loop
    from .geometry import formaldimine_geom

    trav = traverse_loop(
        loop,
        cas=cas,
        cont=cont,
        geom_fn=geom_fn or formaldimine_geom,
        solve_endpoint=solve_endpoint,
        progress=progress,
    )
    return analyse(trav, cont), trav
