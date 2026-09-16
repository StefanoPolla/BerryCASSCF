"""Adaptive step control for walking a loop.

A uniform discretization has to be fine enough for the *hardest* arc of the loop, so it
oversamples everywhere else. The walk already measures its own error at every step --
``m = 1 - |<Psi_prev|Psi_new>|``, the same quantity the continuity check reads -- so the
step size can be steered by it for free.

The controller
--------------
``m`` behaves like ``d^2`` for a smooth transported state, so a step that achieved ``m``
with size ``d`` should use

    d_new = d * sqrt(target / m)

to hit ``target`` next time. That is a proportional controller in ``log d``, clipped to
``[d_min, d_max]`` with a growth cap so it cannot leap into a region it has not probed.
Steering toward a target is deliberately different from rejecting at a floor: rejection
costs a full CASSCF solve and there is no cheap predictor, so the walk is cheapest when it
sits just inside the accept threshold rather than bouncing off it.

Reaching ``d_min`` is a **result, not a safeguard**. If shrinking the step restores
continuity, the problem was discretization. If it does not, the loop is passing through
something -- which is exactly the radius-8 versus radius-6 distinction found by hand in
``docs/findings.md`` §5, and here it is discovered automatically.

Closing the loop
----------------
The walk ends *exactly* on ``t = 1`` (the final step is clipped), so the closing geometry is
identical to the start and the endpoint estimator stays exact. That last point plays the
role of the separate closing solve in :func:`berrycasscf.continuation.traverse_loop`.

The generic seam
----------------
:func:`walk_adaptive` knows nothing about quantum chemistry: it is given a ``solve`` and an
``overlap`` callable. This is not abstraction for its own sake -- it lets the identical
accept/reject/floor logic run on the analytically solvable Jahn-Teller model in
:mod:`berrycasscf.toy` at zero cost, which is the only way to test a controller without
letting a CASSCF failure masquerade as a control failure.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, asdict, field
from typing import Any, Callable

import numpy as np

from .casscf import build_mol
from .config import AdaptiveConfig, CasConfig, ContinuationConfig
from .continuation import LoopTraversal, PointRecord, _solve_point
from .geometry import Loop, formaldimine_geom
from .overlap import cas_overlap


@dataclass
class StepEvent:
    """One attempted step, accepted or not. The rejected ones are the interesting record."""

    t_from: float
    t_to: float
    d: float
    raw_overlap: float | None
    mismatch: float                 # 1 - |overlap|
    accepted: bool
    reason: str                     # "accepted" | "mismatch" | "unconverged" | "step-floor"
    d_next: float
    # Solver work this trial cost, so rejections can be charged rather than hidden.
    n_macro: int = 0
    n_micro: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AdaptiveWalk:
    """The accepted chain plus the full history of what it took to get there."""

    t_values: list[float]
    states: list[Any] = field(repr=False)
    infos: list[dict] = field(repr=False)
    # Every attempted step, rejected ones included: the cost record and the diagnostic.
    events: list[StepEvent] = field(default_factory=list)
    # Only the steps that were kept, one per chain point after the first. Consumers must use
    # THIS to attribute an overlap to a point: indexing `events` by point number silently
    # slips by one for every rejection, which once reported an adjacent overlap of 0.79 on a
    # walk whose accepted steps were all above 0.90.
    accepted_events: list[StepEvent] = field(default_factory=list)
    hit_step_floor: bool = False
    hit_budget: bool = False

    @property
    def closed(self) -> bool:
        """True when the walk actually reached t = 1 rather than stopping early."""
        return bool(self.t_values) and abs(self.t_values[-1] - 1.0) < 1e-9

    @property
    def n_rejected(self) -> int:
        return sum(1 for e in self.events if not e.accepted)

    @property
    def n_accepted(self) -> int:
        return sum(1 for e in self.events if e.accepted)

    @property
    def step_sizes(self) -> np.ndarray:
        return np.array([e.d for e in self.accepted_events])


def control_step(d: float, mismatch: float, cfg: AdaptiveConfig) -> float:
    """Next step size from the achieved mismatch, assuming ``m ~ d^2``.

    Pure and side-effect free so it can be tested on its own.
    """
    m = max(float(mismatch), cfg.mismatch_floor)
    factor = math.sqrt(cfg.target_mismatch / m)
    factor = min(max(factor, cfg.max_shrink), cfg.growth_cap)
    return float(min(max(d * factor, cfg.d_min), cfg.d_max))


def walk_adaptive(
    solve: Callable[[float, Any | None], tuple[Any, dict]],
    overlap: Callable[[Any, Any], float],
    cfg: AdaptiveConfig | None = None,
    negate: Callable[[Any], Any] | None = None,
    progress: Callable[[str], None] | None = None,
) -> AdaptiveWalk:
    """Walk t from 0 to 1, choosing every step size from the measured continuity.

    ``solve(t, previous)`` returns ``(state, info)``; ``overlap(a, b)`` returns the *signed*
    overlap; ``negate`` applies the gauge flip (defaults to ``-state``).

    **Acceptance depends on continuity alone, never on whether the solver reported
    convergence.** Those are different failures needing different responses, and conflating
    them is a trap: a smaller step makes the warm start *better*, so a point whose solver is
    struggling is no likelier to converge once the step shrinks -- the walk just shrinks
    forever. Measured on formaldimine CAS(2,2): steps with a mismatch of 0.011, far inside
    any threshold, were rejected for non-convergence, and one probe ran for over ten minutes
    without finishing. Convergence is handled by the solver fallback ladder and then by the
    ``all_points_converged`` check in :func:`berrycasscf.berry.analyse`, which refuses the
    run without the controller having to guess.

    On rejection the previous accepted state is restored untouched -- the trial is discarded
    entirely, never carried forward. Getting this wrong is silent, so it is done in one place.
    """
    cfg = cfg or AdaptiveConfig()
    say = progress or (lambda _m: None)
    flip = negate or (lambda s: -s)

    state, info = solve(0.0, None)
    walk = AdaptiveWalk(t_values=[0.0], states=[state], infos=[info])

    t = 0.0
    d = cfg.d_max
    while 1.0 - t > 1e-12:
        if len(walk.t_values) >= cfg.max_points or walk.n_rejected >= cfg.max_rejections:
            walk.hit_budget = True
            say(f"  [adaptive] budget exhausted at t={t:.6f} "
                f"({len(walk.t_values)} points, {walk.n_rejected} rejections)")
            break

        d_try = min(d, 1.0 - t)          # clip so the walk lands exactly on t = 1
        t_new = t + d_try
        trial, trial_info = solve(t_new, walk.states[-1])

        raw = float(overlap(walk.states[-1], trial))
        mismatch = 1.0 - abs(raw)
        d_next = control_step(d_try, mismatch, cfg)
        cost = (int(trial_info.get("n_macro", 0)), int(trial_info.get("n_micro", 0)))

        if mismatch <= cfg.accept_mismatch:
            if raw < 0.0:
                trial = flip(trial)      # sequential gauge fixing, as in the uniform walk
            walk.t_values.append(t_new)
            walk.states.append(trial)
            walk.infos.append(trial_info)
            ev = StepEvent(t, t_new, d_try, raw, mismatch, True, "accepted", d_next, *cost)
            walk.events.append(ev)
            walk.accepted_events.append(ev)
            t, d = t_new, d_next
            say(f"  [adaptive] accept t={t:.6f}  d={d_try:.5f}  "
                f"overlap={raw:+.6f}  next d={d:.5f}")
            continue

        # --- rejected: discard the trial, shrink, retry from the same accepted state ---
        # `trial` is simply dropped here. The previous accepted state is still walk.states[-1]
        # and was never modified, which is the whole correctness requirement of a rejection.
        reason = "mismatch"
        at_floor = d_try <= cfg.d_min * (1.0 + 1e-9)
        walk.events.append(StepEvent(t, t_new, d_try, raw, mismatch, False,
                                     "step-floor" if at_floor else reason, d_next, *cost))
        say(f"  [adaptive] REJECT t={t:.6f}->{t_new:.6f}  d={d_try:.5f}  "
            f"overlap={raw:+.6f} ({reason})  retry with d={d_next:.5f}")
        if at_floor:
            # Shrinking did not restore continuity. This is the informative failure.
            walk.hit_step_floor = True
            say(f"  [adaptive] step floor d_min={cfg.d_min} reached without recovering "
                f"continuity (mismatch {mismatch:.4f}): the loop is passing through "
                f"something, not merely undersampled")
            break
        d = d_next

    return walk


def traverse_loop_adaptive(
    loop: Loop,
    cas: CasConfig | None = None,
    cont: ContinuationConfig | None = None,
    adaptive: AdaptiveConfig | None = None,
    geom_fn: Callable[[float, float], str] = formaldimine_geom,
    progress: Callable[[str], None] | None = None,
) -> LoopTraversal:
    """Adaptive counterpart of :func:`berrycasscf.continuation.traverse_loop`.

    Returns the same :class:`LoopTraversal` record, so :func:`berrycasscf.berry.analyse`
    and every downstream plot work unchanged. ``loop.n_points`` is ignored -- the walk
    chooses its own points.
    """
    cas = cas or CasConfig()
    cont = cont or ContinuationConfig()
    adaptive = adaptive or AdaptiveConfig()
    say = progress or (lambda _m: None)
    t_start = time.time()

    timings: list[float] = []

    def solve(t: float, previous):
        alpha, phi = loop.point_at(t)
        tic = time.time()
        wfn, info = _solve_point(geom_fn(alpha, phi), cas, previous,
                                 label=f"{loop.name}[t={t:.4f}]",
                                 use_fallback=cont.use_fallback)
        info = dict(info)
        info["converged"] = bool(wfn.converged)
        info["alpha"], info["phi"] = float(alpha), float(phi)
        info["wall_time"] = time.time() - tic
        timings.append(info["wall_time"])
        return wfn, info

    walk = walk_adaptive(solve, lambda a, b: float(cas_overlap(a, b)),
                         cfg=adaptive, negate=lambda w: w.scaled(-1.0), progress=say)

    # The point at t = 1 is the closing geometry: it plays the role of the separate closing
    # solve in the uniform walk, so it is the endpoint estimator rather than a chain point.
    chain_t = walk.t_values
    endpoint_overlap = None
    endpoint_converged = None
    if walk.closed:
        chain_t = walk.t_values[:-1]
        endpoint_overlap = float(cas_overlap(walk.states[0], walk.states[-1]))
        endpoint_converged = bool(walk.infos[-1].get("converged", True))

    records: list[PointRecord] = []
    for k, t in enumerate(chain_t):
        info = walk.infos[k]
        wfn = walk.states[k]
        ev = walk.accepted_events[k - 1] if k > 0 else None
        gap = None
        if cont.compute_casci_gap:
            from .casscf import casci_roots
            roots = casci_roots(wfn, nroots=2)
            if roots.size >= 2:
                gap = float(roots[1] - roots[0])
        records.append(PointRecord(
            index=k, t=float(t), alpha=info["alpha"], phi=info["phi"],
            energy=wfn.energy, converged=bool(wfn.converged),
            raw_overlap_with_prev=None if ev is None else ev.raw_overlap,
            gauge_sign=1 if ev is None or (ev.raw_overlap or 0.0) >= 0 else -1,
            abs_overlap_with_prev=None if ev is None else abs(ev.raw_overlap),
            casci_gap=gap,
            guess_orthonormality_error=info.get("guess_orthonormality_error"),
            max_mo_change=info.get("max_mo_change"),
            strategy=info.get("strategy", "initial"),
            n_macro=int(info.get("n_macro", 0)), n_micro=int(info.get("n_micro", 0)),
            wall_time=float(info.get("wall_time", 0.0)),
        ))

    closing_overlap = float(cas_overlap(walk.states[len(chain_t) - 1], walk.states[0]))

    warnings: list[str] = []
    if not walk.closed:
        warnings.append(
            "the adaptive walk did not reach t = 1, so the loop was never closed and no "
            "Berry phase is defined for it"
        )
    if walk.hit_step_floor:
        warnings.append(
            f"step floor reached: shrinking to d_min={adaptive.d_min} did not restore "
            "continuity, so the problem is not discretization. Either the loop passes "
            "through or very near a degeneracy, or the solver is switching between nearby "
            "stationary solutions; the walk cannot tell these apart and does not guess"
        )
    if walk.hit_budget:
        warnings.append(
            f"adaptive budget exhausted ({len(walk.t_values)} points, "
            f"{walk.n_rejected} rejections)"
        )
    degraded = [r.index for r in records if r.strategy in ("warm-mo", "warm-mo-loose")]
    if degraded:
        warnings.append(
            f"fell back to a weaker warm start at point(s) {degraded} "
            "(orbitals still continued, CI guess dropped or gradient threshold relaxed)"
        )

    # Rejected trials are real work and are charged to the run: an adaptive walk that
    # rejected ten steps did ten extra solves, and a cost comparison that hid them would
    # flatter the method.
    rej = [e for e in walk.events if not e.accepted]
    rejected_wall = sum(timings) - sum(r.wall_time for r in records) - (
        walk.infos[-1].get("wall_time", 0.0) if walk.closed else 0.0)

    trav = LoopTraversal(
        loop_name=loop.name,
        loop_centre=tuple(float(x) for x in loop.centre),
        loop_radius=tuple(float(x) for x in loop.radius),
        n_points=len(records),
        cas_label=cas.cas_label,
        basis=cas.basis,
        points=records,
        closing_overlap=closing_overlap,
        endpoint_overlap=endpoint_overlap,
        endpoint_converged=endpoint_converged,
        wall_time=time.time() - t_start,
        warnings=warnings,
        rejected_macro=sum(e.n_macro for e in rej),
        rejected_micro=sum(e.n_micro for e in rej),
        adaptive={
            "config": adaptive.to_dict(),
            "closed": walk.closed,
            "hit_step_floor": walk.hit_step_floor,
            "hit_budget": walk.hit_budget,
            "n_rejected": len(rej),
            "rejected_wall_time": float(max(rejected_wall, 0.0)),
            "events": [e.to_dict() for e in walk.events],
            "t_values": [float(t) for t in walk.t_values],
        },
    )
    return trav
