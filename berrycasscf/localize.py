"""Locating a degeneracy with loop transport alone, by shrinking and moving the loop.

The gap scan finds an intersection by *looking* for it: it evaluates the state-averaged
S1-S0 gap over a region and reports where it is smallest. Loop transport never evaluates a
gap. It answers one bit -- "does this loop enclose a degeneracy?" -- so a position has to be
built out of bits.

Bisection on the radius
-----------------------
Scale a loop about a fixed centre ``c`` by a factor ``s``. Writing the loop's shape as
``(Ra, Rp)``, the point ``x`` lies inside the scaled loop exactly when

    rho(x; c) = sqrt( ((x0-c0)/Ra)^2 + ((x1-c1)/Rp)^2 )  <  s

so the Berry phase, as a function of ``s``, is 0 below ``rho`` and pi above it. **The
transition radius is the elliptical distance from the centre to the degeneracy.** Bisecting
on ``s`` measures it. One centre therefore does not give a position -- it gives a curve (an
ellipse) that the degeneracy lies on.

Triangulation from several centres
----------------------------------
Repeat from a second centre and intersect the two ellipses. Rescaling the coordinates by
``(Ra, Rp)`` turns both into circles, so this is the classic circle-circle intersection and
generically has **two** solutions, mirror images in the line joining the centres. A third
centre off that line picks one; with two centres the caller must break the tie by symmetry
or by prior knowledge, and both candidates are returned so the choice is explicit.

Why the answer comes back as a bracket
--------------------------------------
Exactly at ``s = rho`` the loop runs through the degeneracy, and no discretization can
transport a state through it. Adaptive stepping turns that into an honest report rather than
a wrong sign: the step size collapses to its floor and the run is refused
(:mod:`berrycasscf.adaptive`). So bisection converges to a bracket -- the largest radius
that says 0 and the smallest that says pi -- separated by a band of refusals whose width is
the method's resolution. That band is a feature: its half-width is the uncertainty on the
position, measured rather than assumed.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, asdict, field, fields as dataclass_fields
from typing import Callable, Sequence

import numpy as np

from .adaptive import traverse_loop_adaptive
from .berry import NONTRIVIAL, TRIVIAL, analyse
from .config import AdaptiveConfig, CasConfig, ContinuationConfig
from .geometry import Loop, formaldimine_geom

PI = "pi"
ZERO = "zero"
UNDETERMINED = "undetermined"

# Two adaptive settings that must agree before a radius is believed. The stability criterion
# of the uniform study was "at least two discretizations, all passing, same phase"; for an
# adaptive walk the analogue is two different step-control settings, which produce genuinely
# different point sets rather than two rescalings of the same one.
DEFAULT_SETTINGS: tuple[dict, ...] = (
    {"d_max": 0.10, "target_mismatch": 0.02},
    {"d_max": 0.07, "target_mismatch": 0.008},
)


def elliptical_radius(point, centre, shape) -> float:
    """``rho(x; c)``: distance to ``centre`` in units of the loop's own semi-axes."""
    return float(math.hypot((point[0] - centre[0]) / shape[0],
                            (point[1] - centre[1]) / shape[1]))


@dataclass
class RadiusProbe:
    """One radius, evaluated with every required setting."""

    scale: float
    radius: tuple[float, float]
    verdict: str
    reason: str
    runs: list[dict] = field(default_factory=list)
    cost_micro: int = 0
    wall_time: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class BisectionResult:
    """The transition radius for one centre, as a bracket."""

    centre: tuple[float, float]
    shape: tuple[float, float]
    cas_label: str
    probes: list[RadiusProbe]
    lo: float | None                 # largest scale reporting 0
    hi: float | None                 # smallest scale reporting pi
    undetermined: list[float] = field(default_factory=list)
    total_micro: int = 0
    wall_time: float = 0.0

    @property
    def rho(self) -> float | None:
        """Midpoint of the bracket: the measured elliptical distance to the degeneracy."""
        if self.lo is None or self.hi is None:
            return None
        return 0.5 * (self.lo + self.hi)

    @property
    def rho_uncertainty(self) -> float | None:
        if self.lo is None or self.hi is None:
            return None
        return 0.5 * (self.hi - self.lo)

    @property
    def bracketed(self) -> bool:
        return self.lo is not None and self.hi is not None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["rho"] = self.rho
        d["rho_uncertainty"] = self.rho_uncertainty
        return d

    def summary(self) -> str:
        if not self.bracketed:
            return (f"centre {self.centre}: NOT bracketed "
                    f"(lo={self.lo}, hi={self.hi}, {len(self.probes)} probes)")
        return (f"centre ({self.centre[0]:.3f}, {self.centre[1]:.3f}): "
                f"rho = {self.rho:.4f} +- {self.rho_uncertainty:.4f}  "
                f"[0 at {self.lo:.4f}, pi at {self.hi:.4f}]  "
                f"{len(self.probes)} probes, {self.total_micro} micro")


def restore_bisection(record: dict) -> BisectionResult:
    """Rebuild a finished bisection from the dict :meth:`BisectionResult.to_dict` wrote.

    ``rho`` and ``rho_uncertainty`` are properties rather than fields and are dropped; the
    probes stay as plain dicts, which is all that triangulation and reporting read.
    """
    known = {f.name for f in dataclass_fields(BisectionResult)}
    return BisectionResult(**{k: v for k, v in record.items() if k in known})


def resume_bisections(saved: Sequence[dict], centres) -> list[BisectionResult]:
    """The leading centres of ``centres`` that ``saved`` already covers.

    A bisection is hours per centre and days at a large active space, so a run that is
    killed part-way must not start over. Saved centres are matched against the requested
    list *position by position*, and the first mismatch ends the match: a record made with
    different centres is not partially reusable, because the ellipses it measured belong to
    a different construction.
    """
    done: list[BisectionResult] = []
    for centre, record in zip(centres, saved):
        stored = record.get("centre")
        if stored is None or not np.allclose(stored, centre, atol=1e-6):
            break
        done.append(restore_bisection(record))
    return done


def evaluate_radius(
    centre,
    shape,
    scale: float,
    cas: CasConfig | None = None,
    cont: ContinuationConfig | None = None,
    settings: Sequence[dict] = DEFAULT_SETTINGS,
    geom_fn: Callable[[float, float], str] = formaldimine_geom,
    name: str = "L",
    progress: Callable[[str], None] | None = None,
) -> RadiusProbe:
    """Berry phase of the loop at ``scale``, believed only if every setting agrees."""
    say = progress or (lambda _m: None)
    tic = time.time()
    radius = (shape[0] * scale, shape[1] * scale)
    loop = Loop(f"{name}@{scale:.4f}", tuple(centre), radius)

    runs, phases, micro = [], [], 0
    for setting in settings:
        trav = traverse_loop_adaptive(loop, cas=cas, cont=cont,
                                      adaptive=AdaptiveConfig(**setting),
                                      geom_fn=geom_fn)
        res = analyse(trav, cont)
        ad = trav.adaptive or {}
        micro += trav.total_micro
        runs.append({
            "setting": setting, "status": res.status, "phase": res.berry_phase,
            "min_overlap": float(res.min_abs_adjacent_overlap),
            "n_points": trav.n_points, "n_rejected": ad.get("n_rejected", 0),
            "closed": ad.get("closed", True), "hit_step_floor": ad.get("hit_step_floor", False),
            "micro": trav.total_micro, "wall": trav.wall_time,
            "product": float(res.product_estimator), "endpoint": res.endpoint_estimator,
            "messages": res.messages,
        })
        phases.append(res.berry_phase if res.status == "OK" else None)

    if any(p is None for p in phases):
        floored = any(r["hit_step_floor"] for r in runs)
        verdict, reason = UNDETERMINED, (
            "step floor reached: shrinking the step did not restore continuity, so the loop "
            "passes very near something the walk cannot transport through"
            if floored else "at least one setting failed its checks")
    elif len(set(phases)) > 1:
        verdict, reason = UNDETERMINED, "the settings disagree on the phase"
    elif phases[0] == NONTRIVIAL:
        verdict, reason = PI, "all settings agree: pi"
    else:
        verdict, reason = ZERO, "all settings agree: 0"

    probe = RadiusProbe(scale=float(scale), radius=radius, verdict=verdict, reason=reason,
                        runs=runs, cost_micro=micro, wall_time=time.time() - tic)
    say(f"  scale {scale:.4f} (radius {radius[0]:.3f}, {radius[1]:.3f}): "
        f"{verdict.upper():13s} {reason}  [{micro} micro, {probe.wall_time:.1f}s]")
    return probe


def bisect_radius(
    centre,
    shape,
    cas: CasConfig | None = None,
    cont: ContinuationConfig | None = None,
    settings: Sequence[dict] = DEFAULT_SETTINGS,
    geom_fn: Callable[[float, float], str] = formaldimine_geom,
    scale_hi: float = 1.0,
    scale_lo: float = 0.05,
    tol: float = 0.02,
    max_probes: int = 12,
    name: str = "L",
    progress: Callable[[str], None] | None = None,
) -> BisectionResult:
    """Shrink the loop about ``centre`` until the phase turns over.

    ``scale_hi`` must enclose the degeneracy and ``scale_lo`` must not; both are verified
    rather than assumed, and a failure of either is reported instead of being bisected
    through. Bisection is **geometric** (the midpoint of the logarithms), because the
    quantity sought is a distance spanning possibly an order of magnitude and relative
    precision is what matters.
    """
    say = progress or (lambda _m: None)
    t0 = time.time()
    probes: list[RadiusProbe] = []

    def probe(scale: float) -> RadiusProbe:
        p = evaluate_radius(centre, shape, scale, cas=cas, cont=cont, settings=settings,
                            geom_fn=geom_fn, name=name, progress=say)
        probes.append(p)
        return p

    say(f"bisecting about ({centre[0]:.3f}, {centre[1]:.3f}) with shape {tuple(shape)}")
    top, bottom = probe(scale_hi), probe(scale_lo)
    lo = scale_lo if bottom.verdict == ZERO else None
    hi = scale_hi if top.verdict == PI else None
    if hi is None:
        say(f"  outer loop does not report pi ({top.verdict}); nothing to bracket")
    if lo is None and bottom.verdict == PI:
        say(f"  inner loop already reports pi: the degeneracy is inside scale {scale_lo}")

    # Plain bisection gives up the moment a midpoint comes back undetermined, and that is
    # exactly what the first midpoint tends to do: it lands in the refusal band around the
    # transition. So the search tracks THREE things -- the largest scale reporting 0, the
    # smallest reporting pi, and the band of refusals in between -- and narrows the two
    # usable gaps alternately:
    #
    #     lo ....... [ undetermined band ] ....... hi
    #        ^ raise this                ^ lower this
    #
    # Every probe therefore does something: it either classifies (moving lo up or hi down) or
    # widens the refusal band, which shrinks the interval the next probe will target. The
    # answer is the bracket (lo, hi); the band's width is the method's resolution, measured
    # rather than assumed.
    raise_lo = True
    while lo is not None and hi is not None and len(probes) < max_probes:
        if (hi - lo) / hi <= tol:
            break
        band = [p.scale for p in probes if p.verdict == UNDETERMINED and lo < p.scale < hi]
        gap_below = (lo, min(band) if band else hi)
        gap_above = (max(band) if band else lo, hi)
        widths = {"below": gap_below[1] / gap_below[0], "above": gap_above[1] / gap_above[0]}
        # Alternate, but skip a side that is already tighter than the tolerance.
        side = "below" if raise_lo else "above"
        if widths[side] - 1.0 <= tol:
            side = "above" if side == "below" else "below"
            if widths[side] - 1.0 <= tol:
                say("  both sides of the refusal band are inside tolerance; stopping")
                break
        raise_lo = not raise_lo

        a, b = gap_below if side == "below" else gap_above
        mid = math.sqrt(a * b)
        p = probe(mid)
        if p.verdict == PI:
            hi = min(hi, mid)
        elif p.verdict == ZERO:
            lo = max(lo, mid)
        else:
            say("  refused: the band widens and the next probe targets the other side")

    result = BisectionResult(
        centre=tuple(float(c) for c in centre),
        shape=tuple(float(s) for s in shape),
        cas_label=(cas or CasConfig()).cas_label,
        probes=probes, lo=lo, hi=hi,
        undetermined=[p.scale for p in probes if p.verdict == UNDETERMINED],
        total_micro=sum(p.cost_micro for p in probes),
        wall_time=time.time() - t0,
    )
    say("  " + result.summary())
    return result


def _circle_intersections(c1, r1, c2, r2):
    """Both intersection points of two circles, or [] if they do not meet."""
    (x1, y1), (x2, y2) = c1, c2
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy)
    if d < 1e-12 or d > r1 + r2 or d < abs(r1 - r2):
        return []
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    h2 = r1 * r1 - a * a
    h = math.sqrt(max(h2, 0.0))
    xm, ym = x1 + a * dx / d, y1 + a * dy / d
    return [(xm + h * dy / d, ym - h * dx / d), (xm - h * dy / d, ym + h * dx / d)]


@dataclass
class TriangulationResult:
    candidates: list[tuple[float, float]]
    residual: float
    shape: tuple[float, float]
    centres: list[tuple[float, float]]
    rhos: list[float]
    chosen: tuple[float, float] | None = None
    note: str = ""
    # RMS of the bisection bracket half-widths: the precision the residual has to be judged
    # against. Comparing raw residuals between runs is misleading, because a run whose
    # brackets are loose can misfit by a lot and still be perfectly consistent.
    uncertainty: float = float("nan")

    @property
    def residual_over_uncertainty(self) -> float:
        """Misfit in units of the measurement precision. Above ~1 the data are inconsistent.

        This, not the raw residual, is the consistency criterion. Measured on formaldimine:
        CAS(2,2) 1.57 (inconsistent -- no single degeneracy explains three centres),
        CAS(4,4) 0.29, CAS(6,6) 0.76. The raw residuals are 0.197, 0.0088 and 0.0574, which
        would wrongly suggest CAS(6,6) is 6.5x worse than CAS(4,4) rather than equally
        consistent but less precisely measured.
        """
        if not np.isfinite(self.uncertainty) or self.uncertainty <= 0:
            return float("nan")
        return float(self.residual / self.uncertainty)

    @property
    def consistent(self) -> bool | None:
        r = self.residual_over_uncertainty
        return None if not np.isfinite(r) else bool(r <= 1.0)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["residual_over_uncertainty"] = self.residual_over_uncertainty
        d["consistent"] = self.consistent
        return d


def triangulate(
    results: Sequence[BisectionResult],
    prefer: tuple[float, float] | None = None,
) -> TriangulationResult:
    """Intersect the ellipses from several bisections into a position.

    Rescaling by the loop's semi-axes turns each ellipse into a circle, so two centres give
    two mirror-image candidates and three or more generically give one. ``prefer`` picks the
    candidate nearest a point the caller has independent reason to expect -- a symmetry line,
    say -- and the choice is recorded in ``note`` rather than hidden.
    """
    usable = [r for r in results if r.bracketed]
    if len(usable) < 2:
        raise ValueError(f"triangulation needs at least two bracketed bisections, "
                         f"got {len(usable)}")
    shape = usable[0].shape
    if any(tuple(r.shape) != tuple(shape) for r in usable):
        raise ValueError("all bisections must use the same loop shape to be intersected")

    # Work in units of the semi-axes, where every ellipse is a unit-shaped circle.
    scaled = [(r.centre[0] / shape[0], r.centre[1] / shape[1]) for r in usable]
    rhos = [r.rho for r in usable]

    cands = _circle_intersections(scaled[0], rhos[0], scaled[1], rhos[1])
    if not cands:
        return TriangulationResult(candidates=[], residual=float("nan"), shape=shape,
                                   centres=[r.centre for r in usable], rhos=rhos,
                                   note="the measured ellipses do not intersect: the "
                                        "bisections are inconsistent with a single "
                                        "degeneracy")

    def residual(u) -> float:
        return float(np.sqrt(np.mean([
            (math.hypot(u[0] - c[0], u[1] - c[1]) - rho) ** 2
            for c, rho in zip(scaled, rhos)
        ])))

    # With three or more centres the extra ellipses select among the two candidates, and
    # the residual reports how well a single point explains all of them.
    best = min(cands, key=residual)
    cands_real = [(u[0] * shape[0], u[1] * shape[1]) for u in cands]
    chosen = (best[0] * shape[0], best[1] * shape[1])
    note = ""
    if len(usable) == 2:
        if prefer is not None:
            chosen = min(cands_real,
                         key=lambda p: math.hypot(p[0] - prefer[0], p[1] - prefer[1]))
            note = (f"two centres give two mirror-image candidates; chose the one nearer "
                    f"{prefer} on external grounds")
        else:
            note = ("two centres give two mirror-image candidates and nothing here "
                    "distinguishes them; add a third centre off the line joining these two")
    else:
        note = f"{len(usable)} centres; residual is the misfit of a single position"

    half_widths = [r.rho_uncertainty for r in usable if r.rho_uncertainty is not None]
    unc = float(np.sqrt(np.mean(np.square(half_widths)))) if half_widths else float("nan")
    return TriangulationResult(candidates=cands_real, residual=residual(best), shape=shape,
                               centres=[r.centre for r in usable], rhos=rhos,
                               chosen=chosen, note=note, uncertainty=unc)
