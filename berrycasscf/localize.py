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


def restore_probe(record: dict) -> RadiusProbe:
    """Rebuild one evaluated radius from its saved form."""
    known = {f.name for f in dataclass_fields(RadiusProbe)}
    return RadiusProbe(**{k: v for k, v in record.items() if k in known})


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


def merge_centre_records(records: Sequence[dict], centres) -> list[BisectionResult]:
    """Assemble per-centre records into the bisections of one localization.

    Bisections about different centres are independent -- they share no state and no
    intermediate -- so at a large active space, where one centre is a day or more, they are
    run as separate jobs and merged afterwards. Each record must hold exactly one *finished*
    bisection, and the centres must match the requested list in order, so that a merge cannot
    silently combine ellipses from two different constructions, nor a checkpoint written
    part-way through one.
    """
    merged: list[BisectionResult] = []
    for i, (centre, record) in enumerate(zip(centres, records)):
        # A record is also written after every probe, so that a killed run resumes. Those
        # carry complete=False and a bisection with no bracket yet; merging one would quietly
        # contribute a centre that had not finished measuring.
        if record.get("complete", True) is False:
            raise ValueError(f"record {i} is an unfinished checkpoint, not a result")
        bisections = record.get("bisections", [])
        if len(bisections) != 1:
            raise ValueError(f"record {i} holds {len(bisections)} bisections, expected 1")
        stored = bisections[0].get("centre")
        if stored is None or not np.allclose(stored, centre, atol=1e-6):
            raise ValueError(f"record {i} is centred on {stored}, expected {centre}")
        merged.append(restore_bisection(bisections[0]))
    return merged


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
    anchor_factor: float = 2.0,
    max_anchor_moves: int = 3,
    scale_ceiling: float = 4.0,
    scale_floor: float = 1e-3,
    name: str = "L",
    progress: Callable[[str], None] | None = None,
    cached_probes: Sequence[dict] = (),
    on_probe: Callable[[list[RadiusProbe]], None] | None = None,
) -> BisectionResult:
    """Shrink the loop about ``centre`` until the phase turns over.

    ``scale_hi`` must enclose the degeneracy and ``scale_lo`` must not; both are verified
    rather than assumed, and a failure of either is reported instead of being bisected
    through. Bisection is **geometric** (the midpoint of the logarithms), because the
    quantity sought is a distance spanning possibly an order of magnitude and relative
    precision is what matters.

    **Anchors move.** ``scale_hi`` must enclose and ``scale_lo`` must not, but neither is known
    in advance: the first is grown (up to ``scale_ceiling``) while it reports anything but pi,
    the second shrunk (down to ``scale_floor``) while it reports anything but zero, at most
    ``max_anchor_moves`` times each. A refused anchor is usually a grazing loop, and moving it
    away from the seam is exactly what it needs -- the alternative, which this replaced, was to
    abandon the centre with nothing measured.

    **Resuming.** One bisection is days at a large active space, which is longer than a
    walltime, so probes can be carried across runs: ``on_probe`` is called with every probe
    completed so far, and ``cached_probes`` supplies what an earlier run saved. Replay is
    exact rather than approximate — the control flow is a pure function of the verdicts, so
    re-running it with the saved verdicts reproduces the same sequence of radii and then
    continues from where it stopped. A cached probe costs nothing and is not re-solved.
    """
    say = progress or (lambda _m: None)
    done = {round(float(p["scale"]), 9): restore_probe(p) for p in cached_probes}
    # Start the clock behind by whatever the cached probes already spent, so that a resumed
    # bisection reports the cost of the whole measurement rather than of its last session.
    t0 = time.time() - sum(p.wall_time for p in done.values())
    probes: list[RadiusProbe] = []
    if done:
        say(f"  resuming with {len(done)} probes already evaluated: "
            + ", ".join(f"{k:.4f}->{v.verdict}" for k, v in sorted(done.items())))

    def probe(scale: float) -> RadiusProbe:
        key = round(float(scale), 9)
        cached = key in done
        p = done[key] if cached else evaluate_radius(
            centre, shape, scale, cas=cas, cont=cont, settings=settings,
            geom_fn=geom_fn, name=name, progress=say)
        probes.append(p)
        if not cached and on_probe is not None:
            on_probe(probes)
        return p

    say(f"bisecting about ({centre[0]:.3f}, {centre[1]:.3f}) with shape {tuple(shape)}")

    # Move the two anchors until they say what a bracket needs, instead of giving up the
    # moment they do not. Both moves are monotone in the right direction:
    #
    #   outer anchor, wants pi.  A zero means the loop is too small to enclose anything, and a
    #                            refusal usually means it is grazing -- growing fixes both.
    #   inner anchor, wants 0.   A pi means even this loop encloses the object, and a refusal
    #                            again means grazing -- shrinking fixes both.
    #
    # Without this a single refused anchor ends the bisection with nothing measured, which is
    # how ethylene CAS(6,6) spent 18 minutes and returned "not bracketed" while its outer loop
    # had cleanly reported pi.
    top = probe(scale_hi)
    for _ in range(max_anchor_moves):
        if top.verdict == PI or scale_hi * anchor_factor > scale_ceiling:
            break
        scale_hi *= anchor_factor
        say(f"  outer anchor was {top.verdict}; growing it to {scale_hi:.4f}")
        top = probe(scale_hi)

    bottom = probe(scale_lo)
    for _ in range(max_anchor_moves):
        if bottom.verdict == ZERO or scale_lo / anchor_factor < scale_floor:
            break
        scale_lo /= anchor_factor
        say(f"  inner anchor was {bottom.verdict}; shrinking it to {scale_lo:.4f}")
        bottom = probe(scale_lo)

    lo = scale_lo if bottom.verdict == ZERO else None
    hi = scale_hi if top.verdict == PI else None
    if hi is None:
        say(f"  outer loop still does not report pi ({top.verdict}) at scale {scale_hi:.4f}; "
            f"nothing to bracket")
    if lo is None and bottom.verdict == PI:
        say(f"  inner loop still reports pi at scale {scale_lo:.4f}: the degeneracy is closer "
            f"to the centre than any loop tried")

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
