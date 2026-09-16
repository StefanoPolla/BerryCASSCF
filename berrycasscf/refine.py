"""Locating a conical intersection from gap data along a line.

Why this module exists
----------------------
Earlier versions took the apex of a parabola fitted to the *gap*. That is the wrong model.
Near a conical intersection the gap is **linear** in the branching-plane coordinates, so a cut
through the region is a V, and a parabola fitted to three points spanning a V is systematically
dragged toward the middle point.

The right model is available in closed form. For an ideal cone, a straight cut passing at
perpendicular offset ``d`` from the apex gives

    gap(x)^2 = a^2 (x - x0)^2 + b^2,        b = (slope perpendicular) * d

so **the square of the gap is exactly a parabola in x**, whatever the offset. Fitting a parabola
to ``gap**2`` therefore recovers the apex position ``x0`` without bias, and its vertex value
``b`` is the closest approach of the cut to the seam -- a quantity the raw grid minimum only
bounds from above.

The two estimators are both provided so the difference can be reported rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ApexFit:
    """Result of locating a cone apex along a line."""

    position: float           # x0, the apex coordinate
    closest_approach: float   # b, the gap at the apex of the fitted cut (same units as gap)
    slope: float              # a, gap per unit x far from the apex
    n_points: int             # how many points entered the fit
    residual: float           # RMS relative residual of the gap**2 fit
    method: str

    def __repr__(self) -> str:  # pragma: no cover - convenience only
        return (f"ApexFit({self.method}, position={self.position:.3f}, "
                f"closest_approach={self.closest_approach:.4g}, "
                f"n={self.n_points}, residual={self.residual:.2e})")


def parabolic_apex(x, gap) -> ApexFit:
    """Old estimator: parabola through the grid minimum and its two neighbours, on the gap.

    Kept so the bias of the wrong model can be quantified instead of argued about.
    """
    x = np.asarray(x, dtype=float)
    g = np.asarray(gap, dtype=float)
    j = int(np.nanargmin(g))
    if j == 0 or j == len(g) - 1:
        return ApexFit(float(x[j]), float(g[j]), float("nan"), 1, float("nan"), "parabolic-edge")
    y0, y1, y2 = g[j - 1], g[j], g[j + 1]
    denom = y0 - 2 * y1 + y2
    step = float(x[j + 1] - x[j])
    pos = float(x[j]) if abs(denom) < 1e-18 else float(x[j] + 0.5 * (y0 - y2) / denom * step)
    return ApexFit(pos, float(y1), float("nan"), 3, float("nan"), "parabolic")


def cone_apex(x, gap, window: int | None = None, min_points: int = 5) -> ApexFit:
    """Locate the apex by fitting ``gap**2`` with a parabola, which is exact for a cone.

    ``window`` restricts the fit to that many grid steps either side of the grid minimum;
    ``None`` uses every point. A window is worth setting when the cut extends far enough that
    the surfaces stop being conical, which shows up as a large ``residual``.

    Falls back to :func:`parabolic_apex` when there are too few points or the fitted parabola
    opens downward (no apex).
    """
    x = np.asarray(x, dtype=float)
    g = np.asarray(gap, dtype=float)
    ok = np.isfinite(x) & np.isfinite(g)
    x, g = x[ok], g[ok]
    if x.size < min_points:
        return parabolic_apex(x, g)

    j = int(np.argmin(g))
    if window is not None:
        lo, hi = max(0, j - window), min(len(x), j + window + 1)
        xf, gf = x[lo:hi], g[lo:hi]
        if xf.size < min_points:
            xf, gf = x, g
    else:
        xf, gf = x, g

    coef = np.polyfit(xf, gf ** 2, 2)
    a2, b1, c0 = coef
    if a2 <= 0:                       # opens downward: not a cone cut
        return parabolic_apex(x, g)

    x0 = -b1 / (2.0 * a2)
    vertex = c0 - b1 ** 2 / (4.0 * a2)
    closest = float(np.sqrt(max(vertex, 0.0)))
    pred = np.polyval(coef, xf)
    scale = float(np.mean(gf ** 2))
    residual = float(np.sqrt(np.mean((pred - gf ** 2) ** 2)) / scale) if scale > 0 else np.nan
    return ApexFit(float(x0), closest, float(np.sqrt(a2)), int(xf.size), residual, "cone")


def refine_row(scan, axis: str = "phi", window: int | None = 3) -> ApexFit:
    """Apply :func:`cone_apex` along the row of a :class:`~berrycasscf.scan.ScanResult`
    that holds its two-dimensional minimum."""
    gap = scan.gap
    i, _ = np.unravel_index(np.nanargmin(gap), gap.shape)
    if axis == "phi":
        return cone_apex(scan.phis, gap[i], window=window)
    j = int(np.nanargmin(gap[i]))
    return cone_apex(scan.alphas, gap[:, j], window=window)
