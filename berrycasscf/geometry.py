"""Molecular geometries and nuclear-coordinate loops.

Kept free of any electronic-structure dependency so it can be tested instantly.

The formaldimine parameterization and the CI-enclosing loop are reproduced verbatim
from the ``auto_oo`` Berry-phase tutorial; see ``docs/provenance.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Frozen internal coordinates of formaldimine: N-C, C-H, N-H distances (Angstrom)
# and the HCN angle (degrees). Verbatim from the auto_oo tutorial, cell 6.
FORMALDIMINE_FROZEN = (1.498047, 1.066797, 0.987109, 118.359375)


def formaldimine_geom(alpha: float, phi: float) -> str:
    """Z-matrix for formaldimine as a function of the two loop coordinates.

    Parameters
    ----------
    alpha : HNC bending angle, degrees.
    phi : HNCH dihedral angle, degrees.

    All other internal coordinates are frozen at ``FORMALDIMINE_FROZEN``.
    Distances are Angstrom, angles degrees.
    """
    variables = list(FORMALDIMINE_FROZEN) + [alpha, phi]
    return """
                    N
                    C 1 {0}
                    H 2 {1}  1 {3}
                    H 2 {1}  1 {3} 3 180
                    H 1 {2}  2 {4} 3 {5}
                    """.format(*variables)


@dataclass(frozen=True)
class Loop:
    """A closed ellipse in a two-dimensional plane of nuclear coordinates.

    Traversed once for t in [0, 1). The two coordinates are generic: for formaldimine they
    are (alpha, phi) in degrees, for fulvene (exocyclic bond length, methylene torsion).

    ``alpha(t) = centre[0] + radius[0] * cos(2*pi*t + phase)``
    ``phi(t)   = centre[1] + radius[1] * sin(2*pi*t + phase)``

    This is the auto_oo tutorial's parameterization (cell 8). ``n_points`` counts
    *distinct* geometries: the closing point t=1 coincides with t=0 and is not stored
    in :meth:`points`, but is returned by :meth:`closing_point`.
    """

    name: str
    centre: tuple[float, float]
    radius: tuple[float, float] = (10.0, 10.0)
    phase: float = float(np.pi / 20)
    n_points: int = 25

    def point_at(self, t: float) -> tuple[float, float]:
        """The (coordinate 1, coordinate 2) pair at loop parameter ``t``."""
        ang = 2.0 * np.pi * t + self.phase
        return (
            self.centre[0] + self.radius[0] * np.cos(ang),
            self.centre[1] + self.radius[1] * np.sin(ang),
        )

    @property
    def t_values(self) -> np.ndarray:
        """``n_points`` distinct parameter values, t = k / n_points."""
        return np.arange(self.n_points, dtype=float) / self.n_points

    def points(self) -> list[tuple[float, float]]:
        """The distinct (alpha, phi) points around the loop."""
        return [self.point_at(t) for t in self.t_values]

    def closing_point(self) -> tuple[float, float]:
        """The t=1 point. Numerically identical to ``points()[0]`` by construction."""
        return self.point_at(1.0)

    def encloses(self, alpha: float, phi: float) -> bool:
        """Whether (alpha, phi) lies strictly inside the ellipse."""
        da = (alpha - self.centre[0]) / self.radius[0]
        dp = (phi - self.centre[1]) / self.radius[1]
        return bool(da * da + dp * dp < 1.0)

    def bounding_box(self, margin: float = 0.0) -> tuple[float, float, float, float]:
        """(alpha_min, alpha_max, phi_min, phi_max), optionally padded by ``margin`` degrees."""
        return (
            self.centre[0] - self.radius[0] - margin,
            self.centre[0] + self.radius[0] + margin,
            self.centre[1] - self.radius[1] - margin,
            self.centre[1] + self.radius[1] + margin,
        )

    def with_n_points(self, n: int) -> "Loop":
        return Loop(self.name, self.centre, self.radius, self.phase, n)

    def geometries(self) -> list[str]:
        return [formaldimine_geom(a, p) for a, p in self.points()]


# --- The three benchmark loops (see docs/provenance.md) ----------------------------
# Centres are (alpha, phi) in degrees. The CI-enclosing loop is verbatim from the
# notebook; the two control loops follow arXiv:2304.06070 Sec. 6.3 and Fig. 1a, with
# the notebook's phi = 89.9 centre convention retained.

LOOP_CI = Loop("C_x", centre=(130.0, 89.9))       # encloses the S0/S1 CI  -> expect pi
LOOP_CONTROL_LOWER = Loop("C_1", centre=(110.0, 89.9))   # trivial -> expect 0
LOOP_CONTROL_UPPER = Loop("C_2", centre=(150.0, 89.9))   # trivial -> expect 0

BENCHMARK_LOOPS = {
    loop.name: loop for loop in (LOOP_CI, LOOP_CONTROL_LOWER, LOOP_CONTROL_UPPER)
}

# Reference CI location from the FCI gap map of arXiv:2304.06070 Sec. 6.3.
REFERENCE_CI_ALPHA_PHI = (132.0, 90.0)

# A scan *region*, not a traversal path: its bounding box covers all three benchmark loops,
# for the overview gap map that reproduces Fig. 1a of arXiv:2304.06070. Only its
# ``bounding_box`` and ``name`` are used; it is never traversed.
OVERVIEW_REGION = Loop("overview", centre=(130.0, 89.9), radius=(31.0, 11.0))

SCAN_REGIONS = {**BENCHMARK_LOOPS, "overview": OVERVIEW_REGION}
