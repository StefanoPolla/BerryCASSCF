"""Fulvene (C6H6): the follow-up system.

Why fulvene
-----------
It is a real step up from formaldimine without being intractable:

* the full pi manifold is **6 electrons in 6 orbitals**, so CAS(6,6) is the smallest
  chemically honest active space -- three times the orbitals of the formaldimine benchmark,
  and 400 determinants instead of 4;
* S1 is a genuine pi->pi* state, so the state-averaged comparator has real work to do;
* the S1/S0 photochemistry is well studied, and the two coordinates that drive it -- the
  exocyclic C=C bond length and the methylene torsion -- are the standard ones in that
  literature, which makes the scan plane defensible rather than arbitrary;
* it is a closed-shell hydrocarbon with no near-degenerate sigma manifold, so the active
  space stays well defined around a loop, which matters for continuation.

Coordinates
-----------
Two parameters distort a fixed reference geometry:

``r``      length of the exocyclic C1=C6 bond, Angstrom (reference 1.3293)
``theta``  rotation of the methylene (C6) hydrogens about the C1-C6 axis, degrees
           (reference 0, i.e. planar)

``r`` translates C6 and its two hydrogens rigidly along the C1->C6 axis; ``theta`` rotates
only the two methylene hydrogens about that axis. Every other internal coordinate is frozen,
exactly as ``alpha``/``phi`` leave formaldimine's other coordinates frozen.

Status
------
The geometry builder and the loop machinery are tested (``tests/test_fulvene.py``), but **no
production Berry-phase or scan result exists**: per the project brief, anything heavier than
the formaldimine benchmark is prepared for the cluster rather than run locally. See
``docs/followup.md`` for the two-stage plan and ``slurm/`` for the batch templates.
"""

from __future__ import annotations

import numpy as np

from .geometry import Loop

# Reference S0 geometry, RHF/6-31G optimized, planar (C2v), Angstrom.
# Atom order: C1 (ring, bears the exocyclic carbon), C2..C5 (ring), C6 (exocyclic),
# H on C2..C5, then the two methylene hydrogens on C6.
REFERENCE_LABELS = ("C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H")
REFERENCE_COORDS = np.array(
    [
        [0.000007,  1.220542, 0.0],   # C1  ring, attached to C6
        [-1.176007, 0.323692, 0.0],   # C2
        [-0.739563, -0.942701, 0.0],  # C3
        [0.739554, -0.942691, 0.0],   # C4
        [1.176007,  0.323698, 0.0],   # C5
        [0.000005,  2.549844, 0.0],   # C6  exocyclic
        [-2.190885, 0.661821, 0.0],   # H on C2
        [-1.342744, -1.826540, 0.0],  # H on C3
        [1.342740, -1.826527, 0.0],   # H on C4
        [2.190890,  0.661813, 0.0],   # H on C5
        [0.912993,  3.114446, 0.0],   # H on C6 (methylene)
        [-0.912997, 3.114424, 0.0],   # H on C6 (methylene)
    ]
)

IDX_C1, IDX_C6 = 0, 5
IDX_METHYLENE_H = (10, 11)

REFERENCE_R_EXO = float(np.linalg.norm(REFERENCE_COORDS[IDX_C6] - REFERENCE_COORDS[IDX_C1]))
REFERENCE_THETA = 0.0


def _rotate_about_axis(points: np.ndarray, origin: np.ndarray, axis: np.ndarray,
                       degrees: float) -> np.ndarray:
    """Rodrigues rotation of ``points`` about the line through ``origin`` along ``axis``."""
    k = axis / np.linalg.norm(axis)
    ang = np.deg2rad(degrees)
    v = points - origin
    return (
        origin
        + v * np.cos(ang)
        + np.cross(k, v) * np.sin(ang)
        + k * (v @ k)[:, None] * (1.0 - np.cos(ang))
    )


def fulvene_coords(r: float, theta: float) -> np.ndarray:
    """Cartesian coordinates (Angstrom) for the two-parameter fulvene family."""
    xyz = REFERENCE_COORDS.copy()
    c1, c6 = xyz[IDX_C1], xyz[IDX_C6]
    axis = c6 - c1
    axis_hat = axis / np.linalg.norm(axis)

    # 1. set the exocyclic bond length: move C6 and its hydrogens rigidly along the axis
    shift = (r - REFERENCE_R_EXO) * axis_hat
    movers = [IDX_C6, *IDX_METHYLENE_H]
    xyz[movers] += shift

    # 2. twist the methylene hydrogens about the (unchanged) C1-C6 axis
    if theta != 0.0:
        h_idx = list(IDX_METHYLENE_H)
        xyz[h_idx] = _rotate_about_axis(xyz[h_idx], xyz[IDX_C6], axis_hat, theta)
    return xyz


def fulvene_geom(r: float, theta: float) -> str:
    """PySCF geometry string for the two-parameter fulvene family.

    Signature matches :func:`berrycasscf.geometry.formaldimine_geom`, so the same
    continuation and scan drivers work unchanged via their ``geom_fn`` argument.
    """
    xyz = fulvene_coords(r, theta)
    return "\n".join(
        f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(REFERENCE_LABELS, xyz)
    )


# --- scan region for stage 1 -------------------------------------------------------
# A generous region of the (r, theta) plane in which to look for the S1/S0 seam. The CI
# position is NOT assumed: stage 1 locates it, stage 2 builds loops around whatever is found.
# Bounding box: r in [1.20, 1.60] Angstrom, theta in [0, 90] degrees.
SEARCH_REGION = Loop("fulvene_search", centre=(1.40, 45.0), radius=(0.20, 45.0))

# Default active space: the full pi manifold.
DEFAULT_NCAS = 6
DEFAULT_NELECAS = 6


def loops_around(centre: tuple[float, float],
                 radius: tuple[float, float] = (0.06, 15.0),
                 n_points: int = 25) -> dict[str, Loop]:
    """Build a CI-enclosing loop and two control loops around a located intersection.

    Mirrors the formaldimine benchmark's layout: one loop centred on the CI, and two displaced
    by twice the radius along the first coordinate so that they enclose nothing.
    """
    r0, t0 = centre
    dr = 2.0 * radius[0]
    return {
        "F_x": Loop("F_x", (r0, t0), radius, n_points=n_points),
        "F_1": Loop("F_1", (r0 - dr, t0), radius, n_points=n_points),
        "F_2": Loop("F_2", (r0 + dr, t0), radius, n_points=n_points),
    }
