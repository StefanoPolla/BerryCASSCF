"""Ethylene (C2H4): the textbook twisted-pyramidalized S1/S0 conical intersection.

Why ethylene
------------
Its S1/S0 intersection is the canonical example in nonadiabatic photochemistry and, unlike
fulvene's, it is reached by two *rigid* coordinates that leave every bond length untouched:

``tau``    torsion of one CH2 group about the C=C axis, degrees (0 = planar, 90 = twisted)
``phi``    umbrella pyramidalization of that same CH2, degrees (0 = planar sp2)

The minimal pi active space is CAS(2,2). What makes ethylene interesting for an
active-space study is that its S1 near the twisted geometry has substantial zwitterionic
character, which a two-orbital space describes badly -- so the two workflows are expected to
have genuinely different active-space requirements here.

STO-3G has 14 AOs and 16 electrons, of which 2 orbitals are C 1s core. **CAS(12,12) is
therefore the full valence space**, i.e. FCI in all but the two core orbitals. That gives a
near-exact in-basis reference without the cost of true FCI, which is what makes ethylene a good
system for asking "how large an active space is actually needed?".

Geometry construction
---------------------
C1 at the origin, C2 along +z. Group 1 (on C1) is fixed. Group 2 (on C2) is rotated about the
C=C axis by ``tau`` and then tilted as a unit by ``phi`` about its own H-H direction. Both
operations are rigid rotations, so every bond length and both HCH angles are preserved exactly.
"""

from __future__ import annotations

import numpy as np

from .geometry import Loop

# Equilibrium internal coordinates, RHF/6-31G* optimized (Angstrom / degrees).
R_CC = 1.3160
R_CH = 1.0759
ANGLE_HCC = 121.84

LABELS = ("C", "C", "H", "H", "H", "H")


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _rotate(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, degrees: float) -> np.ndarray:
    """Rodrigues rotation of ``points`` about the line through ``origin`` along ``axis``."""
    k = _unit(axis)
    a = np.deg2rad(degrees)
    v = points - origin
    return (
        origin
        + v * np.cos(a)
        + np.cross(k, v) * np.sin(a)
        + k * (v @ k)[:, None] * (1.0 - np.cos(a))
    )


def ethylene_coords(tau: float, phi: float) -> np.ndarray:
    """Cartesian coordinates (Angstrom) for the two-parameter ethylene family."""
    c1 = np.array([0.0, 0.0, 0.0])
    c2 = np.array([0.0, 0.0, R_CC])
    a = np.deg2rad(ANGLE_HCC)

    # Group 1 on C1: fixed, in the xz-plane, pointing away from C2 (-z component).
    h1a = c1 + R_CH * np.array([np.sin(a), 0.0, np.cos(a)])
    h1b = c1 + R_CH * np.array([-np.sin(a), 0.0, np.cos(a)])

    # Group 2 on C2: planar reference, mirror image about the C=C midpoint.
    h2a = c2 + R_CH * np.array([np.sin(a), 0.0, -np.cos(a)])
    h2b = c2 + R_CH * np.array([-np.sin(a), 0.0, -np.cos(a)])
    group2 = np.array([h2a, h2b])

    # 1. torsion about the C=C axis
    if tau != 0.0:
        group2 = _rotate(group2, c2, np.array([0.0, 0.0, 1.0]), tau)

    # 2. umbrella pyramidalization: tilt the group as a unit about its own H-H direction,
    #    which after the torsion lies along (cos tau, sin tau, 0).
    if phi != 0.0:
        hh_axis = group2[0] - group2[1]
        group2 = _rotate(group2, c2, hh_axis, phi)

    return np.array([c1, c2, h1a, h1b, group2[0], group2[1]])


def ethylene_geom(tau: float, phi: float) -> str:
    """PySCF geometry string; signature matches the other systems' ``geom_fn``."""
    xyz = ethylene_coords(tau, phi)
    return "\n".join(
        f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(LABELS, xyz)
    )


# Coarse search region: the full torsion and umbrella ranges, tau and phi both in [0, 180].
# Only the bounding box is used; this is never traversed. The first coarse scan was run over
# phi in [0, 90] and found nothing -- the intersection turned out to sit near phi = 110, so
# the umbrella coordinate has to be followed past the planar-inverted point.
SEARCH_REGION = Loop("ethylene_search", centre=(90.0, 90.0), radius=(90.0, 90.0))

# Focused region around the twisted-pyramidalized intersection, located by a coarse scan at
# 6-31G*: tau in [70, 110], phi in [90, 130].
CI_REGION = Loop("ethylene_ci", centre=(90.0, 110.0), radius=(20.0, 20.0))

# STO-3G has no polarization functions and never closes the S1/S0 gap along these coordinates
# (minimum ~10 mHa); 6-31G* does. The basis, not the active space, is the binding constraint
# there, so all ethylene work uses 6-31G*.
DEFAULT_BASIS = "6-31g*"

DEFAULT_NCAS = 2
DEFAULT_NELECAS = 2
# Full valence space in STO-3G (only the two C 1s orbitals remain core).
FULL_VALENCE = (12, 12)
