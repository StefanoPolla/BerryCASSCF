"""1,3-butadiene (C4H6): a system where the active space is expected to matter more.

Why butadiene
-------------
Its low-lying excited states are the standard hard case for active-space truncation. The
1(1)Bu state is a singly-excited pi->pi* state, but the 2(1)Ag state carries a large
*doubly*-excited component, and the ordering and relative position of the two is notoriously
sensitive to both the active space and dynamic correlation. A conical intersection whose upper
state has that character is exactly the situation in which a small CAS should fail -- in contrast
to ethylene, where the minimal CAS(2,2) turned out to be accidentally excellent.

Unlike formaldimine (13 orbitals) and ethylene (12 valence orbitals in 6-31G*), butadiene's full
valence space is CAS(22,22), far beyond reach. There is therefore **no exact in-basis reference
here**: the largest affordable active space is the best available, and if the answer is still
moving at the top of the ladder that is itself the finding.

Coordinates
-----------
Four rigid distortions of a fixed RHF/6-31G* optimized planar s-trans reference. Each is a pure
rotation or a rigid translation of a fragment, so bond lengths and angles inside each fragment
are preserved exactly.

``tw``    twist of the C1 methylene about the C1=C2 axis (deg)
``pyr``   umbrella pyramidalization of that same methylene (deg)
``tc``    torsion about the central C2-C3 bond, i.e. s-trans (0) to s-cis (180) (deg)
``bend``  in-plane C1-C2-C3 angle change (deg)

The search phase (``examples/search_butadiene_ci.py``) scans pairs of these to locate an
intersection before any ladder is run; nothing here assumes where it is.

Symmetry
--------
Ethylene's gap maps could be validated against ``tau -> 180 - tau``, because its two methylene
groups are equivalent. **Butadiene has no such symmetry**: the two hydrogens on C1 are
inequivalent (one cis, one trans to the C3=C4 unit), so ``tw`` and ``180 - tw`` give genuinely
different geometries -- their interatomic-distance spectra differ by ~3e-03 Angstrom, which is
small but real, and the corresponding energies differ accordingly.

The exact symmetry is reflection through the molecular plane,

    (tw, pyr) -> (-tw, -pyr)     at tc = 0

which is isometric to numerical precision and is what
:func:`berrycasscf.butadiene.mirror_partner` returns. Any path-independence check on this system
must use that pair, not the ethylene one. (With ``tc != 0`` the molecule is non-planar and even
this symmetry is gone.)
"""

from __future__ import annotations

import numpy as np

from .geometry import Loop

# RHF/6-31G* optimized, planar s-trans (C2h). Angstrom.
# Order: C1 C2 C3 C4, then H,H on C1, H on C2, H on C3, H,H on C4.
LABELS = ("C", "C", "C", "C", "H", "H", "H", "H", "H", "H")
REFERENCE_COORDS = np.array([
    [-1.829923,  0.169004, 0.0],   # 0  C1  terminal, bears the twisted methylene
    [-0.626842, -0.381170, 0.0],   # 1  C2
    [ 0.626909,  0.381255, 0.0],   # 2  C3
    [ 1.830038, -0.168842, 0.0],   # 3  C4  terminal
    [-2.723104, -0.428649, 0.0],   # 4  H on C1
    [-1.964532,  1.237086, 0.0],   # 5  H on C1
    [-0.532999, -1.455057, 0.0],   # 6  H on C2
    [ 0.533061,  1.455137, 0.0],   # 7  H on C3
    [ 1.963789, -1.237016, 0.0],   # 8  H on C4
    [ 2.723603,  0.428251, 0.0],   # 9  H on C4
])

IDX_C1, IDX_C2, IDX_C3, IDX_C4 = 0, 1, 2, 3
IDX_METHYLENE_H = (4, 5)                 # the two hydrogens twisted/pyramidalized
# Everything on the C2 side of the central bond, excluding C2 itself: rotated by ``tc``.
IDX_LEFT_FRAGMENT = (0, 4, 5, 6)

DEFAULT_BASIS = "6-31g*"
# Full pi manifold: 4 electrons in 4 pi orbitals.
DEFAULT_NCAS = 4
DEFAULT_NELECAS = 4


def _rotate(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, degrees: float) -> np.ndarray:
    """Rodrigues rotation of ``points`` about the line through ``origin`` along ``axis``."""
    k = axis / np.linalg.norm(axis)
    a = np.deg2rad(degrees)
    v = points - origin
    return (origin + v * np.cos(a) + np.cross(k, v) * np.sin(a)
            + k * (v @ k)[:, None] * (1.0 - np.cos(a)))


def butadiene_coords(tw: float = 0.0, pyr: float = 0.0, tc: float = 0.0,
                     bend: float = 0.0) -> np.ndarray:
    """Cartesian coordinates (Angstrom) for the four-parameter butadiene family."""
    xyz = REFERENCE_COORDS.copy()

    # 1. in-plane C1-C2-C3 bend: rotate the C1 methylene group about the axis through C2
    #    perpendicular to the molecular plane.
    if bend != 0.0:
        grp = [IDX_C1, *IDX_METHYLENE_H]
        xyz[grp] = _rotate(xyz[grp], xyz[IDX_C2], np.array([0.0, 0.0, 1.0]), bend)

    # 2. methylene twist about the (possibly bent) C1=C2 axis
    if tw != 0.0:
        axis = xyz[IDX_C2] - xyz[IDX_C1]
        h = list(IDX_METHYLENE_H)
        xyz[h] = _rotate(xyz[h], xyz[IDX_C1], axis, tw)

    # 3. umbrella pyramidalization of that methylene, about its own H-H direction
    if pyr != 0.0:
        h = list(IDX_METHYLENE_H)
        xyz[h] = _rotate(xyz[h], xyz[IDX_C1], xyz[h[0]] - xyz[h[1]], pyr)

    # 4. central torsion: rotate the whole left fragment about the C2-C3 axis
    if tc != 0.0:
        axis = xyz[IDX_C3] - xyz[IDX_C2]
        grp = list(IDX_LEFT_FRAGMENT)
        xyz[grp] = _rotate(xyz[grp], xyz[IDX_C2], axis, tc)

    return xyz


def butadiene_geom(tw: float, pyr: float) -> str:
    """PySCF geometry string in the (twist, pyramidalization) plane.

    Signature matches the other systems' ``geom_fn`` so the same drivers apply. Planes using
    other coordinate pairs are built with :func:`plane_geom_fn`.
    """
    return _to_string(butadiene_coords(tw=tw, pyr=pyr))


def plane_geom_fn(coord_x: str, coord_y: str, **fixed):
    """Build a ``geom_fn(x, y)`` for any pair of the four coordinates.

    Used by the search phase to explore several planes without duplicating code::

        fn = plane_geom_fn("tw", "tc", pyr=0.0, bend=0.0)
    """
    names = {"tw", "pyr", "tc", "bend"}
    if not {coord_x, coord_y} <= names or coord_x == coord_y:
        raise ValueError(f"coordinates must be two distinct names from {sorted(names)}")

    def geom_fn(x: float, y: float) -> str:
        kwargs = dict(fixed)
        kwargs[coord_x] = float(x)
        kwargs[coord_y] = float(y)
        return _to_string(butadiene_coords(**kwargs))

    geom_fn.__name__ = f"butadiene_{coord_x}_{coord_y}"
    return geom_fn


def mirror_partner(tw: float, pyr: float) -> tuple[float, float]:
    """The geometry related to ``(tw, pyr)`` by reflection through the molecular plane.

    Exact only at ``tc = 0``. Energies at a point and its partner must agree; any difference is
    a solver inconsistency, not physics.
    """
    return (-float(tw), -float(pyr))


def _to_string(xyz: np.ndarray) -> str:
    return "\n".join(f"{s} {x:.10f} {y:.10f} {z:.10f}"
                     for s, (x, y, z) in zip(LABELS, xyz))


# Candidate planes for the search phase. Only bounding boxes are used; none is traversed.
SEARCH_PLANES = {
    # terminal twist + pyramidalization, the motif that hosts ethylene's intersection
    "tw_pyr": (("tw", "pyr"), {"tc": 0.0, "bend": 0.0},
               Loop("tw_pyr", centre=(90.0, 90.0), radius=(90.0, 90.0))),
    # terminal twist + central torsion (s-trans to s-cis)
    "tw_tc": (("tw", "tc"), {"pyr": 0.0, "bend": 0.0},
              Loop("tw_tc", centre=(90.0, 90.0), radius=(90.0, 90.0))),
    # terminal twist + skeletal bend
    "tw_bend": (("tw", "bend"), {"pyr": 0.0, "tc": 0.0},
                Loop("tw_bend", centre=(90.0, 0.0), radius=(90.0, 40.0))),
}
