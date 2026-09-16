"""Butadiene geometry construction: four rigid coordinates on a ten-atom molecule."""

import numpy as np
import pytest

from berrycasscf.butadiene import (
    IDX_C1,
    IDX_C2,
    IDX_C3,
    IDX_C4,
    IDX_LEFT_FRAGMENT,
    IDX_METHYLENE_H,
    LABELS,
    REFERENCE_COORDS,
    SEARCH_PLANES,
    butadiene_coords,
    butadiene_geom,
    plane_geom_fn,
)

BONDS = [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5), (1, 6), (2, 7), (3, 8), (3, 9)]


def _lengths(xyz):
    return np.array([np.linalg.norm(xyz[i] - xyz[j]) for i, j in BONDS])


def test_reference_is_planar_s_trans_with_bond_alternation():
    xyz = REFERENCE_COORDS
    assert xyz.shape == (10, 3)
    assert LABELS.count("C") == 4 and LABELS.count("H") == 6
    assert np.abs(xyz[:, 2]).max() < 1e-9                       # planar
    d = lambda i, j: np.linalg.norm(xyz[i] - xyz[j])
    assert d(IDX_C1, IDX_C2) == pytest.approx(1.3229, abs=2e-3)  # double
    assert d(IDX_C3, IDX_C4) == pytest.approx(1.3229, abs=2e-3)  # double
    assert d(IDX_C2, IDX_C3) == pytest.approx(1.4674, abs=2e-3)  # single, longer
    assert d(IDX_C2, IDX_C3) > d(IDX_C1, IDX_C2) + 0.1


def test_defaults_reproduce_the_reference():
    assert np.allclose(butadiene_coords(), REFERENCE_COORDS, atol=1e-12)


@pytest.mark.parametrize("kw", [
    {"tw": 90.0}, {"pyr": 60.0}, {"tc": 180.0}, {"bend": 20.0},
    {"tw": 90.0, "pyr": 110.0}, {"tw": 45.0, "pyr": 30.0, "tc": 90.0, "bend": 15.0},
])
def test_every_coordinate_is_rigid(kw):
    """No distortion may change any bond length: all four are pure rotations."""
    assert np.allclose(_lengths(butadiene_coords(**kw)), _lengths(REFERENCE_COORDS), atol=1e-10)


@pytest.mark.parametrize("kw", [
    {"tw": 90.0}, {"pyr": 90.0}, {"tc": 90.0}, {"bend": 30.0},
    {"tw": 90.0, "pyr": 110.0, "tc": 90.0, "bend": 20.0},
])
def test_no_atomic_collisions(kw):
    xyz = butadiene_coords(**kw)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    assert d.min() > 0.9


def test_twist_moves_only_the_terminal_methylene_hydrogens():
    ref, tw = butadiene_coords(), butadiene_coords(tw=90.0)
    moved = {i for i in range(10) if not np.allclose(ref[i], tw[i], atol=1e-9)}
    assert moved == set(IDX_METHYLENE_H)
    assert np.abs(tw[list(IDX_METHYLENE_H), 2]).max() > 0.8      # now out of plane


def test_pyramidalization_moves_both_hydrogens_together():
    flat, pyr = butadiene_coords(tw=90.0), butadiene_coords(tw=90.0, pyr=60.0)
    delta = pyr[list(IDX_METHYLENE_H)] - flat[list(IDX_METHYLENE_H)]
    assert np.abs(delta).max() > 0.1
    # umbrella, not a shear: both hydrogens move the same way along the C1->C2 axis
    axis = flat[IDX_C2] - flat[IDX_C1]
    axis /= np.linalg.norm(axis)
    proj = delta @ axis
    assert np.sign(proj[0]) == np.sign(proj[1])


def test_central_torsion_moves_only_the_left_fragment():
    ref, rot = butadiene_coords(), butadiene_coords(tc=90.0)
    moved = {i for i in range(10) if not np.allclose(ref[i], rot[i], atol=1e-9)}
    assert moved == set(IDX_LEFT_FRAGMENT)
    assert np.allclose(rot[[IDX_C2, IDX_C3, IDX_C4]], ref[[IDX_C2, IDX_C3, IDX_C4]], atol=1e-12)


def _dihedral(x, a, b, c, d):
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


@pytest.mark.parametrize("tc,expected", [(0.0, 180.0), (45.0, 135.0), (90.0, 90.0),
                                         (180.0, 0.0)])
def test_central_torsion_sets_the_ccc_c_dihedral(tc, expected):
    """tc = 0 is s-trans (dihedral 180) and tc = 180 is s-cis (dihedral 0), linearly between."""
    x = butadiene_coords(tc=tc)
    assert abs(_dihedral(x, IDX_C1, IDX_C2, IDX_C3, IDX_C4)) == pytest.approx(expected, abs=1e-6)


def test_s_cis_brings_the_terminal_carbons_closer():
    d = lambda x: np.linalg.norm(x[IDX_C1] - x[IDX_C4])
    assert d(butadiene_coords(tc=180.0)) < d(butadiene_coords(tc=0.0)) - 0.5


def test_bend_changes_the_ccc_angle_in_plane():
    ref, bent = butadiene_coords(), butadiene_coords(bend=20.0)
    def ccc(x):
        u = (x[IDX_C1] - x[IDX_C2]) / np.linalg.norm(x[IDX_C1] - x[IDX_C2])
        v = (x[IDX_C3] - x[IDX_C2]) / np.linalg.norm(x[IDX_C3] - x[IDX_C2])
        return np.degrees(np.arccos(np.clip(u @ v, -1, 1)))
    assert abs(ccc(bent) - ccc(ref)) == pytest.approx(20.0, abs=1e-6)
    assert np.abs(bent[:, 2]).max() < 1e-9                       # stays planar


def test_geometry_string_is_well_formed():
    lines = butadiene_geom(90.0, 60.0).strip().splitlines()
    assert len(lines) == 10
    assert [ln.split()[0] for ln in lines] == list(LABELS)
    assert all(len(ln.split()) == 4 for ln in lines)


def test_plane_geom_fn_selects_the_right_two_coordinates():
    fn = plane_geom_fn("tw", "tc", pyr=0.0, bend=0.0)
    got = fn(90.0, 45.0)
    want = butadiene_geom.__globals__["_to_string"](
        butadiene_coords(tw=90.0, tc=45.0, pyr=0.0, bend=0.0)
    )
    assert got == want


def test_plane_geom_fn_rejects_bad_coordinate_names():
    with pytest.raises(ValueError, match="two distinct names"):
        plane_geom_fn("tw", "tw")
    with pytest.raises(ValueError, match="two distinct names"):
        plane_geom_fn("tw", "stretch")


def test_search_planes_are_well_formed():
    for name, ((cx, cy), fixed, region) in SEARCH_PLANES.items():
        assert cx != cy
        assert set(fixed) | {cx, cy} == {"tw", "pyr", "tc", "bend"}
        lo_x, hi_x, lo_y, hi_y = region.bounding_box()
        assert hi_x > lo_x and hi_y > lo_y


def _spectrum(xyz):
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    return np.sort(d.ravel())


@pytest.mark.parametrize("tw,pyr", [(90.0, 110.0), (75.0, 60.0), (45.0, 135.0), (30.0, 0.0)])
def test_reflection_through_the_molecular_plane_is_an_exact_symmetry(tw, pyr):
    """(tw, pyr) -> (-tw, -pyr) is isometric at tc = 0. This is the check butadiene supports."""
    from berrycasscf.butadiene import mirror_partner

    a = butadiene_coords(tw=tw, pyr=pyr)
    b = butadiene_coords(tw=mirror_partner(tw, pyr)[0], pyr=mirror_partner(tw, pyr)[1])
    assert np.allclose(_spectrum(a), _spectrum(b), atol=1e-10)


@pytest.mark.parametrize("tw,pyr", [(75.0, 0.0), (75.0, 105.0), (45.0, 90.0)])
def test_ethylene_style_mirror_is_NOT_a_symmetry_here(tw, pyr):
    """Guards against reusing ethylene's check.

    Butadiene's two methylene hydrogens are inequivalent (one cis, one trans to C3=C4), so
    tw and 180 - tw are genuinely different geometries. Applying ethylene's validation to this
    system would silently compare unrelated points.
    """
    a = butadiene_coords(tw=tw, pyr=pyr)
    b = butadiene_coords(tw=180.0 - tw, pyr=pyr)
    assert not np.allclose(_spectrum(a), _spectrum(b), atol=1e-6)


def test_reflection_symmetry_is_lost_once_the_molecule_is_non_planar():
    a = butadiene_coords(tw=75.0, pyr=60.0, tc=90.0)
    b = butadiene_coords(tw=-75.0, pyr=-60.0, tc=90.0)
    assert not np.allclose(_spectrum(a), _spectrum(b), atol=1e-6)
