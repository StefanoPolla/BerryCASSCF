"""Ethylene geometry construction and the symmetry that validates its gap scans."""

import numpy as np
import pytest

from berrycasscf.ethylene import (
    ANGLE_HCC,
    CI_REGION,
    LABELS,
    R_CC,
    R_CH,
    SEARCH_REGION,
    ethylene_coords,
    ethylene_geom,
)


def _distance_spectrum(xyz):
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    return np.sort(d.ravel())


def test_planar_reference_geometry():
    xyz = ethylene_coords(0.0, 0.0)
    assert xyz.shape == (6, 3)
    assert LABELS == ("C", "C", "H", "H", "H", "H")
    assert np.abs(xyz[:, 1]).max() < 1e-12            # planar, lies in the xz-plane
    assert np.linalg.norm(xyz[1] - xyz[0]) == pytest.approx(R_CC, abs=1e-12)


@pytest.mark.parametrize("tau,phi", [(0, 0), (45, 0), (90, 0), (90, 60), (137, 113), (180, 90)])
def test_both_coordinates_are_rigid(tau, phi):
    """Neither coordinate may change a bond length or either HCH angle."""
    xyz = ethylene_coords(float(tau), float(phi))
    assert np.linalg.norm(xyz[1] - xyz[0]) == pytest.approx(R_CC, abs=1e-10)
    for c, h in [(0, 2), (0, 3), (1, 4), (1, 5)]:
        assert np.linalg.norm(xyz[h] - xyz[c]) == pytest.approx(R_CH, abs=1e-10)
    for c, ha, hb in [(0, 2, 3), (1, 4, 5)]:
        u = (xyz[ha] - xyz[c]) / R_CH
        v = (xyz[hb] - xyz[c]) / R_CH
        expected = 2 * (180.0 - ANGLE_HCC)
        assert np.degrees(np.arccos(np.clip(u @ v, -1, 1))) == pytest.approx(expected, abs=1e-8)


@pytest.mark.parametrize("phi", [0.0, 60.0, 110.0, 130.0])
@pytest.mark.parametrize("dtau", [5.0, 20.0, 45.0])
def test_mirror_symmetry_in_torsion(phi, dtau):
    """tau and 180-tau are exact mirror images.

    This is the symmetry the scans are checked against: any asymmetry in a computed gap map
    is active-space drift along the scan path, never physics.
    """
    a = ethylene_coords(90.0 - dtau, phi)
    b = ethylene_coords(90.0 + dtau, phi)
    assert np.allclose(_distance_spectrum(a), _distance_spectrum(b), atol=1e-10)


def test_torsion_twists_the_second_group_only():
    ref = ethylene_coords(0.0, 0.0)
    tw = ethylene_coords(90.0, 0.0)
    assert np.allclose(ref[:4], tw[:4], atol=1e-12)      # C1, C2 and group 1 unmoved
    assert np.abs(tw[4:, 1]).max() > 0.8                 # group 2 rotated out of the xz-plane


def test_pyramidalization_moves_both_hydrogens_the_same_way():
    """The umbrella mode: both group-2 hydrogens tilt together, not in opposite senses."""
    flat = ethylene_coords(90.0, 0.0)
    pyr = ethylene_coords(90.0, 60.0)
    dz = pyr[4:, 2] - flat[4:, 2]
    assert np.sign(dz[0]) == np.sign(dz[1])
    assert np.abs(dz).min() > 0.1


def test_torsion_is_periodic():
    assert np.allclose(ethylene_coords(30.0, 40.0), ethylene_coords(390.0, 40.0), atol=1e-9)


def test_no_atomic_collisions_across_the_search_region():
    lo_t, hi_t, lo_p, hi_p = SEARCH_REGION.bounding_box()
    for tau in np.linspace(lo_t, hi_t, 7):
        for phi in np.linspace(lo_p, hi_p, 7):
            xyz = ethylene_coords(tau, phi)
            d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
            np.fill_diagonal(d, np.inf)
            assert d.min() > 0.8, f"atoms too close at tau={tau}, phi={phi}"


def test_ci_region_sits_inside_the_coarse_search_region():
    lo_t, hi_t, lo_p, hi_p = CI_REGION.bounding_box()
    Lo_t, Hi_t, Lo_p, Hi_p = SEARCH_REGION.bounding_box()
    assert Lo_t <= lo_t and hi_t <= Hi_t
    assert Lo_p <= lo_p and hi_p <= Hi_p
    assert CI_REGION.encloses(90.0, 110.0)


def test_geometry_string_is_well_formed():
    lines = ethylene_geom(90.0, 110.0).strip().splitlines()
    assert len(lines) == 6
    assert [ln.split()[0] for ln in lines] == list(LABELS)
    assert all(len(ln.split()) == 4 for ln in lines)


@pytest.mark.slow
def test_mirror_symmetry_is_reproduced_by_a_cold_started_casscf():
    """The computational counterpart of the geometric symmetry test above."""
    from pyscf import mcscf

    from berrycasscf.casscf import apply_singlet_constraint, build_mol, run_rhf

    energies = []
    for tau in (85.0, 95.0):
        mol = build_mol(ethylene_geom(tau, 110.0), "sto-3g")
        mc = mcscf.CASSCF(run_rhf(mol), 2, 2)
        apply_singlet_constraint(mc)
        mc = mc.state_average_([0.5, 0.5])
        mc.verbose = 0
        mc.kernel()
        energies.append(np.asarray(mc.e_states[:2]))
    assert np.allclose(energies[0], energies[1], atol=1e-9)


@pytest.mark.slow
@pytest.mark.parametrize("strategy", ["cold", "anchor"])
def test_gap_scan_is_path_independent(strategy):
    """Regression test for active-space drift in the scan.

    A warm-started sweep made the gap map depend on the route taken to each geometry, so
    mirror-image points disagreed by up to 9.8 mHa. The path-independent strategies must
    reproduce the geometric mirror symmetry to numerical precision.
    """
    from berrycasscf.config import CasConfig, ScanConfig
    from berrycasscf.geometry import Loop
    from berrycasscf.scan import scan_gap

    # A small odd grid symmetric about tau = 90, so mirror pairs are actually sampled.
    region = Loop("test", centre=(90.0, 110.0), radius=(10.0, 8.0))
    res = scan_gap(
        region,
        cas=CasConfig(basis="sto-3g", ncas=2, nelecas=2),
        scan=ScanConfig(n_alpha=5, n_phi=3, margin=0.0, strategy=strategy),
        geom_fn=ethylene_geom,
    )
    gap = res.gap
    assert np.isfinite(gap).all()
    assert np.allclose(gap, gap[::-1, :], atol=1e-8), (
        f"strategy={strategy} produced a path-dependent map: "
        f"max asymmetry {np.abs(gap - gap[::-1, :]).max():.2e} Ha"
    )
