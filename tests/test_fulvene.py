"""Fulvene geometry construction and loop placement.

Geometry tests are instant. One slow test confirms the CASSCF/overlap machinery accepts the
larger molecule; it is a smoke test of the code, not the production experiment (which is
prepared for the cluster -- see docs/followup.md).
"""

import numpy as np
import pytest

from berrycasscf.fulvene import (
    IDX_C1,
    IDX_C6,
    IDX_METHYLENE_H,
    REFERENCE_COORDS,
    REFERENCE_LABELS,
    REFERENCE_R_EXO,
    SEARCH_REGION,
    fulvene_coords,
    fulvene_geom,
    loops_around,
    DEFAULT_NCAS,
    DEFAULT_NELECAS,
)


def test_reference_geometry_is_planar_fulvene():
    assert REFERENCE_COORDS.shape == (12, 3)
    assert REFERENCE_LABELS.count("C") == 6 and REFERENCE_LABELS.count("H") == 6
    assert np.abs(REFERENCE_COORDS[:, 2]).max() < 1e-9          # planar
    assert REFERENCE_R_EXO == pytest.approx(1.3293, abs=1e-3)


def test_reference_parameters_reproduce_the_reference_geometry():
    assert np.allclose(fulvene_coords(REFERENCE_R_EXO, 0.0), REFERENCE_COORDS, atol=1e-12)


@pytest.mark.parametrize("r", [1.25, 1.3293, 1.45, 1.60])
def test_bond_length_coordinate_does_what_it_says(r):
    xyz = fulvene_coords(r, 0.0)
    assert np.linalg.norm(xyz[IDX_C6] - xyz[IDX_C1]) == pytest.approx(r, abs=1e-10)


@pytest.mark.parametrize("theta", [0.0, 30.0, 90.0, 180.0])
def test_torsion_preserves_every_bond_length(theta):
    """Rotating the methylene about the C1-C6 axis is rigid: no bond may change."""
    xyz = fulvene_coords(1.40, theta)
    for h in IDX_METHYLENE_H:
        assert np.linalg.norm(xyz[h] - xyz[IDX_C6]) == pytest.approx(
            np.linalg.norm(REFERENCE_COORDS[h] - REFERENCE_COORDS[IDX_C6]), abs=1e-10
        )
    # the C1-C6 axis itself is untouched by the torsion
    assert np.linalg.norm(xyz[IDX_C6] - xyz[IDX_C1]) == pytest.approx(1.40, abs=1e-10)


def test_ring_is_frozen_by_both_coordinates():
    xyz = fulvene_coords(1.55, 75.0)
    assert np.allclose(xyz[:5], REFERENCE_COORDS[:5], atol=1e-12)
    assert np.allclose(xyz[6:10], REFERENCE_COORDS[6:10], atol=1e-12)


def test_torsion_of_90_degrees_takes_the_methylene_out_of_plane():
    xyz = fulvene_coords(1.40, 90.0)
    assert np.abs(xyz[list(IDX_METHYLENE_H), 2]).max() > 0.8


def test_torsion_is_periodic():
    assert np.allclose(fulvene_coords(1.40, 0.0), fulvene_coords(1.40, 360.0), atol=1e-9)


def test_no_atoms_collide_anywhere_in_the_search_region():
    lo_r, hi_r, lo_t, hi_t = SEARCH_REGION.bounding_box()
    for r in np.linspace(lo_r, hi_r, 5):
        for t in np.linspace(lo_t, hi_t, 5):
            xyz = fulvene_coords(r, t)
            d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
            np.fill_diagonal(d, np.inf)
            assert d.min() > 0.7, f"atoms too close at r={r}, theta={t}"


def test_geometry_string_is_parseable_and_complete():
    text = fulvene_geom(1.42, 30.0)
    lines = text.strip().splitlines()
    assert len(lines) == 12
    assert [ln.split()[0] for ln in lines] == list(REFERENCE_LABELS)
    assert all(len(ln.split()) == 4 for ln in lines)


def test_loops_around_places_one_enclosing_and_two_control_loops():
    centre = (1.42, 30.0)
    loops = loops_around(centre, radius=(0.06, 15.0))
    assert set(loops) == {"F_x", "F_1", "F_2"}
    assert loops["F_x"].encloses(*centre)
    assert not loops["F_1"].encloses(*centre)
    assert not loops["F_2"].encloses(*centre)


def test_search_region_brackets_chemically_sensible_bond_lengths():
    lo_r, hi_r, lo_t, hi_t = SEARCH_REGION.bounding_box()
    assert lo_r < REFERENCE_R_EXO < hi_r
    assert (lo_t, hi_t) == (0.0, 90.0)


@pytest.mark.slow
def test_casscf_and_overlap_machinery_accept_fulvene():
    """Smoke test: the workflow runs on the larger molecule and the overlap stays exact."""
    from berrycasscf.casscf import build_mol, run_casscf, transfer_mo, orthonormality_error
    from berrycasscf.overlap import cas_overlap, brute_force_cas_overlap

    mol_a = build_mol(fulvene_geom(1.38, 20.0), "sto-3g")
    assert mol_a.nelectron == 42
    wfn_a = run_casscf(mol_a, DEFAULT_NCAS, DEFAULT_NELECAS)
    assert wfn_a.converged
    assert wfn_a.ncore == 18
    assert wfn_a.ci.shape == (20, 20)
    assert cas_overlap(wfn_a, wfn_a) == pytest.approx(1.0, abs=1e-10)

    mol_b = build_mol(fulvene_geom(1.40, 25.0), "sto-3g")
    guess = transfer_mo(mol_a, wfn_a.mo_coeff, mol_b)
    assert orthonormality_error(mol_b, guess) < 1e-10
    wfn_b = run_casscf(mol_b, DEFAULT_NCAS, DEFAULT_NELECAS,
                       mo_guess=guess, ci0=np.asarray(wfn_a.ci))
    assert wfn_b.converged
    # the factorized overlap must still match brute force with 18 core orbitals
    assert cas_overlap(wfn_a, wfn_b) == pytest.approx(
        brute_force_cas_overlap(wfn_a, wfn_b), abs=1e-10
    )
