"""The nonorthogonal CAS-CI overlap: the technical core of the package.

These run real CASSCF calculations but on a minimal basis and active space, so the whole
file takes a few seconds.
"""

import numpy as np
import pytest
import scipy.linalg as sla
from pyscf import fci

from berrycasscf.casscf import (
    build_mol,
    run_casscf,
    transfer_mo,
    orthonormality_error,
    core_count,
)
from berrycasscf.geometry import formaldimine_geom
from berrycasscf.overlap import cas_overlap, brute_force_cas_overlap, cross_mo_overlap

BASIS = "sto-3g"
GEOM_A = (130.0, 89.9)
GEOM_B = (133.0, 95.0)


@pytest.fixture(scope="module")
def wfn_a():
    return run_casscf(build_mol(formaldimine_geom(*GEOM_A), BASIS), 2, 2)


@pytest.fixture(scope="module")
def wfn_b():
    return run_casscf(build_mol(formaldimine_geom(*GEOM_B), BASIS), 2, 2)


def test_self_overlap_is_one(wfn_a):
    """Required check 1 from the brief: a wavefunction overlaps itself to 1."""
    assert cas_overlap(wfn_a, wfn_a) == pytest.approx(1.0, abs=1e-10)


def test_factorized_matches_brute_force_across_geometries(wfn_a, wfn_b):
    """The Schur-complement core elimination must be exact, not approximate."""
    fast = cas_overlap(wfn_a, wfn_b)
    slow = brute_force_cas_overlap(wfn_a, wfn_b)
    assert fast == pytest.approx(slow, abs=1e-12)
    # and the two geometries really are different, so this exercised nonorthogonality
    assert abs(fast) < 0.999


def test_overlap_is_symmetric_for_real_wavefunctions(wfn_a, wfn_b):
    assert cas_overlap(wfn_a, wfn_b) == pytest.approx(cas_overlap(wfn_b, wfn_a), abs=1e-12)


def test_independent_solves_at_the_same_geometry_give_plus_or_minus_one(wfn_a):
    """Required check 2 from the brief.

    Also demonstrates the arbitrary-sign problem the gauge fixing exists to handle.
    """
    mol = wfn_a.mol
    rng = np.random.default_rng(7)
    k = rng.normal(scale=0.05, size=(wfn_a.mo_coeff.shape[1],) * 2)
    mo_perturbed = np.asarray(wfn_a.mo_coeff) @ sla.expm(k - k.T)
    other = run_casscf(mol, 2, 2, mo_guess=mo_perturbed)

    assert other.energy == pytest.approx(wfn_a.energy, abs=1e-8)
    assert abs(cas_overlap(wfn_a, other)) == pytest.approx(1.0, abs=1e-6)


def test_invariant_under_active_space_rotation_with_compensating_ci_transform(wfn_a):
    """The total wavefunction is unchanged by an active-space rotation.

    Here ``transform_ci_for_orbital_rotation`` is legitimate: one geometry, one AO basis,
    a genuine unitary within the active space. It is *not* legitimate between geometries,
    which is exactly why this package computes the nonorthogonal overlap instead.
    """
    ncore, ncas = wfn_a.ncore, wfn_a.ncas
    rng = np.random.default_rng(3)
    u, _ = np.linalg.qr(rng.normal(size=(ncas, ncas)))     # random orthogonal

    mo_rot = np.array(wfn_a.mo_coeff)
    mo_rot[:, ncore:ncore + ncas] = mo_rot[:, ncore:ncore + ncas] @ u
    ci_rot = fci.addons.transform_ci_for_orbital_rotation(
        np.asarray(wfn_a.ci), ncas, wfn_a.nelecas, u
    )
    rotated = wfn_a.scaled(1.0)
    rotated.mo_coeff, rotated.ci = mo_rot, ci_rot

    assert abs(cas_overlap(wfn_a, rotated)) == pytest.approx(1.0, abs=1e-10)
    assert cas_overlap(rotated, rotated) == pytest.approx(1.0, abs=1e-10)


def test_core_rotation_does_not_change_the_wavefunction_sign(wfn_a):
    """Core orbitals are doubly occupied, so even a det=-1 core rotation leaves the sign."""
    ncore = wfn_a.ncore
    rng = np.random.default_rng(11)
    q, _ = np.linalg.qr(rng.normal(size=(ncore, ncore)))
    q[:, 0] *= -1.0                                   # force det(q) = -1
    assert np.linalg.det(q) < 0
    mo_rot = np.array(wfn_a.mo_coeff)
    mo_rot[:, :ncore] = mo_rot[:, :ncore] @ q
    rotated = wfn_a.scaled(1.0)
    rotated.mo_coeff = mo_rot
    assert cas_overlap(wfn_a, rotated) == pytest.approx(1.0, abs=1e-10)


def test_sign_flip_is_detected(wfn_a):
    assert cas_overlap(wfn_a, wfn_a.scaled(-1.0)) == pytest.approx(-1.0, abs=1e-10)


def test_cross_geometry_ao_overlap_is_not_symmetric(wfn_a, wfn_b):
    """Guards the assumption that makes this a genuinely nonorthogonal problem."""
    s = cross_mo_overlap(wfn_a, wfn_b)
    assert np.abs(s - s.T).max() > 1e-6


def test_mismatched_cas_definitions_raise(wfn_a):
    other = run_casscf(build_mol(formaldimine_geom(*GEOM_A), BASIS), 4, 4)
    with pytest.raises(ValueError, match="CAS definitions differ"):
        cas_overlap(wfn_a, other)


def test_oao_transfer_gives_orthonormal_orbitals_naive_reuse_does_not(wfn_a, wfn_b):
    """Why the OAO frame is used to carry orbitals between geometries."""
    transferred = transfer_mo(wfn_a.mol, wfn_a.mo_coeff, wfn_b.mol)
    assert orthonormality_error(wfn_b.mol, transferred) < 1e-10
    assert orthonormality_error(wfn_b.mol, np.asarray(wfn_a.mo_coeff)) > 1e-3


def test_transfer_between_identical_geometries_is_the_identity(wfn_a):
    same = build_mol(formaldimine_geom(*GEOM_A), BASIS)
    transferred = transfer_mo(wfn_a.mol, wfn_a.mo_coeff, same)
    assert np.abs(transferred - np.asarray(wfn_a.mo_coeff)).max() < 1e-10


def test_core_count():
    mol = build_mol(formaldimine_geom(*GEOM_A), BASIS)
    assert core_count(mol, 2) == 7          # 16 electrons, 2 active -> 7 doubly occupied
    assert core_count(mol, (2, 2)) == 6
    with pytest.raises(ValueError):
        core_count(mol, 3)                  # odd core electron count
