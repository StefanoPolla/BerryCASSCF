"""End-to-end smoke tests for the two workflows.

Short but real: a coarse loop and a tiny scan grid on formaldimine/STO-3G. These are the
reproducible short-run checks required by the brief; the full benchmark lives in
``examples/``.
"""

import numpy as np
import pytest

from berrycasscf import (
    CasConfig,
    ContinuationConfig,
    ScanConfig,
    LOOP_CI,
    LOOP_CONTROL_UPPER,
    REFERENCE_CI_ALPHA_PHI,
    run_loop,
    scan_gap,
)
from berrycasscf.berry import NONTRIVIAL, TRIVIAL
from berrycasscf.continuation import traverse_loop
from berrycasscf.overlap import cas_overlap

CAS22 = CasConfig(ncas=2, nelecas=2)

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def ci_loop_run():
    return run_loop(LOOP_CI.with_n_points(9), cas=CAS22)


def test_ci_enclosing_loop_gives_a_nontrivial_berry_phase(ci_loop_run):
    res, _ = ci_loop_run
    assert res.status == "OK"
    assert res.berry_phase == NONTRIVIAL
    assert res.product_estimator < 0
    assert res.endpoint_estimator < 0


def test_control_loop_gives_a_trivial_berry_phase():
    res, _ = run_loop(LOOP_CONTROL_UPPER.with_n_points(9), cas=CAS22)
    assert res.status == "OK"
    assert res.berry_phase == TRIVIAL
    assert res.product_estimator > 0


def test_endpoint_overlap_is_essentially_plus_or_minus_one(ci_loop_run):
    """The closing geometry is identical to the start, so the overlap must be unimodular."""
    res, _ = ci_loop_run
    assert abs(abs(res.endpoint_estimator) - 1.0) < 1e-6


def test_the_two_estimators_agree(ci_loop_run):
    res, _ = ci_loop_run
    assert np.sign(res.product_estimator) == np.sign(res.endpoint_estimator)


def test_diagnostics_are_recorded_for_every_point(ci_loop_run):
    _, trav = ci_loop_run
    assert len(trav.points) == 9
    assert all(p.converged for p in trav.points)
    # the first point has no predecessor; every later one has a real overlap
    assert trav.points[0].abs_overlap_with_prev is None
    assert all(p.abs_overlap_with_prev is not None for p in trav.points[1:])
    # orbitals transferred into each geometry stay orthonormal there
    errs = [p.guess_orthonormality_error for p in trav.points[1:]]
    assert max(errs) < 1e-10


def test_gauge_fixing_actually_fires(ci_loop_run):
    """Independent solves do return arbitrary signs; the run must be correcting some."""
    _, trav = ci_loop_run
    flips = [p.gauge_sign for p in trav.points if p.gauge_sign == -1]
    assert flips, "expected at least one sign flip to be corrected along the loop"


def test_result_is_stable_against_refining_the_discretization():
    """The Z2 answer must not depend on N; only |Pi| should approach 1."""
    coarse, _ = run_loop(LOOP_CI.with_n_points(9), cas=CAS22)
    fine, _ = run_loop(LOOP_CI.with_n_points(17), cas=CAS22)
    assert coarse.berry_phase == fine.berry_phase == NONTRIVIAL
    assert abs(fine.product_estimator) > abs(coarse.product_estimator)


def test_energies_are_continuous_around_the_loop(ci_loop_run):
    _, trav = ci_loop_run
    e = trav.energies
    assert np.abs(np.diff(e)).max() < 0.05      # Hartree, no discontinuous jumps


def test_traversal_without_endpoint_solve_still_yields_the_product_estimator():
    trav = traverse_loop(LOOP_CI.with_n_points(9), cas=CAS22, solve_endpoint=False)
    assert trav.endpoint_overlap is None
    from berrycasscf.berry import analyse
    res = analyse(trav)
    assert res.berry_phase == NONTRIVIAL


def test_sa_scan_finds_the_conical_intersection_inside_the_loop():
    """CAS(4,4) is needed here: CAS(2,2) cannot resolve S1 (see docs/results.md)."""
    res = scan_gap(
        LOOP_CI,
        cas=CasConfig(ncas=4, nelecas=4),
        scan=ScanConfig(n_alpha=7, n_phi=7, margin=0.0),
    )
    alpha, phi, gap = res.min_gap_point()
    assert res.converged.all()
    assert gap < 0.01                                   # Hartree
    assert LOOP_CI.encloses(alpha, phi)
    assert abs(alpha - REFERENCE_CI_ALPHA_PHI[0]) < 6.0
    assert abs(phi - REFERENCE_CI_ALPHA_PHI[1]) < 6.0
