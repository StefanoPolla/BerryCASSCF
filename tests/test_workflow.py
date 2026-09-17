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


def test_the_initial_point_escalates_its_iteration_budget_when_it_does_not_converge():
    """The first point is the one every other point warm-starts from.

    It was also the only rung of the solver ladder with no recourse: continued points escalate
    through SOLVE_STRATEGIES while point 0 got a single attempt at a fixed budget. Measured on
    formaldimine CAS(6,6) on a loop of radius 1.06 deg, point 0 exhausted exactly its 200 macro
    iterations while every other point converged, the loop closed and the step floor was never
    reached -- so the whole run was refused because of the first solve's budget.

    This geometry is the cheap reproduction of the same thing: PySCF reaches |grad[o]| = 8.8e-06,
    inside the 1e-5 threshold, while dE = 4.4e-10 will not fall below conv_tol = 1e-10 within 200
    macro iterations. It converges at 239, to the same energy.
    """
    from berrycasscf.casscf import build_mol, run_casscf, run_rhf
    from berrycasscf.continuation import _solve_point
    from berrycasscf.geometry import formaldimine_geom

    geom = formaldimine_geom(130.8597, 91.9642)
    cas = CasConfig(ncas=2, nelecas=2)

    # the mechanism: 200 macro iterations is not enough here, 600 is
    mol = build_mol(geom, cas.basis)
    mf = run_rhf(mol)
    short = run_casscf(mol, 2, 2, conv_tol=cas.conv_tol, conv_tol_grad=cas.conv_tol_grad,
                       max_cycle_macro=200, mf=mf)
    assert not short.converged, "this geometry is supposed to be the hard one"

    wfn, info = _solve_point(geom, cas, previous=None, label="t")
    assert wfn.converged
    assert info["strategy"] == "initial-long"
    assert abs(wfn.energy - short.energy) < 1e-7      # same solution, just finished
    # the retry's work is charged, not hidden
    assert info["n_macro"] > 200
