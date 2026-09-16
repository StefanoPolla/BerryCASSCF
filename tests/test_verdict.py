"""The pass/fail verdict logic, exercised without any electronic structure.

Synthetic traversals only: these run instantly and pin down exactly which diagnostic
causes a refusal, which the end-to-end tests cannot do because they need a real failure to
occur first.
"""

import numpy as np
import pytest

from berrycasscf.berry import analyse, NONTRIVIAL, TRIVIAL, UNDETERMINED
from berrycasscf.config import ContinuationConfig
from berrycasscf.continuation import LoopTraversal, PointRecord


def _point(index, strategy="warm-mo+warm-ci", overlap=0.99, converged=True):
    return PointRecord(
        index=index, t=index / 8.0, alpha=130.0, phi=90.0, energy=-93.0,
        converged=converged,
        raw_overlap_with_prev=None if index == 0 else overlap,
        gauge_sign=1,
        abs_overlap_with_prev=None if index == 0 else abs(overlap),
        casci_gap=0.3, guess_orthonormality_error=1e-14, max_mo_change=0.1,
        strategy="initial" if index == 0 else strategy,
        n_macro=3, n_micro=9, wall_time=0.1,
    )


def _traversal(points, closing=-0.98, endpoint=-1.0):
    return LoopTraversal(
        loop_name="T", loop_centre=(130.0, 90.0), loop_radius=(10.0, 10.0),
        n_points=len(points), cas_label="CAS(2,2)", basis="sto-3g",
        points=points, closing_overlap=closing, endpoint_overlap=endpoint,
        endpoint_converged=True, wall_time=1.0, warnings=[],
    )


def test_a_clean_traversal_passes_and_reports_pi():
    res = analyse(_traversal([_point(k) for k in range(8)]))
    assert res.status == "OK"
    assert res.berry_phase == NONTRIVIAL
    assert all(res.checks.values())


def test_a_positive_closing_overlap_reports_the_trivial_phase():
    res = analyse(_traversal([_point(k) for k in range(8)], closing=0.98, endpoint=1.0))
    assert res.status == "OK"
    assert res.berry_phase == TRIVIAL


def test_a_cold_solve_after_the_first_point_breaks_the_chain_and_is_refused():
    """The butadiene radius-6 failure mode: a run that passes every *other* check.

    That run (docs/findings.md §5) had min overlap 0.844 and endpoint 0.984 -- comfortably
    inside both thresholds -- and reported a phase its neighbours contradict. The only
    evidence against it was a cold solve mid-loop, which was reported and not acted on.
    """
    pts = [_point(k) for k in range(8)]
    pts[6] = _point(6, strategy="cold", overlap=0.844)
    res = analyse(_traversal(pts, closing=-0.55, endpoint=-0.984))

    assert res.checks["continuity"] is True          # 0.844 clears the 0.80 threshold
    assert res.checks["endpoint_is_same_state"] is True
    assert res.checks["continuation_chain_intact"] is False
    assert res.status == "FAILED"
    assert res.berry_phase == UNDETERMINED
    assert any("chain broken" in m for m in res.messages)


def test_the_first_point_is_cold_by_definition_and_does_not_break_the_chain():
    pts = [_point(k) for k in range(8)]
    pts[0] = _point(0, strategy="initial")
    assert analyse(_traversal(pts)).checks["continuation_chain_intact"] is True


def test_a_point_that_converged_nowhere_also_breaks_the_chain():
    pts = [_point(k) for k in range(8)]
    pts[3] = _point(3, strategy="none-converged", converged=False)
    res = analyse(_traversal(pts))
    assert res.checks["continuation_chain_intact"] is False


def test_the_chain_check_can_be_switched_off_for_studies_that_need_it():
    pts = [_point(k) for k in range(8)]
    pts[6] = _point(6, strategy="cold")
    cont = ContinuationConfig(require_unbroken_chain=False)
    res = analyse(_traversal(pts), cont)
    assert res.checks["continuation_chain_intact"] is True
    assert res.status == "OK"


def test_a_degraded_warm_start_is_not_a_broken_chain():
    """warm-mo and warm-mo-loose still continue the orbitals, so the chain holds."""
    for strategy in ("warm-mo", "warm-mo-loose"):
        pts = [_point(k) for k in range(8)]
        pts[4] = _point(4, strategy=strategy)
        res = analyse(_traversal(pts))
        assert res.checks["continuation_chain_intact"] is True, strategy
        assert res.status == "OK"


def test_estimator_disagreement_is_a_failure_not_a_result():
    res = analyse(_traversal([_point(k) for k in range(8)], closing=-0.98, endpoint=+1.0))
    assert res.checks["estimators_agree"] is False
    assert res.berry_phase == UNDETERMINED


def test_lost_continuity_is_refused():
    pts = [_point(k) for k in range(8)]
    pts[2] = _point(2, overlap=0.42)
    res = analyse(_traversal(pts))
    assert res.checks["continuity"] is False
    assert res.berry_phase == UNDETERMINED


# --- a walk that never closed must not yield a phase -----------------------------------

def _unclosed(points, closing, endpoint=None):
    trav = _traversal(points, closing=closing, endpoint=endpoint)
    trav.endpoint_converged = None
    trav.adaptive = {"closed": False, "hit_step_floor": True, "n_rejected": 6}
    return trav


def test_a_walk_that_never_closed_the_loop_is_refused():
    """A real failure: an adaptive walk stopped a fifth of the way round and passed.

    The endpoint checks are skipped when there is no endpoint, every accepted step was
    continuous, the chain was intact -- so the run came back OK with phase "trivial (0)",
    read off a closing overlap between two points that are nowhere near each other.
    """
    res = analyse(_unclosed([_point(k) for k in range(5)], closing=+0.93))
    assert res.checks["loop_closed"] is False
    assert res.status == "FAILED"
    assert res.berry_phase == UNDETERMINED
    assert any("not a closed loop" in m for m in res.messages)


def test_a_discontinuous_closing_step_is_refused_even_if_the_walk_claims_to_be_closed():
    """Defence in depth: |<Psi_last|Psi_0>| alone carries the sign of the product estimator."""
    trav = _traversal([_point(k) for k in range(8)], closing=-0.31, endpoint=-1.0)
    res = analyse(trav)
    assert res.checks["closing_step_continuous"] is False
    assert res.status == "FAILED"
    assert res.berry_phase == UNDETERMINED


def test_a_uniform_traversal_closes_by_construction_and_pays_nothing():
    res = analyse(_traversal([_point(k) for k in range(8)]))
    assert res.checks["loop_closed"] is True
    assert res.checks["closing_step_continuous"] is True
    assert res.status == "OK"
