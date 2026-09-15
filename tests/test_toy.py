"""Jahn-Teller toy model: validates the sign convention and both estimators. Instant."""

import numpy as np
import pytest

from berrycasscf.toy import jt_loop_berry_phase, jt_ground_state, jt_gap


@pytest.mark.parametrize("n_points", [5, 8, 13, 24, 64])
def test_enclosing_loop_is_nontrivial(n_points):
    res = jt_loop_berry_phase(centre=(0.0, 0.0), radius=1.0, n_points=n_points)
    assert res.product_estimator < 0
    assert res.endpoint_estimator < 0
    assert res.agrees_with_theory


@pytest.mark.parametrize("centre", [(2.0, 0.0), (0.0, 5.0), (-3.0, 3.0)])
def test_non_enclosing_loop_is_trivial(centre):
    res = jt_loop_berry_phase(centre=centre, radius=1.0, n_points=24)
    assert res.product_estimator > 0
    assert res.endpoint_estimator > 0
    assert res.agrees_with_theory


def test_the_two_estimators_always_agree_in_sign():
    for centre in [(0, 0), (0.5, 0), (2, 0), (0, 0.9), (0, 1.1)]:
        res = jt_loop_berry_phase(centre=centre, radius=1.0, n_points=32)
        assert np.sign(res.product_estimator) == np.sign(res.endpoint_estimator)


def test_result_is_independent_of_arbitrary_state_signs():
    """The estimators must be gauge-invariant: random +-1 per point changes nothing."""
    for centre, expect in [((0.0, 0.0), -1.0), ((3.0, 0.0), 1.0)]:
        signs = {
            np.sign(
                jt_loop_berry_phase(centre=centre, radius=1.0, n_points=11, seed=s).product_estimator
            )
            for s in range(50)
        }
        assert signs == {expect}


def test_just_inside_versus_just_outside_the_seam():
    """The delicate case: loops that barely do or do not enclose the degeneracy."""
    inside = jt_loop_berry_phase(centre=(0.99, 0.0), radius=1.0, n_points=64)
    outside = jt_loop_berry_phase(centre=(1.01, 0.0), radius=1.0, n_points=64)
    assert inside.product_estimator < 0 < outside.product_estimator


def test_model_basics():
    assert jt_gap(3.0, 4.0) == pytest.approx(10.0)
    v = jt_ground_state(1.0, 0.0)
    assert np.linalg.norm(v) == pytest.approx(1.0)
