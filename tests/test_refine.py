"""Locating a cone apex from gap data. Instant: pure numerics, no electronic structure."""

import numpy as np
import pytest

from berrycasscf.refine import ApexFit, cone_apex, parabolic_apex


def ideal_cone(x, x0, slope, offset):
    """Gap along a straight cut passing at perpendicular offset from a conical intersection."""
    return 2.0 * np.sqrt(slope ** 2 * (np.asarray(x) - x0) ** 2 + offset ** 2)


@pytest.mark.parametrize("x0", [100.0, 110.37, 123.9])
@pytest.mark.parametrize("offset", [0.0, 0.3, 1.5])
def test_cone_fit_is_exact_for_an_ideal_cone(x0, offset):
    """The model is exact, so the apex must be recovered regardless of grid alignment."""
    x = np.arange(80.0, 141.0, 5.0)
    g = ideal_cone(x, x0, 1.8, offset)
    fit = cone_apex(x, g, window=3)
    assert fit.method == "cone"
    assert fit.position == pytest.approx(x0, abs=1e-6)
    assert fit.closest_approach == pytest.approx(2.0 * offset, abs=1e-6)


def test_cone_fit_beats_the_parabolic_one_off_grid():
    """A parabola fitted to the gap is dragged toward the grid minimum; the cone fit is not."""
    x0 = 110.37
    x = np.arange(80.0, 141.0, 5.0)
    g = ideal_cone(x, x0, 1.8, 0.9)
    cone_err = abs(cone_apex(x, g, window=3).position - x0)
    para_err = abs(parabolic_apex(x, g).position - x0)
    assert cone_err < 1e-6
    assert para_err > 0.1
    assert cone_err < para_err


def test_recovers_the_closest_approach_the_grid_only_bounds():
    """The grid minimum overestimates how close the cut passes; the fit does not."""
    x = np.arange(80.0, 141.0, 5.0)
    g = ideal_cone(x, 110.37, 1.8, 1.0)
    fit = cone_apex(x, g, window=3)
    assert fit.closest_approach == pytest.approx(2.0, abs=1e-6)
    assert g.min() > fit.closest_approach          # the grid never samples the apex


def test_slope_is_recovered():
    x = np.arange(80.0, 141.0, 2.0)
    fit = cone_apex(x, ideal_cone(x, 110.0, 1.8, 0.5), window=6)
    assert fit.slope == pytest.approx(2.0 * 1.8, rel=1e-6)


def test_falls_back_when_there_are_too_few_points():
    x = np.array([100.0, 105.0, 110.0])
    fit = cone_apex(x, ideal_cone(x, 104.0, 1.0, 0.2))
    assert fit.method.startswith("parabolic")


def test_falls_back_when_the_data_do_not_look_like_a_cone():
    """An inverted (downward) cut has no apex; the estimator must not invent one."""
    x = np.arange(0.0, 11.0, 1.0)
    g = 10.0 - (x - 5.0) ** 2 * 0.1        # a maximum, not a minimum
    fit = cone_apex(x, g)
    assert fit.method.startswith("parabolic")


def test_edge_minimum_is_reported_not_extrapolated():
    x = np.arange(0.0, 6.0, 1.0)
    g = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])     # minimum at the left edge
    fit = parabolic_apex(x, g)
    assert fit.method == "parabolic-edge"
    assert fit.position == x[0]


def test_residual_flags_a_cut_that_is_not_conical():
    """A quartic cut is not a cone; the fit should still return, but with a large residual."""
    x = np.arange(80.0, 141.0, 5.0)
    good = cone_apex(x, ideal_cone(x, 110.0, 1.8, 0.5), window=3)
    bad = cone_apex(x, 1.0 + ((x - 110.0) / 10.0) ** 4, window=3)
    assert good.residual < 1e-10
    assert bad.residual > good.residual


def test_nan_points_are_ignored():
    x = np.arange(80.0, 141.0, 5.0)
    g = ideal_cone(x, 110.37, 1.8, 0.9)
    g[0] = np.nan
    assert cone_apex(x, g, window=3).position == pytest.approx(110.37, abs=1e-6)


def test_apexfit_is_reportable():
    fit = cone_apex(np.arange(0.0, 11.0), ideal_cone(np.arange(0.0, 11.0), 5.0, 1.0, 0.1))
    assert isinstance(fit, ApexFit)
    assert "cone" in repr(fit)
