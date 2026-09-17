"""Bisection and triangulation geometry, tested without any electronic structure.

The expensive part of locating a degeneracy this way is the loop transport; the part that
is easy to get subtly wrong is the geometry that turns transition radii into a position.
These tests pin the geometry down on synthetic inputs where the answer is exact.
"""

import math

import pytest

from berrycasscf.localize import (
    BisectionResult,
    _circle_intersections,
    elliptical_radius,
    triangulate,
)


def _bisection(centre, rho, shape=(10.0, 10.0), half_width=1e-4):
    """A bracketed result with a known rho, for exercising triangulate()."""
    return BisectionResult(centre=centre, shape=shape, cas_label="CAS(2,2)", probes=[],
                           lo=rho - half_width, hi=rho + half_width)


# --- elliptical radius ----------------------------------------------------------------

def test_elliptical_radius_is_one_on_the_loop_boundary():
    shape = (12.0, 18.0)
    centre = (90.0, 100.0)
    for t in (0.0, 0.3, 0.77):
        pt = (centre[0] + shape[0] * math.cos(2 * math.pi * t),
              centre[1] + shape[1] * math.sin(2 * math.pi * t))
        assert elliptical_radius(pt, centre, shape) == pytest.approx(1.0)


def test_elliptical_radius_matches_the_enclosure_test_of_the_loop_itself():
    from berrycasscf.geometry import Loop
    loop = Loop("L", (90.0, 100.0), (12.0, 18.0))
    for pt in [(95.0, 104.0), (90.0, 130.0), (101.0, 100.0), (90.0, 100.1)]:
        inside = elliptical_radius(pt, loop.centre, loop.radius) < 1.0
        assert inside == loop.encloses(*pt), pt


# --- circle intersection --------------------------------------------------------------

def test_two_circles_meet_in_a_mirror_pair():
    pts = _circle_intersections((0.0, 0.0), 1.0, (1.0, 0.0), 1.0)
    assert len(pts) == 2
    assert pts[0][0] == pytest.approx(0.5) and pts[1][0] == pytest.approx(0.5)
    assert pts[0][1] == pytest.approx(-pts[1][1])


def test_circles_that_cannot_meet_return_nothing():
    assert _circle_intersections((0.0, 0.0), 1.0, (5.0, 0.0), 1.0) == []   # too far
    assert _circle_intersections((0.0, 0.0), 5.0, (0.1, 0.0), 1.0) == []   # nested


# --- triangulation --------------------------------------------------------------------

def test_two_centres_recover_a_known_position_up_to_the_mirror_ambiguity():
    truth, shape = (132.6, 90.0), (10.0, 10.0)
    centres = [(130.0, 89.9), (136.0, 89.9)]
    res = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape)
                       for c in centres])
    assert len(res.candidates) == 2
    assert min(math.hypot(p[0] - truth[0], p[1] - truth[1])
               for p in res.candidates) < 1e-3
    # the pair is mirrored in the line joining the centres (here phi = 89.9)
    assert res.candidates[0][0] == pytest.approx(res.candidates[1][0])


def test_a_third_centre_off_the_line_removes_the_ambiguity():
    truth, shape = (132.6, 90.0), (10.0, 10.0)
    centres = [(130.0, 89.9), (136.0, 89.9), (133.0, 97.0)]
    res = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape)
                       for c in centres])
    assert res.chosen is not None
    assert math.hypot(res.chosen[0] - truth[0], res.chosen[1] - truth[1]) < 1e-3
    assert res.residual < 1e-6


def test_prefer_breaks_the_two_centre_tie_and_says_so():
    truth, shape = (132.6, 90.0), (10.0, 10.0)
    centres = [(130.0, 89.9), (136.0, 89.9)]
    res = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape)
                       for c in centres], prefer=(132.0, 90.0))
    assert math.hypot(res.chosen[0] - truth[0], res.chosen[1] - truth[1]) < 1e-3
    assert "chose" in res.note


def test_an_anisotropic_loop_shape_is_handled_in_scaled_coordinates():
    """The loops are ellipses, not circles: tw and pyr are different coordinates."""
    truth, shape = (89.0, 110.0), (12.0, 18.0)
    centres = [(90.0, 101.85), (90.0, 120.0), (96.0, 110.0)]
    res = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape)
                       for c in centres])
    assert math.hypot(res.chosen[0] - truth[0], res.chosen[1] - truth[1]) < 1e-3


def test_inconsistent_measurements_are_reported_not_averaged():
    shape = (10.0, 10.0)
    res = triangulate([_bisection((130.0, 90.0), 0.1, shape),
                       _bisection((160.0, 90.0), 0.1, shape)])
    assert res.candidates == []
    assert "do not intersect" in res.note


def test_triangulation_refuses_unbracketed_input():
    shape = (10.0, 10.0)
    unbracketed = BisectionResult(centre=(130.0, 90.0), shape=shape, cas_label="CAS(2,2)",
                                  probes=[], lo=None, hi=0.5)
    with pytest.raises(ValueError):
        triangulate([unbracketed, _bisection((136.0, 89.9), 0.3, shape)])


# --- the residual has to be judged against the measurement precision -------------------

def test_the_residual_is_normalised_by_the_bracket_widths():
    """Raw residuals are not comparable between runs of different precision.

    Measured on formaldimine: raw residuals 0.197 (CAS(2,2)), 0.0088 (CAS(4,4)) and 0.0574
    (CAS(6,6)) would suggest CAS(6,6) is 6.5x worse than CAS(4,4). Normalised by the RMS
    bracket half-width they are 1.57, 0.29 and 0.76 -- CAS(2,2) is the only one whose misfit
    exceeds its own precision, and it is the only one independently known to be pathological.
    """
    truth, shape = (132.6, 90.0), (10.0, 10.0)
    centres = [(130.0, 89.9), (137.0, 89.9), (133.0, 97.0)]

    tight = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape,
                                    half_width=1e-3) for c in centres])
    loose = triangulate([_bisection(c, elliptical_radius(truth, c, shape), shape,
                                    half_width=0.1) for c in centres])

    # same (near-zero) misfit, very different precision
    assert tight.residual == pytest.approx(loose.residual, abs=1e-9)
    assert loose.uncertainty > 50 * tight.uncertainty
    assert tight.consistent and loose.consistent


def test_a_misfit_larger_than_the_precision_is_reported_inconsistent():
    shape = (10.0, 10.0)
    centres = [(130.0, 89.9), (137.0, 89.9), (133.0, 97.0)]
    truth = (132.6, 90.0)
    # third centre's radius deliberately wrong by far more than its stated precision
    results = [_bisection(c, elliptical_radius(truth, c, shape), shape, half_width=1e-3)
               for c in centres[:2]]
    results.append(_bisection(centres[2], elliptical_radius(truth, centres[2], shape) + 0.3,
                              shape, half_width=1e-3))
    tri = triangulate(results)
    assert tri.residual_over_uncertainty > 1.0
    assert tri.consistent is False
