"""Loop parameterization and geometry construction. No electronic structure: instant."""

import numpy as np
import pytest

from berrycasscf.geometry import (
    Loop,
    formaldimine_geom,
    LOOP_CI,
    LOOP_CONTROL_LOWER,
    LOOP_CONTROL_UPPER,
    REFERENCE_CI_ALPHA_PHI,
    FORMALDIMINE_FROZEN,
)


def test_geometry_matches_brief_verbatim():
    """The Z-matrix must reproduce the parameterization given in the project brief."""
    geom = formaldimine_geom(130.0, 90.0)
    assert "1.498047" in geom and "1.066797" in geom and "0.987109" in geom
    assert "118.359375" in geom
    lines = [ln.strip() for ln in geom.strip().splitlines()]
    assert [ln.split()[0] for ln in lines] == ["N", "C", "H", "H", "H"]
    # Last line is "H 1 <r_NH> 2 <alpha> 3 <phi>": alpha is the HNC angle (token 4),
    # phi the HNCH dihedral (token 6).
    tok = lines[-1].split()
    assert tok[0] == "H" and tok[1] == "1" and tok[3] == "2" and tok[5] == "3"
    assert (float(tok[4]), float(tok[6])) == (130.0, 90.0)


def test_loop_closes():
    for loop in (LOOP_CI, LOOP_CONTROL_LOWER, LOOP_CONTROL_UPPER):
        assert np.allclose(loop.closing_point(), loop.points()[0])


def test_loop_points_are_distinct_and_on_the_ellipse():
    loop = LOOP_CI.with_n_points(12)
    pts = np.array(loop.points())
    assert len(pts) == 12
    assert len({tuple(np.round(p, 9)) for p in pts}) == 12
    da = (pts[:, 0] - loop.centre[0]) / loop.radius[0]
    dp = (pts[:, 1] - loop.centre[1]) / loop.radius[1]
    assert np.allclose(da**2 + dp**2, 1.0)


def test_enclosure_of_the_reference_conical_intersection():
    """The benchmark loops must bracket the literature CI location as intended."""
    assert LOOP_CI.encloses(*REFERENCE_CI_ALPHA_PHI)
    assert not LOOP_CONTROL_LOWER.encloses(*REFERENCE_CI_ALPHA_PHI)
    assert not LOOP_CONTROL_UPPER.encloses(*REFERENCE_CI_ALPHA_PHI)


def test_loop_settings_match_the_auto_oo_tutorial():
    """Guard against silent drift from the documented upstream values."""
    assert LOOP_CI.centre == (130.0, 89.9)
    assert LOOP_CI.radius == (10.0, 10.0)
    assert LOOP_CI.phase == pytest.approx(np.pi / 20)
    assert LOOP_CONTROL_LOWER.centre == (110.0, 89.9)
    assert LOOP_CONTROL_UPPER.centre == (150.0, 89.9)
    assert FORMALDIMINE_FROZEN == (1.498047, 1.066797, 0.987109, 118.359375)


def test_bounding_box():
    lo_a, hi_a, lo_p, hi_p = LOOP_CI.bounding_box(margin=2.0)
    assert (lo_a, hi_a) == (118.0, 142.0)
    assert (lo_p, hi_p) == pytest.approx((77.9, 101.9))
