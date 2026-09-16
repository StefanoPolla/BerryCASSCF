"""Direct gap minimization. Uses analytic surfaces, so it is instant."""

import numpy as np
import pytest

from berrycasscf.minsearch import SearchResult, minimize_gap


class FakeCone:
    """Stands in for a geometry builder plus a CASSCF solve: gap = 2|displacement| from apex."""

    def __init__(self, x0, y0, slope=1.5, offset=0.0):
        self.x0, self.y0, self.slope, self.offset = x0, y0, slope, offset
        self.calls = 0

    def geom(self, x, y):
        return (float(x), float(y))

    def gap(self, x, y, geom_fn, cas, weights=(0.5, 0.5)):
        self.calls += 1
        gx, gy = geom_fn(x, y)
        return 2.0 * np.sqrt(self.slope ** 2 * ((gx - self.x0) ** 2 + (gy - self.y0) ** 2)
                             + self.offset ** 2)


@pytest.fixture
def patched(monkeypatch):
    def _make(x0, y0, **kw):
        cone = FakeCone(x0, y0, **kw)
        monkeypatch.setattr("berrycasscf.minsearch.gap_at", cone.gap)
        return cone
    return _make


def test_finds_a_cone_apex_from_a_displaced_start(patched):
    cone = patched(110.0, 90.0)
    res = minimize_gap(cone.geom, cas=None, start=(107.0, 93.0), step=(1.0, 1.0),
                       max_evaluations=200, tol=1e-6)
    assert res.x == pytest.approx(110.0, abs=1e-2)
    assert res.y == pytest.approx(90.0, abs=1e-2)
    assert res.gap < res.start_gap


def test_respects_the_evaluation_budget(patched):
    cone = patched(110.0, 90.0)
    res = minimize_gap(cone.geom, cas=None, start=(100.0, 80.0), max_evaluations=12)
    assert res.n_evaluations <= 13          # budget plus the initial probe
    assert cone.calls <= 13


def test_caches_repeated_points(patched):
    """The simplex revisits points; each geometry must only be solved once."""
    cone = patched(110.0, 90.0)
    res = minimize_gap(cone.geom, cas=None, start=(109.0, 91.0), max_evaluations=60)
    assert cone.calls == res.n_evaluations   # no evaluation repeated


def test_reports_how_far_it_moved(patched):
    cone = patched(112.0, 88.0)
    res = minimize_gap(cone.geom, cas=None, start=(110.0, 90.0), max_evaluations=200, tol=1e-6)
    assert res.moved == pytest.approx(np.hypot(2.0, 2.0), abs=5e-2)
    assert "->" in res.summary()


def test_handles_an_offset_cone_whose_minimum_is_not_zero(patched):
    """A cut that misses the seam has a nonzero floor; the search must find it, not diverge."""
    cone = patched(110.0, 90.0, offset=0.8)
    res = minimize_gap(cone.geom, cas=None, start=(107.0, 92.0), max_evaluations=200, tol=1e-6)
    assert res.gap == pytest.approx(1.6, abs=1e-2)


def test_a_failed_evaluation_does_not_kill_the_search(monkeypatch):
    calls = {"n": 0}

    def flaky(x, y, geom_fn, cas, weights=(0.5, 0.5)):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("SCF blew up")
        return float(abs(x - 110.0) + abs(y - 90.0))

    monkeypatch.setattr("berrycasscf.minsearch.gap_at", flaky)
    res = minimize_gap(lambda a, b: (a, b), cas=None, start=(108.0, 92.0),
                       max_evaluations=40)
    assert isinstance(res, SearchResult)
    assert np.isfinite(res.gap)


def test_trace_records_every_distinct_evaluation(patched):
    cone = patched(110.0, 90.0)
    res = minimize_gap(cone.geom, cas=None, start=(108.0, 88.0), max_evaluations=25)
    assert len(res.trace) == res.n_evaluations
    assert all(len(t) == 3 for t in res.trace)
