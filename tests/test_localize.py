"""Bisection and triangulation geometry, tested without any electronic structure.

The expensive part of locating a degeneracy this way is the loop transport; the part that
is easy to get subtly wrong is the geometry that turns transition radii into a position.
These tests pin the geometry down on synthetic inputs where the answer is exact.
"""

import json
import math

import pytest

from berrycasscf.localize import (
    BisectionResult,
    RadiusProbe,
    bisect_radius,
    _circle_intersections,
    elliptical_radius,
    merge_centre_records,
    restore_bisection,
    resume_bisections,
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


# --- resuming a killed run ------------------------------------------------------------
#
# A localization at a large active space runs for days, and the record used to be written
# only after the last centre, so a walltime kill lost all of it. These pin the round trip
# that makes a resubmission continue instead of starting over.

def test_a_saved_bisection_round_trips():
    original = _bisection((90.0, 101.85), rho=0.5385, shape=(12.0, 18.0), half_width=0.0843)
    restored = restore_bisection(original.to_dict())
    # In memory the centre is a tuple; through JSON it comes back a list. Both must work,
    # since a resumed cluster job always takes the JSON path.
    assert tuple(restored.centre) == (90.0, 101.85)
    assert tuple(restore_bisection(json.loads(json.dumps(original.to_dict()))).centre) \
        == (90.0, 101.85)
    assert restored.rho == pytest.approx(0.5385)
    assert restored.rho_uncertainty == pytest.approx(0.0843)
    assert restored.bracketed
    # The derived keys to_dict() adds must not be mistaken for constructor fields.
    assert restore_bisection({**original.to_dict(), "rho": 999.0}).rho == pytest.approx(0.5385)


def test_resume_takes_the_leading_centres_that_match():
    centres = [(90.0, 101.85), (90.0, 90.0), (99.0, 110.0)]
    saved = [_bisection(c, rho=0.5).to_dict() for c in centres[:2]]
    assert len(resume_bisections(saved, centres)) == 2
    assert len(resume_bisections([], centres)) == 0


def test_resume_stops_at_the_first_changed_centre():
    """A record made with different centres measured different ellipses; it is not reusable."""
    centres = [(90.0, 101.85), (90.0, 90.0), (99.0, 110.0)]
    saved = [_bisection(centres[0], rho=0.5).to_dict(),
             _bisection((85.0, 95.0), rho=0.5).to_dict()]
    assert len(resume_bisections(saved, centres)) == 1


def test_resume_ignores_extra_saved_centres():
    """Asking for fewer centres than are saved uses only the ones asked for."""
    centres = [(90.0, 101.85), (90.0, 90.0)]
    saved = [_bisection(c, rho=0.5).to_dict()
             for c in centres + [(99.0, 110.0)]]
    assert len(resume_bisections(saved, centres[:1])) == 1


# --- merging centres that ran as separate jobs ----------------------------------------

def _record(bisections):
    return {"system": "butadiene", "bisections": [b.to_dict() for b in bisections]}


def test_merge_combines_one_bisection_per_record():
    centres = [(90.0, 101.85), (90.0, 90.0), (99.0, 110.0)]
    records = [_record([_bisection(c, rho=0.5 + i * 0.1)]) for i, c in enumerate(centres)]
    merged = merge_centre_records(records, centres)
    assert [tuple(r.centre) for r in merged] == centres
    assert [r.rho for r in merged] == pytest.approx([0.5, 0.6, 0.7])


def test_merge_refuses_a_record_centred_somewhere_else():
    """Ellipses from a different construction must not be silently combined."""
    centres = [(90.0, 101.85), (90.0, 90.0)]
    records = [_record([_bisection(centres[0], rho=0.5)]),
               _record([_bisection((12.0, 34.0), rho=0.5)])]
    with pytest.raises(ValueError, match="centred on"):
        merge_centre_records(records, centres)


def test_merge_refuses_a_record_holding_more_than_one_bisection():
    centres = [(90.0, 101.85)]
    records = [_record([_bisection(centres[0], rho=0.5), _bisection((90.0, 90.0), rho=0.6)])]
    with pytest.raises(ValueError, match="expected 1"):
        merge_centre_records(records, centres)


# --- resuming one bisection, probe by probe -------------------------------------------
#
# A bisection at a large active space is days long -- one centre is one cluster job -- so a
# walltime kill must not cost the whole measurement. Replay is only trustworthy if the
# control flow is a pure function of the verdicts, and that is what these check, with a
# synthetic oracle in place of the CASSCF so the logic is tested rather than the chemistry.

def _oracle(rho_true=0.5):
    """A degeneracy at elliptical distance rho_true, and a counter of real evaluations."""
    calls = []

    def evaluate(centre, shape, scale, **kw):
        calls.append(scale)
        return RadiusProbe(scale=scale, radius=(shape[0] * scale, shape[1] * scale),
                           verdict=("pi" if scale > rho_true else "zero"),
                           reason="synthetic", cost_micro=1, wall_time=1.0)

    return evaluate, calls


def test_bisection_brackets_a_synthetic_degeneracy(monkeypatch):
    evaluate, calls = _oracle(0.5)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    res = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0, tol=0.02)
    assert res.bracketed
    assert res.lo <= 0.5 <= res.hi
    assert res.rho == pytest.approx(0.5, abs=0.05)
    assert len(calls) > 2


def test_a_resumed_bisection_reaches_the_same_bracket_without_redoing_work(monkeypatch):
    evaluate, calls = _oracle(0.5)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    saved: list[list[dict]] = []
    full = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0, tol=0.02,
                         on_probe=lambda ps: saved.append([p.to_dict() for p in ps]))
    assert len(saved) == len(calls)                    # every probe was checkpointed

    # Kill it after three probes and resume from what was written.
    partial = saved[2]
    calls.clear()
    resumed = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0,
                            tol=0.02, cached_probes=partial)
    assert (resumed.lo, resumed.hi) == (full.lo, full.hi)
    assert [p.scale for p in resumed.probes] == [p.scale for p in full.probes]
    assert len(calls) == len(full.probes) - 3          # the first three were not re-solved


def test_a_fully_cached_bisection_solves_nothing(monkeypatch):
    evaluate, calls = _oracle(0.5)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    saved: list[list[dict]] = []
    full = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0, tol=0.02,
                         on_probe=lambda ps: saved.append([p.to_dict() for p in ps]))
    calls.clear()
    again = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0, tol=0.02,
                          cached_probes=saved[-1])
    assert (again.lo, again.hi) == (full.lo, full.hi)
    assert calls == []
    # The cost carried in the record is the whole measurement's, not the last session's.
    assert again.total_micro == full.total_micro


def test_merge_refuses_an_unfinished_checkpoint():
    """A checkpoint looks like a result: one bisection, right centre, no bracket yet."""
    centres = [(90.0, 101.85)]
    partial = _record([_bisection(centres[0], rho=0.5)])
    partial["complete"] = False
    with pytest.raises(ValueError, match="unfinished checkpoint"):
        merge_centre_records([partial], centres)


# --- anchors that move -----------------------------------------------------------------
#
# A bisection needs an outer loop that encloses and an inner one that does not, and neither
# is known in advance. Before these, a single anchor that came back refused ended the whole
# bisection with nothing measured -- which is how ethylene CAS(6,6) burned 18 minutes while
# its outer loop had cleanly reported pi.

def _band_oracle(rho_true=0.5, band=0.0):
    """Verdicts for a degeneracy at rho_true, refused within `band` of it (the grazing zone)."""
    calls = []

    def evaluate(centre, shape, scale, **kw):
        calls.append(scale)
        if abs(scale - rho_true) <= band:
            verdict = "undetermined"
        else:
            verdict = "pi" if scale > rho_true else "zero"
        return RadiusProbe(scale=scale, radius=(shape[0] * scale, shape[1] * scale),
                           verdict=verdict, reason="synthetic", cost_micro=1, wall_time=1.0)

    return evaluate, calls


def test_inner_anchor_shrinks_when_the_degeneracy_is_inside_it(monkeypatch):
    """A degeneracy closer than scale_lo used to be reported as 'not bracketed'."""
    evaluate, calls = _band_oracle(0.05)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    res = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.08, scale_hi=1.0, tol=0.05)
    assert res.bracketed
    assert res.lo < 0.05 < res.hi
    assert min(calls) < 0.08                       # it really did shrink the anchor


def test_outer_anchor_grows_when_the_degeneracy_is_outside_it(monkeypatch):
    evaluate, calls = _band_oracle(1.5)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    res = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.05, scale_hi=1.0, tol=0.05)
    assert res.bracketed
    assert res.lo < 1.5 < res.hi
    assert max(calls) > 1.0


def test_a_refused_anchor_is_moved_rather_than_fatal(monkeypatch):
    """The ethylene CAS(6,6) case: the inner anchor lands in the refusal band."""
    evaluate, _ = _band_oracle(0.09, band=0.02)     # 0.08 falls inside the band
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", evaluate)
    res = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.08, scale_hi=1.0, tol=0.05)
    assert res.bracketed
    assert res.lo < 0.09 < res.hi


def test_anchor_moves_are_bounded(monkeypatch):
    """A loop that can never be walked must stop, not escalate forever."""
    def always_refused(centre, shape, scale, **kw):
        return RadiusProbe(scale=scale, radius=(shape[0] * scale, shape[1] * scale),
                           verdict="undetermined", reason="synthetic", cost_micro=1,
                           wall_time=1.0)
    monkeypatch.setattr("berrycasscf.localize.evaluate_radius", always_refused)
    res = bisect_radius((90.0, 100.0), (12.0, 18.0), scale_lo=0.08, scale_hi=1.0,
                        max_anchor_moves=3)
    assert not res.bracketed
    assert len(res.probes) <= 2 + 2 * 3            # both anchors, each moved at most 3 times
