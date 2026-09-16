"""Adaptive step control, tested where the answer is known in closed form.

Everything here runs on the 2x2 Jahn-Teller model or on the controller alone, so the whole
file is instant and no CASSCF failure can masquerade as a control failure.
"""

import numpy as np
import pytest

from berrycasscf.adaptive import control_step, walk_adaptive
from berrycasscf.config import AdaptiveConfig
from berrycasscf.toy import jt_loop_adaptive, jt_loop_berry_phase


# --- the controller on its own --------------------------------------------------------

def test_the_controller_holds_its_step_when_the_target_is_met():
    cfg = AdaptiveConfig()
    assert control_step(0.05, cfg.target_mismatch, cfg) == pytest.approx(0.05)


def test_the_controller_shrinks_above_target_and_grows_below():
    cfg = AdaptiveConfig()
    assert control_step(0.05, 4 * cfg.target_mismatch, cfg) == pytest.approx(0.025)
    assert control_step(0.02, cfg.target_mismatch / 4, cfg) == pytest.approx(0.04)


def test_the_controller_respects_every_clip():
    cfg = AdaptiveConfig(d_max=0.1, d_min=1e-3, growth_cap=2.0, max_shrink=0.1)
    assert control_step(0.09, 0.0, cfg) == 0.1                 # growth cap then d_max
    assert control_step(0.05, 1e-12, cfg) == 0.1               # never exceeds d_max
    assert control_step(0.05, 1e9, cfg) == pytest.approx(0.005)  # max_shrink caps the shrink
    assert control_step(1e-3, 1e9, cfg) == pytest.approx(1e-3)  # never below d_min


# --- the walk on the Jahn-Teller model ------------------------------------------------

def test_an_enclosing_loop_gives_pi_and_a_control_loop_gives_zero():
    inside, _ = jt_loop_adaptive(centre=(0.0, 0.0), radius=1.0)
    outside, _ = jt_loop_adaptive(centre=(3.0, 0.0), radius=1.0)
    assert inside.agrees_with_theory and inside.is_nontrivial
    assert outside.agrees_with_theory and not outside.is_nontrivial


def test_the_walk_lands_exactly_on_t_equals_one():
    """The endpoint estimator is only exact if the closing geometry is identical."""
    _, walk = jt_loop_adaptive(radius=1.0)
    assert walk.closed
    assert walk.t_values[-1] == pytest.approx(1.0, abs=1e-12)


def test_the_answer_is_insensitive_to_the_arbitrary_sign_of_each_solve():
    signs = {jt_loop_adaptive(radius=1.0, seed=s)[0].is_nontrivial for s in range(6)}
    assert signs == {True}


def test_every_accepted_step_is_inside_the_accept_threshold():
    cfg = AdaptiveConfig()
    _, walk = jt_loop_adaptive(centre=(0.0, 0.25), radius=1.0, cfg=cfg)
    accepted = [e for e in walk.events if e.accepted]
    assert accepted
    assert max(e.mismatch for e in accepted) <= cfg.accept_mismatch


def test_a_rejection_does_not_advance_the_walk():
    """The previous accepted state must be restored untouched -- a silent failure mode."""
    cfg = AdaptiveConfig(d_max=0.5, target_mismatch=0.02, accept_mismatch=0.05)
    _, walk = jt_loop_adaptive(centre=(0.0, 0.0), radius=1.0, cfg=cfg)
    assert walk.n_rejected > 0, "a huge d_max must force at least one rejection"
    # t advances only on accepted events, and strictly.
    ts = walk.t_values
    assert all(b > a for a, b in zip(ts, ts[1:]))
    assert len(ts) == 1 + sum(1 for e in walk.events if e.accepted)


def test_accepted_events_line_up_one_per_chain_point():
    """Guards a real bug: indexing `events` by point slips by one for every rejection.

    A walk with one rejection reported an adjacent overlap of 0.79 for a chain whose
    accepted steps were all above 0.90 -- the rejected trial's overlap, attributed to a
    point that never used it.
    """
    cfg = AdaptiveConfig(d_max=0.5, accept_mismatch=0.05)
    _, walk = jt_loop_adaptive(centre=(0.0, 0.0), radius=1.0, cfg=cfg)
    assert walk.n_rejected > 0
    assert len(walk.accepted_events) == len(walk.t_values) - 1
    assert all(e.accepted for e in walk.accepted_events)
    # every accepted event must join consecutive accepted t values
    for ev, (a, b) in zip(walk.accepted_events, zip(walk.t_values, walk.t_values[1:])):
        assert ev.t_from == pytest.approx(a) and ev.t_to == pytest.approx(b)


def test_reported_min_overlap_never_undercuts_the_accept_threshold():
    cfg = AdaptiveConfig(d_max=0.5, accept_mismatch=0.05)
    res, walk = jt_loop_adaptive(centre=(0.0, 0.0), radius=1.0, cfg=cfg)
    assert walk.n_rejected > 0
    assert res.min_abs_adjacent_overlap >= 1.0 - cfg.accept_mismatch - 1e-12


def test_the_step_size_adapts_rather_than_staying_uniform():
    """A loop far off-centre is easy on one arc and hard on the other."""
    _, walk = jt_loop_adaptive(centre=(0.0, 0.9), radius=1.0)
    d = walk.step_sizes
    assert d.max() / d.min() > 3.0


def test_a_loop_through_the_degeneracy_hits_the_step_floor_instead_of_lying():
    """Shrinking cannot restore continuity when the loop passes through the degeneracy.

    This is the informative failure: the walk reports that the loop is passing through
    something rather than returning a phase that happens to have a sign.
    """
    _, walk = jt_loop_adaptive(centre=(1.0, 0.0), radius=1.0)   # passes exactly through (0,0)
    assert walk.hit_step_floor
    assert not walk.closed


def test_adaptive_resolves_a_near_degenerate_loop_that_uniform_stepping_fails():
    """The result that makes bisection possible: shrink the step, recover the answer.

    A loop passing 0.02 from the degeneracy defeats a uniform walk at N=24 -- an adjacent
    overlap collapses to 0.69 and the run is refused. The adaptive walk spends its solves
    where they are needed, closes the loop with a worst overlap of 0.92, and does it in
    about the same number of points.
    """
    centre, radius = (0.0, 0.98), 1.0      # inner edge passes 0.02 from the degeneracy
    uniform = jt_loop_berry_phase(centre=centre, radius=radius, n_points=24)
    assert uniform.min_abs_adjacent_overlap < 0.80    # a refused run

    res, walk = jt_loop_adaptive(centre=centre, radius=radius)
    assert walk.closed and not walk.hit_step_floor
    assert res.min_abs_adjacent_overlap >= 0.90
    assert res.agrees_with_theory and res.is_nontrivial
    assert res.n_points < 2 * uniform.n_points        # not bought with brute force


def test_uniform_stepping_can_be_accidentally_right_while_failing_its_own_check():
    """Keep 'gets the right answer' separate from 'can be trusted without knowing it'.

    At every closest approach tested, uniform N=24 returns the correct sign. From 0.1
    inwards it also fails its own continuity threshold, so a user with no reference would
    have to refuse it. Being right is not the same as being trustworthy.
    """
    res = jt_loop_berry_phase(centre=(0.0, 0.9), radius=1.0, n_points=24)
    assert res.is_nontrivial                          # the right answer ...
    assert res.min_abs_adjacent_overlap < 0.80        # ... from a run that must be refused


def test_the_step_floor_obeys_its_predicted_resolution_limit():
    """How close a loop may pass before the walk gives up is predictable, not mysterious.

    A loop of radius R passing at distance eps sweeps the state's angle at d(theta)/dt
    ~ 2*pi*R/eps near closest approach. Holding the per-step angle below
    ``dtheta_max = 2*arccos(1 - accept_mismatch)`` therefore needs a step of about
    ``eps*dtheta_max/(2*pi*R)``, and the walk gives up when that falls below ``d_min``:

        eps_min ~ 2*pi*R*d_min/dtheta_max

    With the defaults this is 0.0070, i.e. 0.7% of the loop radius. Measured: the walk
    closes at eps = 0.007 and hits the floor at eps = 0.005.
    """
    cfg = AdaptiveConfig()
    dtheta_max = 2.0 * np.arccos(1.0 - cfg.accept_mismatch)
    eps_min = 2.0 * np.pi * 1.0 * cfg.d_min / dtheta_max
    assert eps_min == pytest.approx(0.0070, abs=5e-4)

    _, ok = jt_loop_adaptive(centre=(0.0, 1.0 - 0.007), radius=1.0)
    _, floored = jt_loop_adaptive(centre=(0.0, 1.0 - 0.005), radius=1.0)
    assert ok.closed and not ok.hit_step_floor
    assert floored.hit_step_floor
