# Planned work

Things that are designed and justified but not yet implemented. Each entry says what it is, why
it is worth doing, and what would make it fail. The level of detail is deliberate: enough to fix
the intent and the failure mode, not enough to prescribe the implementation.

---

## Standing requirement: every new technique gets explained, with visuals

**This applies to every item below and to anything added later.** A technique or finding is not
finished when it runs and produces a number in `results/`. It is finished when a human reader who
has never seen it can follow *what it does and why it works* from a notebook — either its own
notebook or a clearly marked section of a more general one.

That means, for each new method:

* a short statement of the question it answers, before any code;
* **enough figures to make the mechanism visible**, not only the result. For bisection (§2), that
  means drawing the shrinking loops on the plane with their measured phases, so the reader sees the
  turnover happen; for triangulation, the two circles and their intersection; for adaptive stepping
  (§1), the achieved step sizes around the loop against the overlap that drove them;
* the failure cases shown, not just described — a plot of the thing going wrong is worth more than
  a sentence saying it can;
* numbers traceable to a file in `results/`, so the notebook reports rather than recomputes.

**Why it is a requirement and not a nicety.** These methods are the contribution; the numbers are
evidence for them. A method nobody can follow cannot be assessed, reused, or published, and the
notebooks are the only place where the reasoning and the figures sit together. The current gap in
§3 is exactly what happens when this slips.

---

## 1. Adaptive step size around the loop

**What.** Replace uniform discretization with step control driven by the adjacent overlap:

1. start at `d = dmax` (say 1/10 of the loop);
2. step by `d`, solve, measure `s = |<Psi_prev|Psi_new>|`;
3. if `s` is below the accept threshold, **reject**: restore the previous accepted state,
   shrink `d -> alpha*d`, retry;
4. otherwise accept and grow `d -> beta*d`, clipped at `dmax`.

**Why.** The per-step difficulty is very non-uniform, so uniform stepping sizes the whole loop for
its worst arc. Measured spread of adjacent overlaps, and the resulting estimate of how many fewer
points an equal-margin schedule would need:

| loop | min overlap | median | est. fewer points |
|---|---|---|---|
| formaldimine C_x CAS(2,2) | 0.9863 | 0.9910 | 1.3x |
| formaldimine C_x CAS(6,6) | 0.9782 | 0.9920 | 1.6x |
| butadiene B_x CAS(8,8) | 0.8897 | 0.9882 | 3.1x |
| butadiene B_x CAS(12,12) | 0.8807 | 0.9896 | 3.2x |
| ethylene E_x CAS(8,8) | 0.8550 | 0.9936 | **4.7x** |

It is also **complementary to single-update stepping** (`docs/findings.md` §7): adaptivity reduces
the *number of points*, which attacks the fixed per-point overhead that single-stepping cannot
touch, while single-stepping reduces updates *per* point. On the harder loops the two should
multiply to roughly 6x in wall time.

Two further reasons it suits this problem: the error indicator (`|<Psi_prev|Psi_new>|`) already
exists and costs ~0.1 ms against a solve of 30 ms to 3 min; and unequal spacing costs nothing,
because a topological readout needs the chain to be continuous and closed, not evenly sampled.

**Design points that matter.**

* **Control toward a target, do not merely reject at a floor.** A hard accept/reject at 0.80
  oscillates and wastes solves. Steer toward ~0.98 using `d_new = d * sqrt(m_target/m_achieved)`
  with `m = 1 - s` (since `1 - s ~ d^2`), clipped to `[dmin, dmax]` with a growth cap
  (`beta <= 2`). A rejection costs a full solve and there is no cheap predictor, so feedback
  control is what keeps rejections rare.
* **Land exactly on `t = 1`.** Clip the final step so the closing geometry is identical to the
  start; otherwise the endpoint estimator stops being exact, which is its whole value.
* **`dmin` plus a hard failure at it.** This is a feature, not a safeguard: if shrinking restores
  continuity the problem was discretization; if it does not, the loop is passing through something.
  That is exactly the radius-8 versus radius-6 distinction found by hand in
  `docs/findings.md` §5, and an adaptive scheme would discover it automatically.
* **Rejection must restore the previous accepted `(mo_coeff, ci)`**, not the trial. Easy to get
  subtly wrong, and the failure mode is silent.

**Reporting consequence.** Adaptive point sets make "at fixed N" comparisons meaningless; the
metric becomes cost-to-target (updates or wall time to reach a given `1 - |Pi|`), which is what
`notebooks/stepping_comparison.ipynb` already uses.

---

## 2. Bisection and triangulation to locate what loop transport encircles

**What.** Two related experiments, neither yet implemented:

* **Bisection on radius.** Shrink the loop about a fixed centre until the Berry phase flips, and
  bisect to converge on the radius at which it does. That radius locates the enclosed degeneracy
  along the centre-to-edge direction.
* **Triangulation from two centres.** Repeat from a second centre displaced along `pyr`. Each
  bisection gives a circle on which the degeneracy lies; two circles intersect in points, which
  pins the position independently of any gap scan. For butadiene the intersection is expected on
  the `tw = 90` line, so two centres suffice to reduce it to a point ("biangulation").

**Why it is now more valuable than when first discussed.** The direct gap search places the
**state-averaged** CAS(8,8) intersection at `pyr = 121.08` (gap 0.0026 mHa), while the radius
sweep bounds whatever **loop transport** encircles to `pyr < 109.85`. Those differ by ~11 deg.
Bisection plus triangulation would measure the second position directly, instead of bounding it,
and so test the central claim of `docs/findings.md` §3 — that the two methods are sensing
different objects — rather than inferring it.

**What would make it fail.** The radius sweep already showed that a loop close to the degeneracy
can be non-monotonic in `N` (radius 6 flipped sign between N = 31 and N = 61). A bisection that
trusts a single run per radius would converge on noise. Each bisection step must therefore require
the phase to be stable across at least two discretizations, which makes the procedure several
times more expensive than a naive bisection and must be budgeted for.

---

## 3. Reporting gaps in the notebooks

Results that exist in `results/` and are written up in `docs/`, but are not shown in any notebook:

* the **direct gap-minimum search** (`results/butadiene/gap_minimum_search.json`) — including the
  finding that the fitted closest approach overestimates by ~100x, and that CAS(6,6)'s
  intersection sits at `tw = 89.11`, off the line the cuts were taken along;
* the **CAS(12,12) exception** to that search (searched gap 1.496 mHa against ~0.005 at every
  other rung) together with the **5 x 5 loop-area scan** that rules out an intersection hiding off
  the sampled cross — this is the evidence behind the project's sharpest open question and appears
  in no notebook at all;
* the **localisation conclusion** from the radius sweep — the sweep itself appears, but framed only
  as a failure-mode study, not as evidence about where the encircled degeneracy is.

See §4 for the structural fix.

---

## 4. Split `active_space_study.ipynb`

At 36 cells covering two systems and several sub-studies it is doing too much, and the reporting
gaps in §3 are a symptom. Proposed split, with the per-system notebooks each self-contained:

| notebook | contents |
|---|---|
| `formaldimine_benchmark.ipynb` | unchanged — the primary benchmark |
| `ethylene_ladder.ipynb` | active-space convergence *with* an exact reference |
| `butadiene_ladder.ipynb` | no reference; drift, Berry ladder, gap search, localisation, assumption checks |
| `stepping_comparison.ipynb` | unchanged — single update vs converged |
| `summary.ipynb` (new, short) | the cross-system verdict table and the two headline figures |

---

## 5. Make the continuity diagnostics binding

**What.** Two changes to how loop transport decides a run is trustworthy:

* **Promote "continuation chain broken" from a message to a check.** When the solver's fallback
  ladder falls through to a cold start, that point was not reached by continuation at all. It is
  currently reported and then ignored by the pass/fail logic in `berry.analyse`.
* **Revisit the 0.80 adjacent-overlap threshold.** It is too permissive: the clearest suspect run
  in the whole study passed at 0.844.

**Why.** `docs/findings.md` §5 documents a run (butadiene radius 6, N=31) that passed every per-run
check and reported a phase its neighbours contradict — and it carried *both* warning signs. Gauge
fixing repairs a cold start's arbitrary sign but not a change of branch, so a broken chain is a
concrete mechanism for a wrong answer, and 0.844 is exactly the regime where a branch change hides.
Either change alone would have refused that run on its own evidence, without appeal to the
multi-discretization criterion. Layered checks are good, but a layer that never fires is not one.

**What would make it fail.** Both changes trade false positives for false negatives. Runs that are
in fact correct will start being refused — the CAS(12,12) loops report weaker warm starts at
several points and may not survive a stricter rule, which would cost the project its reference
rung. The threshold must therefore be *calibrated*, not guessed: sweep it against the runs already
in `results/` and find where it separates the stable-across-N answers from the unstable ones. If no
such value exists, that is itself worth knowing, and adaptive stepping (§1) becomes the fix instead
— it removes the failures rather than detecting them.

---

## 6. A criterion for when a gap-scan minimum counts as an intersection

**What.** Decide, explicitly and in code, when a state-averaged minimum gap is small enough to be
called a conical intersection. Some combination of an absolute tolerance and the local cone slopes
(a floor of `g` with slopes `a` means the cut misses the apex by roughly `g/a`), reported alongside
every located position instead of left to the reader.

**Why.** The project currently has no such criterion, and the absence has already produced a wrong
statement twice. CAS(12,12) butadiene reaches 1.496 mHa under direct search where every other rung
reaches ~0.005 — a factor of ~300, which is real and not a fitting, sampling or convergence
artifact. It was written up as "no intersection in this plane", which over-reads it: 1.5 mHa is
0.04 eV, inside the error of the model itself. Without a threshold, neither "there is one" nor
"there is not" is sayable, and that is precisely the quantity the comparison with loop transport
turns on.

**What would make it fail.** Any fixed tolerance is arbitrary and will be system- and
basis-dependent; a threshold chosen after seeing which answer it gives is worthless. Safer to
report the *ingredients* (floor, slopes, implied miss distance, convergence noise) and treat a
single number as a convenience, stating it once and never tuning it. The honest outcome may be that
this question cannot be settled by the gap scan at all — in which case §2 is the answer, since
bisection measures the enclosed degeneracy without needing a gap threshold.

---

## 7. Solver robustness, identified but not implemented

Three items noticed during the study and deferred. None blocks a current result; each would remove
a known caveat.

* **Deterministic multi-start for the state-averaged solve.** Cold start from RHF is
  path-independent and respects the molecular symmetry, which is why it is the default, but it is
  not guaranteed to find the *global* SA-CASSCF minimum. A measured case: CAS(4,4) ethylene cold
  gives -77.79438 against -77.79469 from an anchored start. It did not reach the observable there
  (same position, same minimum gap), but "did not this time" is not "cannot". A fixed, seeded set
  of orbital perturbations per point would make the claim checkable; the cost is a multiple of the
  scan, so it belongs on a sub-grid near the located intersection rather than everywhere.
* **Outlier repair for isolated cold-start failures.** A single point settling on a different
  solution shows up as a spike against its neighbours. Detect it from the neighbour gaps and re-solve
  that point from a neighbour's orbitals, keeping the result only if it *lowers* the energy — which
  preserves determinism in the sense that matters, since a lower stationary point is a better answer
  regardless of how it was reached.
* **Orbital-character tracking across a scan.** Nothing currently verifies that the active space at
  one grid point contains the same orbitals as at the next; active-space drift is inferred
  indirectly, from mirror asymmetry. Projecting each point's active orbitals onto the previous
  point's and reporting the smallest singular value would make drift directly visible, and is cheap.

**What would make these fail.** All three add machinery whose own failures are silent. Multi-start
that always picks the lowest energy can *lose* state-averaged solutions that are the physically
meaningful branch; outlier repair can smooth over a genuine discontinuity, which is exactly the
CAS(6,6) solution switch that is currently a *result*. Each must therefore report every time it
intervenes, and the notebooks must show those interventions rather than presenting a repaired
surface as if it came out that way.

---

## 8. A larger butadiene reference on the cluster — worth queueing, not worth waiting for

**What.** A CAS(16,16) (or larger) run for butadiene, as a SLURM job with the usual
`#### SITE ####` placeholders, covering at minimum the `tw = 90` row and a direct gap search at its
minimum, so the ladder gains a rung above CAS(12,12).

**Why, and why it is last.** Butadiene has no exact in-basis reference — full valence is CAS(22,22)
— and the ladder is still moving at the top, so every statement about "the" intersection position
is relative to a rung that is not itself converged. A higher rung would show whether the position
is settling or still wandering, and whether the 1.5 mHa floor at CAS(12,12) persists.

**But it is not on the critical path, and should not become one.** The absence of a converged
reference is itself one of the project's clearer results: it is the sharpest available answer to
"is there a system needing more than CAS(4,4)?". And the open question in `docs/findings.md` §3 is
about *the two methods disagreeing on a fixed loop*, which one more rung does not settle — §2 does,
by measuring rather than inferring. CAS(16,16) is ~5 M determinants against 853 k, so several hours
per point even on a cluster; it buys one more data point in a sequence that has already refused to
converge over six.

**Practical form.** Queue it when the cluster is free, take whatever it returns, and let the
conclusions stand or fall on the local work. If it arrives, it goes into `butadiene_ladder.ipynb`
as an extra rung; if it does not, nothing in the write-up depends on it.
