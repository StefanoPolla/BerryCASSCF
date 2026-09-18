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

## Done — §1 to §5, implemented 2026-09-17

Kept here in one paragraph each so the plan reads as a record rather than a wish list. Detail is in
`docs/progress.md` (2026-09-17) and `docs/findings.md` §8.

* **§1 adaptive step size** — `berrycasscf/adaptive.py`, explained in
  `notebooks/adaptive_stepping.ipynb`. The predicted cost saving did **not** materialise: measured
  at matched quality it is 0.60x–1.23x, a wash. Its real value is that it walks loops uniform
  discretization cannot walk at any affordable N, with a resolution limit that is predicted and
  then measured.
* **§2 bisection and triangulation** — `berrycasscf/localize.py`, explained in
  `notebooks/locating_intersections.ipynb`. Works; the residual from three or more centres is a
  built-in consistency check and it correctly refuses to produce a position at formaldimine
  CAS(2,2).
* **§3 reporting gaps** — the direct gap search is now in `butadiene_ladder.ipynb`. The radius
  sweep is reframed there as evidence about the checks rather than only as a failure-mode study.
* **§4 notebook split** — `active_space_study.ipynb` became `ethylene_ladder.ipynb` and
  `butadiene_ladder.ipynb`; `summary.ipynb` added.
* **§5 binding checks** — `require_unbroken_chain` plus two checks the session's own bugs revealed
  (`loop_closed`, `closing_step_continuous`). Calibrated against all 116 saved runs. **The
  suspicion that the 0.80 overlap threshold was too permissive is withdrawn**: with the chain check
  on, no threshold from 0.70 to 0.92 admits a contradiction or a lone dissenter.

---

## 6. A criterion for when a gap-scan minimum counts as an intersection — done 2026-09-18

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

**Done**, and no threshold was invented. `examples/report_gap_criterion.py` converts each rung's
floor into the distance a cut would have to pass from a cone's apex to show it, which is a
quantity with a natural comparison — the loop radius — rather than an arbitrary one. Butadiene
CAS(12,12)'s 1.46 mHa becomes 0.83 deg against an 18 deg semi-axis, so the floor never supported
the claim that the loop encloses nothing (`docs/findings.md` §3). The third caveat there is the
live one: the conversion assumes a cone, and an avoided crossing has no apex to miss. Only the §9
measurement separates those.

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

## 8. A larger butadiene reference on the cluster — sizing measured 2026-09-18

**What.** A rung above CAS(12,12) for butadiene, so the ladder gains a point and the question
"has the position stopped moving?" gains an answer. `slurm/rung14.job` runs CAS(14,14): the
`tw = 90` gap-scan row plus the fine cut in one task, loop transport on all three loops in the
other.

**Why, and why it is last.** Butadiene has no exact in-basis reference — full valence is
CAS(22,22) — and the ladder is still moving at the top, so every statement about "the"
intersection position is relative to a rung that is not itself converged. A higher rung would
show whether the position is settling. It does **not** settle the §3 disagreement, which is
about the two methods disagreeing on a fixed loop; §9 does that, by measuring.

**Status.** Written and validated, **not submitted**. It is sized by assuming CASSCF cost tracks
determinant count (CAS(14,14) is 11.8 M against CAS(12,12)'s 853 k), which is the weakest number
in `docs/compute.md`. `slurm/bench.job` measures one solve at that active space; submit `rung14`
only once that measurement exists, and record the measurement either way.

**What would make it fail.** If a CAS(14,14) point costs hours rather than tens of minutes, the
rung is not reachable by this route and that should be *recorded as the finding* rather than left
as an unsized job in `slurm/`. Nothing in the write-up depends on it.

---

## 9. Butadiene localization — the ladder is running, CAS(12,12) will not finish in a session

**Done at CAS(2,2)** (2026-09-17): the bisection brackets rho = 0.5385 ± 0.0843 about the `B_x`
centre while the gap-scan intersection implies rho = 0.2121 — decisively outside. A second centre
reports 0 cleanly, leaving the object at pyr in (105.6, 111.5), tw in (84, 96); a third was refused
at full size, so no triangulation and no position is claimed. Cost 2 h at the cheapest rung.

**Running since 2026-09-18** on ALICE (`docs/progress.md`): `slurm/localize.job` covers CAS(4,4),
(6,6), (8,8) and (10,10) as twelve tasks — one centre each, merged per rung with `--merge` —
and `slurm/localize_cas12.job` covers the reference rung the same way. Expect 2-7 h per centre on
the ladder and 40-60 h per centre at CAS(12,12), measured from a node that is 2.7x slower per
point than the laptop.

**What still needs deciding, and it is a choice of centres.** The ladder tasks use the **default**
third centre (99, 101.85), the one that was refused at CAS(2,2) because its full-size loop grazes
the object. CAS(12,12) uses (99, 110) instead, chosen from the CAS(2,2) bracket. Neither choice can
be right for every rung: the gap-scan intersections span pyr 102-121, and a centre that brackets
cleanly at one rung grazes at another. The correct order is to pick each rung's third centre from
*its own* first-centre measurement, once that exists — so if a third centre comes back refused,
the follow-up is a single extra task with a centre chosen from that rung's bracket, not a re-run.

**Why it matters.** `docs/findings.md` §3 is the sharpest open question in the project. Every
objection to reading it as a disagreement has been tested and answered; what remains is to
**measure** the position of whatever loop transport encircles at the rung where the gap scan finds
nothing, which is what bisection does and the only way to interrogate the state-specific object
directly.

**What would make it fail.** §10 below: if a loop encloses an even number of degeneracies the
phase is 0, and bisection would measure the wrong boundary.

---

## 10. Mirror pairs: an even number of enclosed degeneracies reads as zero

**What.** Handle the case where a symmetry plane forces degeneracies to come in pairs, so a loop
enclosing both reports 0 and is indistinguishable from a loop enclosing none.

**Why.** This is not hypothetical. A direct map of the state-specific in-CAS gap for formaldimine
is **exactly symmetric about phi = 90**, so any degeneracy off that line has a mirror partner. A
loop centred on the line grows to enclose both at nearly the same radius, and its Berry phase never
flips — which would make bisection report "nothing here" for a region containing two intersections.
It is a plausible contributor to the inconsistent triangulation measured at formaldimine CAS(2,2).

**What would make it work.** Centres deliberately placed *off* the symmetry line break the tie,
because the two partners are then enclosed at different radii and the phase flips twice. That is
cheap to add — it is a choice of centres, not new machinery — but it needs the bisection to look
for *two* transitions rather than stopping at the first.

---

## 11. Smaller items the session surfaced

* **`d_min` is a resolution knob, not a safety knob.** It sets how close a loop may pass
  (`eps_min ~ 2*pi*R*d_min/dtheta_max`). Bisection's precision is therefore bounded by it, and it
  should probably be tightened for localization runs specifically, at a cost in rejected solves.
  Worth measuring rather than guessing.
* **The formaldimine FCI reference is a one-dimensional cut** (17 points in alpha at fixed phi), so
  it fixes `alpha_x` along that line and does *not* establish that the minimum over `phi` is at 90.
  Any claim comparing a triangulated position to it in two dimensions is weaker than it looks.
* **Charge the rejected solves everywhere.** `LoopTraversal.total_micro` now includes rejected
  adaptive trials. Any future cost comparison must keep doing so; hiding them flatters the method.
* **Drivers must merge, not overwrite.** One single-rung re-run destroyed a five-rung result file,
  recoverable only from git. `run_gap_minimum_search.py` is fixed; the other drivers that write a
  single summary file should be audited for the same pattern.
