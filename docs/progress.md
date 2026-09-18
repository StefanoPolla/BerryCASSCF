# Progress log

Append-only. Each entry: what was implemented or tested, key decisions, current results,
unresolved issues, next step.

---

## 2026-09-15 — Orientation, literature, backend choice

**Implemented / tested.** Repository skeleton, project-local `.venv` (Python 3.12.10,
PySCF 2.14.0, NumPy 2.5.3, SciPy 1.18.1).

**Literature.** Archived in `docs/literature/`:
- arXiv:2304.06070 (the Berry-phase algorithm paper) — PDF;
- arXiv:2009.11417 — PDF;
- `three_loops_FCI.png` and `formaldimine.png` from the `auto_oo` repository.

**Key decisions.**
- The brief's link `EmielKoridon/auto_oo` 404s; the package now lives at **`Emieeel/auto_oo`**.
  The tutorial notebook was retrieved from there and is the source of the loop definitions.
- **Backend: PySCF**, chosen after a capability check against Psi4 and OpenMolcas. The
  deciding factor is that PySCF is the only one exposing all of: CI vectors in a documented
  determinant-string basis, MO coefficients as plain arrays, and — critically —
  `gto.intor_cross`, which gives AO overlap integrals **between two different geometries**.
  Full reasoning in `docs/plan.md` §1.
- The notebook defines only the CI-enclosing loop. The two control loops come from
  arXiv:2304.06070 §6.3 and Fig. 1a. All provenance recorded in `docs/provenance.md`.

**Next step.** Implement the nonorthogonal overlap.

---

## 2026-09-15 — Core implementation and validation

**Implemented.** `geometry`, `overlap`, `casscf`, `config`, `continuation`, `berry`, `scan`,
`toy`, `store`, `report`.

**Key decisions.**
- **Nonorthogonal overlap by Schur-complement core elimination.** The common doubly-occupied
  core is eliminated analytically, leaving `det(S_cc)^2` times a Frobenius contraction over
  minors of a single `ncas x ncas` matrix. Exact, and cheap enough that overlaps are never
  the bottleneck. PySCF's `transform_ci_for_orbital_rotation` is *not* used for cross-geometry
  overlaps anywhere, as the brief requires.
- **Orbital continuity through the OAO frame**, taken from the auto_oo notebook. Measured:
  a transferred guess is orthonormal at the new geometry to `6e-15`, whereas naively re-using
  the previous geometry's AO coefficients gives an orthonormality error of `3.6e-2`.
- **Two estimators** (cyclic product, gauge-invariant by construction; and the exact
  same-geometry endpoint overlap), reported together and cross-checked.

**Validation results.**
- `<A|A> = 1` to `1e-13`.
- Factorized overlap vs brute-force determinant-pair sum at two *different* geometries:
  agree to `4.4e-16`.
- Two independent CASSCF solves at the same geometry: overlap `-1.000000000` — the arbitrary
  sign problem is real and is detected cleanly. This is what the gauge fixing exists for.
- Invariance under an active-space rotation with a compensating CI transform: `1.0` to `1e-10`.
- Jahn-Teller toy model: Berry phase `pi` inside / `0` outside for every loop tested,
  including the delicate pair centred at `(0.99, 0)` and `(1.01, 0)` with radius 1. Verified
  gauge-invariant over 200 random sign realizations. (The brief made this optional; it was
  cheap and is kept as a permanent regression test.)

**Next step.** Run the benchmark sweep.

---

## 2026-09-15 — Two real bugs found and fixed

**1. Triplet contamination in the state-averaged comparator.** PySCF's default
`direct_spin1` solver returns the Ms=0 **triplet** as a root. The first SA-CASSCF scan was
therefore reporting an S0/T1 gap, not S0/S1: at every geometry tested the lowest "state" had
`<S^2> = 2`. The resulting gap map was uniformly small (≤0.008 Ha) with no cone structure.

Fixed by spin-adapting **every** CAS solve in the package
(`berrycasscf.casscf.apply_singlet_constraint`, applied before the state-average wrapper).
After the fix all roots have `<S^2> = 0`.

**2. Over-tight `conv_tol_grad` broke warm starts.** With `conv_tol_grad = 1e-6`, a CASSCF
solve that was warm-started in *both* orbitals and CI vector failed its convergence test
while sitting on the correct energy to `8e-11`. Isolated by varying the warm start:

| start | converged | dE from reference |
|---|---|---|
| warm MO + warm CI | **no** | `+8.2e-11` |
| warm MO only | yes | `+1.5e-10` |
| cold MO + warm CI | yes | `-1e-13` |
| cold | yes | — |

The predictor lands so close to the solution that the macro iterations rattle around below
the energy threshold without meeting an over-tight gradient threshold. Default changed to
`conv_tol_grad = 1e-5` (PySCF's own default is `sqrt(conv_tol)`); the affected solve then
converges in 0.10 s instead of failing after 0.76 s.

**Also added.** A four-rung solver fallback ladder (warm MO+CI → warm MO → warm MO with a
relaxed gradient threshold → cold). The rung used is recorded per point, and any mid-loop
fall-through to a cold start raises a warning, because it breaks the continuation chain.
Every rung is still overlap-checked, so a fallback that lands on a different branch is caught
by the continuity test rather than silently accepted.

**Next step.** Complete the sweep and write up.

---

## 2026-09-15 — Formaldimine benchmark complete

**Implemented / run.** Full sweep: 3 loops × 3 active spaces × 4 discretizations = 36
Berry-phase traversals (~4 min total), plus SA-CASSCF gap scans at CAS(2,2)/(4,4)/(6,6) and an
FCI reference.

**Results.** All 36 runs succeeded, and all nine (loop, CAS) combinations are **stable** by the
stated criterion: `C_x` → π, `C_1`/`C_2` → 0, unchanged across N = 9, 17, 25, 41. The endpoint
estimator came out at exactly ±1.000000 in every run. Details in `docs/results.md`.

**Key finding — the two methods disagree at CAS(2,2), and the comparator is wrong.** The
SA-CASSCF(2,2) gap scan reports a 0.34 mHa near-degeneracy *inside the trivial loop* `C_2` and
finds nothing inside the true CI loop `C_x` — two of three loops backwards. The Berry phase at
the same CAS(2,2) gets all three right.

Diagnosed and understood: recomputing the gap at the spurious minimum with CAS(4,4), CAS(6,6)
and CAS(8,8) gives 0.038, 0.038 and 0.040 Ha respectively — three successively larger spaces
agree there is no degeneracy there. The CAS(2,2) result is a converged artifact. The reason the
two workflows have different requirements is that the Berry phase is state-specific (it only
tracks S0) while the comparator must place S1, which CAS(2,2) cannot represent.

**Minimal reliable CAS: CAS(4,4)**, being the smallest tested at which both methods are stable
and mutually consistent on all three loops. CAS(6,6) changes no conclusion, which is the
evidence that CAS(4,4) is converged rather than lucky.

**FCI calibration (new).** STO-3G is small enough for exact FCI (1.66 M determinants, ~28 s per
point), so `examples/run_fci_reference.py` was added and run along φ = 89.9°. The gap is linear
on both sides of a minimum at **α× = 132.6°**, matching the paper's ≈132°. This resolves a
previously open question: the CAS(4,4)→CAS(6,6)→FCI sequence for the CI position
(130.0° → 130.9° → 132.6°) shows the ~2° offset is active-space truncation converging
monotonically. The Berry phase is unaffected — a 2° error is irrelevant to a loop of radius 10°.

**Also observed.** `1 − |Π|` scales as `N^-0.89` (fitted over all nine series), consistent with
the expected `O(1/N)`: adjacent states differ by an angle of order 1/N, so each overlap is
`1 − O(1/N²)` and the product of N of them is `1 − O(1/N)`.

**Deliverables.** `notebooks/formaldimine_benchmark.ipynb` (executed, with figures), `docs/results.md`,
`docs/limitations.md`, `docs/compute.md` with measured timings, SLURM templates.

**Next step.** The follow-up system.

---

## 2026-09-15 — Follow-up system prepared (fulvene), not run

**Decision: fulvene (C6H6).** Chosen over benzene because its two driving coordinates — the
exocyclic C=C stretch and the methylene torsion — form a natural two-dimensional plane, whereas
benzene's S1/S0 seam needs out-of-plane ring puckering and a defensible 2D plane is harder to
justify. The full π manifold is CAS(6,6): 400 determinants against formaldimine's 4, with a
genuine π→π\* S1 for the comparator to resolve. Reasoning in `docs/followup.md`.

**Implemented.** `berrycasscf/fulvene.py` (RHF/6-31G reference geometry, two rigid distortion
coordinates), two-stage drivers, `slurm/fulvene_scan.sbatch` and `slurm/fulvene_berry.sbatch`,
18 tests.

**Design decision — the CI position is not assumed.** Stage 1 locates the intersection by
SA-CASSCF scan; stage 2 places the loops around whatever stage 1 found. This avoids committing
to a literature CI geometry that would have to be taken on trust, and matches how the
calculation would actually be done.

**Verified locally (machinery, not the experiment).** Geometry builder exact to 1e-10 with the
ring frozen and no atomic collisions across the search region; CASSCF(6,6)/STO-3G converges in
0.8 s with `ncore = 18`; OAO transfer orthonormal to 1.7e-14; and the factorized nonorthogonal
overlap still matches the brute-force determinant-pair reference with 18 core orbitals — the
check that the Schur-complement core elimination scales past formaldimine's 7 core orbitals.
The overlap between two neighbouring fulvene geometries came out at −0.974, so the arbitrary
CASSCF sign appears here too.

**Not run.** Per the brief, the production fulvene calculation is prepared for cluster
submission rather than executed locally. Two jobs, ~4–11 h and ~1–3 h at CAS(6,6)/6-31G\*.
**Awaiting submission by the user.**

**Unresolved.** See `docs/limitations.md`. The most interesting open question is whether the
CAS(2,2) Berry phase is right by construction or by luck — a direct check would be to compute
the overlap between the CAS(2,2) and CAS(6,6) S0 wavefunctions around the loop, which needs a
generalization of the overlap code to two different active-space definitions.

---

## 2026-09-16 — Active-space convergence study (ethylene)

**Question.** Is there a system for which the workflows need larger active spaces than
formaldimine's CAS(2,2) (Berry) and CAS(4,4) (state-averaged comparator)?

**Fulvene was tried first and abandoned.** Its S1/S0 intersection is not reachable with the
rigid two-coordinate model in `berrycasscf/fulvene.py`: scanning the whole (exocyclic bond
length, methylene torsion) plane at SA-CASSCF(6,6)/STO-3G gives a gap that never falls below
59 mHa and decreases monotonically toward the corner, and extending the stretch to 2.5 A at full
twist does not close it either. The real fulvene intersection needs ring deformation, which that
model freezes. Recorded rather than forced; the fulvene module and its cluster jobs are unchanged
and still valid for the two-stage scan-then-loop protocol in `docs/followup.md`.

**Ethylene** was used instead: the twisted-pyramidalized S1/S0 intersection, two rigid
coordinates (CH2 torsion `tau`, umbrella pyramidalization `phi`). **6-31G\* is required** —
STO-3G has no intersection along these coordinates at all (gap bottoms out near 10 mHa), because
twisted-ethylene S1 is the zwitterionic V state. CAS(12,12) is the full valence space and serves
as the in-basis reference.

**A real defect in the scan was found and fixed.** Ethylene's mirror symmetry (`tau` and
`180 - tau` geometries are exactly isometric) means a correct gap map must be symmetric about
`tau = 90`. Warm-started CAS(8,8) returned **19.1** and **29.4 mHa** at two mirror-image
geometries, against **26.9 mHa** cold: the active space drifts along the scan path, so the answer
depended on the route taken to a geometry. SA-CASSCF turns out to have more than one stationary
point here, and a warm sweep reaches the lower one only from some directions. Added a
path-independent `anchor` strategy (now the `ScanConfig` default) and made mirror asymmetry a
reported reliability criterion. **The formaldimine results are unaffected** — recomputing its
CAS(4,4) scan reproduces the committed one to 2.5e-06 Ha across the whole grid; drift needs an
active space large enough to have somewhere to drift to.

**Result — the intersection position across the ladder** (reference CAS(12,12): `phi` = 110.00):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| `phi` | 110.02 | 114.85 | 102.09 | 103.11 | 109.10 | 110.00 |
| error | +0.02 | +4.85 | −7.91 | −6.90 | −0.90 | ref |

Convergence is **not monotonic**: the minimal pi space is essentially exact, the intermediate
spaces are displaced by 5–8 degrees in both directions, and the reference is recovered only from
CAS(10,10). Every rung passes the mirror-symmetry test (1e-11 to 1e-5 mHa), so this is the
physics of truncation, not noise.

**Answer.** For the comparator, **yes** — but the useful statement distinguishes two questions.
The smallest active space that *gets the right answer* is CAS(2,2); the smallest that could be
*trusted without already knowing the answer* is CAS(10,10), since only there does enlarging the
space stop moving the result — which is the test anyone would actually apply. Formaldimine fails
the opposite way (CAS(2,2) qualitatively wrong, CAS(4,4) already converged), which is the
argument for checking convergence per system rather than reusing an active-space recipe.

For the **Berry phase, no**: CAS(2,2) suffices on both systems, and it is insensitive to the
8-degree misplacement that defeats the comparator, because a loop of radius 12 degrees encloses
the intersection either way. Its cost with active space is *discretization* — min adjacent
overlap on the CI-enclosing loop at N=13 falls from 0.895 at CAS(2,2) to ~0.67 at CAS(8,8), so
the larger spaces fail the continuity test at N=13 and pass at N=21. That failure is detected and
reported, never silently absorbed.

**Deliverables.** `berrycasscf/ethylene.py`, `examples/run_active_space_study.py`,
`notebooks/active_space_study.ipynb`, `docs/active_space.md`, 28 new tests including a regression test
for the path-independence fix.

**Next step.** A system whose S1 is genuinely doubly excited (a polyene 2Ag state, e.g.
butadiene or hexatriene) would test whether a case exists where *no* small active space is
accurate, as opposed to ethylene's accidental-but-unverifiable CAS(2,2). That needs a CI search
in a new coordinate plane and is cluster work.


---

## 2026-09-16 — Cold standard applied to all results; butadiene added

**Audit (task 1).** All ten warm-started formaldimine SA scans were recomputed cold and compared
against the committed ones. **Every minimum position and every inside/outside conclusion is
unchanged**, so no documented conclusion moved. The grids do differ away from the minima, and in
every case *warm* found the lower energy — at one `C_2` CAS(2,2) corner point a cold start
converges **0.115 Ha above** the warm solution. That is the measured price of choosing
reproducibility; it changes nothing here but could matter elsewhere. Berry records were untouched:
they come from `traverse_loop`, whose warm start is the method rather than an optimization.

**Butadiene (task 2).** Chosen because its 2¹Ag state carries a large doubly-excited component.
Its full valence space is CAS(22,22), so unlike the first two systems **there is no exact
reference**.

The plane was searched for, not assumed: three candidate planes of four rigid coordinates, only
`tw_pyr` containing an intersection (0.73 mHa; `tw_tc` bottoms out at 99.5 mHa, `tw_bend` at
36.2). The search reports strain because the plane's *unrestricted* minimum (0.19 mHa) sits where
the methylene has folded onto its own C–C bond, 266 kcal/mol up.

**Butadiene has no `tw` → `180 − tw` symmetry** — its two methylene hydrogens are inequivalent.
Reusing ethylene's validation would have compared unrelated geometries. The exact symmetry is
`(tw, pyr)` → `(−tw, −pyr)`; every rung passes it at exactly 0.0 mHa. Both facts are locked in by
tests.

**Result — the comparator does not converge anywhere on the ladder:**

| CAS | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|
| `pyr` | 109.83 | 114.81 | 121.05 | 105.04 | 101.85 |

Spread 19.2°, top two rungs still 3.2° apart. This is the clearest answer to whether some system
needs more than CAS(4,4): here even CAS(12,12) is not demonstrably enough.

**A prediction that failed, instructively.** CAS(8,8) places its intersection 19.35° from the
loop centre, outside the 18° radius, so it was predicted to report a trivial phase. **It reports
π, as do all four rungs** (controls trivial throughout, all 24 runs `OK`). The assumption behind
the prediction — that the state-averaged gap minimum is the point the Berry phase encircles — is
false: the Berry phase transports the *state-specific* ground state, whose degeneracy need not
coincide with a state-averaged gap minimum when the active space is truncated. Shrinking the loop
bounds the degeneracy it actually encircles to `pyr < 113.9`, not 121.2.

**The inversion.** The Berry phase gives the same answer at every rung while the comparator's CI
estimate scatters over 19° and never settles — the topological method is more stable than the
quantity meant to validate it. That does not make it right by default (with no reference, "π at
every rung" could be four consistent errors), but it is consistent, passes every internal check,
and costs a fraction of the comparator.

**Next step.** Localize the degeneracy the state-specific transport senses, by loop-shrinking with
enough discretization to keep continuity, and compare it against the state-averaged seam. That
would turn the `pyr < 113.9` bound into a measurement and test directly whether state-specific
transport is less truncation-sensitive than the state-averaged surface.


---

## 2026-09-16 — Refinement fixed, CAS(2,2) added, assumptions attacked

**The apex estimator was wrong and is fixed.** It fitted a parabola to the gap, but near a
conical intersection the gap is linear, so a cut is a V and the fit is dragged toward the grid
minimum. For an ideal cone a cut at perpendicular offset satisfies `gap² = a²(x−x₀)² + b²` —
exactly a parabola in `gap²` — so fitting that recovers the apex without bias and also returns the
cut's closest approach. Exact on synthetics; shifts real positions by up to ~1.8 deg, **including
the ethylene reference itself** (110.00 → 110.90). Every verdict is unchanged, which is worth
recording: the correction matters for the numbers, not for the conclusions.

**CAS(2,2) on butadiene** — deliberately below the π manifold — returns π on the enclosing loop and
0 on both controls, with *higher* overlaps than any larger rung. Loop transport now succeeds at the
smallest active space on all three systems.

**Assumptions attacked, three claims overturned:**

1. **The loop-shrinking bound.** Withdrawn: it rested on one passing run at radius 6, and N=61
   flipped its sign.
2. **A per-run check passed a probably-wrong answer** (radius 6, N=31: overlap 0.844, endpoint
   0.984, reported π, while N=15 and N=61 both say 0). Only the multi-N stability criterion
   rejects it. The protocol is safe because the checks are layered, not because any one is sound.
3. **"The seam misses the plane at CAS(12,12)."** Withdrawn as stated, twice over. The
   observation was that loop transport returns π there (stable at N=13 and N=21) while the gap
   scan finds no degeneracy. Two objections were raised and both have since been *tested and
   answered*: the 1.7 mHa floor is not a line-fitting artifact (a free 2D search reaches only
   1.496 mHa, against ~0.005 at every other rung), and it is not an intersection hiding off the
   sampled cross (a 5×5 grid over the loop's area has its minimum at the centre). What remains
   withdrawn is the *interpretation*: 1.5 mHa is 0.04 eV, inside this model's own error, so it
   cannot be called an absence of an intersection without a threshold nobody has set. The claim
   that survives is relative — the search behaves differently at this rung by ~300x — and the open
   question is now sharper rather than closed.

**Assumptions that survived:** the intersection sits at tw = 90 for every rung tested (89.982,
89.992, 89.958), with no symmetry forcing it.

**Also found:** the CAS(6,6) cut is a *solution discontinuity*, not a cone (0.85 mHa then a jump to
4.1 mHa over 0.25 deg), so its position is not a fit; and the fitted "closest approach" on a 1-deg
grid is a resolution artifact — sub-degree cuts drop CAS(4,4) from 0.92 to 0.098 mHa.

**Infrastructure.** Runs now log to `logs/<job>.log` in the repository, one line per completed
point with elapsed time and ETA, so a slow job is distinguishable from a dead one; previously the
scan logged per grid row, which meant half an hour of silence at CAS(12,12). Notebooks renamed to
`formaldimine_benchmark.ipynb` and `active_space_study.ipynb`.

**Next step.** A genuine 2D scan at CAS(12,12) over the loop interior (~1 h for a coarse 5×5) to
settle whether the intersection is off the sampled cross. Until then that rung's position is a
point on a cross, not a located intersection.

---

## 2026-09-17 — Adaptive stepping, loop localization, and three checks that should have existed

Implements `docs/todo.md` §1–§5. Two of this session's most useful results are bugs the new code
found in the old code, and one is a negative result about the new code itself.

### §5 — the continuation-chain check, calibrated rather than guessed

`examples/calibrate_thresholds.py` re-decides all **116** saved Berry records under candidate
thresholds and scores them by *self-consistency*: runs of the same loop and active space at
different N form a "question", and physics forbids them to disagree. Two failure modes are counted
— a **contradiction** (two passing runs disagree) and a **lone dissenter** (one run passes while a
refused run of the same question says the opposite, so the multi-N criterion cannot fire).

| chain check | min overlap | decidable | contradictions | lone dissenters |
|---|---|---|---|---|
| off | 0.70 | 42 | **1** | 0 |
| off | 0.80 | 40 | 0 | **1** |
| off | 0.86 | 35 | 0 | 0 |
| **on** | **0.80** | **39** | **0** | **0** |

**The overlap threshold was never the danger.** With the chain check on, *no* threshold from 0.70
to 0.92 admits either failure mode. Without it, 0.80 has to rise to 0.86 to be safe, costing five
decidable questions against the chain check's one. The earlier suspicion that 0.80 was too
permissive (`docs/findings.md` §5) is therefore **withdrawn**: 0.80 stays, and
`require_unbroken_chain` defaults to True.

The one run in the corpus with a broken chain that the check now refuses is butadiene `B_x`
CAS(4,4) N=21, which reports π consistently with its ladder — so the cost is real and is a *lost
confirmation*, not a lost answer.

### §1 — adaptive step control (`berrycasscf/adaptive.py`)

Step size steered by the measured continuity, `d_new = d*sqrt(m_target/m)`, clipped, with the final
step clipped so the walk lands exactly on `t = 1` and the endpoint estimator stays exact. Validated
first on the Jahn–Teller model, where a **resolution limit can be predicted and then measured**:

    eps_min ~ 2*pi*R*d_min/dtheta_max  =  0.0070 of the loop radius

The walk closes at eps = 0.007 and hits the floor at 0.005. Uniform N=24 returns the *correct sign*
at every eps tested and fails its own continuity check from eps = 0.1 inwards — right, but not
trustworthy.

**Negative result, reported as one: adaptive stepping never beat uniform on cost.** At matched
quality: 0.60x–1.23x on formaldimine, 0.67x–1.15x on butadiene CAS(2,2), **0.75x–1.04x on ethylene
CAS(8,8)** — the loop predicted to save 4.7x, the largest prediction in the project. All three
predictions (1.3x, 3.1x, 4.7x) are **retracted**.

The cause was measured rather than guessed. Cost is not proportional to point count: hard points
are also expensive points (r between per-point micro-iterations and adjacent overlap is −0.40 at
formaldimine CAS(2,2), −0.67 at CAS(6,6)). On ethylene `E_x` CAS(8,8), comparing the zero-rejection
adaptive run against uniform at equal worst-case overlap — so rejection cost is excluded — adaptive
uses **0.79x the points at 1.47x the cost each**, netting 1.16x. Rejected trials make it worse
again where they occur.

Where it *does* pay is not cost: it walks loops uniform discretization cannot walk at any
affordable N (an order of magnitude closer to a degeneracy), removes the need to guess N, and
diagnoses whether a failure was undersampling or proximity.

### §2 — bisection and triangulation (`berrycasscf/localize.py`)

Shrinking a loop about a fixed centre measures the **elliptical distance** to the degeneracy;
intersecting the ellipses from several centres gives a position. Two centres leave a mirror pair,
three resolve it, and with three or more the construction is over-determined so its **residual is a
built-in consistency check**. The geometry is validated exactly on synthetic input.

Plain bisection was useless here: the first midpoint lands in the refusal band around the
transition and the search gives up. It now tracks `lo`, `hi` *and* the refusal band, narrowing the
two usable gaps alternately, so every probe either classifies or shrinks the next target interval.

**The formaldimine ladder, comparing both methods at the same active space** so that only the
state-specific/state-averaged difference is in play:

| CAS | loop transport (SS) | gap scan (SA) | SS − SA | SA − FCI | residual | residual/precision |
|---|---|---|---|---|---|---|
| (2,2) | (131.67, 86.30) | 151.61 | −19.94 | +19.00 | 0.197 | **1.57 inconsistent** |
| (4,4) | (128.99, 90.24) | 130.43 | **−1.44** | −2.18 | 0.0088 | 0.29 |
| (6,6) | (129.73, 88.49) | 131.35 | **−1.63** | −1.26 | 0.0574 | 0.76 |

**The state-specific degeneracy sits consistently below the state-averaged one, by ~1.5 deg** at
both rungs where the measurement is self-consistent — the same size as the gap scan's own error
against FCI. This turns `docs/findings.md` §3 from an inference into a measurement.

**The consistency check had to be normalised to mean anything.** Raw residuals of 0.197, 0.0088 and
0.0574 read as CAS(6,6) being six times worse than CAS(4,4). Divided by each run's own precision
(RMS bracket half-width) they are 1.57, 0.29 and 0.76: CAS(6,6) is less *precisely* measured, not
less consistent, and CAS(2,2) is the only rung whose misfit exceeds its own uncertainty. An earlier
revision of this log compared the raw numbers.

CAS(2,2) failing the check is the expected outcome, not a surprise: the gap scan there finds no
minimum in the region at all, and a direct map of the state-specific in-CAS S1/S0 gap over the
region **never falls below 321 mHa**. That map is also exactly symmetric about phi = 90, so any
off-axis degeneracy has a mirror partner and a loop centred on the line would enclose both
(`docs/todo.md` §10).

**Butadiene CAS(2,2) — the experiment the open question needs.** Bisection about the `B_x` centre
brackets the transition at **rho = 0.5385 ± 0.0843** from a verdict sequence monotone in radius,
while the gap scan's intersection for the same active space implies **rho = 0.2121** — a factor 2.5
inside the bracket and decisively outside it. A second centre reports 0 *cleanly* at full size,
excluding 63% of the measured circle and leaving pyr in (105.6, 111.5), tw in (84, 96). A third was
refused at full size (step floor), so **no triangulation and no position is claimed**: two
constraints confine it to an arc. Cost: 50 557 micro-iterations, 2 h, at the cheapest active space.

### Bugs found, all of them ours

1. **Fail-loudly hole.** An adaptive walk that stopped a fifth of the way round reported
   `status OK`, phase `trivial (0)`. The endpoint checks are skipped when there is no endpoint,
   every *accepted* step was continuous, the chain was intact — and `|<Psi_last|Psi_0>|`, the
   single factor carrying the sign of the product estimator, **was never checked at all**. Now two
   checks (`loop_closed`, `closing_step_continuous`), regression-tested. A uniform walk closes by
   construction and pays nothing. Found by looking at a figure, not by reading the code.
2. **The controller pulled the wrong lever.** It rejected a step when CASSCF failed to converge and
   then shrank — but a smaller step makes the warm start *better*, not the solver happier, so the
   walk shrank forever. Measured: steps with a mismatch of 0.011 rejected as unconverged, one probe
   running over ten minutes. Acceptance now depends on continuity alone; convergence is the
   fallback ladder's job and then `analyse`'s.
3. **The fallback ladder had a blind spot.** Every rung relaxed the *gradient* threshold; none
   raised the *iteration budget*. Measured at formaldimine CAS(2,2) (130.86, 91.96): `|grad[o]|` =
   8.8e-06 is inside the 1e-5 threshold, but `dE` = 4.4e-10 will not reach `conv_tol` = 1e-10 in
   200 macro iterations. It converges at 239, to the same energy to 1e-8. New rung
   `warm-mo+warm-ci-long`.
4. **Adaptive event indexing** slipped by one per rejection, once attributing an overlap of 0.79 to
   a chain whose accepted steps were all above 0.90.
5. **`run_gap_minimum_search.py` overwrote its result file** instead of merging, so a CAS(12,12)
   re-run silently destroyed five rungs; they survived only in git. Now merges; rungs restored.

### A claim of ours corrected

"The CAS(12,12) gap search converged in 29 evaluations" was **wrong**. *No* rung's search
converged — every one stopped on its evaluation budget (40 for the five light rungs, 30 for
CAS(12,12)). The comparison is between comparable *efforts*, not converged minima. A
60-evaluation CAS(12,12) run is under way to equalise the budgets.

### Notebooks (§3, §4)

`active_space_study.ipynb` split into `ethylene_ladder.ipynb` (with an exact reference) and
`butadiene_ladder.ipynb` (without one). New: `adaptive_stepping.ipynb`,
`locating_intersections.ipynb`, `summary.ipynb`. The butadiene notebook gains the direct
gap-minimum search, which existed in `results/` and `docs/` but in no notebook.

### One more solver fix, found by the bisection

The initial point of a walk was the only rung of the solver ladder with **no recourse at all**:
continued points escalate through `SOLVE_STRATEGIES` while point 0 got a single attempt at a fixed
budget. Diagnosed on formaldimine CAS(6,6), a loop of radius 1.06 deg: the *only* unconverged point
was point 0, which ran exactly its 200 macro iterations, while every other point converged, the
loop closed and the step floor was never reached — the probe was refused entirely because of the
first solve's budget. With a retry at 3x the budget (convergence criteria untouched) that probe
returns OK/π. It matters most for bisection, which probes small loops by construction, and small
loops mean near-perfect warm starts and exactly the flat-direction crawl that budget was cutting
short. The CAS(6,6) localization above predates this fix and would be tighter with it.

### The CAS(12,12) gap floor, settled

The earlier search used 30 evaluations against the other rungs' 40, leaving the contrast open to
being a budget artifact. Re-run with 60: **1.4595 mHa in 58 evaluations**, against 0.003–0.006 mHa
in 40 for the five lighter rungs (a fall of 47x–163x, against 1.2x here). Doubling the effort moved
the answer by 2.5%. The floor is not an artifact of fitting along a line, of an intersection off the
sampled cross, or of stopping early — and it still does not license "no intersection", because
1.46 mHa is 0.04 eV and no criterion exists for when a state-averaged minimum gap is compatible
with a true crossing.

### Next step

`docs/todo.md` §9: butadiene localization at the rungs where the disagreement is sharpest,
especially CAS(12,12). That is a cluster job — CAS(2,2) alone took two hours. §10 (mirror pairs)
and §6 (a gap threshold) are the two open method questions that bear on interpreting it.

---

## 2026-09-18 — The cluster campaign: ALICE jobs submitted

**Implemented.** `slurm/` became five runnable ALICE job scripts plus a bench job, replacing
the `#### SITE ####` templates for every workload that actually needs a cluster. Site-specific
lines are marked `# SITE:`. The cluster clone is `~/git_repos/BerryCASSCF` on ALICE (Leiden),
built on Python 3.11.5 / PySCF 2.14.0 / NumPy 2.4.6 / SciPy 1.17.1.

**Three things the cluster taught us before any science came back.**

1. **The node is 2.7x slower per point than the laptop**, consistently across active spaces:
   butadiene CAS(2,2) 9.34 s against 3.65 s, CAS(12,12) 349 s (cold) against 127 s. Threading
   does not rescue it. Every runtime estimate in `docs/compute.md` that had been scaled from
   laptop timings was therefore 2.7x optimistic, and the walltimes are now sized from measured
   cluster figures.
2. **Memory, not CPU, is the scarce resource on `cpu_lorentz`.** All three nodes sat at
   253 952 of 254 897 MB allocated with 174 of 384 cores idle, because `DefMemPerCPU` is
   4027 MB and wide jobs reserve hundreds of GB. The first submission asked for 16-32 GB per
   task against a **measured peak of 593 MB** at CAS(12,12), and that over-request blocked this
   account's own flagship job (PENDING, reason `Resources`). Re-submitted at 4 GB, all 18 tasks
   started within seconds. The jobs were 15 minutes old, so the fix cost nothing.
3. **A job that writes only at the end cannot be run for days.** `run_localization.py` saved
   its record after the last centre, so a walltime kill lost everything. It now saves after
   every centre and resumes from the first unfinished one, and — because bisections about
   different centres share no state — it can run **one centre per task** (`--only-centre`) and
   merge afterwards (`--merge`). That turns the CAS(12,12) localization from ~6 days in series
   into ~2 days in parallel. The merge refuses records whose centres do not match the requested
   list, so ellipses from different constructions cannot be combined by accident.

**Submitted** (all on `cpu_lorentz`, 8 CPUs and 4 GB per task):

| job | tasks | what it answers | expected |
|---|---|---|---|
| `berry_cas12.job` | 3 | butadiene CAS(12,12) loop transport on **all three loops** at N = 13, 21, 31 | ~5 h |
| `localize.job` | 12 | bisection ladder: CAS(4,4), (6,6), (8,8), (10,10) x 3 centres | 2-7 h per centre |
| `localize_cas12.job` | 3 | the same at CAS(12,12), one centre per task | ~40-60 h per centre |
| `bench.job` | 1 | one CAS(14,14) solve, to size `rung14.job` from measurement | < 4 h or it TIMEOUTs |

**Why these, in this order.** `docs/findings.md` §3 — the sharpest open question in the project
— rests on the enclosing loop returning pi at CAS(12,12) where the gap scan finds nothing below
1.5 mHa inside it. Two gaps in that evidence are closable cheaply and were closed first: the
rung had **no control loops at all** (so its pi had nothing to be contrasted against) and only
two discretizations (the bare minimum the stability criterion accepts). The localization ladder
then measures, rather than infers, where the state-specific degeneracy sits at each rung; at
CAS(2,2) that measurement already contradicts the gap scan, and one rung is an observation
rather than a result. `rung14.job` (`docs/todo.md` §8) stays last and unsubmitted until the
bench says what a CAS(14,14) point costs.

### When each job returns

* **`berry_cas12.job`.** Read the three verdicts together. If `B_1` and `B_2` return 0 while
  `B_x` returns pi at all three N, the CAS(12,12) pi survives its first real control and
  `docs/findings.md` §3 can drop the "loop transport is wrong here" reading to third place. If a
  control returns pi, that reading becomes the leading one and the section needs rewriting, not
  amending. Either way N=31 gives the stability criterion a third discretization.
* **`localize.job`.** Merge each rung (`--merge`), then compare its measured rho against the
  gap-scan intersection for the *same* rung, exactly as the formaldimine table in §3 does. The
  question is whether the CAS(2,2) offset is a property of that rung or of the method.
* **`localize_cas12.job`.** Same, at the rung where the methods disagree. This is the one that
  can separate the three readings in §3, and it will **not** finish inside a working session —
  plan on two days.
* **`bench.job`.** If a CAS(14,14) point is under ~25 min, `rung14.job` is viable as written;
  if it TIMEOUTs at 4 h, the rung is not reachable this way and §8 should record that rather
  than leaving an unsized job in `slurm/`.

**Unresolved.** Everything the jobs are for. Nothing here changes a scientific conclusion yet.

---

## 2026-09-18 — What the cluster returned, and one correction

Results from the campaign submitted earlier the same day. The two multi-day jobs are still
running; everything below is finished work.

### 1. The CAS(12,12) pi is selective (butadiene, `berry_cas12.job`)

The rung that carries `docs/findings.md` §3 had never been run on a control loop. Both are now
done, at three discretizations each:

| loop | N=13 | N=21 | N=31 | worst adjacent overlap |
|---|---|---|---|---|
| `B_x` encloses | **pi** | **pi** | **pi** | 0.83, 0.88, 0.92 |
| `B_1` control | 0 | 0 | 0 | 0.986, 0.995, 0.997 |
| `B_2` control | 0 | 0 | 0 | 0.986, 0.994, 0.997 |

Every run `OK`, every endpoint estimator ±1.000000 exactly, and the two controls — mirror-image
loops traversed independently — agree to 1e-4. The enclosing loop's worst overlap *improves* with
refinement (0.83 → 0.88 → 0.92) as |Pi| grows (0.49 → 0.63 → 0.72): a converging discretization,
not a marginal one. The stability criterion is now met at this rung with three discretizations
and both controls, where before it had two and none. The enclosing loop is also the *hard* one by the
method's own measure (0.83 against 0.99), which is what passing near a degeneracy looks like.
"The loop-transport result is simply wrong at this rung" now needs a failure selective enough to
spare both controls at three discretizations each.

### 2. The gap floor was never evidence about enclosure (`docs/todo.md` §6, now done)

No threshold was invented. Near a conical intersection the gap is linear in the branching-plane
coordinates, so a floor `g` on a cut of local slope `a` is what a cut passing `g/a` from the apex
would show. `examples/report_gap_criterion.py` reports that conversion per rung. Butadiene
CAS(12,12): floor 1.4595 mHa, slope 1.749 mHa/deg, **miss distance 0.83 deg — 4.6% of the loop's
18 deg semi-axis, and therefore inside the loop.** Every other rung lands at 0.002-0.005 deg. So
the gap scan never said the loop encloses nothing; `docs/findings.md` said it did, and that
framing is withdrawn. The 300x contrast between this rung and the others is real and unexplained,
but it is not about enclosure.

### 3. CAS(14,14) is out of reach, measured (`docs/todo.md` §8, closed)

`bench.job` timed one state-specific solve: **13 123 s — 3.6 h — returning `converged=False`**.
That is 37x the CAS(12,12) cost where determinant counting predicts 13.8x. The butadiene ladder
therefore stops at CAS(12,12) **by measurement, not by choice**, and going higher needs a
different solver (DMRG, selected CI), not a longer walltime. `slurm/rung14.job` is kept, marked
not-viable at the top, as the record.

### 4. Ethylene joins the localization experiment, and immediately disagrees

`docs/findings.md` §4 said the open question — is loop transport right, or consistently wrong? —
needs a system with an exact reference whose ladder still misbehaves, and concluded none of the
three provides one. **Ethylene does**: CAS(12,12) is its full valence space, so the intersection
position is exact in this basis, and its ladder is the badly behaved one.

At CAS(2,2), the bisection about (90, 98) brackets **rho = 0.8216 ± 0.0324**, while the same
rung's gap scan implies 0.7267 and the exact full-valence reference implies 0.7167. Both lie
**outside** the measured bracket. The state-specific object sits about 1.7 deg beyond the
state-averaged one — the same magnitude as formaldimine's 1.44 and 1.63 deg, now measured against
a reference that is exact rather than against another approximation.

### 5. The correction: a refusal is not an exclusion

A first reading of ethylene's other two centres treated their non-bracketing as "the loop
encloses nothing" and built an inconsistency on it. Checking the probe records, **three of the
four non-bracketing loops across CAS(2,2) and CAS(4,4) were refused at the adaptive step floor**,
which says the loop could not be walked, not that it is empty. Only (80, 110.9) at CAS(2,2)
returned a clean 0.

What survives: at CAS(2,2), that one clean 0 plus the ladder's `E_x` pi leaves an allowed region
of `tw` in (92.0, 99.5), `pyr` in (103.0, 113.0) — **off the mirror line**, where every degeneracy
must be paired and every `tw`-symmetric loop must read even. That is a real tension and
`slurm/mirror_test.job` is testing it with two deliberately off-line centres. What does not
survive: the same claim at CAS(4,4), where both other loops were refused, nothing excludes the
line, and the allowed region straddles it.

**The mirror test has since returned, and its answer is a third one I had not enumerated: every
full-size loop in that region is refused.** Three of its four tasks are done and none produced a
verdict, each failing differently and each diagnosed:

| centre | rung | full-size verdict | why |
|---|---|---|---|
| (100, 106) | CAS(2,2) | refused | the adaptive walk never reached `t = 1`; the loop was never closed |
| (80, 106) | CAS(2,2) | refused | endpoint `\|<Psi_0\|Psi_N>\|` = 0.457 — the loop did not return to the same state |
| (100, 106) | CAS(4,4) | refused | CASSCF did not converge at one point |

So the experiment designed to separate "an off-line mirror pair" from "the large loops are not
measuring enclosure" could not be run — and that failure is itself the evidence for the second
option. **Seven of the eight full-size loops tried in this region of ethylene's plane fail their
checks**, in four distinct ways. The one that passes, (80, 110.9) at CAS(2,2), is the sole basis
for the CAS(2,2) tension reported above, and a single passing verdict among siblings that fail
that consistently is not a sound foundation for a parity claim.

**So the tension is withdrawn as a finding and kept as a caution.** What is established is
narrower and more useful: at full loop size — a (12, 18) ellipse reaching `pyr` = 124 and `tw` =
68 — ethylene's plane is not traversable by this continuation. Localization there needs loops
that stay in the region where transport works: centres nearer the object so that `scale = 1` is
already small, or a smaller shape. The bracketing centre (90, 98) works precisely because its
transition radius, 0.82, keeps every probe inside that region.

### 6. A reported adaptive/uniform conflict — withdrawn, it was a misattribution

An earlier version of this entry reported that uniform transport returned pi on butadiene's
`B_x` loop at CAS(8,8) while the adaptive walk returned 0 on the identical loop, and flagged it
as blocking both the ladder result and every bisection bracketed from a full-size loop. **That
was wrong, and the error was mine in reading, not the code's in running.**

The `ZERO` verdict at `scale 1.0000` came from `butadiene_cas8-8_centre1` — the loop about
**(90, 90)**, not the `B_x` centre (90, 101.853). Three concurrent array tasks were appending to
one shared log file, because they started before the per-task log naming was added, and the line
carries no centre. The saved records settle it: the probe with 8117 micro-iterations at scale 1.0
belongs to centre 1, and a loop about (90, 90) reporting 0 at full size is the expected result —
CAS(4,4) does the same.

Checked directly rather than argued, by running the adaptive walk on the actual `B_x` loop
locally at that active space:

```python
loop = Loop("B_x", (90.0, 101.85321091497578), (12.0, 18.0))
traverse_loop_adaptive(loop, cas=CasConfig(basis="6-31g*", ncas=8, nelecas=8), ...)
```

| step control | phase | Π | endpoint | worst overlap | points |
|---|---|---|---|---|---|
| `d_max` 0.10, target 0.020 | **pi** | −0.6583 | −1.0000 | 0.929 | 19 |
| `d_max` 0.07, target 0.008 | **pi** | −0.7581 | −1.0000 | 0.955 | 29 |

Both agree with the uniform runs (pi at N=13 and 21). **There is no conflict**, nothing is
blocked, and the two transports agree wherever they have been compared.

Two things worth keeping from the episode. First, the per-task log naming added earlier today was
not cosmetic: a shared log across concurrent tasks produced a false finding within hours of being
introduced. Second, the correct reflex was the cheap one — butadiene CAS(8,8) is ~5 s per point,
so the "unresolvable until a cluster task lands" claim was itself wrong: the check took four
minutes on the laptop.

A genuine observation does come out of those records, unrelated to the false alarm: the (90, 90)
loop at full size reports **0** at CAS(4,4) and CAS(8,8) but **pi** at CAS(6,6). What that loop
encloses is rung-dependent, which is worth following up when the ladder completes.

### Operational notes worth keeping

* The node is **2.7x slower per point** than the laptop, consistently across active spaces.
* `cpu_lorentz` saturates on **memory**, not CPUs (`DefMemPerCPU` 4027 MB). A 16 GB request
  against a measured 593 MB peak blocked this account's own flagship job; 4 GB everywhere now.
* Threading is worth 2.5x at CAS(12,12) (349 s on 8 threads, 876 s on one), so 8 CPUs is justified
  there and 4 is enough for the small rungs.
* `results/` is tracked in git *and* written by the cluster, so committing a cluster-produced file
  breaks the next `git pull` there. The check-then-drop recipe is in `docs/compute.md`.

---

## Resume here — jobs still running on ALICE as of 2026-09-18 05:13 CEST

Eighteen tasks were still running when this session ended (`berry_cas12.job` finished during
the write-up; its N=31 result is in the table above). They need no attention while they
run; the drivers skip finished work, so **re-submitting any of these jobs is safe** and picks up
only what is missing.

| SLURM id | job | what it is doing | left |
|---|---|---|---|
| `5028447_[0-2]` | `localize_cas12.job` | butadiene CAS(12,12) bisection, one centre per task | ~2 d 20 h |
| `5028404_[0,3-7,9,11]` | `localize.job` | butadiene ladder, the centres still unfinished | ~19 h |
| `5028474_[0-3]` | `localize_ethylene.job` | ethylene CAS(8,8) all centres, CAS(6,6) centre 0 | ~20 h |
| `5029232_[0-2]` | `mirror_test.job` | the off-line mirror centres (one task already finished) | ~7 h |

Check with `squeue -u pollas1`, or
`~/.claude/skills/alice-hpc/scripts/alice-jobs.sh triage <jobid>` for a finished one.

### What to do when they land, in order

1. **Pull and commit.** `rsync -az alice:git_repos/BerryCASSCF/results/ results/`, commit from
   the laptop, then clear the cluster's now-duplicate untracked copies with the recipe in
   `docs/compute.md` §Pulling results back — otherwise the next `git pull` there aborts.
2. **Merge each rung whose three centres are done**:
   `python examples/run_localization.py <system> --cas <ne> <ncas> --merge`. It refuses
   unfinished checkpoints, so a merge that errors means a centre is still running, not that
   something is broken. For the CAS(12,12) run add the centre list:
   `--centre-xy 90.0,101.85321091497578 90.0,90.0 99.0,110.0`.
3. **Rebuild the notebook**: `python notebooks/build_locating_intersections_notebook.py` then
   `jupyter nbconvert --to notebook --execute --inplace notebooks/locating_intersections.ipynb`.
   Its ladder table and constraint-region figures pick up new rungs automatically.
4. **Look at what the (90, 90) loop encloses per rung.** It reports 0 at full size at CAS(4,4)
   and CAS(8,8) but pi at CAS(6,6) (item 6). If that survives the remaining rungs it is a
   rung-dependent statement about enclosure, which is the same phenomenon the ladder measures
   at `B_x` and worth reporting alongside it.
5. **The CAS(12,12) localization is the one to wait for.** When its three centres merge, compare
   the measured rho against the gap scan's 0.83 deg miss distance (item 2). If they agree, the
   two methods are seeing the same object at different resolution and §3 closes; if the measured
   object is somewhere else, the disagreement is real and about position, not enclosure.

### Two things deliberately not done

* **The ethylene third centres were not re-chosen per rung.** `docs/todo.md` §9 says each rung's
  third centre should come from *its own* first-centre bracket; these were all placed from the
  CAS(12,12) reference, and three of four were refused. The fix is one extra task per rung, not
  a re-run, and it should wait until the first centres have reported. The mirror-test outcome
  sharpens what that replacement must satisfy: a centre is only usable if its **full-size** loop
  stays inside the traversable region, which for ethylene means a transition radius well under 1
  — so centres should be placed *closer* to the expected object, not further, and the loop shape
  shrunk with them.
* **`rung14.job` was not submitted** and should not be, at 3.6 h per non-converged point.
