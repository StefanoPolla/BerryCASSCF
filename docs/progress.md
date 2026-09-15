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

**Deliverables.** `notebooks/results.ipynb` (executed, with figures), `docs/results.md`,
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

**Result — the intersection position across the ladder** (reference CAS(12,12): `phi` = 110.01):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| `phi` | 110.02 | 114.85 | 102.09 | 103.11 | 109.10 | 110.01 |
| error | +0.01 | +4.84 | −7.92 | −6.90 | −0.91 | ref |

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
`notebooks/active_space.ipynb`, `docs/active_space.md`, 28 new tests including a regression test
for the path-independence fix.

**Next step.** A system whose S1 is genuinely doubly excited (a polyene 2Ag state, e.g.
butadiene or hexatriene) would test whether a case exists where *no* small active space is
accurate, as opposed to ethylene's accidental-but-unverifiable CAS(2,2). That needs a CI search
in a new coordinate plane and is cluster work.
