# How large an active space does each method need?

The formaldimine benchmark settled on CAS(2,2) for the Berry phase and CAS(4,4) for the
state-averaged comparator (`docs/results.md`). Those are small, and the obvious question is
whether that reflects the methods or just an easy molecule. This study runs the same two
workflows over a ladder of active spaces on a second system.

Narrative version with figures: `notebooks/active_space_study.ipynb`.
Reproduce with `python examples/run_active_space_study.py scan` then `... berry`.

## System and design

**Ethylene at its twisted-pyramidalized S1/S0 conical intersection** — the textbook case in
nonadiabatic photochemistry. Two rigid coordinates, both pure rotations, so every bond length
and both HCH angles are preserved exactly:

| | |
|---|---|
| `tau` | torsion of one CH2 about the C=C axis (deg) |
| `phi` | umbrella pyramidalization of that same CH2 (deg) |

**Basis: 6-31G\*.** STO-3G was tried first and has no intersection along these coordinates at
all — the gap bottoms out around 10 mHa. The twisted-ethylene S1 is the zwitterionic V state,
which a minimal basis describes badly. Here the *basis*, not the active space, is the binding
constraint, which is worth separating from the question being asked.

**Ladder:** CAS(2,2), (4,4), (6,6), (8,8), (10,10), (12,12). CAS(12,12) spans the full valence
space — only the two carbon 1s orbitals stay frozen — and is the in-basis reference.

**The loops are fixed across the whole ladder**, centred on the intersection located by the
reference active space, with radius 12° in both coordinates. `E_x` encloses it; `E_1` and `E_2`
are displaced by ±30° in `tau` and enclose nothing. Letting the loops track each active space's
own CI estimate would mean every rung was answering a different question.

## Two methodological findings, before any result

### 1. Active-space drift makes a warm-started scan path-dependent

The gap scan warm-starts each grid point from its neighbour, for speed. For small active spaces
that is harmless; for larger ones the active space itself drifts along the scan path, so the
answer depends on the route taken to a geometry.

Ethylene gives a clean way to catch this: geometries at `tau` and `180 - tau` are exact mirror
images (verified by comparing full interatomic-distance spectra), so **a correct gap map must be
symmetric about `tau = 90`**. Measured for CAS(8,8)/6-31G\*, naive warm starting returned
**19.1 mHa** and **29.4 mHa** at two geometries that are mirror images; a cold start at either
gives **26.9 mHa**.

### 2. The obvious fixes trade determinism against solution quality

For CAS(4,4)/6-31G\* on an 11×11 grid:

| strategy | mirror asymmetry | mean SA energy | comment |
|---|---|---|---|
| `warm` | 9.8 mHa | −77.79478813 | finds the *lowest* energy, but path-dependent |
| `cold` | 1.1e-07 mHa | −77.79438077 | deterministic, but misses the better solution |
| `best` (warm+cold, lower wins) | 7.9 mHa | — | inherits warm's path dependence |
| **`anchor`** | **4.2e-04 mHa** | **−77.79468688** | deterministic *and* finds most of the improvement |

So SA-CASSCF here has **more than one stationary point**, and a warm sweep reaches the lower one
only from some directions. `anchor` solves once, cold, at the
centre of the region, then gives every grid point the same two path-independent guesses: a cold
start and that one anchor transferred in.

**The ladder uses `cold`, which is also the `ScanConfig` default.** Both are path-independent, but this region's
anchor sits at the intersection, and transferring those orbitals outward is a poor guess:
CAS(6,6) ran at **8.6 s/point** with `anchor` against **0.65 s** cold, which would have put the
full ladder beyond the compute budget. Cold's cost is that it can settle on a slightly higher
SA-CASSCF solution than the best reachable — for CAS(4,4), −77.79438 against −77.79469. That
difference does not reach the observable: on the same grid the two strategies give an identical
intersection position *and* an identical minimum gap, **1.362 mHa at (90, 114)** for both. Using
one strategy for every rung also keeps the rungs directly comparable.

**Mirror asymmetry is reported for every scan** and is treated as a first-class reliability
criterion: a map that fails the symmetry test is path-dependent and not trustworthy, whatever
number it reports.

### What this does and does not fix

None of this *solves* the multiple-minimum problem — nothing here searches for the global
SA-CASSCF solution. What it does is make the ambiguity **detectable** (the symmetry check) and
**reproducible** (no guess is chained along the scan path). Those are weaker claims than
correctness, and the difference is measurable: on the `tau = 90` row at CAS(4,4) the three
strategies reach mean state-averaged energies spanning 1.3e-03 Ha, with `warm` lowest. Cold is
provably *not* always the best solution available; it is merely always the same one.

The reason the conclusions survive is that the ambiguity does not reach the observable. On the
`tau = 90` line — where the intersection sits and from which `phi` is measured:

| CAS | strategy | refined phi | min gap | mean SA energy |
|---|---|---|---|---|
| (4,4) | cold / anchor / warm | 114.85 / 114.85 / 114.85 | 1.362 / 1.362 / 1.362 mHa | −77.79058 / −77.79101 / −77.79191 |
| (6,6) | cold / anchor / warm | 102.09 / 102.09 / 102.09 | 0.155 / 0.155 / 0.156 mHa | identical |
| (8,8) | cold / anchor / warm | 103.11 / 103.11 / 103.10 | 2.033 / 2.033 / 2.032 mHa | identical |

All three strategies agree on `phi` to within 0.01 deg and on the gap to within 0.001 mHa, at
every rung tested — at CAS(4,4) even while disagreeing about the energy. At CAS(6,6) and (8,8)
they converge to the *same* solution on this line: there the ambiguity lives off-axis, which is
precisely where the 2D map was corrupted and where the reported observable does not live.

**Remaining limitations.** The symmetry test is necessary, not sufficient: a solver that lands on
the same *wrong* branch at both mirror points passes it. It catches inconsistency, not consistent
error.

This is sharper than it first appears. A cold start should *automatically* respect the symmetry:
mirror-image geometries give integral matrices related by a signed permutation, PySCF's `minao`
initial guess is built from atomic densities and transforms the same way, and everything after
that is deterministic — so the two solves are the same calculation in different coordinates.
The measured cold asymmetries (1e-11 to 3e-05 mHa, i.e. at or just above the 1e-09 Ha
convergence threshold) are consistent with exact equivariance limited only by convergence, not
with any real asymmetry.

So **for a cold-started scan the mirror test is close to vacuous** — it is guaranteed to pass by
construction and mostly re-confirms determinism. Its diagnostic power was against the *warm*
sweep, where chaining guesses along the path breaks the equivariance. Treating
"mirror-symmetric" as evidence of solution *quality* would be a mistake: cold is symmetric and
still lands 1.3e-03 Ha above the best solution found at CAS(4,4). The places it could genuinely
fail are where the argument's premises break — a near-degenerate orbital ordering flipping which
orbitals enter the active space, or a symmetry-broken SCF solution. And it is a symmetry ethylene happens to possess; formaldimine's grid is centred at
`phi = 89.9` and so does not sample exact mirror pairs, which is why that system was checked a
different way (cold against warm at identical geometries, agreeing to 2.5e-06 Ha).

### The price of cold, measured

When the formaldimine scans were redone cold and compared against the committed warm-started
ones, every minimum position and every inside/outside conclusion was unchanged. But the grids
differ elsewhere, and **in every case it is warm that found the lower energy**:

| scan | points differing > 1 mHa | worst | where | lower solution |
|---|---|---|---|---|
| `C_2` CAS(2,2) | 1 / 625 | 279 mHa | (161.0, 89.9), grid corner | warm, by **0.115 Ha** |
| `C_1` CAS(6,6) | 17 / 625 | 3.9 mHa | (113.7, 97.2) | warm |
| `C_x` CAS(6,6) | 24 / 625 | 3.7 mHa | (127.2, 89.0) | warm |

That single corner point is a genuine cold-start failure: a fresh RHF guess there converges to a
solution 0.115 Ha above the one a warm start reaches. It sits far from the intersection and
changes nothing, but it is the concrete cost of choosing reproducibility over energy, and on a
different system such a point could land somewhere that matters.

The targeted fix is cheap and was not implemented: flag grid points whose state-averaged energy
is an outlier against their neighbours, and retry only those from a neighbour's orbitals. That
keeps the bulk of the grid deterministic while repairing isolated failures, and it is a better
use of effort than multi-starting every point.

**What would actually solve it** is a deterministic multi-start: some fixed number of seeded
orbital perturbations at each geometry, keeping the lowest state-averaged energy. That is
path-independent, embarrassingly parallel, and costs a constant factor. It is the obvious next
step and was not done here.

### The formaldimine results are unaffected

Since this defect was found after `docs/results.md` was written, the formaldimine headline scan
was recomputed with `anchor` and compared against the committed warm-started one:

| | min gap | at (alpha, phi) | mean SA energy |
|---|---|---|---|
| committed (warm) | 0.000882 Ha | (130.000, 89.900) | −92.74471004 |
| `anchor` | 0.000881 Ha | (130.000, 89.900) | −92.74471004 |

Maximum gap difference anywhere on the 25x25 grid: **2.5e-06 Ha**. Same intersection position,
same energies. Drift needs a large enough active space to have somewhere to drift *to*;
CAS(4,4) on formaldimine does not. All conclusions in `docs/results.md` stand.

## Criteria

**Berry phase** succeeds at an active space if, for some discretization `N`, all three loops
pass every check (all points converged, continuity held, endpoint returned to the same state,
both estimators agreeing in sign) *and* the phases are correct: non-trivial on `E_x`, trivial on
`E_1` and `E_2`. The smallest such `N` is reported, because it turns out to depend on the active
space.

**SA comparator** succeeds if its gap minimum is (a) inside `E_x`, (b) close to the reference
active space's minimum, (c) low enough to read as an intersection, and (d) on a
mirror-symmetric map.

## The control loops really do enclose nothing

The scan region covers only `tau` in [70, 110], so the controls at `tau = 60` and `120` sit
outside it. Sampling the SA-CASSCF(2,2) gap on and inside each control loop (centre, plus two
rings of 8 points at radius 6 deg and 12 deg) gives a minimum of **81.8 mHa** for both, against
**0.177 mHa** at the intersection: a factor of 460. Nothing is enclosed.

Both controls return 81.80 mHa at mirror-image points, which is the symmetry check again.

## A free end-to-end check

`E_1` and `E_2` are mirror images of one another, and the two loops are traversed completely
independently — different geometries, independently converged orbitals, independently chosen CI
signs. Nothing in the continuation knows they are related. Every quantity computed on one must
therefore equal its counterpart on the other, and agreement is a genuine end-to-end test of the
whole pipeline on a molecule it was not developed on.

Measured: the largest discrepancy in the product estimator across the ladder is **~1e-06**.

## Results

Live tables are generated from the saved records in `notebooks/active_space_study.ipynb`; the numbers
below are the summary.

### Where each active space puts the intersection

The reference is CAS(12,12) (full valence), which along `tau = 90` gives a clean V with its
minimum at **phi = 110.00**, gap 0.189 mHa. Sub-grid positions come from parabolic
interpolation through the grid minimum and its neighbours (grid step 4 deg).

| active space | own minimum gap | phi (refined) | error vs reference | mirror asymmetry |
|---|---|---|---|---|
| CAS(2,2) | 0.175 mHa | 111.08 | **+0.18 deg** | 2.8e-11 |
| CAS(4,4) | 1.362 mHa | 116.41 | **+5.50 deg** | 1.1e-07 |
| CAS(6,6) | 0.155 mHa | 102.64 | **−8.26 deg** | 2.8e-05 |
| CAS(8,8) | 2.033 mHa | 103.90 | **−7.00 deg** | 1.1e-06 |
| CAS(10,10) | 1.850 mHa | 110.93 | +0.03 deg | 3.5e-09 |
| **CAS(12,12)** | **0.189 mHa** | **110.90** | reference | 4.3e-11 |

Positions come from the cone model (a parabola fitted to `gap`&sup2;, exact for any cut through a
cone). An earlier version fitted a parabola to the gap itself, which is a V near an intersection;
that shifted every position by up to ~1.8 deg, the reference included, without changing the
verdict.

Every scan passes the mirror-symmetry test (asymmetries 1e-11 to 1e-5 mHa), so these
differences are the physics of the truncation, not numerical noise.

**Convergence with active space is not monotonic.** The error runs
+0.18, +5.50, −8.26, −7.00, +0.03, 0 degrees. The minimal pi space is essentially exact;
CAS(4,4), (6,6) and (8,8) are displaced by 5–8 degrees *in both directions*; the reference is
recovered only from CAS(10,10). An intermediate active space here is markedly **worse** than the
smallest one.

A sharper way to see the same thing is to ask what each active space believes the gap to be
**at the position the reference puts the intersection**:

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| gap at reference CI | 0.175 | 6.758 | 11.305 | 9.250 | 1.850 | 0.189 mHa |

CAS(6,6) reports 11.3 mHa where the gap is actually ~0.2: it does not see an intersection there
at all. CAS(2,2) reports 0.175 mHa, essentially the exact value.

The likely reason is balance rather than size. CAS(2,2) contains exactly the pi/pi\* pair that
the twisted-ethylene S1/S0 pair is built from, so both states are described equally. Adding
orbitals a few at a time brings in correlation that lowers one state more than the other until
the space is large enough to treat them evenly again.

### Berry phase

Every active space tested gives the correct topology — non-trivial (pi) on `E_x`, trivial on
both controls — provided the loop is discretized finely enough. What changes with the active
space is not the answer but how hard it is to obtain:

| active space | min adjacent overlap on `E_x`, N=13 | smallest working N |
|---|---|---|
| CAS(2,2) | 0.895 | 13 |
| CAS(4,4) | 0.718 | 21 |
| CAS(6,6) | 0.689 | 21 |
| CAS(8,8) | 0.673 | 21 |
| CAS(10,10) | 0.692 | 21 |
| CAS(12,12) | 0.686 | 21 |

A larger active space has more freedom to rearrange, so the wavefunction turns faster around the
same loop and the same `N` yields smaller adjacent overlaps. Below the 0.80 continuity
threshold the run is reported FAILED — correctly, since continuity has not been established —
and refining `N` fixes it. This is a statement about **loop resolution**, not about whether the
active space can describe the physics, and the two must not be conflated.

One further sign of the same effect: at CAS(10,10) the endpoint overlap comes back as
**−0.9860** rather than −1.0000. The closing solve, at a geometry identical to the start, lands
on a very slightly different solution. It still passes the 0.90 threshold, but the margin
shrinks as the active space grows.

## Answering the question

**Is there a system for which larger active spaces are needed?** For the state-averaged
resolver, **yes — ethylene**, and the way it fails is instructive.

* **Smallest active space that gets the right answer: CAS(2,2)** (error 0.18 deg).
* **Smallest active space you could *trust* without already knowing the answer: CAS(10,10)–(12,12)**,
  because only there does enlarging the space stop changing the result.

Those are different questions and for ethylene they differ by four rungs. CAS(2,2) is right, but
nothing available at CAS(2,2) tells you so: the very next rung moves the intersection by 5 deg,
and the one after that by 8 deg in the other direction. The usual practical test — "enlarge the
active space and see whether the answer moves" — would reject CAS(2,2), and would keep rejecting
until CAS(10,10). So an honest protocol needs a much larger space than CAS(4,4) here, even
though a lucky small one exists.

Compare formaldimine, where CAS(2,2) was *qualitatively* wrong (it invented an intersection
inside a trivial loop) and CAS(4,4) was already converged. The two systems fail in opposite ways,
which is the argument for checking convergence per system rather than reusing an active-space
recipe.

**For the Berry phase: no.** CAS(2,2) suffices on both systems tested, and the reason it is
insensitive to the errors that defeat the comparator can be made quantitative. `E_x` is centred
at `phi` = 110 with radius 12 deg, and the largest displacement anywhere on the ladder is
CAS(6,6)'s 7.92 deg. **Every active space tested still places the intersection inside the loop**,
so every one of them returns the same, correct, topological answer — even the rungs whose
geometry is badly wrong. The topological question is
far more forgiving than the geometric one — the same conclusion the formaldimine FCI calibration
reached, now with a second, independent illustration.

The cost of that robustness is discretization, not active space: larger spaces need finer loops
(N = 21 rather than 13 here), which is cheap and, more importantly, *detected* — the continuity
check fails loudly rather than returning a wrong phase.
