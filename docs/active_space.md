# How large an active space does each method need?

The formaldimine benchmark settled on CAS(2,2) for the Berry phase and CAS(4,4) for the
state-averaged comparator (`docs/results.md`). Those are small, and the obvious question is
whether that reflects the methods or just an easy molecule. This study runs the same two
workflows over a ladder of active spaces on a second system.

Narrative version with figures: `notebooks/active_space.ipynb`.
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
only from some directions. `anchor` (now the default in `ScanConfig`) solves once, cold, at the
centre of the region, then gives every grid point the same two path-independent guesses: a cold
start and that one anchor transferred in. Same cost as `best`, without the path dependence.

**Mirror asymmetry is reported for every scan** and is treated as a first-class reliability
criterion: a map that fails the symmetry test is path-dependent and not trustworthy, whatever
number it reports.

## Criteria

**Berry phase** succeeds at an active space if, for some discretization `N`, all three loops
pass every check (all points converged, continuity held, endpoint returned to the same state,
both estimators agreeing in sign) *and* the phases are correct: non-trivial on `E_x`, trivial on
`E_1` and `E_2`. The smallest such `N` is reported, because it turns out to depend on the active
space.

**SA comparator** succeeds if its gap minimum is (a) inside `E_x`, (b) close to the reference
active space's minimum, (c) low enough to read as an intersection, and (d) on a
mirror-symmetric map.

## A free end-to-end check

`E_1` and `E_2` are mirror images of one another, and the two loops are traversed completely
independently — different geometries, independently converged orbitals, independently chosen CI
signs. Nothing in the continuation knows they are related. Every quantity computed on one must
therefore equal its counterpart on the other, and agreement is a genuine end-to-end test of the
whole pipeline on a molecule it was not developed on.

Measured: the largest discrepancy in the product estimator across the ladder is **~1e-06**.

## Results

*(filled in below once the ladder completed — see the tables in
`notebooks/active_space.ipynb`, which are generated from the saved records.)*
