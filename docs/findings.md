# Findings

Consolidated observations across the three systems studied, written to be read on its own.
Numbers are reproducible from the saved records; see `docs/results.md` (formaldimine),
`docs/active_space.md` (ethylene) and `docs/butadiene.md` (butadiene) for the runs behind them,
and the two notebooks for the figures.

## Terminology

Two workflows are compared throughout. On first mention in any document they are named in full;
short forms are used thereafter.

| full name | short form | what it is |
|---|---|---|
| **SS-CASSCF Berry-phase loop transport** | **loop transport** | continuation of a *state-specific* CASSCF ground state around a closed nuclear loop, with the Z2 Berry phase read from the sign of the initial–final nonorthogonal overlap |
| **SA-CASSCF gap scan** | **gap scan** | equal-weight *state-averaged* CASSCF on a grid over the enclosed region, locating the S1–S0 minimum |

"Berry phase" is reserved for the quantity loop transport returns, not for the method itself.

## Systems

| system | basis | plane | exact in-basis reference? |
|---|---|---|---|
| formaldimine H2C=NH | STO-3G | (bend α, dihedral φ) | **yes** — FCI, 1.66 M determinants |
| ethylene C2H4 | 6-31G\* | (twist τ, pyramidalization φ) | **yes** — CAS(12,12) is the full valence space |
| butadiene C4H6 | 6-31G\* | (methylene twist tw, pyramidalization pyr) | **no** — full valence is CAS(22,22) |

---

## 1. The intersection position drifts with active space, non-monotonically

This is the central result, and "it drifts" understates it: on two of three systems, **enlarging
the active space made the answer worse** at intermediate rungs.

Ethylene, error in the intersection position against the CAS(12,12) reference:

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| error (deg) | **+0.02** | +4.85 | −7.91 | −6.90 | −0.90 | ref |

Butadiene, from fine 1-deg cuts refined with the cone model (no reference exists):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| pyr (deg) | 105.86 | 109.68 | 114.46 | 120.87 | 105.07 | 102.17 |
| shift | — | +3.83 | +4.78 | +6.41 | **−15.80** | −2.90 |
| closest approach (mHa) | 0.67 | 0.92 | **3.13** | 0.00 | 0.66 | 2.94 |

Spread **18.7 deg**, and the top two rungs still differ by **2.90 deg** — above the 2-deg
tolerance ethylene met at CAS(10,10), so the sequence has not settled even at the largest
affordable active space.

**CAS(2,2) is included deliberately although it is below the pi space and not a chemically
defensible choice.** It turns out to be well defined here, and it lands 3.7 deg from the largest
rung while CAS(8,8) — three times the orbitals — lands 18.7 deg away. That is the ethylene pattern
again: the minimal space is accidentally the better one, and nothing available at that level would
tell you so.

Formaldimine is the third pattern again: CAS(2,2) is not merely imprecise but *qualitatively*
wrong — it reports a 0.34 mHa near-degeneracy inside a **trivial** control loop and finds nothing
inside the true one. CAS(4,4) is already converged there.

Three systems, three different failure shapes. **No active-space recipe transferred between
them.** What did transfer is the procedure: enlarge until the answer stops moving, and check an
exact symmetry of the system to confirm the solver is not what moved.

### The practical consequence

The heuristic everyone actually uses — bigger active space, more trustworthy — fails here. On
ethylene, comparing CAS(2,2) against CAS(4,4) and trusting the larger would have been wrong: the
minimal space was accurate to 0.02 deg and the larger one to 4.85.

This separates two questions that are usually conflated:

* **accurate** — does this active space land on the right answer?
* **stable** — does *enlarging* it leave the answer unchanged?

Only the second is answerable without already knowing the answer, and for ethylene they differ by
four rungs: accurate from CAS(2,2), stable only from CAS(10,10).

## 2. Loop transport succeeded at the smallest active space tested, on every system

| system | smallest tested | loop transport verdict |
|---|---|---|
| formaldimine | CAS(2,2) | correct (π on the enclosing loop, 0 on both controls) |
| ethylene | CAS(2,2) | correct |
| butadiene | CAS(2,2) | correct at every rung tested |

The cost asymmetry is large: loop transport at the smallest active space is seconds, against hours
for the gap scan ladder.

## 3. The gap-scan minimum and the point loop transport encircles are different objects

At butadiene CAS(8,8) the gap-scan minimum sits at pyr = 121.2, **outside** the enclosing loop
(centre 101.85, radius 18 deg). The natural prediction was that this rung would return a trivial
phase. **It returns π**, as do all other rungs.

The prediction rested on equating the gap-scan minimum with the degeneracy loop transport
encircles. That equation is false. Loop transport follows the *state-specific* ground state,
whose degeneracy lies wherever the state-specific surfaces touch; the gap-scan minimum is computed
with *state-averaged* orbitals, a different set. For a converged active space the two coincide;
for a truncated one they need not.

Shrinking the loop bounds the degeneracy loop transport actually encircles:

| radius (pyr) | reach | N | min adjacent overlap | verdict |
|---|---|---|---|---|
| 18 | 119.9 | 15 | 0.89 | OK, π |
| 12 | 113.9 | 31 | 0.87 | OK, π |
| 8 | 109.9 | 31 | 0.77 | refused (continuity) |
| 8 | 109.9 | 15 | 0.24 | refused (continuity) |

A 12 deg loop reaching only to pyr = 113.9 still returns a clean π, so **what carries the phase
lies below 113.9, not at 121.2**.

## 4. The inversion: the cheap method is the stable one

Across butadiene's ladder, loop transport returns the same answer at every rung while the gap scan
— the method nominally validating it — scatters over 19 deg and never settles, at three orders of
magnitude more cost.

This does **not** make loop transport right by default. With no exact reference for butadiene,
"π at every rung" is equally consistent with "very robust" and "consistently wrong in the same
direction". What can be said is that it is self-consistent, passes every internal check, and is
cheap. Settling it needs a system where an exact reference exists *and* the ladder still
misbehaves — which none of the three provides.

## 5. Failure behaviour: what the checks catch, what they cost, and what they miss

Probed directly at butadiene CAS(8,8) by shrinking the loop toward the degeneracy and varying the
discretization independently (`examples/run_failure_modes.py`). Every run records what it *would*
have reported had the checks not run, so refusals of correct answers count as a cost rather than
a save.

| radius (reach in pyr) | N=15 | N=31 | N=61 |
|---|---|---|---|
| 6 (107.8) | refused, would say **0** | reported π (ovl 0.844) | refused, would say **0** |
| 8 (109.8) | refused, would say π | refused, would say π | **reported π** (ovl 0.821) |
| 12 (113.8) | refused, would say π | reported π (ovl 0.874) | — |
| 18 (119.8) | reported π (ovl 0.801) | reported π (ovl 0.948) | — |
| displaced control | reported 0 (ovl 0.946) | reported 0 (ovl 0.987) | — |

**It catches wrong answers, not just imprecise ones.** At radius 6, N=15 the product estimator was
+0.328 — a confident-looking *trivial* verdict, the opposite of what every well-resolved loop
says. The continuity check refused it (min overlap 0.473).

**It distinguishes "too coarse" from "too close".** Radius 8 is refused at N=15 and N=31 but its
sign is stable (−0.18, −0.52, −0.66) and it passes at N=61: a pure discretization problem,
resolved by refinement. Radius 6 is different — its sign *flips* (+0.33, −0.55, +0.57) and it
never passes twice. Sign stability under refinement is what separates the two.

**The conservatism has a measured price.** Five of twelve runs were refused; three of those would
have given the same sign as the well-resolved ones. The check errs toward refusal.

**It refuses lost continuity, not mere proximity.** The displaced control encloses nothing but
passes near the seam, and is accepted with overlaps of 0.95–0.99, correctly reporting a trivial
phase. Without this control the refusals above could have been the check rejecting anything near
a seam.

### A single passing run is not enough — and this is where a per-run check fails

Radius 6 at N=31 **passed** every per-run check (overlap 0.844, endpoint 0.984) and reported π.
Its neighbours at N=15 and N=61 both say 0. Since radius 8 reports π reliably and radius 6 does
not, the degeneracy most likely sits *between* their reaches — pyr between 107.8 and 109.8 —
which would make **0 the correct answer at radius 6 and the passing run the wrong one**.

So the per-run checks are necessary but not sufficient: one of them accepted an answer that is
probably wrong. What rejects it is the **stability criterion**, which demands at least two
discretizations *all* passing with the same phase. Radius 6 has exactly one passing N, so no phase
is claimed for it. The layering is what makes the protocol safe, not any individual check.

## 6. The required loop discretization grows with the active space

Ethylene, minimum adjacent overlap on the enclosing loop at fixed N = 13:

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| min overlap | 0.895 | 0.718 | 0.689 | 0.673 | 0.692 | 0.686 |

CAS(4,4) and above fall below the 0.80 threshold and are refused at N = 13; all pass at N = 21. A
larger active space has more freedom to rearrange, so the wavefunction turns faster around the
same loop. This is a statement about *loop resolution*, not about whether the active space can
describe the physics, and the two must not be conflated.

## 7. Methodological results that are not about active spaces

**Locating an intersection: fit the square of the gap.** Near a conical intersection the gap is
linear, so a cut through the region is a V and a parabola fitted to the *gap* is biased toward the
grid minimum. For an ideal cone, a straight cut at perpendicular offset satisfies
`gap² = a²(x − x0)² + b²` — exactly a parabola in `gap²`. Fitting that recovers the apex without
bias and additionally returns the cut's closest approach, which the grid minimum only bounds from
above. On a synthetic cone with a 5 deg grid the cone fit is exact where the parabolic fit is off
by 0.16 deg; on the real scans the two differ by 0.5–2.3 deg. See `berrycasscf/refine.py`.

**A warm-started gap scan is path-dependent.** Warm starting each grid point from its neighbour
lets the active space drift along the scan path, so the answer depends on the route taken to a
geometry. Caught by an exact symmetry: ethylene geometries at τ and 180 − τ are isometric, so the
map must be symmetric. Warm starting gave **19.1 and 29.4 mHa at mirror-image points** (cold gives
26.9). All scans now use cold starts, which respect molecular symmetry automatically and are
path-independent by construction; the measured price is that cold occasionally converges to a
worse solution (one formaldimine grid corner, 0.115 Ha high). Details and the alternatives
considered are in `docs/active_space.md`.

**Symmetry checks do not transfer between systems.** Ethylene's τ → 180 − τ check does not apply
to butadiene, whose two methylene hydrogens are inequivalent; its exact symmetry is reflection
through the molecular plane, `(tw, pyr) → (−tw, −pyr)`. Reusing the ethylene check would have
compared unrelated geometries. Both facts are locked in by tests.

---

## Open questions

1. **Is loop transport right on butadiene, or consistently wrong?** Unresolvable without an exact
   reference. The strongest available evidence is self-consistency across rungs plus agreement
   with the two better-converged gap-scan rungs.
2. **Where exactly does the state-specific degeneracy sit, per rung?** At CAS(8,8) it is bounded
   to roughly pyr 107.8–109.8 by the radius sweep above — inside the radius-8 loop, outside or
   marginal at radius 6. That is already ≥11 deg from the gap-scan intersection at 120.87.
   Bisecting per rung, with enough discretization to keep continuity, would turn the bound into a
   measurement for every rung and let the two methods' estimates be plotted against each other.
   An earlier reading of this sweep took a single passing run at radius 6 at face value and
   inferred a tighter bound; refining N contradicted it. The bound above rests only on runs whose
   sign is stable across discretizations.
3. **Solution discontinuities in the gap maps.** The butadiene CAS(8,8) cut jumps 2.84 → 27.35 mHa
   between pyr 125 and 130, and the CAS(12,12) cut drops 28.72 → 7.04 between 130 and 135. These
   are branch changes, not cone structure, and they sit near the region the control loops occupy.
4. **Does the intersection really sit at tw = 90 for every rung?** Assumed from a five-point grid
   with no symmetry argument for butadiene; being tested directly.
