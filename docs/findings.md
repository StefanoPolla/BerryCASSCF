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

Butadiene, shift from the previous rung (no reference exists):

| CAS | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|
| pyr (deg) | 109.83 | 114.81 | 121.05 | 105.04 | 101.85 |
| shift | — | +4.98 | +6.24 | **−16.01** | −3.19 |

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
| butadiene | CAS(4,4) → CAS(2,2) *(in progress)* | correct at every rung tested |

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

## 5. Failures are detected, not silent — and the check is conservative

When a loop passes too close to a degeneracy the adjacent overlaps collapse and the continuity
check **refuses to report a phase**. The clearest case: at CAS(8,8), radius 8, N=15, the product
estimator was **−0.176**. Reported naively that reads as a weak non-trivial phase — a plausible
wrong answer. The check rejected it (min overlap 0.24).

The cost of that conservatism is real and is reported alongside: the same loop at N=31 reaches
min overlap 0.77 and is *still* refused, though its sign agrees with the better-resolved runs. The
check errs toward refusal, which means it will sometimes reject answers that were correct.

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
2. **Where exactly does the state-specific degeneracy sit, per rung?** Currently only bounded
   (< 113.9 deg at CAS(8,8)). Bisecting on loop radius per rung would turn the bound into a
   measurement and let the two methods' intersection estimates be plotted against each other.
3. **Solution discontinuities in the gap maps.** The butadiene CAS(8,8) cut jumps 2.84 → 27.35 mHa
   between pyr 125 and 130, and the CAS(12,12) cut drops 28.72 → 7.04 between 130 and 135. These
   are branch changes, not cone structure, and they sit near the region the control loops occupy.
4. **Does the intersection really sit at tw = 90 for every rung?** Assumed from a five-point grid
   with no symmetry argument for butadiene; being tested directly.
