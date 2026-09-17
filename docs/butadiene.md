# Butadiene: a third system, and the search that chose its plane

`docs/active_space.md` found that ethylene's state-averaged comparator is accurate at CAS(2,2)
but not *stable* until CAS(10,10). Butadiene is the follow-up, chosen because it should fail for
a different and more fundamental reason.

Narrative with figures: `notebooks/active_space_study.ipynb` (final section).
Reproduce: `python examples/search_butadiene_ci.py` then `python examples/run_butadiene_study.py scan` / `... berry`.

## Why butadiene

Its low-lying excited states are the standard hard case for active-space truncation. The
1¹B<sub>u</sub> state is a singly-excited π→π\*, but **2¹A<sub>g</sub> carries a large doubly-excited
component**, and the ordering and relative position of the two is notoriously sensitive to both
the active space and dynamic correlation. That is the situation in which a small CAS should fail
outright — in contrast to ethylene, where CAS(2,2) turned out to be accidentally excellent.

There is one structural difference that changes what can be claimed. Formaldimine's full space is
reachable by FCI (13 orbitals) and ethylene's full valence space is CAS(12,12); butadiene's is
**CAS(22,22)**, far out of reach. So **there is no exact in-basis reference here**. The largest
affordable rung is the best available, and "converged" can only mean "the answer has stopped
moving" — never "the answer is exact". If it is still moving at the top of the ladder, that is
itself the finding.

## Geometry and coordinates

RHF/6-31G\* optimized planar s-trans reference (C1=C2 1.3229 Å, C2–C3 1.4674 Å, C3=C4 1.3229 Å —
the expected bond alternation). Four rigid distortions, all pure rotations, verified to preserve
every bond length to 1e-10:

| | |
|---|---|
| `tw` | twist of the C1 methylene about the C1=C2 axis (deg) |
| `pyr` | umbrella pyramidalization of that same methylene (deg) |
| `tc` | torsion about the central C2–C3 bond: s-trans (0) to s-cis (180) |
| `bend` | in-plane C1–C2–C3 angle change (deg) |

## The search: three planes, only one contains an intersection

Nothing downstream assumes where the intersection is. Three candidate planes were scanned at
SA-CASSCF(4,4)/6-31G\* on 13×13 grids (`examples/search_butadiene_ci.py`):

| plane | minimum gap | at | verdict |
|---|---|---|---|
| `tw_pyr` | **0.73 mHa** | tw=90, pyr≈110 | **contains a conical intersection** |
| `tw_tc` | 99.49 mHa | tw=90, tc=180 | no intersection |
| `tw_bend` | 36.21 mHa | tw=90, bend=−40 (box edge) | no intersection |

So the twist alone does not close the gap, and neither central torsion nor skeletal bending helps;
**pyramidalization of the twisted methylene is what brings the states together** — the same motif
as ethylene, now carrying a vinyl substituent.

### A trap the search caught

The unrestricted minimum of the `tw_pyr` plane is **not** the usable one. At `pyr = 180` the gap
is 0.19 mHa, but the methylene has folded back onto its own C1–C2 bond (the HCH bisector lies
0.1° from the C1→C2 axis) at **266 kcal/mol** of strain: two states degenerate at a chemically
meaningless structure. The search summary therefore reports **strain alongside the gap**, and the
usable candidate is the strain-filtered one at `pyr ≈ 110`, 161 kcal/mol.

Refining along `tw = 90` in 5° steps gives the CAS(4,4) intersection at **pyr = 109.83**, gap
0.729 mHa, with a clean linear cone on both sides (19.8, 12.6, 5.7, **0.7**, 6.4, 11.2 mHa).

## Symmetry: butadiene is not ethylene

Ethylene's gap maps were validated against `tau → 180 − tau`, which is exact there because its two
methylene groups are equivalent. **Butadiene has no such symmetry.** Its two C1 hydrogens are
inequivalent — one cis and one trans to the C3=C4 unit — so `tw` and `180 − tw` are genuinely
different geometries, differing by ~3e-03 Å in their interatomic-distance spectra.

The exact symmetry is reflection through the molecular plane,

    (tw, pyr) -> (-tw, -pyr)      at tc = 0

which is isometric to numerical precision, and which maps *outside* the scanned region — so the
path-independence check is a handful of extra solves rather than a property of the grid
(`symmetry_spot_check` in the study driver). With `tc != 0` even that symmetry is gone.

This is why the ladder grid samples `tw` with five rows rather than assuming the intersection sits
at `tw = 90`: nothing pins it there.

Reusing ethylene's check here would have silently compared unrelated geometries. Both facts are
locked in by tests (`test_reflection_through_the_molecular_plane_is_an_exact_symmetry` and
`test_ethylene_style_mirror_is_NOT_a_symmetry_here`).

## Results: the intersection position across the ladder

Grid 5 x 13 over `tw` in [70, 110] and `pyr` in [80, 140]; sub-grid positions by parabolic
interpolation along the row holding the two-dimensional minimum. Every rung passes the
plane-reflection check at **exactly 0.0 mHa**.

| active space | `pyr` (refined) | grid minimum | shift from previous rung |
|---|---|---|---|
| CAS(4,4) | 109.83 | 0.724 mHa | — |
| CAS(6,6) | 114.81 | 6.162 mHa | +4.98 |
| CAS(8,8) | 121.05 | 0.906 mHa | +6.24 |
| CAS(10,10) | 105.04 | 0.255 mHa | **−16.01** |
| CAS(12,12) | 101.85 | 4.555 mHa | −3.19 |

**Total spread: 19.2 degrees. The top two rungs still differ by 3.2 degrees** — above the 2-degree
tolerance used for ethylene, so the sequence has not settled even at the largest affordable
active space. Every rung sits at `tw = 90` and passes the plane-reflection check at exactly
0.0 mHa, so this is the physics of truncation, not solver noise.

CAS(12,12) is run on the single `tw = 90` row rather than the full grid: a cold solve there
measured at **over 2.4 minutes per point** (853k determinants on 68 basis functions), so a
5 x 13 grid would have taken ~2.6 h for one rung. Every other rung's two-dimensional minimum
lies on `tw = 90`, and the refined `pyr` position is the only quantity compared across the
ladder, so the row carries the comparison. The cost is that this rung's `tw` is assumed rather
than resolved.

The Berry stage stops at CAS(10,10) for the same reason: a state-specific solve at CAS(12,12)
costs minutes, so three loops at two discretizations would run to several hours without changing
a conclusion already established over four rungs.

Each rung produces a clean, well-formed cone — these are genuine intersections, and it is their
*position* that is active-space dependent. The sequence climbs by ~5 degrees a rung and then
falls back by 16, spanning **16 degrees** in total with no sign of settling.

This is a qualitatively worse situation than either earlier system. Formaldimine converged by
CAS(4,4); ethylene was unstable until CAS(10,10) but had a full-valence reference to prove it.
Butadiene is **still moving at the top of the ladder, with no exact reference available to say
which rung (if any) is right**. The honest statement is that the state-averaged comparator is
*not converged anywhere on the affordable ladder for this system* — which is the clearest
possible answer to "is there a system needing more than CAS(4,4)?": here even CAS(12,12) is not
demonstrably enough.

### Direct search: the fitted positions hold, the fitted gaps do not

The positions above come from cone fits along a cut, which assume the cut passes through the apex.
`examples/run_gap_minimum_search.py` tests that by seeding a derivative-free 2D minimization of the
gap from each fitted position (cold SA-CASSCF evaluations, cached, hard budget):

| CAS | fitted (tw, pyr) | fitted gap | searched (tw, pyr) | searched gap | moved | evals |
|---|---|---|---|---|---|---|
| (2,2) | (90.00, 105.86) | 0.276 mHa | (89.97, 105.67) | **0.003** | 0.19 | 40 |
| (4,4) | (90.00, 109.68) | 0.341 | (89.98, 109.41) | **0.005** | 0.27 | 40 |
| (6,6) | (90.00, 114.46) | 0.990 | (**89.11**, 114.61) | **0.006** | **0.91** | 40 |
| (8,8) | (90.00, 120.87) | 0.201 | (89.98, 121.08) | **0.003** | 0.22 | 40 |
| (10,10) | (90.00, 105.07) | 0.303 | (89.96, 104.90) | **0.006** | 0.18 | 40 |
| (12,12) | (90.00, 102.17) | 1.776 | (90.00, 101.96) | **1.496** | 0.21 | 29 |

The positions move by at most 0.91 degrees, so the ladder's 19.2-degree spread is not a fitting
artifact. The *gaps* are another matter: five of six rungs fall by about two orders of magnitude,
confirming that those cuts pass essentially through a real intersection and that a fitted "closest
approach" measured along a line means very little.

**CAS(12,12) does not follow.** The same search improved the gap by only 16% in 29
evaluations. **No rung's search converged** -- every one stopped on its evaluation budget, 40 for
the five lighter rungs and 30 here, because a single CAS(12,12) evaluation costs ~2.2 minutes. The
comparison is therefore between comparable efforts rather than between converged minima. Two
follow-ups were run to see whether the difference is an artifact:

* **Was the minimum off the sampled cross?** A 5 x 5 grid over the loop's bounding box
  (`tw` 78-102, `pyr` 83.9-119.9) puts its minimum at **(90.00, 101.85)** — the centre, on the row
  already sampled — with a clean bowl around it. Nothing is hiding off-axis.
* **Is 1.5 mHa near the model's numerical floor?** No: convergence noise is ~1e-06 mHa, and the
  local cone slopes are ~1.3 mHa/deg.

So at the reference rung the state-averaged gap does not drop below ~1.5 mHa anywhere in the loop,
while every other rung is driven to ~0.005 mHa by the same procedure — a factor of ~300.

**This is deliberately not reported as "CAS(12,12) has no intersection in this plane."** 1.5 mHa is
0.04 eV, well inside the error of the model itself (6-31G\*, truncated active space, no dynamic
correlation), so a surface that genuinely touches could plausibly present a floor that size.
Calling it an absence would require a threshold for when a state-averaged minimum gap is compatible
with a true crossing, and no such criterion has been established here (see `docs/todo.md` §6). What
is established is the contrast between this rung and the others, which is large, reproducible, and
not explained by fitting, sampling or convergence.

It matters because loop transport returns **pi** on that same loop at CAS(12,12), stable at N=13
and N=21 with endpoint overlap -1.000000. The two methods are hardest to reconcile exactly where
the active space is largest.

### Locating what loop transport encircles, by bisection

The ladder above locates the *state-averaged* intersection. Loop transport never evaluates a gap,
so its object has to be located differently: shrink the loop about a fixed centre until the phase
turns over, and the transition radius is the elliptical distance to whatever is enclosed
(`examples/run_localization.py`, method explained in `notebooks/locating_intersections.ipynb`).

At CAS(2,2), about the `B_x` centre (90, 101.85) with the `B_x` shape (12, 18):

| probe scale | 1.00 | 0.73 | 0.62 | 0.53 | 0.49 | 0.45 | 0.39 | 0.28 | 0.08 |
|---|---|---|---|---|---|---|---|---|---|
| verdict | pi | pi | pi | refused | refused | 0 | 0 | 0 | 0 |

**rho = 0.5385 +- 0.0843**, bracketed between 0.4542 and 0.6228, from a sequence that is monotone
in verdict with the two refusals falling inside the bracket exactly as expected -- that band is the
resolution, not a defect.

**The gap scan's intersection for this same active space is at rho = 0.2121**, from the direct 2D
search position (89.97, 105.67). That is a factor 2.5 inside the measured bracket, and decisively
outside it. The two methods do not localize to the same point.

A second centre at (90, 90) reports **0 cleanly** at full size -- both settings agreeing, not a
refusal -- which excludes 63% of the measured circle and leaves the enclosed object at `pyr`
between 105.6 and 111.5 with `tw` between 84 and 96. A third centre at (99, 101.85) hit the step
floor at full size and was refused, so **no triangulation is available and no position is claimed**.
Two constraints confine it to an arc; they do not pin it to a point, and with only one distance
there is no residual to check consistency with.

Cost: 50 557 micro-iterations, 2 hours, at the *cheapest* active space. Repeating this at
CAS(12,12), where the disagreement is sharpest, is a cluster job (`docs/todo.md` §9).

### A consequence for the Berry phase

The loops are centred on the CAS(12,12) intersection at `pyr = 101.85`, with radius 18 degrees in
`pyr`. Measuring each rung's own intersection against that centre:

| CAS | displacement from loop centre | inside the 18 deg loop? |
|---|---|---|
| CAS(4,4) | 7.98 | yes |
| CAS(6,6) | 12.96 | yes |
| CAS(8,8) | **19.20** | **no — just outside** |
| CAS(10,10) | 3.19 | yes |
| CAS(12,12) | 0 | yes |

On the first two systems the Berry phase was immune to the errors that defeated the comparator,
because every active space still enclosed the intersection. CAS(8,8) appears to break that: it
places its intersection *outside* the loop it is being asked about. The prediction was therefore
that CAS(8,8) would report a **trivial** phase while the others report pi.

**That prediction was wrong, and the way it failed is the most informative result in this file.**

## Berry phase: pi at every rung

| loop | CAS(4,4) | CAS(6,6) | CAS(8,8) | CAS(10,10) |
|---|---|---|---|---|
| `B_x` | **pi** | **pi** | **pi** | **pi** |
| `B_1` (control) | 0 | 0 | 0 | 0 |
| `B_2` (control) | 0 | 0 | 0 | 0 |

All 24 runs report status `OK`, with endpoint overlaps of exactly +-1.000000 and minimum adjacent
overlaps between 0.81 and 0.99. The two controls agree with each other to ~3e-04 in the product
estimator. Nothing is marginal.

### Why the prediction failed

It was not a grid artifact. Refining the CAS(8,8) state-averaged minimum on a 2-degree grid
confirms it at `pyr = 121.20` — a displacement of **19.35 degrees** from the loop centre, genuinely
outside the 18-degree radius:

| pyr | 112 | 114 | 116 | 118 | **120** | **122** | 124 | 126 | 128 |
|---|---|---|---|---|---|---|---|---|---|
| SA gap (mHa) | 8.88 | 6.69 | 4.62 | 2.69 | **0.90** | **0.74** | 2.20 | 3.47 | 4.56 |

The prediction rested on an assumption that turns out to be false: **that the state-averaged gap
minimum is the point the Berry phase encircles.** It is not. The Berry phase transports the
*state-specific* CASSCF ground state, whose degeneracy with S1 sits wherever the state-specific
surfaces touch; the scan minimum is computed with *state-averaged* orbitals, which are a
different set. For a converged active space the two coincide. For a truncated one they need not,
and for butadiene they evidently do not.

Shrinking the loop about the same centre at CAS(8,8) bounds where that degeneracy actually is:

| radius in `pyr` | `pyr` range | N | min adjacent overlap | product | status |
|---|---|---|---|---|---|
| 18 | [83.9, 119.9] | 15 | 0.89 | −0.537 | OK, **pi** |
| 12 | [89.9, 113.9] | 31 | 0.87 | −0.654 | OK, **pi** |
| 8 | [93.9, 109.9] | 31 | 0.77 | −0.519 | FAILED (continuity) |
| 8 | [93.9, 109.9] | 15 | 0.24 | −0.176 | FAILED (continuity) |

The 12-degree loop reaches only to `pyr = 113.9` and still returns a clean pi. **Whatever the
Berry phase is encircling therefore lies below 113.9 — it is not the state-averaged minimum at
121.2.** The 8-degree loops trend to the same sign but their adjacent overlaps collapse (0.24 at
N=15, recovering to 0.77 at N=31), which is what a loop passing close to a degeneracy looks like;
the continuity check refuses them rather than reporting a phase, which is the intended behaviour.

Suggestively, the bound `pyr < 113.9` is much closer to where the two largest rungs put the
intersection (105.04 and 101.85) than to CAS(8,8)'s own state-averaged estimate. The
state-specific transport appears less distorted by the truncation than the state-averaged gap
surface computed in the same active space — though with no exact reference this cannot be
settled here.

### The inversion worth noting

The Berry phase returns the same answer — pi on the enclosing loop, 0 on both controls — at every
active space tested, while the state-averaged estimate of *where the intersection is* scatters
over 19 degrees and never settles. **The topological method is more stable across the ladder than
the quantity that was supposed to validate it.** That is the same robustness seen on formaldimine
and ethylene, but here it holds even though the comparator has given up entirely.

It does not make the Berry phase right by default: with no exact reference, "pi at every rung"
could in principle be four consistent errors. But it is consistent, it passes every internal
check, and it costs a fraction of the comparator.
