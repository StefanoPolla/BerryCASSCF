# Butadiene: a third system, and the search that chose its plane

`docs/active_space.md` found that ethylene's state-averaged comparator is accurate at CAS(2,2)
but not *stable* until CAS(10,10). Butadiene is the follow-up, chosen because it should fail for
a different and more fundamental reason.

Narrative with figures: `notebooks/active_space.ipynb` (final section).
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
because every active space still enclosed the intersection. Butadiene is the first case where
that protection is lost: CAS(8,8) places the intersection *outside* the loop it is being asked
about. If the Berry phase is honest, that rung should report a **trivial** phase while the others
report pi — not because the method failed, but because at that active space the loop genuinely
does not enclose anything. The Berry stage tests exactly this.
