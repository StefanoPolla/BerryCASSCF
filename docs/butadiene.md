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

## Results

*(filled in when the ladder completes)*
