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

Ethylene, error in the intersection position against the CAS(12,12) reference (cone model;
the reference itself sits at φ = 110.90):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| φ (deg) | 111.08 | 116.41 | 102.64 | 103.90 | 110.93 | 110.90 |
| error (deg) | **+0.18** | +5.50 | −8.26 | −7.00 | +0.03 | ref |

These are the cone-model values. The earlier parabola-on-gap estimates gave +0.02, +4.85, −7.91,
−6.90, −0.90 against a reference at 110.00 — shifted by up to ~1.8 deg, including the reference
itself. **The verdict is unchanged either way** (accurate from CAS(2,2), stable only from
CAS(10,10)), which is worth stating: the model correction matters for the numbers but not for the
conclusion they support.

Butadiene, from fine 1-deg cuts refined with the cone model (no reference exists):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| pyr (deg) | 105.86 | 109.68 | 114.46 † | 120.87 | 105.07 | 102.17 |
| shift | — | +3.83 | +4.78 | +6.41 | **−15.80** | −2.90 |
| closest approach (mHa), 1-deg cut | 0.67 | 0.92 | 3.13 † | 0.00 | 0.66 | 2.94 |

† The CAS(6,6) cut is **not a V**. Re-sampled at 0.25 deg it reads

| pyr | 114.50 | 114.75 | 115.00 | 115.25 | 115.50 | ... | 116.50 |
|---|---|---|---|---|---|---|---|
| gap (mHa) | 0.955 | **0.847** | 1.706 | 4.049 | 4.118 | | 4.534 |

so there is a genuine near-degeneracy at **pyr ≈ 114.75, gap 0.85 mHa**, and then a hard jump at
115.1 onto a much shallower branch near 4 mHa. That is a state-averaged CASSCF **solution switch**,
not cone structure: the scan is following two different stationary points either side of it. The
cone fit correctly declines to model it (`nan` residual) rather than returning a number. The
position quoted for this rung is the grid minimum, not a fit, and the rung should not be read as
evidence about where a cone sits.

The "closest approach" column must also be read with care: on a 1-deg grid it is dominated by
resolution, not by how near the cut passes. Re-sampling at 0.25 deg drops CAS(4,4) from 0.92 to
**0.098 mHa** and CAS(10,10) from 0.66 to **0.237 mHa**, both with cone-fit residuals near 1e-02 —
so those cuts do pass essentially through a genuine intersection.

### The fitted closest approach is unreliable, and a direct search says so

Seeding a derivative-free minimization of the gap from each fitted position
(`examples/run_gap_minimum_search.py`, 40 cold SA-CASSCF solves per rung) gives:

| CAS | fitted position | fitted closest approach | search found | search gap | moved |
|---|---|---|---|---|---|
| (2,2) | (90.00, 105.86) | 0.276 mHa | (89.97, 105.67) | **0.003 mHa** | 0.19 deg |
| (4,4) | (90.00, 109.68) | 0.341 | (89.98, 109.41) | **0.005** | 0.27 |
| (6,6) | (90.00, 114.46) | 0.990 | (**89.11**, 114.61) | **0.006** | **0.91** |
| (8,8) | (90.00, 120.87) | 0.201 | (89.98, 121.08) | **0.003** | 0.22 |
| (10,10) | (90.00, 105.07) | 0.303 | (89.96, 104.90) | **0.006** | 0.18 |

Three results, one of which corrects this document:

1. **The fitted positions are good to ~0.2 deg.** The cone model does its job.
2. **The fitted closest approach overestimates by about two orders of magnitude — wherever an
   intersection exists.** CAS(2,2) through CAS(10,10) all reach 0.003–0.006 mHa under direct
   search, against fitted values of 0.2–1.0 mHa. So **those active spaces do contain a genuine
   conical intersection in this plane**, and the earlier reading of a large fitted closest approach
   as evidence of an avoided crossing was wrong for them. A fit constrained to a line through an
   assumed position cannot do better than that line allows.

3. **CAS(12,12) behaves differently under the same search, and not because of the fit.** From its
   fitted position the search reaches **1.4595 mHa** in 58 evaluations — an 18% improvement on
   the fitted 1.776, against the 47x–163x collapse at every other rung. **The budget was
   deliberately raised to settle this**: an earlier 30-evaluation run gave 1.4963, so doubling the
   effort moved the answer by 2.5%. The other five rungs reach 0.003–0.006 mHa in 40 evaluations.
   None of these searches converged in the optimizer's sense — each stops on its evaluation budget,
   since a CAS(12,12) evaluation costs ~2.2 minutes — but the budgets are now comparable and the
   contrast is not a budget artifact. **What that floor means is a separate question** — see the
   note on thresholds below.
4. **CAS(6,6)'s intersection is at `tw` = 89.11, not 90.** That is why its `pyr` cut at `tw` = 90
   fitted badly and appeared to jump: the cut simply missed the apex. The one-dimensional `tw`
   check at the fitted `pyr` reported 89.982 and did not catch it, because that cut was taken at
   the wrong `pyr`. Two one-dimensional cuts through a two-dimensional surface are not a
   substitute for a search.

**What counts as "a conical intersection" needs a threshold this work does not have.** It is
tempting to read CAS(12,12)'s 1.5 mHa floor as "no intersection in the plane", and an earlier
revision of this document did. That over-reads it. 1.5 mHa is 0.04 eV — comfortably inside the
error of the model itself (basis set, active space, no dynamic correlation), so a surface that
truly touches could easily present a floor of that size here. The statement that survives is
**relative**: the search drives every other rung to 0.003–0.006 mHa and this one only to 1.5, some
300x higher, so whatever the search is converging onto at CAS(12,12) is a different kind of object
numerically. In absolute terms 1.5 mHa is still a near-degeneracy, and loop transport returning π
there is not in contradiction with it. Deciding the question properly needs a criterion for when a
state-averaged minimum gap is compatible with a true crossing — some combination of an absolute
tolerance and the local cone slopes, which here are ~1.3 mHa/deg. For now the honest description is
"reasonably small, and unusually large relative to the other rungs". See `docs/todo.md` §6.

### How close the plane comes to an intersection is itself active-space dependent

Re-sampling every rung at 0.25 deg gives a sharper and more awkward result:

| CAS | (4,4) | (6,6) | (10,10) | (12,12) |
|---|---|---|---|---|
| minimum gap at 0.25 deg (mHa) | **0.098** | 0.847 † | **0.237** | **1.479** |
| branches extrapolate to | ~0 | discontinuity | ~0 | **~1.7** |

CAS(4,4) and CAS(10,10) close to within 0.1–0.2 mHa: the plane contains a genuine intersection for
them. **CAS(12,12) stops an order of magnitude higher** — its two branches have clean opposite
slopes (−1.31 and +1.35 mHa/deg) meeting at ≈1.7 mHa, far above the 1e-06 mHa convergence noise,
and the free search below gets only to 1.496. Whether that counts as an intersection depends on a
threshold this work has not fixed; what is certain is that it is not resolution or fit error.

Loop transport at CAS(12,12) nonetheless returns **π**, stable across two discretizations
(N=13 and N=21, min overlaps 0.83 and 0.88, endpoint −1.000000 exactly).

### The controls at that rung, which the claim had been missing (2026-09-18)

Until now CAS(12,12) had been run on the **enclosing loop only**. A π with nothing to contrast
it against is weak evidence: a solver in difficulty at a large active space could produce one
anywhere. Both control loops have now been run at the same rung, at three discretizations each
(`slurm/berry_cas12.job`):

| loop | N=13 | N=21 | N=31 | min abs adjacent overlap |
|---|---|---|---|---|
| `B_x` (encloses) | **π**, Π = −0.493 | **π**, Π = −0.626 | **π**, Π = −0.715 | 0.83, 0.88, 0.92 |
| `B_1` (control) | 0, Π = +0.868 | 0, Π = +0.915 | 0, Π = +0.941 | 0.986, 0.995, 0.997 |
| `B_2` (control) | 0, Π = +0.868 | 0, Π = +0.915 | 0, Π = +0.941 | 0.986, 0.994, 0.997 |

Every run reports `OK`, and every endpoint estimator is ±1.000000 exactly. The enclosing loop's
worst overlap improves monotonically with refinement (0.83 → 0.88 → 0.92) while |Π| grows
(0.49 → 0.63 → 0.72), which is what a converging discretization looks like rather than a
marginal one. Three things follow:

* **the π is selective.** The two controls are trivial at every discretization, so whatever
  produces π on `B_x` does not produce it on loops of the same size and shape 30 deg away;
* **the enclosing loop is the hard one**, by the method's own difficulty measure: its worst
  adjacent overlap is 0.83–0.88 where the controls sit at 0.99. That is the signature of a loop
  passing near a degeneracy, and it is *absent* on the controls;
* **the two controls agree with each other to 1e-4** (Π = +0.8679 against +0.8678 at N=13,
  +0.9409 against +0.9409 at N=31). They are mirror-image loops traversed completely
  independently, so this is an end-to-end check of the pipeline, not a shared intermediate.

The third of the three readings above — that the loop-transport result is simply wrong at this
rung — now requires a failure that is selective in exactly the way a real degeneracy would be,
and that leaves the two controls untouched at three discretizations each. It is no longer the
economical explanation.

**Two objections to reading that as a disagreement were raised and tested; both have now been
answered, and the claim still should not be pushed as far as it first was.**

*Objection 1: the floor came from a fit along a line.* Fits of that kind run 47x–163x too high on
every rung where a direct search was tried, so the 1.7 mHa could have been an artifact. A free 2D
minimization from that position reaches **1.4595 mHa** in 58 evaluations, 18% below the fit, where
the same procedure on every other rung falls by two orders of magnitude in 40. Doubling the budget
from 30 to 58 moved it by 2.5%, so the floor is not an artifact of stopping early either.

*Objection 2: only a cross through the loop had been sampled.* The scan covered one row at
`tw` = 90 over `pyr` 80–140 and one column at `pyr` = 102.2 over `tw` 84–96, while the loop spans
`tw` 78–102 and `pyr` 83.9–119.9 — so an intersection off that cross would have produced exactly
what was seen. **A 5×5 grid over the loop's bounding box has since been run.** Its minimum is at
(90.00, 101.85) — the centre, on the cross already sampled — with a clean bowl around it and
nothing hiding off-axis.

So both checks came back in favour of the original observation, which is the opposite of what
happened to several other claims in this document. **But "the plane contains no intersection at
CAS(12,12)" is still more than the data supports**, for the threshold reason given in §1: 1.5 mHa
is 0.04 eV, within the model's own error, and a genuine crossing could present that floor.

What can be stated is the contrast. At CAS(12,12) the gap scan finds nothing below 1.5 mHa
anywhere in the loop's area, where every other rung is driven to ~0.005 mHa by the same search;
and on that same loop, loop transport returns π with endpoint overlap −1.000000 at two
discretizations. Whether that is a real disagreement about *what is enclosed*, or the two methods
describing the same near-degeneracy with different resolution, is **the sharpest open question in
the project**. Three readings are live:

* the state-specific surfaces have a degeneracy the state-averaged ones do not — finding 3 in its
  strongest form;
* both methods see the same crossing, and 1.5 mHa is simply what this model's numerical floor
  looks like at this rung; or
* the loop-transport result is wrong here — both runs report falling back to a weaker warm start at
  several points, which is the failure mode §5 discusses.

### The floor, converted into a distance, is not a disagreement about enclosure

The floor was being compared against a tolerance nobody had fixed. The quantity that can be
fixed is a **distance**: near a conical intersection the gap is linear in the branching-plane
coordinates, so a floor `g` on a cut whose local slope is `a` is exactly what a cut passing
`g / a` from the apex would show. `examples/report_gap_criterion.py` reports that conversion for
every rung (`results/<system>/gap_criterion.json`):

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | **(12,12)** |
|---|---|---|---|---|---|---|
| floor found by direct search (mHa) | 0.0032 | 0.0045 | 0.0061 | 0.0026 | 0.0064 | **1.4595** |
| local slope (mHa/deg) | 1.310 | 1.237 | 1.366 | 0.816 | 1.266 | 1.749 |
| implied miss distance (deg) | 0.002 | 0.004 | 0.004 | 0.003 | 0.005 | **0.835** |
| as a fraction of the loop's 18 deg semi-axis | 0.0001 | 0.0002 | 0.0002 | 0.0002 | 0.0003 | **0.046** |

**CAS(12,12)'s floor is consistent with an intersection 0.83 deg from the searched point —
4.6% of the loop radius, and therefore inside the loop.** So the gap scan does *not* say the
loop encloses nothing; it says that if what it is approaching is a cone, the cut misses the apex
by under a degree. A π from loop transport on an 18 deg loop is not in contradiction with that,
and this document's earlier framing — "the gap scan finds nothing below 1.5 mHa anywhere in the
loop" — invited a stronger reading than the number supports. The 300x contrast between this rung
and the others is real and still unexplained; what it is *not* is evidence about enclosure.

Three caveats, all of which make the conclusion stronger rather than weaker except the last:

* the direct search stops on its evaluation budget, so its floor is an upper bound. A lower true
  floor means a *smaller* miss distance;
* the conversion assumes the perpendicular slope equals the fitted in-cut slope. Anisotropy of
  a factor two moves the miss distance by a factor two — 0.4 or 1.7 deg, both still inside;
* it assumes the surface near the minimum is a cone at all. If it is an avoided crossing with a
  genuine 1.5 mHa gap, there is no apex to miss and the distance is meaningless. Nothing here
  distinguishes those two, which is exactly why the **measurement** in the next section — where
  loop transport is asked to locate what it encircles — is the experiment that settles it.

The reading that survives is therefore the second of the three above: both methods may be seeing
the same near-degeneracy, at different resolution. What remains genuinely open is not enclosure
but **position**: at CAS(2,2) the located state-specific object is a factor 2.5 further out than
the gap-scan minimum, far outside the measurement's own uncertainty, and that is not explained by
a floor-to-distance conversion.

### Asked directly, in one probe: the answer is no (2026-09-18)

Reconstructing a position first was the long way round. A single loop, centred **on** the gap-scan
minimum and small enough to contain nothing else, answers the question by itself:

> butadiene CAS(12,12), circular loop of radius **1.25°** about (89.968, 101.944) — the searched
> minimum itself — returns a clean **0**. Both step-control settings agree, worst adjacent overlap
> 0.978, at a cost of 7 066 micro-iterations (2.9 h on one dedicated cluster core).

So the state-specific transport does **not** encircle the state-averaged minimum at the rung where
the disagreement is sharpest. That is finding 3 in its strong form, and it now rests on one probe
rather than on an inference chain.

The argument closes tightly with §1's floor-to-distance conversion. That conversion says the
1.46 mHa floor is what a cut passing **0.83°** from a cone apex would show — and 0.83° is *inside*
this loop. So even granting the gap scan its own most favourable reading, the object it is
approaching lies within the loop that loop transport says encloses nothing.

Two caveats, neither of which rescues the alternative. A zero verdict means an **even** count, so
two degeneracies inside a 1.25° loop would also read 0 — possible in principle, unmotivated here.
And the same probe at ethylene CAS(2,2) behaves the same way (§2), where the geometry is understood
well enough to show that what loop transport encircles is not a point degeneracy at all. The
pattern across the two systems is that **loop transport and the gap scan are not sensing the same
object**, and at ethylene the state-specific side looks like solution structure rather than
topology.

**And the larger radii measure how far away it actually is.** The same scan, continued on the
cluster:

| radius about the gap-scan minimum | 10° | 5° | 2.5° | 1.25° |
|---|---|---|---|---|
| verdict | **π** | refused | *running* | **0** |
| worst adjacent overlap | 0.908 | 0.980 | | 0.978 |
| cost (micro-iterations) | 8 323 | 11 614 | | 7 066 |

Enclosed at 10°, grazing at 5°, empty at 1.25°: **whatever loop transport encircles sits about 5°
from where the gap scan puts its minimum.** Three probes, against a bisection that ran twelve hours
without completing one.

The small-loop localization running alongside agrees and narrows it further. A (6, 9) loop about
(90, 106) returns π, so the object is inside that ellipse; the same loop about (90, 120) returns 0
at full size and at 0.08, so it is outside that one. Together with the 5° bracket, and taking the
`tw` = 90 line, that puts it near **pyr ≈ 107** against the gap scan's 101.94.

Separating them needs the position of whatever loop transport encircles to be **measured** rather
than inferred, which is what the bisection and triangulation experiment in `docs/todo.md` §2 is
for: shrink the loop until the phase turns over, and repeat from a second centre.

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
| butadiene | CAS(2,2) | correct at every rung, CAS(2,2) through CAS(10,10) |

Butadiene's CAS(2,2) is the sharpest version of this. It is **below the pi space** and not a
chemically defensible active space for the molecule, yet loop transport returns π on the
enclosing loop and 0 on both controls, with *higher* adjacent overlaps (0.92 at N=13, 0.96 at
N=21) than any larger rung. The topological answer survived an active space that cannot represent
the pi system.

That cuts both ways, and the second reading matters: if an active space this poor still gives the
same answer, the answer is either very robust or insensitive to something it ought to be sensitive
to. Nothing in the present data separates those.

The cost asymmetry is large: loop transport at the smallest active space is seconds, against hours
for the gap scan ladder.

### The second reading wins at ethylene CAS(2,2): the π is not the intersection (2026-09-18)

Small-loop probing (`notebooks/probing_by_small_loops.ipynb`) asked the question this table never
did — not "what phase does the loop report?" but "**does the loop encircle the intersection?**" At
ethylene CAS(2,2) the answer is no, and the evidence is cheap and direct:

* **five concentric loops containing the state-averaged intersection** — radii 3°, 2°, 1°, 0.5°,
  0.25° about (90, 110.9) — all report a **trivial** phase. The exact full-valence position is at
  the centre of every one of them;
* shrinking `E_x` gives a monotone sequence — π at 12°, 10°, 8°, 6°, refused at 5° and 4°, zero at
  3° and below — so whatever carries the phase sits **4.5 ± 1.5° from the intersection**;
* nine probes of radius 2° covering the mirror line from `pyr` 99 to 123, the whole span of `E_x`,
  all report zero: it is not on the line either;
* the only probes that react are at (94, 111) and (86, 111), both **refused**. But
  `E(tw = 94) = E(tw = 86)` to 6e-13 — the mirror is exact — so anything off the line comes in
  pairs, a loop centred on the line encloses both or neither, and an **even count cannot give π**.

No arrangement of point degeneracies compatible with the symmetry explains the observations. What
the region does contain is **competing state-specific CASSCF solutions**: the in-CAS gap map is
ragged at exactly `tw` = 87 and 93, and cold-start energies around a 3° circle jump at four of 72
points. A sign picked up crossing between CASSCF solutions is not a Berry phase.

**So the ethylene CAS(2,2) row above says "correct" for what now looks like the wrong reason.** The
verdict π happens to match the topology a converged calculation would give, and it is not obtained
by encircling the intersection. The row stands as a description of what the method *reported*; it
should not be read as evidence that the method *worked* at that rung.

Three limits on this, all real. It is one rung of one system. CAS(2,2) is below the π space and
this document already calls it indefensible for ethylene. And the cold-start ring maps where
independent solves disagree, not what a warm-started continuation does — it establishes that
competing solutions exist there, not that the walk crossed between them.

**The first limit has since weakened: CAS(4,4) and CAS(8,8) behave the same way.** Concentric loops
on the best available position at each rung return **no π at any radius tested**:

| ethylene rung | r = 4° | 2° | 1° | 0.5° | 0.25° | 0.1° |
|---|---|---|---|---|---|---|
| CAS(4,4), about (90.00, 108.40) | refused | refused | refused | **0** | refused | **0** |
| CAS(8,8), about (89.88, 107.79) | refused | refused | *(stopped after two refusals)* | | | |

Formaldimine CAS(4,4) under the identical procedure is monotone and unambiguous — π at 8°, 4°, 2°
and 1°, refusals at 0.5° and 0.25°, a clean 0 at 0.05°.

*(Costs here are micro-iterations, never laptop seconds: the ethylene CAS(8,8) probes took 15 403
and 10 721 micro-iterations — the second is less work than the first — while taking 655 s and
14 192 s of wall time on a shared, sleeping machine. An earlier revision cited those seconds as
the reason for stopping that run; the work counts do not support it, and the real reason is that
two refusals had already answered the question.)* Ethylene's verdicts do not even order: the
only two clean ones say "nothing inside", at radii a factor of five apart, and every intermediate
radius is refused. **At three ethylene rungs there is nothing locally encircle-able where the
methods place the intersection**, which is the same conclusion §2 reached at CAS(2,2) by mapping
rather than by scanning radius.

### The check, run systematically: it agrees where it should (2026-09-18)

`examples/verify_encirclement.py` puts one 2° loop on each rung's gap-scan intersection and asks
whether loop transport encircles it. Formaldimine settles whether the check is trustworthy:

| system | CAS | gap-scan intersection | verdict | encircles it? |
|---|---|---|---|---|
| formaldimine | (2,2) | (141.00, 90.00) | 0 | **no** |
| formaldimine | (4,4) | (130.00, 90.00) | π | **yes** |
| formaldimine | (6,6) | (130.92, 90.07) | π | **yes** |
| butadiene | (2,2) | (89.97, 105.67) | 0 | **no** |
| butadiene | (4,4) | (89.98, 109.41) | 0 | **no** |
| butadiene | (6,6) | (89.11, 114.61) | 0 | **no** |
| butadiene | (8,8) | (89.98, 121.08) | 0 | **no** |
| butadiene | (10,10) | (89.96, 104.90) | refused | *grazing — the two positions are 2° apart here* |
| butadiene | (12,12) | (89.97, 101.94) | 0 | **no** (r = 1.25°) |

The CAS(2,2) row is a control nobody designed: that rung's gap scan is independently known to be
qualitatively wrong — a spurious minimum near α = 141° while the real intersection is at 132.6°
(`docs/results.md` §3) — and the probe rejects it. Both converged formaldimine rungs are accepted,
and at CAS(4,4) the object encircled is the one triangulation found at (128.99, 90.24), 1° inside
the loop.

So the check agrees where the two methods are known to agree and refuses where one is known to be
wrong. **That is what licenses reading butadiene's "no" as evidence** rather than as an artefact of
the procedure — and butadiene says no at every rung that returns a verdict, with the one refusal
falling exactly where the two positions are closest.

## 3. The gap-scan minimum and the point loop transport encircles are different objects

> Two subsections bearing directly on this question sit in §1, where the CAS(12,12) floor is
> discussed: **the controls at that rung** (both trivial at three discretizations, so the π is
> selective) and **the floor converted into a distance** (1.46 mHa is a 0.83 deg miss — inside
> the loop, so the floor was never evidence about enclosure). Read them with this section.

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

### The difference has now been measured, not only bounded

Bisection on the loop radius (`berrycasscf/localize.py`, `notebooks/locating_intersections.ipynb`)
measures the elliptical distance from a centre to whatever a loop encircles, so the state-specific
object can be located rather than bracketed by hand. Two systems have been done.

**Formaldimine, where an FCI reference exists**, comparing at the *same* active space so that only
the state-specific/state-averaged difference is in play:

| CAS | loop transport (state-specific) | gap scan (state-averaged) | difference | gap scan − FCI | residual | residual / precision |
|---|---|---|---|---|---|---|
| (2,2) | (131.67, 86.30) | 151.61 | −19.94 | +19.00 | 0.197 | **1.57 — inconsistent** |
| (4,4) | (128.99, 90.24) | 130.43 | **−1.44** | −2.18 | 0.0088 | 0.29 |
| (6,6) | (129.73, 88.49) | 131.35 | **−1.63** | −1.26 | 0.0574 | 0.76 |

**The offset is stable where the method is applicable: −1.44 deg at CAS(4,4) and −1.63 deg at
CAS(6,6)**, with the state-specific degeneracy consistently *below* the state-averaged one. That is
comparable to the 2.18 and 1.26 deg by which the gap scan itself misses FCI — so the
state-specific/state-averaged difference is not a small correction on top of active-space error, it
is the same size.

**The last column is the consistency criterion, and it had to be normalised to be meaningful.**
Raw residuals of 0.197, 0.0088 and 0.0574 would suggest CAS(6,6) is six times worse than CAS(4,4).
Divided by each run's own precision — the RMS of its bisection bracket half-widths — they become
1.57, 0.29 and 0.76. CAS(6,6) is not less consistent, only less precisely measured (its brackets
are ±0.06–0.09 against CAS(4,4)'s ±0.03), and **CAS(2,2) is the only rung whose misfit exceeds its
own uncertainty**. Comparing raw residuals across runs of different precision is meaningless, and
this document did it in an earlier revision.

The CAS(2,2) row is therefore not a measurement at all: three centres cannot be explained by one
degeneracy there, and the construction declines to report a position (see §8).

**Butadiene CAS(2,2), where the disagreement was first noticed.** The bisection about the loop
centre (90, 101.85) brackets the transition at **rho = 0.5385 ± 0.0843**, from a verdict sequence
monotone in radius (π at 1.00, 0.73, 0.62; 0 at 0.45, 0.39, 0.28, 0.08). The gap scan's
intersection for the same active space, located by direct 2D search at (89.97, 105.67), sits at
**rho = 0.2121** — a factor 2.5 inside the bracket and decisively outside it.

A second centre at (90, 90) reports 0 *cleanly* at full size, which excludes 63% of the measured
circle and leaves the enclosed object at `pyr` between 105.6 and 111.5 with `tw` between 84 and 96.
A third centre was refused at full size (step floor), so **no triangulation is available and no
position is claimed** — two constraints confine it to an arc, they do not pin it to a point.

What the butadiene measurement does establish is the thing §3 asserts: **loop transport encircles
something that is not where the gap scan puts its intersection.** Until now that was inferred from
a loop still returning π when the gap-scan minimum had moved outside it; it is now measured.

The same experiment at CAS(12,12), where the disagreement is sharpest, is a cluster job — this one
took two hours at the cheapest active space (`docs/todo.md` §9).

## 4. The inversion: the cheap method is the stable one

Across butadiene's ladder, loop transport returns the same answer at every rung while the gap scan
— the method nominally validating it — scatters over 19 deg and never settles, at three orders of
magnitude more cost.

This does **not** make loop transport right by default. With no exact reference for butadiene,
"π at every rung" is equally consistent with "very robust" and "consistently wrong in the same
direction". What can be said is that it is self-consistent, passes every internal check, and is
cheap. Settling it needs a system where an exact reference exists *and* the ladder still
misbehaves.

### Ethylene is that system, and the answer is "comparable", not "better" (2026-09-18)

An earlier revision of this section said no system in the project provides the combination. That
was wrong about ethylene: CAS(12,12) is its **full valence space**, so the intersection position
at that rung is exact in this basis, and its ladder is the badly behaved one — accurate at
CAS(2,2), displaced by 5–8 deg in the middle.

Localizing the state-specific object there by bisection, at four rungs, and comparing both methods
against that exact position (`notebooks/locating_intersections.ipynb`):

| CAS | loop transport region (pyr) | its centre | error | gap-scan position | its error |
|---|---|---|---|---|---|
| (2,2) | (110.3, 111.5) | 110.92 ± 0.58 | **0.02** | 111.08 | 0.18 |
| (4,4) | (102.9, 109.1) | 106.00 ± 3.10 | 4.90 | 116.41 | 5.51 |
| (6,6) | (99.4, 109.2) | 104.30 ± 4.90 | 6.60 | 102.64 | 8.26 |
| (8,8) | (99.7, 111.1) | 105.40 ± 5.70 | 5.50 | 103.90 | 7.00 |

**Loop transport is nearer the exact answer at every rung, and the margin is not significant.**
It wins by 0.2–1.7 deg while its own half-width is 0.6–5.7 deg, so at the three upper rungs the
difference is smaller than the uncertainty on it. The one unambiguous statement is CAS(2,2), where
loop transport is exact to 0.02 deg within a 0.6 deg band — at an active space that cannot
represent the π system.

The sharper question is whether the measurement can **tell the two candidate answers apart**, and
mostly it cannot: the region contains both at CAS(2,2) and CAS(8,8), neither at CAS(4,4), and at
CAS(6,6) it contains the gap scan's position while **excluding** the exact one. So this does not
vindicate loop transport; it bounds how well the present bisections resolve anything, which is to
±3–6 deg once the active space stops being trivial.

Two honest caveats, in the direction that weakens the result. Regions are not positions: a wide
region is easier to be "near" an answer with. And CAS(8,8)'s triangulated point, 107.79, sits
2.4 deg from its own region centre of 105.40 — two estimates from the same data at the same rung,
which is the plainest available measure of the precision actually on offer.

### Butadiene, both methods' positions, all six rungs (2026-09-23)

The phase result above is qualitative — same answer at every rung. Localizing the state-specific
object by bisection turns it into positions, and puts the two methods on the same footing
(`notebooks/butadiene_ladder.ipynb`, final section):

| CAS | ρ about the `B_x` centre | implied `pyr` * | triangulated | gap scan | difference |
|---|---|---|---|---|---|
| (2,2) | 0.5385 ± 0.0843 | 111.55 | 1/3 centres | 105.67 | +5.88 |
| (4,4) | 0.7016 ± 0.0277 | 114.48 | 1/3 centres | 109.41 | +5.07 |
| (6,6) | 0.2367 ± 0.0461 | 106.11 | (91.18, 105.73) | 114.61 | −8.50 |
| (8,8) | 0.2367 ± 0.0461 | 106.11 | 1/3 centres | 121.08 | −14.97 |
| (10,10) | 0.2820 ± 0.1058 | 106.93 | (89.87, 96.78) | 104.90 | +2.03 |
| (12,12) | not bracketed † | | | 101.94 | |

\* assuming the object lies on the `tw` = 90 line; the two successful triangulations put it within
1.2° of it, so the assumption is close but not exact.
† its three centres were run twice on the cluster and the bracketing centre never completed a probe
inside three days; the question was answered instead by the cheaper probe in §3, which is the right
lesson about instrument choice.

**Loop transport's positions are the steadier ones, quantitatively.** Over the five rungs both
methods answered, loop transport spans 8.4° and the gap scan 16.2° (19.1° including the CAS(12,12)
row loop transport could not bracket). From CAS(6,6) upward loop transport sits at 106–107° — under
a degree of movement across three rungs — while the gap scan over those same rungs reports 114.6,
121.1 and 104.9.

**And they never coincide.** The difference is +5.9, +5.1, −8.5, −15.0, +2.0: neither the size nor
the sign settles. Since the brackets are ±0.5–1.9° in `pyr`, only the CAS(10,10) row has the two
within reach of each other — and that is exactly the rung where the encirclement probe comes back
*refused* rather than a clean zero, which is what a 2° loop grazing an object 2° away should do.
The internal consistency there is worth noting: two independent constructions agree about which rung
is the near miss.

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

Radius 6 at N=31 **passed** the per-run checks (overlap 0.844, endpoint 0.984) and reported π,
while N=15 and N=61 both give the opposite sign — though both of those were themselves refused, so
their signs carry little weight either.

**The defensible statement is that radius 6 is marginal: the Berry phase for that loop is not
determined by these calculations.** An earlier version of this section went further and called the
passing run wrong, arguing that since radius 8 reports π reliably and radius 6 does not, the
degeneracy sits between their reaches. That inference is reasonable, but it rests on the signs of
two refused runs, which this document elsewhere declines to trust. What can be said is that radius
6 passes close enough to whatever is being encircled that the answer depends on the
discretization.

**A diagnostic was in fact present, and was not made binding.** The N=31 run carries the message
*"continuation chain broken"* — the solver's fallback ladder fell through to a **cold start** at one
point, so that point was not reached by continuation at all. Gauge fixing repairs a cold start's
arbitrary *sign*, but not a change of *branch*, and a modest branch change still clears the 0.80
overlap threshold (this run's worst overlap, 0.844, is barely above it). That is a concrete
mechanism for a spurious sign, and it was reported as a warning rather than a failure. Promoting it
to a check would have refused this run on its own evidence, without appeal to the multi-N
criterion; so would a less permissive overlap threshold. Both are planned — `docs/todo.md` §5.

So the per-run checks are necessary but not sufficient *as currently wired*: the information needed
to refuse this run was there and simply not acted on. What rejects it today is the **stability
criterion**, which demands at least two discretizations *all* passing with the same phase. Radius 6
has exactly one passing N, so no phase is claimed for it. The layering is what makes the protocol
safe, not any individual check.

## 6. The required loop discretization grows with the active space

Ethylene, minimum adjacent overlap on the enclosing loop at fixed N = 13:

| CAS | (2,2) | (4,4) | (6,6) | (8,8) | (10,10) | (12,12) |
|---|---|---|---|---|---|---|
| min overlap | 0.895 | 0.718 | 0.689 | 0.673 | 0.692 | 0.686 |

CAS(4,4) and above fall below the 0.80 threshold and are refused at N = 13; all pass at N = 21. A
larger active space has more freedom to rearrange, so the wavefunction turns faster around the
same loop. This is a statement about *loop resolution*, not about whether the active space can
describe the physics, and the two must not be conflated.

## 7. One update per point versus converging each point

arXiv:2304.06070 takes **one parameter update per loop point** rather than optimizing to
convergence, on the argument that updates are the cost and a fixed budget is better spent on more
points with fewer updates each. This package converges every point. Both were run on the same
loops (`notebooks/stepping_comparison.ipynb`, 96 runs over two loops, two active spaces, eight
discretizations and three update budgets).

**The paper's argument holds in the currency it is stated in.** At a fixed number of parameter
updates, single-stepping buys a finer loop, and the discretization error is what limits accuracy:

| budget (updates) | converged | single-step | advantage |
|---|---|---|---|
| 500 | N=13, 1−\|Π\| = 0.348 | N=97, 1−\|Π\| = 0.056 | **6.3x** |
| 900 | N=25, 1−\|Π\| = 0.199 | N=97, 1−\|Π\| = 0.056 | 3.6x |

**In wall time the advantage is real but roughly halves**, to 1.9–2.3x. The cost decomposition
says why: fitting `wall ≈ a·(points) + b·(updates)` separates the per-point overhead — integrals,
the AO→MO transformation, the mean-field solve — from the work that scales with updates. That
overhead does not shrink when fewer updates are taken per point, so doubling the number of points
doubles it regardless.

| | per-point overhead | per-update | ratio | overhead share of a single-step run |
|---|---|---|---|---|
| CAS(2,2) | 33.6 ms | 2.3 ms | 14.4 | 67% |
| CAS(6,6) | 75.2 ms | 19.4 ms | 3.9 | 26% |

**The advantage grows with the active space.** The overhead-to-update ratio falls from 14.4 to 3.9
between CAS(2,2) and CAS(6,6), because the fixed per-point work is roughly unchanged while each
update gets much more expensive. Extrapolating, in the large-active-space regime the classical
cost model converges toward the quantum one and the paper's reasoning should hold in wall time
too — though two active spaces is not enough to establish that limit, and it is stated here as an
expectation rather than a result.

**The lag vanishes steeply with refinement.** `1 − |ω|`, the failure of the transported state to
return to itself around a closed loop, scales as `N^-6.4` at CAS(2,2) and `N^-4.4` at CAS(6,6).
It is exactly zero for converged continuation by construction. That steepness is why
single-stepping is safe at moderate N: the accumulated drift collapses much faster than the
discretization error it competes with.

**Neither mode ever produced a wrong Z2 answer.** Of 96 runs, 94 passed the checks and all 94 were
correct; the two refusals were single-step at N = 9, where the endpoint lag reached 0.47 and 0.61
and the check refused rather than returning a sign. The cheap mode fails loudly, which is what
makes it safe to use.

**So the methodological difference is understood and quantified.** On hardware, where updates are
the currency, single-stepping is the right choice and the paper is right to make it. Classically
the per-point overhead pulls the optimum toward fewer, better-converged points for small active
spaces, and back toward the paper's choice as the active space grows. Converging every point, as
this package does, costs a factor of a few in wall time and buys an endpoint fidelity of 1e-10
rather than 1e-4 — worth it while the calculations are cheap, and not obviously worth it once they
are not.

## 8. Methodological results that are not about active spaces

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

**Adaptive step control never beat uniform discretization on cost, and the prediction that it
would is retracted.** Steering the step size by the measured continuity
(`berrycasscf/adaptive.py`) was expected to save 1.3x on formaldimine, 3.1x on butadiene `B_x` at
CAS(8,8) and 4.7x on ethylene `E_x` at CAS(8,8), estimated from the spread of adjacent overlaps on
saved runs. Measured at *matched quality*, against a swept family of uniform N:

| system | loops | measured |
|---|---|---|
| formaldimine STO-3G, CAS(2,2) and CAS(6,6) | `C_x`, `C_2` | 0.60x–1.23x |
| butadiene 6-31G\*, CAS(2,2) | `B_x`, `B_2` | 0.67x–1.15x |
| ethylene 6-31G\*, CAS(8,8) | `E_x`, `E_2` | 0.75x–1.04x |

**All three predictions were wrong, and the largest one belongs to the loop that did worst against
it.** The cause is measurable and is not the controller. Cost is not proportional to the number of
points: the points where the state turns fastest are also the ones where CASSCF needs the most
iterations (measured correlation between per-point micro-iterations and adjacent overlap:
r = −0.40 at formaldimine CAS(2,2), −0.67 at CAS(6,6)). Adaptive stepping removes cheap points and
adds expensive ones. Taking ethylene `E_x` at CAS(8,8), comparing the *zero-rejection* adaptive run
against uniform at the same worst-case overlap so that rejection cost is excluded entirely:

| | uniform N=33 | adaptive d_max=0.05 |
|---|---|---|
| points | 33 | 26 (**0.79x**) |
| micro-iterations per point | 80.2 | 117.6 (**1.47x**) |
| total | 2647 | 3058 (1.16x) |

The point-count saving is real and is simply outweighed. Rejected trials — full CASSCF solves that
produce nothing, 1 to 7 of them on the harder loops — make it worse again where they occur.

Its value is elsewhere, and is not a speed-up:

* **it saves 19x–43x on loops that pass close to a degeneracy.** Measured on the Jahn–Teller
  model, the number of points each method needs to hold a worst-case adjacent overlap of 0.90:

  | closest approach (loop radius 1) | uniform N | adaptive points | ratio |
  |---|---|---|---|
  | 0.10 | 100 | 19 | 5x |
  | 0.05 | 200 | 20 | 10x |
  | 0.02 | 400 | 21 | 19x |
  | 0.007 | 1000 | 23 | **43x** |

  This is what reconciles "cost-neutral" with "worth building". The loops in the active-space
  studies are of uniform difficulty and gain nothing; the loops a *search* walks are close to a
  degeneracy by construction, because a search probes near the thing it is looking for. It is also
  the difference between an answer and a refusal, which is what makes bisection possible at all;
* **the resolution limit is predictable rather than empirical**: a loop of radius `R` passing at
  distance `eps` needs a step `~ eps*dtheta_max/(2*pi*R)`, so the walk gives up below
  `eps_min ~ 2*pi*R*d_min/dtheta_max` = 0.0070 of the radius with the defaults. Measured: closes at
  0.007, floors at 0.005;
* **N no longer has to be guessed** before anything is known about the loop;
* **failures are diagnosed**: shrinking the step either restores continuity (it was undersampling)
  or does not (the loop passes through something).

The same experiment restates a distinction this document keeps insisting on. Uniform N=24 returns
the **correct sign at every closest approach tested**, including ones where its worst adjacent
overlap has collapsed to 0.66 — and fails its own continuity check from 0.1 inwards. Being right is
not the same as being trustworthy.

**A method returning one bit per loop can still return a position.** Scaling a loop about a fixed
centre by `s` encloses a point `x` exactly when `s > rho(x; c)`, the distance to `x` in units of
the loop's own semi-axes. Bisecting on `s` therefore *measures* that distance, and intersecting the
ellipses from several centres gives a position — two centres leaving a mirror pair, three resolving
it (`berrycasscf/localize.py`). Three properties earn it trust:

* the answer is a **bracket with a refusal band**, because exactly at `s = rho` the loop runs
  through the degeneracy and cannot be walked at any step size. The band's width is the resolution,
  measured rather than assumed;
* with three or more centres the construction is **over-determined**, so its residual is a
  consistency check that refuses to produce a position from data no single degeneracy explains;
* the geometry is validated **exactly**, on synthetic radii from a known position, separately from
  the physics — so a failure can be attributed to one or the other.

On formaldimine CAS(2,2) the individual brackets are tight (rho = 0.3968 ± 0.0046 from one centre)
but the **residual is 0.197 in units of the semi-axes, about 2 deg**: the three ellipses do not meet
at a point. That is the check firing, and it fires where it should. CAS(2,2) is independently known
to be pathological for formaldimine — the gap scan there finds no minimum in the region at all —
and a direct map of the state-specific in-CAS S1/S0 gap over the region **never falls below 321
mHa**. That map is also exactly symmetric about phi = 90, so any off-axis degeneracy must have a
mirror partner, and a loop centred on the line would enclose both and report 0.

**A check that is never wired to a verdict is not a check.** Two of this project's diagnostics were
printed for most of its life and never acted on: "continuation chain broken", and — worse — nothing
at all tested `|<Psi_last|Psi_0>|`, the single overlap whose sign *is* the Berry phase. An adaptive
walk that stopped a fifth of the way round consequently reported `status OK`, phase `trivial (0)`,
because the endpoint checks are skipped when there is no endpoint and every *accepted* step was
continuous. Both are now checks. The general lesson is cheap to state and was expensive to learn:
**a diagnostic that cannot cause a refusal will eventually be ignored at exactly the moment it
matters.**

**Convergence failures and continuity failures need different responses.** The adaptive controller
originally rejected a step when CASSCF failed to converge, then shrank it — which is the wrong
lever, because a smaller step makes the *warm start better*, not the solver happier. The walk shrank
forever: steps with a mismatch of 0.011, far inside any threshold, rejected as unconverged. The
underlying solver issue was also mis-diagnosed for most of this project as a `conv_tol_grad`
problem. It is not: at formaldimine CAS(2,2) (130.86, 91.96) PySCF reaches `|grad[o]|` = 8.8e-06,
*inside* the 1e-5 threshold, while `dE` = 4.4e-10 refuses to fall below `conv_tol` = 1e-10 within
200 macro iterations — the optimizer crawling along a flat direction. It converges at 239 macro
iterations to the same energy to 1e-8. The fallback ladder could not fix it because every rung
relaxed the *gradient* threshold and none raised the *iteration budget*.

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
4. ~~Does the intersection really sit at tw = 90 for every rung?~~ **Tested and it holds.** Fine
   1-deg cuts in `tw` give 89.982 (CAS(6,6)), 89.992 (CAS(8,8)) and 89.958 (CAS(12,12)) against an
   assumed 90. Butadiene has no symmetry forcing this, so it is an empirical result rather than a
   consequence — and CAS(6,6), which looked like the most likely exception, is not one.
5. **Does the gap actually reach zero in this plane, for every rung?** Sub-degree cuts say yes
   where they have been run: CAS(4,4) reaches 0.098 mHa and CAS(10,10) 0.237 mHa, with cone-fit
   residuals near 1e-02. On a 1-deg grid the same rungs looked like they bottomed out near 1-3 mHa,
   so **the apparent "closest approach" is a resolution artifact** unless the cut is fine enough.
   CAS(6,6) and CAS(12,12) are being checked at sub-degree resolution; until then their positions
   carry that caveat.
