# Results: formaldimine benchmark

All numbers below are reproducible from this repository; see `README.md` for the commands.
Basis STO-3G, C1 symmetry, PySCF 2.14.0. The narrative version with figures is
`notebooks/results.ipynb`.

## The three loops

Circles of radius 10° in the (α, φ) plane, centred at φ = 89.9°:

| Loop | centre α | encloses the CI? |
|---|---|---|
| `C_1` | 110° | no |
| `C_x` | 130° | **yes** |
| `C_2` | 150° | no |

Provenance in `docs/provenance.md`. Reference CI position from FCI in arXiv:2304.06070:
(α, φ) ≈ (132°, 90°).

## 1. Berry phase (state-specific CASSCF continuation)

36 runs: 3 loops × 3 active spaces × 4 discretizations (N = 9, 17, 25, 41).
**Every run succeeded and every run gave the expected answer.**

| Loop | CAS(2,2) | CAS(4,4) | CAS(6,6) |
|---|---|---|---|
| `C_x` | **π** (all N) | **π** (all N) | **π** (all N) |
| `C_1` | 0 (all N) | 0 (all N) | 0 (all N) |
| `C_2` | 0 (all N) | 0 (all N) | 0 (all N) |

Representative numbers (N = 41):

| Loop | CAS | product Π | endpoint ω | min \|adjacent overlap\| |
|---|---|---|---|---|
| `C_x` | CAS(2,2) | −0.8736 | −1.000000 | 0.9949 |
| `C_x` | CAS(4,4) | −0.8675 | −1.000000 | 0.9911 |
| `C_x` | CAS(6,6) | −0.8686 | −1.000000 | 0.9918 |
| `C_1` | CAS(4,4) | +0.9263 | +1.000000 | 0.9916 |
| `C_2` | CAS(4,4) | +0.9777 | +1.000000 | 0.9976 |

The endpoint estimator comes out at exactly ±1 to six decimals in all 36 runs, as it must:
the closing geometry is *identical* to the starting geometry, so that overlap involves no
change of AO basis. Refining N drives |Π| → 1 monotonically but never threatens its sign;
**the Z2 answer is already correct at N = 9**, which costs about one second at CAS(2,2).

### Stability criteria

A (loop, CAS) combination is called **stable** when:

1. every traversal has `status == "OK"` — all CASSCF points converged, continuity held, and the
   endpoint returned to the same physical state;
2. at least two discretizations N were run;
3. all of them give the same Berry phase;
4. the product and endpoint estimators agree in sign in every run.

All nine (loop, CAS) combinations are stable by this criterion. It is implemented in
`berrycasscf.report.stability_verdict` and printed by `examples/summarize.py`.

## 2. State-averaged CASSCF comparator

Equal-weight two-state SA-CASSCF, 25 × 25 grids over each loop's bounding box (+1° margin),
resolution 0.917° in both coordinates. All scans use the cold strategy
(`ScanConfig.strategy = "cold"`); see `docs/active_space.md` for why, and for the audit
confirming that redoing them cold left every minimum position and every conclusion unchanged.

| Loop | CAS | min gap (Ha) | at (α, φ) | inside the loop? | agrees with Berry phase? |
|---|---|---|---|---|---|
| `C_1` | CAS(2,2) | 0.064351 | (121.00, 78.90) | no | yes |
| `C_1` | CAS(4,4) | 0.018371 | (121.00, 89.90) | no | yes |
| `C_1` | CAS(6,6) | 0.019873 | (121.00, 89.90) | no | yes |
| `C_x` | CAS(2,2) | 0.020570 | (141.00, 89.90) | no | **NO** |
| `C_x` | CAS(4,4) | **0.000881** | (130.00, 89.90) | **yes** | yes |
| `C_x` | CAS(6,6) | **0.000905** | (130.92, 89.90) | **yes** | yes |
| `C_2` | CAS(2,2) | **0.000341** | (150.92, 89.90) | **yes** | **NO** |
| `C_2` | CAS(4,4) | 0.016808 | (139.00, 89.90) | no | yes |
| `C_2` | CAS(6,6) | 0.015137 | (139.00, 89.90) | no | yes |

For the control loops at CAS(4,4)/CAS(6,6) the minimum sits on the *box edge nearest the seam*
with a gap roughly twenty times larger than at the CI — the signature of a slope, not a cone.

**CI location.** CAS(4,4) gives (130.00°, 89.90°) and CAS(6,6) gives (130.92°, 89.90°), both
with a residual gap below 1 mHa at 0.917° resolution.

## 2b. Calibration against FCI in the same basis

Because STO-3G is small, the *exact* answer is computable: FCI in the full 13-orbital space
(1.66 M determinants, ~28 s per point, `examples/run_fci_reference.py`). This is the same
quantity the gap map of arXiv:2304.06070 Fig. 1a reports, so it calibrates the truncated
active spaces against the paper directly.

Scanning α at φ = 89.9°:

| α | FCI gap (Ha) |
|---|---|
| 128.750 | 0.007693 |
| 130.000 | 0.005274 |
| 131.250 | 0.002841 |
| **132.500** | **0.000442** |
| 133.750 | 0.002098 |
| 135.000 | 0.004566 |

The gap is linear in α on both sides of the minimum — the signature of a cone, not an avoided
crossing. Parabolic interpolation around the minimum gives **α× = 132.6°**, in agreement with
the paper's α× ≈ 132°.

Convergence of the CI position with active space:

| method | α of the gap minimum | error vs FCI |
|---|---|---|
| SA-CASSCF(2,2) | no minimum in the region (spurious one at 150.9°) | qualitatively wrong |
| SA-CASSCF(4,4) | 130.00° | −2.6° |
| SA-CASSCF(6,6) | 130.92° | −1.7° |
| **FCI** | **132.61°** | — |

So the ~2° offset noted above **is** active-space truncation: the estimate moves monotonically
toward the exact value as the space grows. This settles what was previously an open question.
Note that the Berry phase is unaffected — a 2° error in the CI position is irrelevant to a loop
of radius 10° that encloses it either way, which is precisely the robustness that makes the
topological question easier than the geometric one.

## 3. Where the methods disagree, and why

**They disagree at CAS(2,2), and the comparator is the one that is wrong.** At that active
space the SA scan gets two of the three loops backwards: it reports a 0.34 mHa near-degeneracy
inside the *trivial* loop `C_2`, and finds nothing inside the true CI loop `C_x`.

The cause is understood. Recomputing the gap at the spurious minimum (α, φ) = (150.92, 89.90)
with larger active spaces (`examples/diagnose_cas22_artifact.py`):

| active space | gap (Ha) |
|---|---|
| CAS(2,2) | 0.000345 ← claims a conical intersection |
| CAS(4,4) | 0.038347 |
| CAS(6,6) | 0.037863 |
| CAS(8,8) | 0.039675 ← converged: no degeneracy there at all |

Three successively larger active spaces agree the gap is ≈0.04 Ha. The CAS(2,2) degeneracy is
a converged artifact of the two-orbital active space, not a physical feature.

**Why the two methods have different active-space requirements.** The Berry-phase workflow is
state-specific: it only ever tracks S0, and needs the active space to describe *one* state well
enough to follow it continuously around the loop. The comparator must place S1 correctly, and
S1 at these geometries involves orbitals a two-orbital active space does not contain. Averaging
two states over a CAS that cannot represent one of them yields a compromised orbital set and a
meaningless second root.

This asymmetry — the comparator needing a larger active space than the Berry-phase calculation
— is the behaviour the project brief anticipated.

## 4. Minimal reliable active space

**CAS(4,4)**, by the criterion *"the smallest tested active space at which both methods are
stable and mutually consistent on all three loops"*.

Broken down:

* Berry phase alone: **CAS(2,2)** suffices and is stable across all four discretizations.
* SA comparator alone: **CAS(4,4)** is the smallest that is even qualitatively correct.
* Both, mutually consistent: **CAS(4,4)**.

CAS(6,6) changes no conclusion, which is the evidence that CAS(4,4) is converged for this
question rather than merely lucky.
