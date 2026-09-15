# Provenance of the formaldimine benchmark settings

This file records exactly which numbers were taken from upstream sources, which were
inferred, and which were changed — as required by the project brief.

## Sources

| Tag | Source |
|-----|--------|
| **NB** | `Emieeel/auto_oo`, `examples/Tutorial_Berry_phase.ipynb` (main branch) |
| **FIG** | `Emieeel/auto_oo`, `examples/three_loops_FCI.png` (copy in `docs/literature/auto_oo_three_loops_FCI.png`) |
| **PAP** | arXiv:2304.06070, *A hybrid quantum algorithm to detect conical intersections*, Sec. 6 |

Note: the brief's link points at `EmielKoridon/auto_oo`, which 404s. The package now lives at
`Emieeel/auto_oo`; the notebook path is otherwise identical. Sources archived under
`docs/literature/`.

## Taken verbatim

| Quantity | Value | Source |
|---|---|---|
| Z-matrix geometry function | `[1.498047, 1.066797, 0.987109, 118.359375] + [alpha, phi]` | NB cell 6 (identical to brief) |
| Loop parameterization | `alpha = a0 + r_a cos(2*pi*t + phase)`, `phi = p0 + r_p sin(2*pi*t + phase)` | NB cell 8 |
| CI-enclosing loop centre | `(alpha, phi) = (130, 89.9)` | NB cell 10 |
| Loop radius | `(10, 10)` degrees | NB cell 10 |
| Phase offset | `pi/20` | NB cell 10 |
| Basis | STO-3G | NB cell 14, PAP |
| Starting active space | CAS(2,2) | NB cell 14, PAP |
| Continuity reference frame | orthonormal atomic orbitals (OAO), `C_AO = S^{-1/2} C_OAO` | NB cell 3 |

## Inferred (not present in notebook code)

The notebook defines **only** the CI-enclosing loop. The two control loops appear
in FIG and are described in PAP:

> "These loops are centered around φ = 90° and α = 110° (C1), α = 130° (C×) and
>  α = 150° (C2), and all have a radius of 10°."  — PAP Sec. 6.3

| Quantity | Value | Basis for choice |
|---|---|---|
| Control loop `C1` (lower) | centre `(alpha, phi) = (110, 89.9)`, radius `(10, 10)` | PAP text + FIG |
| Control loop `C2` (upper) | centre `(alpha, phi) = (150, 89.9)`, radius `(10, 10)` | PAP text + FIG |
| Reference CI location | `phi = 90`, `alpha ~ 132` | PAP Sec. 6.3 (from FCI gap map) |

We keep `phi = 89.9` rather than PAP's nominal `90.0` for the control loops, so that all
three loops share the notebook's centre convention and the same `phase` offset. The 0.1°
offset keeps every sampled geometry away from the exactly-planar/symmetric `phi = 90`
configurations; see "Deliberate changes" below.

## Deliberate changes from the notebook

| Change | Reason |
|---|---|
| **Backend**: full CASSCF in PySCF instead of an orbital-optimized PQC (`auto_oo` + PennyLane + torch). | The brief asks for the classical CASSCF variational continuation problem. Removes the PQC ansatz error and the torch/PennyLane dependency stack. |
| **Overlaps**: genuine nonorthogonal CAS-CI overlap using `gto.intor_cross` AO overlaps between the two geometries. | NB approximates the adjacent-point overlap by identifying OAOs at different geometries (`oao_mo_coeff_l[i-1].T @ oao_mo_coeff_l[i]`), which it flags as "an approximation ... as the AO on which the MOs are defined are moving". The brief explicitly requires the genuine nonorthogonal overlap. Our endpoint overlap agrees with NB's exactly, because there the two geometries coincide. |
| **Optimization**: converged CASSCF at every point (warm-started), not a single damped Newton step. | NB deliberately uses one NR update per point to mimic a quantum-device budget. A converged solve is the natural classical analogue, is cheap here, and makes "the tracked state is *the* state-specific CASSCF ground state" well defined. Single-step mode is *not* implemented; see `docs/limitations.md`. |
| **Discretization**: `n_points` *distinct* geometries at `t_k = k/N`, plus one extra continuation step back to `t = 1` (the same geometry as `t = 0`). | NB uses `np.linspace(0, 1, n_points)`, whose first and last entries are the same geometry, so only `n_points - 1` points are distinct. Our convention makes `N` mean "number of distinct points" and still yields an exact same-geometry endpoint overlap. |
| **Estimator**: report both the gauge-fixed endpoint overlap and the gauge-invariant cyclic product of adjacent overlaps. | Required by the brief as a continuity/gauge cross-check. NB reports the endpoint overlap only. |
