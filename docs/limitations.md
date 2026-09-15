# Known limitations and unresolved questions

## What the Berry phase does and does not tell you

* It is a **Z2** quantity under the real-wavefunction assumption: trivial (0) or non-trivial (π).
  A non-trivial phase means the loop encloses an **odd** number of conical intersections. Two CIs
  inside the same loop give a trivial phase and are indistinguishable from none. The comparator
  scan is what disambiguates this, and is the reason both workflows exist.
* It says nothing about **where** inside the loop the CI is, nor how far the loop passes from it.
* `|Π|` and `|ω|` are **continuity diagnostics only**. Magnitude carries no phase information.
  A value of |Π| = 0.5 with a clean sign is a perfectly good result; it means the discretization
  is coarse, not that the answer is half-confident.

## Methodological limitations

* **Real wavefunctions assumed throughout.** Orbitals and CI coefficients are real, and PySCF is
  run in C1 with real integrals. Complex or spin-orbit-coupled cases are out of scope; the Z2
  reduction would not apply.
* **Only the lowest root is tracked.** The workflow follows the state-specific CASSCF ground
  state. Berry phases of excited states are not implemented.
* **The seam is not located by the Berry-phase workflow.** It answers a yes/no question about a
  loop, by design.
* **No minimum-energy CI optimization.** Explicitly out of scope per the brief; the comparator
  is a scan. Consequently the reported CI position is limited by grid resolution (0.917° here)
  and the reported "minimum gap" is a grid minimum, not the true seam minimum.
* **The in-CAS CASCI gap diagnostic is not the physical S0/S1 gap.** It is computed at the
  state-specific orbitals within the tracked active space and exists only to warn about root
  flipping during continuation. In CAS(2,2) it is around 1 Ha because the second root there is a
  doubly-excited configuration. Use the SA scan for physical gaps.
* **Single-Newton-step mode is not implemented.** The auto_oo tutorial deliberately takes one
  damped Newton step per loop point to mimic a quantum-device budget; this package converges
  CASSCF at every point instead. That makes "the tracked state is *the* state-specific CASSCF
  ground state" well defined, but it does not reproduce the paper's convergence analysis, and
  the package cannot currently be used to study the single-step regime.

## Numerical limitations

* **Core-block singularity.** The overlap uses a Schur complement on the core-core block and
  raises `LinAlgError` if it is singular. That can only happen if the two core spaces become
  orthogonal, which would mean the two points are not connected by continuation — but it is a
  hard failure rather than a graceful degradation.
* **Convergence thresholds interact with warm starts.** Tightening `conv_tol_grad` beyond ~1e-5
  makes warm-started CASSCF fail its convergence *test* while sitting on the correct answer to
  1e-10 (see `docs/progress.md`). The default was chosen with this in mind; if you tighten it,
  expect spurious FAILED verdicts.
* **The fallback ladder can mask a problem.** If a point falls through to a cold start, the
  continuation chain is broken there. This is warned about loudly and the overlap checks still
  apply, but a cold-started point that happens to land on the same branch will pass silently.
* **No parallelism.** Points are computed serially. The loop is inherently sequential (each
  point warm-starts from the previous), but the SA scan is embarrassingly parallel and is not
  parallelized.

## Scope limitations

* **Only formaldimine has been run.** The follow-up system (fulvene) is implemented and its
  machinery is tested, but no production result exists; see `docs/followup.md`.
* **STO-3G only.** No basis-set convergence study. Within STO-3G the results are calibrated
  against FCI (`docs/results.md` §2b), but the CI position is basis-dependent and a minimal
  basis is not quantitative in absolute terms.
* **The comparator's state averaging is 2-state equal-weight only.** Cases needing three or more
  states, or unequal weights, are supported by `ScanConfig.weights`/`nroots` but untested.

## Unresolved questions

1. **Does the CAS(2,2) Berry phase remain correct by luck or by construction?** It gets all
   three loops right while the CAS(2,2) SA scan is qualitatively wrong. The plausible
   explanation is that the state-specific S0 surface is well described even when S1 is not, so
   S0 transport is faithful. This has not been proven — a useful check would be to confirm
   that the CAS(2,2) S0 wavefunction has high overlap with the CAS(6,6) S0 wavefunction all
   the way around the loop.
2. **How close can the loop pass to the seam before the answer degrades?** The toy model
   degrades gracefully down to |Π| ≈ 0.5 at a 1% seam offset, but no such study was run at
   CASSCF level.
3. ~~The residual offset between the computed CI position and the literature value~~ —
   **resolved.** FCI in the same STO-3G basis puts the minimum at α = 132.6°, and the
   CAS(4,4) → CAS(6,6) → FCI sequence (130.0° → 130.9° → 132.6°) shows the offset is
   active-space truncation converging monotonically. See `docs/results.md` §2b.
