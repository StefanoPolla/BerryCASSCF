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
* ~~Single-Newton-step mode is not implemented.~~ **Now implemented and measured.** Setting
  `CasConfig.max_cycle_macro` with `ContinuationConfig(require_converged=False,
  use_fallback=False)` reproduces the single-update regime of arXiv:2304.06070, and
  `docs/findings.md` §7 quantifies the trade-off against converged continuation. The remaining
  difference from the paper is that the update is a CASSCF macro-iteration rather than a damped
  Newton step on a combined (theta, kappa) parameter vector, so the *cost per update* is not
  identical even though the regime is.

## Numerical limitations

* **Core-block singularity.** The overlap uses a Schur complement on the core-core block and
  raises `LinAlgError` if it is singular. That can only happen if the two core spaces become
  orthogonal, which would mean the two points are not connected by continuation — but it is a
  hard failure rather than a graceful degradation.
* **Convergence thresholds interact with warm starts, and not only through the gradient.**
  Tightening `conv_tol_grad` beyond ~1e-5 makes warm-started CASSCF fail its convergence *test*
  while sitting on the correct answer to 1e-10. But the gradient is not always the binding
  criterion: measured at formaldimine CAS(2,2) (130.86, 91.96), PySCF reaches `|grad[o]|` = 8.8e-06
  — *inside* the 1e-5 threshold — while the per-iteration `dE` of 4.4e-10 will not fall below
  `conv_tol` = 1e-10 within 200 macro iterations. The optimizer is crawling along a flat direction;
  it converges at 239 iterations to the same energy to 1e-8. The fallback ladder gained a rung that
  raises the iteration budget, but **a point can still be refused for a numerical technicality
  while its wavefunction is correct**, and there is no stationarity test independent of PySCF's own
  flag.
* ~~**The fallback ladder can mask a problem.**~~ **Resolved.** A fall-through to a cold start
  breaks the continuation chain, and gauge fixing repairs an arbitrary sign but not a change of
  branch. This was a warning for most of the project; it is now a binding check
  (`ContinuationConfig.require_unbroken_chain`), calibrated against all 116 saved runs in
  `examples/calibrate_thresholds.py`.
* **Adaptive stepping has a resolution floor, and it is a knob rather than a guarantee.** A loop
  passing closer than roughly `2*pi*R*d_min/dtheta_max` — 0.7% of the loop radius with the defaults
  — cannot be walked, and the run is refused. Tightening `d_min` buys proximity at the cost of
  rejected solves. This directly bounds how precisely bisection can locate a degeneracy.
* **The step floor cannot say *why* it fired.** Shrinking failing to restore continuity means the
  loop passes through something, but "a degeneracy" and "the solver switching between nearby
  stationary solutions" produce the same signature, and the walk does not distinguish them.
* **No parallelism.** Points are computed serially. The loop is inherently sequential (each
  point warm-starts from the previous), but the SA scan is embarrassingly parallel and is not
  parallelized.

## Scope limitations

* **Three systems have been run** — formaldimine (STO-3G), ethylene and butadiene (6-31G\*).
  Fulvene is implemented and its machinery is tested, but its intersection is not reachable in the
  current rigid two-coordinate model; see `docs/followup.md`.
* **No basis-set convergence study.** Each system is calibrated within its own basis where a
  reference exists, but an intersection position is basis-dependent and none of these numbers is
  quantitative in absolute terms.
* **Butadiene has no reachable reference at all** (full valence is CAS(22,22)), so nothing there
  can be checked against truth — only against self-consistency across rungs.
* **The comparator's state averaging is 2-state equal-weight only.** Cases needing three or more
  states, or unequal weights, are supported by `ScanConfig.weights`/`nroots` but untested.

## Unresolved questions

1. **Does the CAS(2,2) Berry phase remain correct by luck or by construction?** It gets all
   three loops right while the CAS(2,2) SA scan is qualitatively wrong. The plausible
   explanation is that the state-specific S0 surface is well described even when S1 is not, so
   S0 transport is faithful. This has not been proven — a useful check would be to confirm
   that the CAS(2,2) S0 wavefunction has high overlap with the CAS(6,6) S0 wavefunction all
   the way around the loop.
2. ~~**How close can the loop pass to the seam before the answer degrades?**~~ — **answered for
   the discretization, still open for the physics.** With adaptive stepping the limit is
   *predicted* and then measured: a loop is walkable down to a closest approach of
   `eps_min ~ 2*pi*R*d_min/dtheta_max` = 0.7% of its radius (measured: closes at 0.007, floors at
   0.005 on the Jahn-Teller model, against a prediction of 0.0070). Below that the run is refused
   rather than wrong. What remains open is whether the *CASSCF solution itself* degrades before
   that limit is reached — the step floor cannot distinguish a real degeneracy from a solver
   switching branches, and at formaldimine CAS(2,2) it fired in a region whose in-CAS gap never
   falls below 321 mHa, which suggests the solver.
4. **Can a triangulated position be trusted when a symmetry forces degeneracies into mirror
   pairs?** A loop enclosing two symmetry-related degeneracies reports 0, indistinguishable from
   one enclosing none. Formaldimine's in-CAS gap map is exactly symmetric about phi = 90, so this
   is a live concern rather than a hypothetical; see `docs/todo.md` §10.
3. ~~The residual offset between the computed CI position and the literature value~~ —
   **resolved.** FCI in the same STO-3G basis puts the minimum at α = 132.6°, and the
   CAS(4,4) → CAS(6,6) → FCI sequence (130.0° → 130.9° → 132.6°) shows the offset is
   active-space truncation converging monotonically. See `docs/results.md` §2b.
