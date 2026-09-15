# Implementation plan and design choices

## 1. Backend feasibility check

Requirement checklist against candidate open-source backends:

| Capability needed | PySCF 2.14 | Psi4 | OpenMolcas |
|---|---|---|---|
| State-specific CASSCF | `mcscf.CASSCF` | yes | yes |
| State-averaged CASSCF | `mcscf.state_average_` | yes | yes |
| Direct access to `mo_coeff` as a NumPy array | yes | partial (wrappers) | file-based |
| Direct access to the CI vector in a documented string basis | `mc.ci`, `pyscf.fci.cistring` | limited | file-based |
| Warm start from arbitrary MO coefficients **and** CI vector | `mc.kernel(mo_coeff, ci0=...)` | partial | awkward |
| AO overlap integrals **between two different geometries** | `gto.intor_cross('int1e_ovlp', molA, molB)` | not exposed | not exposed |
| Scriptable from Python without file I/O | yes | yes | no |

**Decision: PySCF.** It is the only candidate that exposes all of (a) the CI vector in a
documented determinant-string basis, (b) MO coefficients as plain arrays, and (c) cross-geometry
AO overlap integrals. Together these are exactly what a nonorthogonal CI overlap needs.
Verified working in `.venv`: PySCF 2.14.0, NumPy 2.5.3, SciPy 1.18.1.

### What PySCF does *not* provide, and we must write

1. **A nonorthogonal CAS-CI overlap** between two CASSCF wavefunctions defined on different
   geometries (different AO *and* MO bases). PySCF's `fci.addons.transform_ci_for_orbital_rotation`
   and `fci.addons.overlap` assume both orbital sets span the *same* space and are related by a
   unitary; that is false here. These must not be used, and are not.
2. **Continuation / solution tracking** across a nuclear loop, including the OAO-referenced
   transfer of MO coefficients between geometries.
3. **Gauge fixing and Berry-phase extraction.**

## 2. The nonorthogonal overlap (the technical core)

Two CASSCF wavefunctions `|A>` and `|B>` at geometries `Ra`, `Rb`, each with `ncore`
doubly-occupied core orbitals and `ncas` active orbitals holding `nelecas` electrons.

1. Cross-geometry AO overlap `S_AO = <chi_mu(Ra) | chi_nu(Rb)>` via `gto.intor_cross`. Note this
   matrix is **not** symmetric.
2. MO-basis overlap over the occupied-capable space: `S_MO = C_A^T S_AO C_B`, split into core (`c`)
   and active (`A`) blocks.
3. A determinant pair `(I, J)` contributes `det(D_IJ)` per spin, where `D_IJ` is the overlap matrix
   between the occupied orbitals of `I` (from set A) and of `J` (from set B). The core block is
   common to every determinant, so by the Schur determinant identity

   ```
   det(D_IJ) = det(S_cc) * det( T[occ_I, occ_J] ),      T = S_AA - S_Ac S_cc^{-1} S_cA
   ```

   with `T` a single `ncas x ncas` matrix. Only *minors of `T`* depend on the determinant pair.
4. Because PySCF orders a determinant as (alpha orbitals ascending)(beta orbitals ascending), the
   spin blocks factorize with no extra sign, so with
   `M_sigma[I, J] = det(T[occ_I, occ_J])` over the alpha/beta string lists,

   ```
   <A|B> = det(S_cc)^2 * sum_{Ia,Ib} c_A[Ia,Ib] * (M_alpha c_B M_beta^T)[Ia,Ib]
   ```

   i.e. a Frobenius inner product. Cost is dominated by `n_str^2` small `(nelec_sigma x nelec_sigma)`
   determinants — negligible for the active spaces in scope.

This is exact (no orthogonality assumption anywhere), reduces to the trivial answer when
`Ra == Rb` and `C_A == C_B`, and is validated by the tests in §5.

## 3. Continuation around the loop

For each loop point `k = 1..N`:

1. Build `mol_k`; compute `S_k` and `S_k^{-1/2}`.
2. **Transfer the orbitals through the OAO frame** (taken from the auto_oo notebook): carry
   `C_oao = S_{k-1}^{1/2} C_{k-1}` and set the initial guess `C_k^guess = S_k^{-1/2} C_oao`.
   This is orthonormal at `R_k` by construction and is the identity map when `R_k = R_{k-1}`.
   A plain re-use of `C_{k-1}` is *not* orthonormal at `R_k`.
3. Solve state-specific CASSCF warm-started from `(C_k^guess, ci_{k-1})`.
4. Compute the nonorthogonal overlap `<Psi_{k-1} | Psi_k>` with §2.
5. **Gauge fix**: if that overlap is negative, flip the sign of `ci_k`. Record the raw (pre-flip)
   overlap; the flip only defines the gauge for reporting.
6. Diagnostics: `|overlap|`, CASSCF convergence flag, energy, the S0/S1 CASCI gap at the converged
   orbitals (a seam-proximity warning), and the largest change in any single MO.

After `k = N` the geometry equals `R_0` exactly, so the endpoint overlap involves no AO-basis change
and is exact.

## 4. Two Berry-phase estimators

Real orbitals and real CI coefficients are assumed throughout, so the Berry phase is Z2.

* **Endpoint estimator** `omega = <Psi_0 | Psi_N>` after gauge fixing, with `R_N = R_0`.
  Phase is `0` if `omega > 0`, `pi` if `omega < 0`.
* **Cyclic product estimator** `Pi = prod_{k=1..N} <Psi~_{k-1} | Psi~_k>` computed from the **raw**
  wavefunctions and closing the cycle with `<Psi~_{N-1} | Psi~_0>`. Every wavefunction appears once
  as bra and once as ket, so each arbitrary sign enters squared and cancels: this estimator is
  manifestly gauge-invariant and needs no sign fixing.

The two must agree in sign. Disagreement means a continuity or gauge bug, not physics, and is
reported as a hard failure. `|Pi|` and `min_k |overlap_k|` are the continuity diagnostics; magnitude
carries no phase information.

**Failure conditions** (reported, never silently absorbed):
- any CASSCF point not converged;
- `min_k |overlap_k|` below `min_abs_overlap` (default 0.8);
- endpoint and product estimators disagreeing in sign;
- `|omega|` not close to 1 at the endpoint, which means the final solve did not return to the
  same physical state.

## 5. Validation strategy

Cheap tests that do not require a full loop:
1. `<A|A> = 1` for a CASSCF wavefunction against itself.
2. Two independent CASSCF solves at the *same* geometry (different initial guesses) overlap to `+-1`.
3. Overlap is invariant under an active-space orbital rotation with a compensating CI transform.
4. Overlap against a brute-force determinant-by-determinant reference for a tiny system.
5. Symmetry: `<A|B> = <B|A>` for real wavefunctions.
6. The 2x2 linear Jahn-Teller toy model reproduces Berry phase `pi` inside / `0` outside using the
   *same* estimator code paths.

## 6. State-averaged comparator

Equal-weight 2-state SA-CASSCF on a grid in the `(alpha, phi)` plane covering the loop's bounding
box, reporting `E1 - E0`. The minimum of the gap over the grid locates the CI (or the closest
approach). Warm-started along grid rows for speed and smoothness; restartable. Minimum-energy CI
optimization is explicitly out of scope per the brief.

## 7. Module layout

```
berrycasscf/
  geometry.py      formaldimine Z-matrix, loop parameterization  (no PySCF)
  overlap.py       nonorthogonal CAS-CI overlap                  (the technical core)
  casscf.py        thin PySCF wrappers: build mol, solve, OAO transfer
  continuation.py  walk a loop, gauge fix, collect diagnostics
  berry.py         Berry-phase estimators + pass/fail verdict
  scan.py          SA-CASSCF 2D gap scan
  toy.py           2x2 linear Jahn-Teller validation model
  config.py        dataclasses for every scientific choice
  store.py         JSON/NPZ checkpointing for restartability
```
