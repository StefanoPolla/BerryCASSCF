# BerryCASSCF

Detecting conical intersections with a **variational Berry phase at CASSCF level**, with an
independent **state-averaged CASSCF comparator**.

Given a closed loop in nuclear coordinate space, the package answers one question two different
ways:

> *Does this loop enclose a conical intersection?*

1. **Berry phase.** Track the state-specific CASSCF ground state S0 around the loop by
   continuation, and read the Z2 Berry phase off the *sign* of the initial–final wavefunction
   overlap. Non-trivial (π) ⇒ the loop encloses (an odd number of) CIs.
2. **State-averaged comparator.** Resolve S0 and S1 with equal-weight SA-CASSCF on a grid over
   the enclosed region and look for the gap closing.

The benchmark is formaldimine (H₂C=NH) in STO-3G, following arXiv:2304.06070 and the `auto_oo`
Berry-phase tutorial.

**It runs in minutes on a laptop.** The full benchmark sweep — 36 Berry-phase traversals — takes
about four minutes.

## Results in one table

| Loop (radius 10°, φ=89.9°) | Berry phase | SA-CASSCF comparator |
|---|---|---|
| `C_x`, centre α=130° — encloses the CI | **π** at CAS(2,2)/(4,4)/(6,6) | CI found inside (CAS(4,4)+) |
| `C_1`, centre α=110° — control | **0** at CAS(2,2)/(4,4)/(6,6) | no CI inside (CAS(4,4)+) |
| `C_2`, centre α=150° — control | **0** at CAS(2,2)/(4,4)/(6,6) | no CI inside (CAS(4,4)+) |

All 36 Berry-phase runs agree across four discretizations (N = 9, 17, 25, 41), and the endpoint
overlap comes out at exactly ±1.000000 in every one.

The two methods **disagree at CAS(2,2)**, where the comparator invents a conical intersection
inside the trivial loop `C_2` and misses the real one. That is a genuine finding, not a bug; it
is diagnosed in [docs/results.md](docs/results.md) §3. The smallest active space at which both
methods are stable and mutually consistent is **CAS(4,4)**.

Because STO-3G is small, the results are calibrated against **exact FCI in the same basis**,
which puts the intersection at α× = 132.6° — matching the reference value ≈132° — and confirms
that the ~2° offset of the truncated active spaces is truncation error converging monotonically
(CAS(4,4) 130.0° → CAS(6,6) 130.9° → FCI 132.6°). The Berry phase is untouched by that offset:
**the topological question is far more forgiving than the geometric one.**

Full narrative with figures: [notebooks/results.ipynb](notebooks/results.ipynb).

## Installation

Python ≥ 3.10. Use a project-local virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e ".[notebook,dev]"
```

Dependencies are PySCF, NumPy and SciPy (plus matplotlib/jupyter for the notebook).
Developed against PySCF 2.14.0, NumPy 2.5.3, SciPy 1.18.1.

## Backend

**PySCF**, chosen after a capability check against Psi4 and OpenMolcas. It is the only candidate
exposing all three of: CI vectors in a documented determinant-string basis, MO coefficients as
plain arrays, and — decisively — `gto.intor_cross`, which gives AO overlap integrals **between
two different geometries**. Reasoning in [docs/plan.md](docs/plan.md) §1.

## Running the benchmark

```bash
# Berry-phase sweep: 3 loops x 3 active spaces x 4 discretizations   (~4 min)
python examples/run_formaldimine_berry.py

# State-averaged gap scans                                          (~2 min)
python examples/run_formaldimine_scan.py --loops C_x C_1 C_2 --cas 4,4 --grid 25 25

# Summary table and stability verdicts
python examples/summarize.py

# Why CAS(2,2) disagrees                                            (~1 min)
python examples/diagnose_cas22_artifact.py

# Exact FCI reference in the same basis                             (~10 min)
python examples/run_fci_reference.py --loop C_x --line

# Rebuild the results notebook from saved results                   (seconds)
python notebooks/build_notebook.py
jupyter nbconvert --to notebook --execute --inplace notebooks/results.ipynb
```

### Active-space convergence study (ethylene)

How large an active space does each method actually need? Same two workflows, one system, a
ladder of active spaces, and — critically — the **same loops throughout**:

```bash
python examples/run_active_space_study.py scan     # SA-CASSCF gap maps, CAS(2,2) .. (12,12)
python examples/run_active_space_study.py berry    # Berry phase on three fixed loops
python notebooks/build_active_space_notebook.py
jupyter nbconvert --to notebook --execute --inplace notebooks/active_space.ipynb
```

Results and interpretation: [notebooks/active_space.ipynb](notebooks/active_space.ipynb) and
[docs/active_space.md](docs/active_space.md).

Every driver **skips work already saved** under `results/`, so all of them are restartable; the
gap scans additionally checkpoint after each grid row.

### Minimal example

```python
from berrycasscf import LOOP_CI, CasConfig, run_loop

result, traversal = run_loop(LOOP_CI.with_n_points(25), cas=CasConfig(ncas=2, nelecas=2))
print(result.summary())
# Loop C_x  centre=(130.0, 89.9)  radius=(10.0, 10.0)  N=25  CAS(2,2)/sto-3g
#   product estimator   Pi = -0.801136
#   endpoint estimator  w  = -1.000000
#   min |adjacent overlap| = 0.9863
#   min CASCI S0/S1 gap    = 0.3226 Ha
#   BERRY PHASE: non-trivial (pi)   [OK]
```

The last digits of `Pi` wander by ~1e-5 between runs, since each point is an independently
converged CASSCF solve. The sign — the entire physical content — does not.

## Expected output

* `results/berry/<loop>_cas<ne>-<ncas>_N<n>.json` — Berry-phase verdict plus per-point
  diagnostics (energy, convergence, raw overlap, gauge sign, solver strategy, timings).
* `results/scan/<loop>_cas<ne>-<ncas>_<na>x<np>.npz` — gap grids, via `ScanResult.load`.
* `results/diagnostics/cas22_artifact.json` — the active-space diagnosis.

## How it works

Two things make this more than a sequence of independent CASSCF calculations.

**1. Genuine nonorthogonal overlaps.** Neighbouring loop points have different AO bases (the
nuclei moved) and different MO bases (the orbitals were re-optimized), so their overlap is a
nonorthogonal CI overlap. PySCF's `fci.addons.transform_ci_for_orbital_rotation` and friends
assume a common orbital space and are **not** valid here; they are not used for this purpose.
Instead the common doubly-occupied core is eliminated analytically with a Schur complement,
leaving `det(S_cc)²` times a Frobenius contraction over minors of a single `ncas × ncas` matrix.
Exact, and ~0.1 ms per overlap — three to four orders of magnitude cheaper than the CASSCF solve
it compares. See [berrycasscf/overlap.py](berrycasscf/overlap.py).

**2. Continuation with explicit gauge fixing.** Orbitals are carried between geometries through
the orthonormal-atomic-orbital frame, where they stay orthonormal (transferred guesses are
orthonormal to ~1e-15; naive re-use of the previous AO coefficients is off by ~4e-2). CASSCF
returns an **arbitrary overall sign** at each point — verified: two independent solves at the
same geometry overlap to exactly −1 — so the sign is fixed sequentially against the previous
point rather than trusted.

Two estimators are then reported and cross-checked: the gauge-invariant cyclic product of
adjacent overlaps, and the exact same-geometry endpoint overlap. They must agree in sign; a
disagreement is reported as a failure, never as a result.

## Package layout

```
berrycasscf/
  geometry.py      Loop parameterization, formaldimine Z-matrix, the benchmark loops
  overlap.py       nonorthogonal CAS-CI overlap  (the technical core)
  casscf.py        PySCF wrappers: build, solve, OAO orbital transfer, spin adaptation
  continuation.py  walk a loop, gauge fix, fallback ladder, collect diagnostics
  berry.py         the two estimators and the pass/fail verdict
  scan.py          SA-CASSCF and FCI gap scans
  ethylene.py      twisted-pyramidalized CI; the active-space convergence study
  fulvene.py       prepared for the cluster: reference geometry and its two coordinates
  toy.py           2x2 linear Jahn-Teller validation model
  config.py        every scientific choice, as dataclass fields
  store.py         JSON/NPZ records
  report.py        result tables and stability verdicts
examples/          experiment drivers (see "Running the benchmark")
notebooks/         results.ipynb (formaldimine), active_space.ipynb (ethylene ladder)
slurm/             batch templates for the cluster
docs/              plan, provenance, results, limitations, compute, progress, follow-up
tests/             pytest suite
```

## Tests

```bash
pytest -q            # ~15 s
pytest -q -m "not slow"   # instant: geometry + toy model only
```

58 tests, ~11 s, covering the overlap engine (including agreement with a brute-force
determinant-pair reference to 4e-16, for both molecules), the Jahn-Teller model, the loop
parameterization, the fulvene geometry builder, and end-to-end smoke runs of both workflows.

## Compute and the follow-up system

The formaldimine benchmark is a laptop calculation; see [docs/compute.md](docs/compute.md) for
measured timings. Anything heavier is meant for a cluster, with templates in [slurm/](slurm/)
whose site-specific fields are marked `#### SITE ####`.

The follow-up system is **fulvene**, implemented and machine-tested but **not yet run** — per
the brief it is prepared for cluster submission rather than executed locally. It is a two-stage
job: stage 1 locates the S1/S0 intersection in the (exocyclic bond length, methylene torsion)
plane, stage 2 places loops around whatever was found and runs the continuation. The CI position
is not assumed anywhere. See [docs/followup.md](docs/followup.md).

```bash
mkdir -p logs
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",GRID="21 21" slurm/fulvene_scan.sbatch
# once that finishes:
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",NPOINTS="17 25" slurm/fulvene_berry.sbatch
```

## Documentation

| File | Contents |
|---|---|
| [docs/plan.md](docs/plan.md) | backend feasibility check, design of the overlap and continuation |
| [docs/provenance.md](docs/provenance.md) | which settings came from where, and what was changed |
| [docs/results.md](docs/results.md) | the formaldimine numbers and their interpretation |
| [docs/limitations.md](docs/limitations.md) | what this does not do, and open questions |
| [docs/compute.md](docs/compute.md) | measured costs, cluster jobs, scaling guidance |
| [docs/active_space.md](docs/active_space.md) | how large an active space each method needs |
| [docs/followup.md](docs/followup.md) | fulvene: status and how to run it |
| [docs/progress.md](docs/progress.md) | dated progress log |

## References

- Y. Kim *et al.* / E. Koridon *et al.*, *A hybrid quantum algorithm to detect conical
  intersections*, [arXiv:2304.06070](https://arxiv.org/abs/2304.06070), Quantum 8, 1259 (2024).
- `auto_oo`: <https://github.com/Emieeel/auto_oo> — `examples/Tutorial_Berry_phase.ipynb` is the
  source of the geometry parameterization and the CI-enclosing loop.

Archived copies of the open-access sources are in [docs/literature/](docs/literature/).
