# Follow-up system: fulvene

**Status: prepared and machine-tested, not yet run.** Per the project brief, anything heavier
than the formaldimine benchmark is set up for cluster submission rather than executed locally.
The two jobs below are ready to submit; **no production fulvene result exists yet.**

## Why fulvene

Candidates considered were benzene, fulvene, a protonated Schiff-base model, and a nucleobase
model. Fulvene was chosen because:

* **The active space is a real step up but stays honest.** The full π manifold is 6 electrons
  in 6 orbitals, so CAS(6,6) is the smallest chemically defensible choice — three times the
  orbitals of the formaldimine benchmark, and 400 determinants instead of 4. Formaldimine's
  CAS(2,2) is small enough that the Berry phase is almost too easy; CAS(6,6) is where
  active-space drift and root flipping around a loop become genuine risks.
* **S1 is a real π→π\* state.** The comparator has to do actual work, which is exactly the
  aspect that broke at CAS(2,2) for formaldimine.
* **The driving coordinates are documented.** The exocyclic C=C stretch and the methylene
  torsion are the standard coordinates in the fulvene photochemistry literature, so the scan
  plane is defensible rather than arbitrary.
* **It stays well conditioned.** A closed-shell hydrocarbon with no low-lying σ manifold, so the
  active space remains well defined around a loop — a precondition for continuation.

Benzene was the main alternative. Its S1/S0 seam (the "prefulvene" channel) is equally well
studied, but reaching it requires out-of-plane ring puckering, so a *two*-parameter plane
containing the intersection is harder to justify. Fulvene's two coordinates are natural.

## Geometry and coordinates

Reference S0 geometry: RHF/6-31G optimized, planar (C2v), stored as Cartesians in
`berrycasscf/fulvene.py`. Exocyclic C1=C6 bond 1.3293 Å.

Two parameters distort it; everything else is frozen:

| Parameter | Meaning | Reference | Search range |
|---|---|---|---|
| `r` | exocyclic C1=C6 bond length (Å) | 1.3293 | 1.20 – 1.60 |
| `theta` | methylene torsion about the C1–C6 axis (deg) | 0 (planar) | 0 – 90 |

`r` translates C6 and its two hydrogens rigidly along the C1→C6 axis; `theta` rotates only the
two methylene hydrogens about that axis. Both operations are rigid — verified in
`tests/test_fulvene.py` that no bond length changes and the ring never moves.

## Two-stage plan

**The CI position is not assumed anywhere.** Stage 1 finds it; stage 2 builds loops around
whatever was found. This is why there are two jobs rather than one.

### Stage 1 — locate the intersection

```bash
mkdir -p logs
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",GRID="21 21" slurm/fulvene_scan.sbatch
```

Equal-weight two-state SA-CASSCF(6,6)/6-31G\* on a 21×21 grid over the (`r`, `theta`) plane
(441 points).

*Expected runtime.* SA-CASSCF(6,6)/6-31G\* for fulvene (42 electrons, 111 basis functions) is
roughly 30–90 s per point warm-started, so **4–11 h**. The 12 h limit in the template has
margin; the scan checkpoints after every grid row, so re-submitting the identical job resumes.

*Output.* `results/fulvene/fulvene_scan_cas6-6_6-31gs_21x21.npz`, plus the CI candidate
printed at the end of the job log.

*What to check before proceeding.* The reported minimum gap should be ≲ 5 mHa. If it is much
larger, the seam is probably outside the search box or outside this plane — widen the region
(`SEARCH_REGION` in `berrycasscf/fulvene.py`) or enlarge the active space before running stage 2.
The stage-2 driver prints this warning itself.

### Stage 2 — Berry phase on loops around it

```bash
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",NPOINTS="17 25" slurm/fulvene_berry.sbatch
```

Picks up the newest stage-1 scan automatically (or pass `SCAN=path/to/....npz`), places three
loops — `F_x` centred on the located CI, `F_1` and `F_2` displaced by ±2 radii along `r` so they
enclose nothing — and runs the continuation on each at N = 17 and 25.

Default loop radii are 0.06 Å and 15°; adjust with `RADIUS="dr dtheta"`. They should be small
enough to stay inside the region scanned in stage 1 and large enough that the loop is not
dominated by the seam's immediate neighbourhood.

*Expected runtime.* 3 loops × 2 discretizations × ~21 points × 30–90 s ≈ **1–3 h**.

*Output.* `results/fulvene/fulvene_<loop>_cas6-6_6-31gs_N<n>.json`, in the same format as the
formaldimine records, readable by `berrycasscf.report.load_berry_records`.

### What success looks like

`F_x` gives a non-trivial (π) Berry phase; `F_1` and `F_2` give 0; all runs report status `OK`;
and the phases do not change between N = 17 and N = 25. That is the same stability criterion
applied to formaldimine in `docs/results.md`.

## What has been verified locally

The machinery, not the experiment:

* the geometry builder is exact — bond lengths preserved to 1e-10 under both coordinates, ring
  frozen, no atomic collisions anywhere in the search region (18 tests in `tests/test_fulvene.py`);
* a CASSCF(6,6)/STO-3G point converges (0.8 s), `ncore = 18`, CI vector 20×20;
* the OAO orbital transfer between two fulvene geometries is orthonormal to 1.7e-14;
* a warm-started second point converges, and the **factorized nonorthogonal overlap agrees with
  the brute-force determinant-pair reference** with 18 core orbitals — the important check that
  the Schur-complement core elimination scales beyond formaldimine's 7 core orbitals;
* the overlap between two neighbouring fulvene geometries came out at −0.974, i.e. the arbitrary
  CASSCF sign appears here too and the gauge fixing is doing real work.

## A cheaper option

If you would rather see a fulvene result before committing cluster time, STO-3G/CAS(6,6) is
about 0.8 s per point, which puts a full stage-1 + stage-2 pilot at roughly 10–15 minutes on a
laptop:

```bash
python examples/run_fulvene_scan.py --basis sto-3g --cas 6,6 --grid 15 15
python examples/run_fulvene_berry.py --scan results/fulvene/fulvene_scan_cas6-6_sto-3g_15x15.npz \
       --basis sto-3g --cas 6,6 --npoints 17
```

STO-3G is not adequate for a quantitative S1/S0 intersection in a π system, so treat this as a
pipeline shakedown rather than a result. It was not run, because the brief directs the follow-up
system to the cluster.
