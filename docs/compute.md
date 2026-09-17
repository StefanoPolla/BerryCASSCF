# Compute budget and cluster jobs

## What runs locally

The **entire formaldimine benchmark runs locally** on a laptop in a few minutes. Nothing in
Experiment A needs a cluster. Measured on an Apple M-series laptop, single process,
STO-3G, PySCF 2.14:

| Workload | Cost |
|---|---|
| One CASSCF point, CAS(2,2) | 0.061 s (mean over 123 points) |
| One CASSCF point, CAS(4,4) | 0.101 s (mean over 123 points) |
| One CASSCF point, CAS(6,6) | 0.65 s (mean over 82 points) |
| One nonorthogonal overlap, CAS(2,2) / CAS(6,6) | 0.07 ms / 0.20 ms |
| Full Berry loop, CAS(2,2), N=41 | 2-3 s |
| Full Berry loop, CAS(4,4), N=41 | 2-6 s |
| Full Berry loop, CAS(6,6), N=41 | 13-30 s |
| SA-CASSCF scan, CAS(4,4), 25x25 = 625 points | 36-48 s |
| **Whole Berry sweep** (3 loops x 3 CAS x 4 discretizations = 36 runs) | **~4 min** |
| FCI reference point (STO-3G, 1.66M determinants) | ~28 s |
| One CASSCF point, fulvene CAS(6,6)/STO-3G | 0.8 s |
| One SA-CASSCF point, ethylene CAS(8,8)/6-31G\* | 0.9 s warm, ~4 s with `strategy="best"` |
| One SA-CASSCF point, ethylene CAS(10,10)/6-31G\* | 4.5 s |
| One SA-CASSCF point, ethylene CAS(12,12)/6-31G\* | 28 s |
| One SA-CASSCF point, butadiene CAS(4,4)/6-31G\* | ~2 s |
| One SA-CASSCF point, butadiene CAS(10,10)/6-31G\* | ~14 s |
| One SA-CASSCF point, butadiene CAS(12,12)/6-31G\* | **>144 s** (853k determinants, 68 AOs) |
| Test suite (`pytest -q`) | ~20 s (185 tests) |
| Adaptive loop walk, formaldimine CAS(2,2) | 1.4 s (17 points) |
| Adaptive loop walk, ethylene CAS(8,8)/6-31G\* | ~330-710 s (23-35 points) |
| Adaptive/uniform stepping study, formaldimine (40 runs) | ~4 min |
| Adaptive/uniform stepping study, butadiene CAS(2,2) (20 runs) | ~45 min |
| Adaptive/uniform stepping study, ethylene CAS(8,8) (20 runs) | ~75 min |
| One bisection probe (2 adaptive walks), formaldimine CAS(2,2) | 1.5-55 s; slowest near the transition |
| Localization, formaldimine CAS(2,2), 3 centres | 6 min, 127k micro-iterations |
| Localization, formaldimine CAS(4,4), 3 centres | 39 min, 307k micro-iterations |
| Localization, formaldimine CAS(6,6), 3 centres | 2 h 46 min, 245k micro-iterations |
| Localization, butadiene CAS(2,2), 3 centres | 2 h 1 min, 51k micro-iterations |
| Gap-minimum search, butadiene CAS(12,12), 58 evaluations | 2 h 20 min (~145 s/evaluation) |
| Threshold calibration (reads saved runs only) | < 1 s |

Overlaps are three to four orders of magnitude cheaper than the CASSCF solve they compare,
so the exact nonorthogonal treatment costs essentially nothing.

Reproduce with:

```bash
python examples/run_formaldimine_berry.py       # ~4 min
python examples/run_formaldimine_scan.py --loops C_x C_1 C_2 --cas 4,4 --grid 25 25   # ~2 min
python examples/summarize.py
```

## What needs the cluster

Anything beyond the STO-3G/CAS(6,6) formaldimine benchmark. In particular the follow-up
system (Experiment B) and any larger-basis formaldimine check. **These have not been run**;
the templates below are prepared for submission.

Templates live in `slurm/`. Every site-specific field is marked `#### SITE ####`:
partition, account, time limit, `--cpus-per-task`, `--mem`, module loads, and the
environment-activation line. Nothing else needs editing.

### `slurm/berry_loop.sbatch` — Berry-phase continuation

*What it does.* Runs `examples/run_formaldimine_berry.py` for the requested loops, active
spaces and discretizations.

*Submit.*
```bash
mkdir -p logs
sbatch --export=ALL,LOOPS="C_x C_1 C_2",CAS="6,6",NPOINTS="25 41" slurm/berry_loop.sbatch
```

*Resources.* Single node, single task, threaded BLAS. Memory scales with the CI vector and
integral storage; `32G` is ample for anything up to about CAS(10,10) in a double-zeta basis.
Wall time: the dominant cost is `n_loops x n_CAS x sum(N) x t_point`. With a 20-minute
CAS(12,12)/cc-pVDZ point and N=25, a single loop is roughly 8 h — set `--time` accordingly.

*Expected output.* One JSON per `(loop, CAS, N)` in `results/berry/`, plus a summary table on
stdout. Existing records are skipped, so a re-submission resumes rather than recomputes.

*Where results go.* `results/berry/<loop>_cas<ne>-<ncas>_N<n>.json`, inside the repository.

### `slurm/sa_scan.sbatch` — state-averaged gap scan

*What it does.* Runs `examples/run_formaldimine_scan.py` over a 2D grid.

*Submit.*
```bash
sbatch --export=ALL,LOOPS="C_x",CAS="6,6",GRID="41 41" slurm/sa_scan.sbatch
```

*Resources.* Cost is `n_alpha x n_phi x t_SA_point`, and SA-CASSCF is roughly 2-4x a
state-specific point. A 25x25 CAS(4,4)/STO-3G grid measures 36-48 s; a 41x41 CAS(6,6)/STO-3G grid is ~1681
points and takes roughly 15-30 min. The same grid
in cc-pVDZ with CAS(8,8) is closer to 24 h; request accordingly.

*Restartability.* The scan **checkpoints after every grid row**. Re-submitting the identical
job resumes from the checkpoint, so it is safe to run under a short time limit and requeue.

*Expected output.* `results/scan/<loop>_cas<ne>-<ncas>_<na>x<np>.npz`, loadable with
`berrycasscf.scan.ScanResult.load`.

### `slurm/fulvene_scan.sbatch` and `slurm/fulvene_berry.sbatch` — the follow-up system

Two stages, submitted in order. Stage 1 locates the S1/S0 intersection in the
(exocyclic bond length, methylene torsion) plane; stage 2 places loops around whatever it found
and runs the Berry-phase continuation. Neither has been run.

```bash
mkdir -p logs
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",GRID="21 21" slurm/fulvene_scan.sbatch
# then, once that finishes:
sbatch --export=ALL,BASIS="6-31g*",CAS="6,6",NPOINTS="17 25" slurm/fulvene_berry.sbatch
```

*Estimates.* SA-CASSCF(6,6)/6-31G\* for fulvene (42 electrons, 111 basis functions) is roughly
30–90 s per point warm-started. Stage 1 (441 points) is **4–11 h**; stage 2 (~126 points)
is **1–3 h**. Stage 1 checkpoints per grid row and stage 2 skips completed records, so both
resume on re-submission.

*Outputs.* `results/fulvene/fulvene_scan_*.npz` and `results/fulvene/fulvene_<loop>_*.json`,
both inside the repository. Full detail in `docs/followup.md`.

### The third system: butadiene

All local, in three stages (`docs/butadiene.md`):

```bash
python examples/search_butadiene_ci.py              # 3 candidate planes, ~45 min
python examples/run_butadiene_study.py scan         # CAS(4,4)..(12,12), ~1.5 h
python examples/run_butadiene_study.py berry        # CAS(4,4)..(10,10), ~1 h
```

Two deliberate cost caps, both documented where they bite:

* **CAS(12,12) is scanned on one row** (`tw = 90`, 13 points) rather than the 5x13 grid. At
  >2.4 min per cold point a full grid is ~2.6 h for a single rung, and every other rung's
  two-dimensional minimum lies on that row anyway. The cost is that its `tw` is assumed.
* **The Berry ladder stops at CAS(10,10).** A state-specific CAS(12,12) solve costs minutes, so
  three loops at two discretizations would run to several hours without changing a conclusion
  already established over four rungs.

### Watching a run

Every driver writes a timestamped log to **`logs/<job>.log`** inside the repository (gitignored),
one line per completed point, carrying elapsed time, seconds per step and a projected finish:

```
18:12:41 [7/14] ( 90.000, 113.850)  gap = 0.004182 Ha  [0:21:03 elapsed, 180.4 s/step, ~0:21:03 left, ETA 18:33]
```

So a run can be followed with `tail -f logs/butadiene_scan_cas12-12.log`, and a slow job is
distinguishable from a dead one. Earlier versions logged per *row*, which meant a CAS(12,12) scan
could go silent for half an hour at a time; progress is now emitted per point.

To check a job is alive: `pgrep -fl run_butadiene_study` .

### Threading note

Both templates set `OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`OPENBLAS_NUM_THREADS` from
`SLURM_CPUS_PER_TASK`. PySCF parallelizes through threaded BLAS and its own OpenMP kernels;
oversubscription makes CASSCF slower, so do not leave these unset.

There is no MPI in this package. Do not request more than one task.

### Locating what loop transport encircles (`slurm/localize.sbatch`)

**What it does.** Bisects the loop radius about several centres until the Berry phase turns over,
then intersects the resulting ellipses. The transition radius is the elliptical distance to
whatever the loop encircles, so this *measures* the state-specific degeneracy rather than inferring
it. Method and figures: `notebooks/locating_intersections.ipynb`; motivation: `docs/findings.md` §3.

**Why it needs the cluster.** Each probe is two adaptive loop walks (two step-control settings that
must agree), each walk is 15-35 CASSCF points, and a bisection uses up to 11 probes per centre.
Measured locally at the cheapest active space, butadiene CAS(2,2): **2 h 1 min** for three centres.
A CAS(12,12) point costs >144 s against ~4 s at CAS(2,2), so the reference rung is far out of
laptop range.

**Estimated runtime**, butadiene CAS(12,12), 3 centres:

| | |
|---|---|
| points per adaptive walk | ~20 (more near the transition) |
| walks per probe | 2 |
| probes per centre | up to 11 |
| points per centre | ~440 |
| seconds per point | ~150 (state-specific CAS(12,12)/6-31G\*) |
| **per centre** | **~18 h** |
| **three centres** | **~55 h** |

That exceeds a typical 48 h wall limit, so **submit one centre per job** using `CENTRES=1` and
raising `--centres` as results accumulate, or split by hand. The driver skips a completed result
file, so re-submission is safe and restartable.

**How to submit.**

```bash
sbatch --export=ALL,SYSTEM=butadiene,CAS="12 12",CENTRES=1 slurm/localize.sbatch
```

**Expected output.** `results/localize/butadiene_cas12-12.json` — every probe with its verdict and
cost, the bracket per centre, and the triangulation with its uncertainty-normalised residual.
`logs/localize_butadiene_cas12-12.log` carries one line per probe with an ETA.

**What to look for.** Whether the measured rho agrees with the gap-scan intersection for that rung.
At CAS(2,2) it does not — 0.5385 ± 0.0843 measured against 0.2121 implied — and CAS(12,12) is where
the two methods disagree most sharply.

**Known failure modes**, both seen at CAS(2,2): a centre whose full-size loop hits the step floor
contributes nothing (choose centres so the loop does not graze the degeneracy), and a centre whose
full-size loop does not enclose the target cannot be bracketed at all. Both are reported, not
hidden.

## Scaling guidance for the follow-up system

The cost driver is the CASSCF point, and within it the active space: the CI vector scales as
`C(ncas, n_alpha) * C(ncas, n_beta)`. Practical guidance:

| Active space | Determinants | Comment |
|---|---|---|
| CAS(6,6) | 400 | trivial |
| CAS(8,8) | 4 900 | trivial |
| CAS(10,10) | 63 504 | comfortable |
| CAS(12,12) | 853 776 | fine, minutes per point |
| CAS(14,14) | 11 778 624 | needs real memory; budget carefully |

The overlap machinery costs `n_str^2` small determinants and is never the bottleneck at
these sizes.
