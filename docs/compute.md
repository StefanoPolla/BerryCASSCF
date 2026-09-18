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
| One CASSCF point, butadiene CAS(12,12)/6-31G\* — **ALICE node, 8 threads** | 349 s |
| the same, 1 thread | 876 s (so threading is worth 2.5x) |
| One SA-CASSCF point, butadiene CAS(12,12)/6-31G\* — ALICE node | 144 s |
| One CASSCF point, butadiene CAS(14,14)/6-31G\* — ALICE node | **13 123 s, and not converged** |
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

## What needs the cluster

Anything beyond the STO-3G/CAS(6,6) formaldimine benchmark: the large butadiene active
spaces, and above all the localization runs, which are hours per centre at a cheap active
space and days at the reference one.

`slurm/` holds five **runnable** job scripts with concrete values for ALICE (Leiden
University), plus four older generic templates that still carry `#### SITE ####`
placeholders for workloads that have never needed a cluster (the formaldimine workflows and
the abandoned fulvene follow-up). In the runnable ones, every site-specific line is marked
`# SITE:` — partition, account, module stack, venv path — and those are the only lines
another cluster needs changed.

### Setting up on ALICE

```bash
git clone git@github.com:StefanoPolla/BerryCASSCF.git ~/git_repos/BerryCASSCF
cd ~/git_repos/BerryCASSCF
ml purge && ml load ALICE/default && ml load Python/3.11.5-GCCcore-13.2.0
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"
cd slurm && sbatch --test-only smoke.job && sbatch smoke.job
```

Thereafter: `git pull --ff-only` on the cluster, then `cd slurm && sbatch <job>`. Code goes
up through git, never rsync; results come back down with rsync (below). Submit from inside
`slurm/`, because `--output=./out/...` is relative to the submission directory.

The cluster stack is Python 3.11.5 / PySCF 2.14.0 / NumPy 2.4.6 / SciPy 1.17.1 against the
laptop's 3.12.10 / 2.14.0 / 2.5.3 / 1.18.1. The **only** difference the two have shown is a
fitted parameter whose true value is zero coming back as 8.5e-06 instead of 0.0 at a kink
(`tests/test_refine.py`); fitted positions agree exactly. `smoke.job` is what establishes
that on a new machine, and it runs the test suite before it measures anything.

### Partition choice

`cpu_lorentz` for everything real: it is private, allows 7 days, and was carrying 5 queued
jobs against `cpu-skylake`'s 227 and `cpu-zen4`'s 199 when these runs were submitted, so
eight tasks all started within seconds. It needs **both** `--partition=cpu_lorentz` and
`--account=cpu_lorentz`; ALICE rejects the submission with only the first. `cpu-short` (4 h
cap, its own small queue) is right for `smoke.job` and for anything else under four hours.
`--time` is mandatory on ALICE — there is no default and a job without it is rejected.

## The ALICE jobs

Everything below lands in the repository working directory on the cluster
(`~/git_repos/BerryCASSCF`): results under `results/<system>/`, live progress under
`logs/<name>.log`, SLURM's own streams under `slurm/out/` and `slurm/err/`. Nothing writes
outside the repo, so pulling results back is one rsync.

**All five jobs skip work already saved**, printing `SKIP: <path>` and exiting 0, so a
resubmission runs only what is missing and a job that finds everything present costs
seconds. That also means `sacct` reporting `FAILED` means a real failure.

### `smoke.job` — validate and measure (cpu-short, 40 min)

Runs the fast test suite, then times one state-specific and one state-averaged CASSCF solve
at CAS(2,2) and CAS(12,12), threaded and single-threaded (`examples/bench_point.py`). Run it
first on any new machine: every estimate below is a per-point cost times a point count, and
the per-point cost is what changes between machines.

### `berry_cas12.job` — the missing controls at CAS(12,12) (array of 3)

*What.* Butadiene loop transport at CAS(12,12) on all three loops, at N = 13, 21 and 31.

*Why.* `docs/findings.md` §3 rests on the enclosing loop `B_x` returning pi at this rung
while the gap scan finds nothing below 1.5 mHa inside it. Two things that claim needed:

* **the control loops at the same rung.** Every other rung has them; CAS(12,12) had only
  `B_x`, so its pi had nothing to be contrasted against. A control returning pi would make
  the result a solver artifact rather than a topological statement.
* **a third discretization.** The stability criterion is "at least two discretizations, all
  passing, same phase"; `B_x` had exactly two, which is the minimum the criterion accepts.

*Cost.* 97.6 s per point measured over the saved N=21 run, so 13+21+31 = 65 points is ~1.8 h
per loop; `B_x` only needs N=31 and is ~0.9 h. Asked for 16 h.

*Output.* `results/butadiene/butadiene_{B_x,B_1,B_2}_cas12-12_N{13,21,31}.json`.

### `localize.job` and `localize_cas12.job` — where the degeneracy actually is

*What.* Bisects the loop radius about several centres until the Berry phase turns over, then
intersects the resulting ellipses. The transition radius is the elliptical distance to
whatever the loop encircles, so this **measures** the state-specific degeneracy rather than
inferring it. Method and figures: `notebooks/locating_intersections.ipynb`; motivation:
`docs/findings.md` §3.

*Why a ladder.* The state-specific/state-averaged offset has been measured for butadiene at
CAS(2,2) only (rho = 0.5385 +/- 0.0843, against the gap scan's 0.2121 — decisively outside).
One rung is an observation. `localize.job` runs CAS(4,4), (6,6), (8,8) and (10,10), one per
array task; `localize_cas12.job` runs the reference rung, where the two methods disagree
most sharply and which is far out of laptop range.

*Why it is expensive.* Each probe is two adaptive loop walks (two step-control settings that
must agree), each walk is 15-35 CASSCF points, and a bisection uses up to 9 probes per
centre. Measured at CAS(2,2) on a laptop: **2 h 1 min** for three centres, 50 557
micro-iterations.

| rung | s/point (laptop, state-specific) | scaled from the CAS(2,2) run | asked |
|---|---|---|---|
| CAS(4,4) | 14.0 | ~8 h | 36 h |
| CAS(6,6) | 3.8 | ~2 h | 36 h |
| CAS(8,8) | 5.3 | ~3 h | 36 h |
| CAS(10,10) | 12.4 | ~7 h | 36 h |
| CAS(12,12) | 97.6 | ~57 h | 5 days |

The requests are deliberately generous. A probe that lands *on* the degeneracy costs several
times one that does not — the CAS(2,2) run's slowest probe took 40 minutes against a 3-minute
median — and walltime on a 7-day partition costs only a little backfill priority.

*Restartability.* The record is written **after every centre**, and a resubmission resumes
from the first unfinished one (`resume_bisections`). Before that, the file was written only
at the end, so a walltime kill lost days; the older advice in this file to submit one centre
per job is obsolete. A saved centre is reused only while it matches the requested centre list
position by position — changing a centre invalidates it and everything after it, because the
ellipses it measured belong to a different construction.

*Centres for the reference rung.* `localize_cas12.job` uses (90, 101.85), (90, 90) and
(99, 110) rather than the defaults. At CAS(2,2) the third default centre was refused at full
size — its loop grazes the object — so no triangulation was possible. That run confined the
enclosed object to `pyr` in (105.6, 111.5), `tw` in (84, 96); (99, 110) sits rho ~ 0.75 from
the middle of that box against the old centre's ~0.88, far enough inside to bracket rather
than graze. It is a choice of centre, not an assumed answer: if the degeneracy is elsewhere,
the bisection returns a different radius or refuses.

*Output.* `results/localize/butadiene_cas<ne>-<ncas>.json` — every probe with its verdict and
cost, the bracket per centre, the triangulation and its uncertainty-normalised residual, and
a `complete` flag distinguishing a finished record from a resumable one.

*Known failure modes*, both seen at CAS(2,2) and both reported rather than hidden: a centre
whose full-size loop hits the step floor contributes nothing, and a centre whose full-size
loop does not enclose the target cannot be bracketed at all.

### `rung14.job` — one more rung, queued and not waited for (array of 2)

*What.* CAS(14,14) butadiene: task 0 runs the `tw = 90` gap-scan row and the fine `pyr` cut,
task 1 runs loop transport on all three loops at N = 13 and 21. The two are independent,
because the loops are fixed across the ladder by construction.

*Why last.* Butadiene has no exact in-basis reference (full valence is CAS(22,22)) and the
ladder has not settled over six rungs — the top two differ by 2.9 deg. One more rung shows
whether the position is settling or still wandering; it does not settle the §3 disagreement,
which is what the localization jobs are for. `docs/todo.md` §8 is explicit that this is not
on the critical path.

*Cost, and the weakest number in this file.* CAS(14,14) is 11.8 M determinants against
CAS(12,12)'s 853 k. If cost tracked determinant count that would be ~23 min per
state-specific point and ~50 min per state-averaged one — ~13 h per task. That extrapolation
has not been checked; `python examples/bench_point.py butadiene --cas 14 14` replaces it with
a measurement, and should be run before trusting it. Three days asked, so that a factor of
two is survivable.

*Output.* `results/butadiene/butadiene_scan_cas14-14_1x13.npz`,
`butadiene_refinepyr_cas14-14.npz`, and `butadiene_{B_x,B_1,B_2}_cas14-14_N{13,21}.json`.

### Pulling results back

Laptop-initiated, additive, and safe to repeat:

```bash
rsync -avz --exclude='.venv/' --exclude='__pycache__/' \
      alice:~/git_repos/BerryCASSCF/results/ results/
```

Then commit from the laptop, where the notebooks are rebuilt. Do not commit on the cluster:
the cluster clone is a worker, and the record of what ran is the `.job` file plus the commit
hash each job echoes into its log.

**One trap, and it will happen on the second sync.** `results/` is tracked in git *and* is
where the cluster writes. Committing a cluster-produced file from the laptop makes the next
`git pull` on the cluster abort:

```
error: The following untracked working tree files would be overwritten by merge:
        results/localize/butadiene_cas4-4_centre1.json
```

The file is not in danger — it is the same file, now committed — so the fix is to drop the
cluster's untracked copy, but only after checking it really is identical:

```bash
git fetch -q origin
git status --porcelain | awk '$1=="??"{print $2}' | while read -r f; do
    git cat-file -e "origin/main:$f" 2>/dev/null || continue
    if [ "$(md5sum < "$f" | cut -d' ' -f1)" = "$(git show "origin/main:$f" | md5sum | cut -d' ' -f1)" ]; then
        rm "$f"                      # already in origin/main, byte for byte
    else
        echo "KEPT $f (differs from origin/main)"
    fi
done
git pull --ff-only
```

A file that differs is kept and reported, because that means the cluster recomputed something
the laptop also has — which is a result to look at, not a conflict to clear.

**The same trap has a second form, and `rm` is the wrong fix for it.** A localization writes its
record after every probe, so a sync taken mid-run picks up an *incomplete* checkpoint. Commit
that, let the job finish, and the cluster's copy is now a **tracked modification** rather than an
untracked file:

```
error: Your local changes to the following files would be overwritten by merge:
        results/localize/ethylene_cas6-6_centre0.json
```

The cluster's version is the finished one and the committed version is the stale checkpoint, so
the order matters: **pull it to the laptop first, commit it there, then** `git checkout -- <file>`
on the cluster and pull. Doing the checkout first would discard the finished record in favour of
the checkpoint. Better still, do not commit records whose `complete` flag is false — they are
resumption state, not results.

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

## Generic templates, not site-specific

These predate the ALICE setup and still carry `#### SITE ####` placeholders, because none
of them has ever needed a cluster: the formaldimine workflows run locally in minutes, and
the fulvene follow-up was abandoned when its intersection turned out to be unreachable in
a rigid two-coordinate model (`docs/followup.md`). They are kept as starting points, not
as things to submit.

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
