# Cluster jobs

These are **runnable**, not templates: they carry concrete values for ALICE (Leiden
University), where the production runs of this project are submitted. Every line that is
specific to that site is marked `# SITE:` — partition, account, module stack, venv path —
and those are the only lines another cluster needs changed.

```bash
# on the cluster, once
git clone git@github.com:StefanoPolla/BerryCASSCF.git ~/git_repos/BerryCASSCF
cd ~/git_repos/BerryCASSCF
ml purge && ml load ALICE/default && ml load Python/3.11.5-GCCcore-13.2.0
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"

# thereafter
cd ~/git_repos/BerryCASSCF && git pull --ff-only
cd slurm && sbatch --test-only smoke.job && sbatch smoke.job
```

Submit **from inside `slurm/`**: the `--output=./out/...` paths are relative to the
submission directory, and each script's first action is `cd ..` to reach the repository.

| job | what it runs | shape | walltime asked |
|---|---|---|---|
| `smoke.job` | fast tests, then per-point CASSCF cost at CAS(2,2) and CAS(12,12), threaded and single-threaded | one task | 40 min |
| `berry_cas12.job` | butadiene loop transport at CAS(12,12): the two missing control loops and a third discretization | array, one loop per task | 16 h |
| `localize.job` | butadiene bisection/triangulation at CAS(4,4), (6,6), (8,8), (10,10) | array, one active space per task | 36 h |
| `localize_cas12.job` | the same at CAS(12,12), where the two methods disagree | one task | 5 days |
| `rung14.job` | CAS(14,14): gap scan + fine cut (task 0), loop transport (task 1) | array of 2 | 3 days |

Three generic templates remain, still carrying `#### SITE ####` placeholders and
deliberately not runnable as they stand: `berry_loop.sbatch` and `sa_scan.sbatch` (the
formaldimine workflows, which run locally in minutes and have never needed a cluster), and
`fulvene_berry.sbatch` / `fulvene_scan.sbatch` (the abandoned follow-up system — see
`docs/followup.md`). `localize.sbatch` was superseded by the two localization jobs above and
is deleted; it is in the git history if a placeholder version is ever wanted again.

**Every driver skips work already saved under `results/`**, printing a line beginning
`SKIP:` and exiting 0. Re-submitting a job therefore runs only what is missing, and a job
that finds everything present costs seconds. Pass `--force` to redo.

`out/` and `err/` are gitignored; the `.job` files are not, because they are the record of
what was actually run. Results land in `results/<system>/` and progress is written live to
`logs/<name>.log` inside the repository — `tail -f` that to tell a slow job from a dead one.

Sizing, resources and what to do when a job comes back: `docs/compute.md`.
