# SLURM templates

Every line marked `#### SITE ####` must be filled in before submission. The templates are
deliberately not runnable as-is, so nothing is submitted with a guessed account or partition.

See `docs/compute.md` for what each job does, expected runtimes, and where output lands.

| template | job |
|---|---|
| `berry_loop.sbatch` | Berry-phase loops for formaldimine at a chosen active space |
| `sa_scan.sbatch` | SA-CASSCF gap scan over a region |
| `fulvene_berry.sbatch`, `fulvene_scan.sbatch` | the follow-up system |
| `localize.sbatch` | bisection + triangulation, to measure where the state-specific degeneracy is (`docs/todo.md` §9) |
