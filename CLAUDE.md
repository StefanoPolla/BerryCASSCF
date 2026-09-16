# CLAUDE.md

Guidance for Claude Code working in this repository.

## What this is

A small research package for detecting conical intersections two independent ways — a
**variational Berry phase at CASSCF level** (Z2, from the sign of the initial–final overlap of
the state-specific S0 wavefunction transported around a nuclear loop) and a **state-averaged
CASSCF gap-scan comparator** over the enclosed region. The point is methodological: establish
when each method answers "does this loop enclose a CI?" correctly, and how large an active space
each one needs.

The original brief is `prompt.md`. Read it before any substantial change — it is the spec, and
the constraints below are drawn from it.

**Living state of the project:** `docs/progress.md` (append-only log, newest entry last) and
`README.md`. Those are authoritative; this file describes *how to work here*, not what the
current results are.

## Environment

Project-local venv, Python 3.12, PySCF 2.14 / NumPy 2.5 / SciPy 1.18:

```bash
source .venv/bin/activate      # always; do not use a system Python
pytest -q                      # full suite (~15 s)
pytest -q -m "not slow"        # no electronic structure, instant
```

Dependencies stay minimal: PySCF, NumPy, SciPy, plus matplotlib/jupyter for notebooks. Do not
add a dependency without a reason recorded in `docs/plan.md` or the progress log.

## Non-negotiable scientific constraints

These come from the brief. Violating one invalidates the result, so they are not judgement calls.

1. **Overlaps between adjacent loop points are genuine nonorthogonal CI overlaps.** The AO basis
   moves with the nuclei and the MOs are re-optimized. PySCF's
   `fci.addons.transform_ci_for_orbital_rotation` and any helper assuming a shared orbital space
   are **invalid** for this and must not be used for cross-geometry overlaps. The implementation
   is `berrycasscf/overlap.py` (Schur-complement core elimination), validated against a
   brute-force determinant-pair sum.
2. **Real wavefunctions; the phase is Z2.** Only the *sign* of an overlap carries phase
   information. Magnitudes (`|Π|`, `|ω|`) are continuity diagnostics — a clean sign at
   `|Π| = 0.5` is a good result, not a half-confident one. Never present a magnitude as
   confidence.
3. **Two estimators, always both.** The gauge-invariant cyclic product of adjacent overlaps and
   the gauge-fixed same-geometry endpoint overlap. They must agree in sign; disagreement is
   reported as a **failure**, never as a Berry phase.
4. **Continuation, not a series of independent solves.** Each point warm-starts from the previous
   one through the OAO frame, and its sign is fixed against the previous point (CASSCF returns an
   arbitrary overall sign — this is measured, not hypothetical).
5. **A run that fails a continuity criterion reports FAILED.** The thresholds live in
   `ContinuationConfig` and are the documented stability criteria. Do not widen a threshold to
   make a run pass; report the failure and diagnose it.
6. **C1 symmetry, and every CAS solve is spin-adapted.** PySCF's default solver returns the Ms=0
   triplet as a root; `apply_singlet_constraint` must be applied to every CAS/SA-CASSCF object.
   A regression of this once produced an entire S0/T1 gap map.
7. **The comparator is a 2D gap scan.** Minimum-energy CI optimization is explicitly out of
   scope. Keep basis, active space and coordinate conventions comparable to the Berry run.
8. **Do not force agreement between the two methods.** Disagreements are findings. Diagnose,
   document, and say whether the cause is understood (see the CAS(2,2) formaldimine case).

## Engineering conventions

- **Transparent research code, not a framework.** Flat module layout, no class hierarchies, no
  premature abstraction. Readability for a reviewing scientist beats generality.
- **Every scientific choice is a named field in `berrycasscf/config.py`** (`CasConfig`,
  `ContinuationConfig`, `ScanConfig`), with a comment explaining *why* that default. Nothing
  scientific belongs buried in code. Numerical defaults that exist because of a measured failure
  carry that measurement in the comment.
- **Reusable logic in `berrycasscf/`, experiments in `examples/`.** System modules
  (`geometry.py` for formaldimine, `ethylene.py`, `fulvene.py`) each expose a
  `<name>_geom(x, y) -> str` with the same signature, so drivers stay system-agnostic via their
  `geom_fn` argument. Follow that contract when adding a system, along with a `Loop` region and
  cheap geometry tests.
- **Drivers are restartable.** Skip work already saved under `results/`, support `--force`, and
  checkpoint long scans row by row. Result paths are content-identifying:
  `results/berry/<loop>_cas<ne>-<ncas>_N<n>.json`,
  `results/scan/<loop>_cas<ne>-<ncas>_<na>x<np>.npz`.
- **Notebooks are generated, never hand-edited.** `notebooks/build_*.py` writes the `.ipynb`
  from saved results only (no electronic structure); execution is a separate
  `jupyter nbconvert --execute --inplace` step. Edit the builder, rebuild, re-execute.
  Keep notebooks per-topic and readable end to end without rerunning; start a new one rather
  than growing an existing one past a comfortable length.
- **Tests:** cheap ones (geometry, overlap algebra, toy model) must run without any CASSCF;
  anything needing an SCF gets `@pytest.mark.slow`. Every bug that is fixed gets a regression
  test. Prefer tests that check an invariant the physics guarantees (self-overlap = 1, mirror
  symmetry of a gap map, gauge invariance under random sign flips) over golden numbers.

## Compute budget

- Local runs are fine up to a few hours total, and may be run in parallel or in sequence. The
  formaldimine benchmark must stay a minutes-long laptop calculation.
- Anything heavier is **not** to be run locally. Estimate the runtime, add or update a SLURM
  template in `slurm/` with every site-specific field marked `#### SITE ####` (partition,
  account, time, modules, env activation), document the job in `docs/compute.md` (what it does,
  how to submit, expected runtime and resources, expected outputs, where results land), and ask
  the user to submit it.
- Keep the measured-cost table in `docs/compute.md` current; timings there are measurements, not
  estimates.

## Documentation and record-keeping

`docs/` is part of the deliverable, not an afterthought:

| File | Rule |
|---|---|
| `docs/progress.md` | Append an entry at every meaningful milestone: what was implemented/tested, key decisions, results, unresolved issues, next step. Append-only, dated. Another agent must be able to resume from it. |
| `docs/provenance.md` | Every setting taken from an upstream source, inferred, or deliberately changed — with the reason. Extend it whenever a new system or setting arrives. |
| `docs/plan.md` | Design decisions, including backend capability reasoning. |
| `docs/results.md`, `docs/active_space.md` | Numbers and their interpretation. |
| `docs/limitations.md` | What the workflows do not do, and open questions. Resolved questions are struck through and pointed at the resolution rather than deleted. |
| `docs/compute.md`, `docs/followup.md` | Cluster jobs; follow-up-system status. |
| `docs/literature/` | Archived open-access sources. Everything stays in-repo. |

Nothing important lives only in chat.

## Writing style for results and commits

- **Quantify every claim.** "Agrees to 4.4e-16", "0.65 s/point over 82 points" — not "works
  well". State the sample the number came from.
- **Report negative and null results.** Fulvene's intersection being unreachable in its rigid
  two-coordinate model is recorded in the progress log, not erased. A method that fails is a
  finding.
- **Distinguish "gets the right answer" from "can be trusted without knowing the answer".**
  That distinction is the core of the active-space conclusions; keep it explicit.
- **Say what a test does *not* establish.** Several commits exist purely to narrow an earlier
  overclaim; that is the expected standard, not an exception.
- Commits: imperative subject line, body giving the numbers and the reasoning, including what
  remains unfixed. End with the attribution line the session requires.

## Known traps

- `conv_tol_grad` tighter than ~1e-5 makes *warm-started* CASSCF fail its convergence test while
  sitting on the correct energy to 1e-10. Expect spurious FAILED verdicts if you tighten it.
- SA-CASSCF has multiple stationary points. A warm-swept scan is path-dependent (measured: two
  mirror-image ethylene geometries gave 19.1 vs 29.4 mHa). `ScanConfig.strategy` defaults to
  `anchor` for path independence; `cold` is used for the ethylene ladder for cost. The
  mirror-symmetry check is necessary but not sufficient — a consistently wrong branch passes it.
- Ethylene needs 6-31G*; STO-3G has no intersection along its two coordinates at all.
- The in-CAS CASCI gap recorded along a loop is a root-flipping indicator, **not** the physical
  S0/S1 gap. Use the SA scan for physical gaps.
- `auto_oo` lives at `Emieeel/auto_oo`; the link in `prompt.md` 404s.

## Current scope (snapshot — check `docs/progress.md` for the truth)

Formaldimine/STO-3G is the benchmark (both workflows, FCI-calibrated). Ethylene/6-31G* is the
active-space ladder, CAS(2,2)–(12,12), in a separate notebook. Fulvene is implemented and
prepared for the cluster but its intersection is not reachable in the current rigid coordinate
model. Goals may change; the constraints and conventions above are meant to outlive them.
