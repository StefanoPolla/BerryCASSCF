# Task: reusable CASSCF implementation of a variational Berry-phase CI detector

Work in this repository to develop a small, readable Python research package and reproducible examples for detecting conical intersections (CIs) with a variational Berry-phase method based on CASSCF. The project should also implement a state-averaged CASSCF comparison workflow.

The goal is methodological and exploratory: establish whether the Berry-phase approach can robustly answer the question **“does this nuclear-coordinate loop enclose a conical intersection?”**, and compare that answer with a conventional state-averaged CASSCF resolution of the CI in the enclosed two-dimensional region.

Do not aim for production software or an elaborate abstraction layer. Prefer transparent, well-documented research code whose assumptions, numerical choices, and failure modes can be inspected and modified.

## First: orient yourself and make an implementation plan

1. Download and read the minimal relevant literature, beginning with:
   - arXiv:2304.06070, *A hybrid quantum algorithm to detect conical intersections*;
   - arXiv:2009.11417 and the relevant cited work on state-averaged CASSCF CI resolution.
2. Save open-access copies or stable bibliographic records in the repository for future consultation. Download further sources only when they are useful for a technical decision or benchmark.
3. The implementation details for the formaldimine benchmark are documented in the `auto_oo` repository, specifically:

   https://github.com/EmielKoridon/auto_oo/blob/main/examples/Tutorial_Berry_phase.ipynb

   Consult this notebook explicitly for the geometry parameterization, the loop definitions (both the CI-enclosing loop and the control loop), and the numerical settings used. Record in the repository which choices were taken from it and which were changed, and why. Do not reconstruct these from the paper text if the notebook specifies them.
4. Before committing to an electronic-structure backend, perform a short feasibility check. Choose a popular open-source backend that can support the required CASSCF calculations and the necessary access to orbitals, CI vectors, overlaps, and state averaging. PySCF is a likely candidate, but make the choice based on practical suitability rather than assuming it.
5. Write a concise implementation plan and record the main design choices, including what the backend can and cannot provide directly.

Use a Python virtual environment installed within the project directory, preferably `.venv`.

## Scientific scope

### 1. Variational Berry-phase workflow at full CASSCF level

Implement a full CASSCF version of the variational Berry-phase idea in arXiv:2304.06070. The variational variables include both orbital rotations and CI coefficients.

The tracked state is the state-specific CASSCF ground state (S0); the conical intersection of interest is S0/S1.

**Real-wavefunction assumption and estimator.** Assume throughout that orbitals and CI coefficients are real. The Berry phase is then a Z2 quantity: trivial (0) or non-trivial (pi), read off from the *sign* of the initial-final overlap. Magnitude carries no phase information and serves only as a continuity diagnostic. Note that independently converged solves may return an arbitrary overall sign even under warm start, so the endpoint overlap is only meaningful when every intermediate overlap is close to +1 in a consistent gauge. As cheap insurance, also accumulate the running product of consecutive overlaps around the loop: under the stated assumptions it must equal the endpoint overlap, so any disagreement indicates a gauge or continuity bug rather than physics. Report both.

**Nonorthogonal overlaps.** Wavefunctions at adjacent geometries live in different AO and MO bases, so the overlap between neighbouring points is a genuine nonorthogonal CI overlap. Routines that assume the two orbital sets span the same space -- including PySCF's `fci.addons.transform_ci_for_orbital_rotation` and similar helpers -- are **not valid** here and must not be used for this purpose. The choice of method is left to the implementer, but validate it with cheap tests: the overlap of a wavefunction with itself must be 1, and two independent solves at the same geometry must give +/-1.

The central requirement is **continuous tracking of one intended CASSCF stationary solution around a closed nuclear loop**. This is not merely a series of independent CASSCF calculations. Design the workflow to use the preceding solution as the initial guess and to preserve the intended state/solution by suitable continuation and overlap-based diagnostics.

Pay particular attention to:

- continuity of orbitals, CI vectors, and the total CASSCF wavefunction;
- phase/sign conventions and the orbital gauge, including the fact that individually optimized wavefunctions may acquire arbitrary signs;
- avoiding accidental root/solution switching or convergence to a different local stationary point;
- checking that the endpoint represents the same physical state as the starting point before interpreting its signed overlap;
- recording intermediate overlaps and diagnostics so that loss of continuity is detectable rather than silently producing a Berry phase;
- allowing practical predictor–corrector, augmented-Hessian/Newton-like, or warm-started optimization strategies as appropriate for the chosen backend.

The initial implementation does not need to reproduce every quantum-device-specific ingredient of the paper. It should faithfully implement the corresponding **classical CASSCF variational continuation problem** and infer the Berry phase from the consistently gauged endpoint overlap. Use the literature to determine appropriate diagnostics and conventions.

A successful calculation should report, at minimum:

- the loop and discretization used;
- CASSCF state energies and convergence status along the loop;
- continuation/overlap diagnostics between adjacent points;
- the magnitude and sign of the initial–final wavefunction overlap;
- the inferred trivial or non-trivial Berry phase;
- warnings or a clear failure status when continuity or endpoint equivalence is not established.

### 2. State-averaged CASSCF comparator

Implement a state-averaged CASSCF workflow that can independently investigate the same question for a loop.

The comparison target is deliberately different: use equal-weight or otherwise justified state-averaged CASSCF to resolve the relevant two-state near-degeneracy within the two-dimensional region enclosed by the loop. A two-dimensional gap scan over the plane enclosed by the loop is sufficient; minimum-energy CI optimization is explicitly **not** required. An adaptive or refined search may be added if it is cheap and helpful.

The output should make it possible to determine whether a same-symmetry CI lies in the chosen plane and inside the loop, and to compare this conclusion with the Berry-phase result. Report the estimated CI location or closest-approach region, energy gap, search/scan resolution, and limitations. Keep the active space, basis, molecular Hamiltonian, and nuclear-coordinate conventions comparable to the Berry-phase calculation whenever that is scientifically reasonable.

Do not force artificial agreement: document discrepancies and investigate obvious numerical or physical explanations, such as an insufficient active space, an incomplete search, a different tracked state, or a path too close to the seam.

## Required experiments

### 0. Optional toy-model validation

Optional, at your discretion: before running any quantum chemistry, validate the sign/continuation machinery on a 2x2 linear Jahn-Teller (E x e) model Hamiltonian, whose Berry phase is known analytically (pi for loops enclosing the intersection, 0 otherwise). This is cheap and catches sign and continuation bugs before CASSCF cost obscures them. Skip it if you judge the effort better spent elsewhere, but say so in the progress log.

### A. Formaldimine benchmark

Use the known formaldimine example from arXiv:2304.06070 as the first benchmark, following the implementation details in the `auto_oo` Berry-phase tutorial notebook linked above.

**Fixed inputs.** Use the following parameterized geometry verbatim. All other internal coordinates are frozen; the two loop coordinates are `alpha` (the HNC angle) and `phi` (the HNCH dihedral). Distances are in Angstrom, angles in degrees. Run in C1 (no point-group symmetry), since the loop breaks symmetry.

```python
def formaldimine_geom(alpha, phi):
    variables = [1.498047, 1.066797, 0.987109, 118.359375] + [alpha, phi]
    geom = """
                    N
                    C 1 {0}
                    H 2 {1}  1 {3}
                    H 2 {1}  1 {3} 3 180
                    H 1 {2}  2 {4} 3 {5}
                    """.format(*variables)
    return geom
```

Use the STO-3G basis for formaldimine. Start from CAS(2,2); note that the state-averaged comparator may require a larger active space than the Berry-phase calculation to resolve the two-state near-degeneracy, and this asymmetry is acceptable provided it is documented.

- Take the CI-enclosing loop and the control loop from the tutorial notebook where it defines them; otherwise construct both in the (`alpha`, `phi`) plane and state the centre and radius used.
- Run both the Berry-phase continuation and state-averaged CASSCF comparator.
- Explore active-space sizes systematically, starting from the smallest chemically sensible choice and increasing only where justified.
- Identify the smallest tested CAS for which **both** approaches give a stable, mutually consistent conclusion. “Stable” should include reasonable robustness to loop discretization and to relevant continuation/search settings, not only a single successful run. Define the concrete stability criteria yourself, state them explicitly in the documentation, and apply them consistently.
- If the methods disagree, do not hide the result. Diagnose it and state whether the cause is understood.

### B. A more demanding follow-up system

Only after the formaldimine workflows are working and tested, apply them to one more complicated system. Reconsider the choice after the benchmark, based on backend capabilities and available literature/geometries.

Possible candidates include:

- **Benzene**, where the full π manifold provides a clear step beyond a two-orbital active space and the S1/S0 CI seam is well studied.
- **Fulvene**, which offers a modestly larger π-space problem with documented conical-intersection chemistry.
- A small **protonated Schiff-base/retinal model**, if a documented, computationally tractable geometry and reference protocol can be found.
- A carefully selected **nucleobase model**, if the workflow can reliably handle the relevant lone-pair and π/π* state character.

Choose a system for which a meaningful loop and a reasonable reference comparison can be established. The candidate list is guidance, not a commitment. Record why the final system was chosen and what active-space considerations make it more demanding than formaldimine.

## Engineering and reproducibility requirements

- Keep all code, scripts, inputs, downloaded open-access reference material, generated configurations, and progress notes inside this repository.
- Write Python. Use a project-local virtual environment.
- Keep the code human-readable and modest in scope. Avoid premature general frameworks, excessive class hierarchies, or unnecessary dependencies.
- Separate reusable calculation logic from example-specific molecular geometries and experiment drivers where this improves clarity.
- Include concise documentation explaining installation, backend choice, how to run each example, expected outputs, and known limitations.
- Provide lightweight automated tests for components that can be tested cheaply without long electronic-structure calculations. Include at least smoke tests or reproducible short-run checks for the full workflows where feasible.
- Make calculations restartable where practical. Save enough intermediate data to resume a loop or scan rather than rerunning successful points unnecessarily.
- Periodically append a brief progress record in the repository. At each meaningful milestone, state: what was implemented or tested, key decisions, current results, unresolved issues, and the next recommended step. This should allow another agent or researcher to pause and resume the project with minimal rediscovery.
- Present the results in a narrative form that can be read without rerunning anything -- a Jupyter notebook is the expected default, but choose whatever format explains the results best. At minimum it should show energies along each loop, the overlap diagnostics along each loop, and the two-dimensional gap scan.
- **Compute budget and cluster jobs.** The formaldimine benchmark must run locally in reasonable wall-clock time; keep runtimes modest and say so in the documentation. For anything heavier -- in particular the follow-up system -- do not attempt to run it locally. Instead prepare SLURM batch script templates with clearly marked placeholders for the site-specific fields (partition, account, time limit, modules, environment activation), and ask me to submit them. Record in a dedicated file (e.g. `docs/compute.md`) what each job does, how to submit it, the expected runtime and resource needs, the expected outputs, and where results should be written on return. This information must live in the repository, not only in the chat.
- Use explicit configuration files or clearly visible parameters for molecular geometry, basis set, active space, state averaging, loop discretization, convergence thresholds, and scan/search settings. Do not bury scientific choices in opaque code.

## Working style

Make reasonable implementation choices without requesting confirmation for ordinary engineering details. If a scientific ambiguity materially affects interpretation—for example, the definition of the benchmark loop, a state-label ambiguity, or an active-space choice that changes the physical target—record the alternatives, select a defensible initial choice, and flag it clearly in the progress log and documentation.

Prioritize a correct, inspectable formaldimine benchmark over broad feature coverage. Once that is trustworthy, extend carefully to the second system.

## Final handoff

At the end, leave the repository in a runnable state and provide a short final summary containing:

1. backend and package structure;
2. how to reproduce the formaldimine results;
3. Berry-phase and state-averaged-CASSCF conclusions for each tested loop and CAS;
4. the identified minimal reliable CAS for formaldimine, with the criterion used;
5. results and limitations for the follow-up system, if completed;
6. known limitations, unresolved questions, and the highest-value next steps.
