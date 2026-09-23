# The notebooks

Eight notebooks. Each is **generated** by its `build_*.py` script from saved records in
`results/`, so none of them runs any electronic structure — they report, they do not recompute.
To refresh one after new results land:

```bash
python notebooks/build_<name>_notebook.py
jupyter nbconvert --to notebook --execute --inplace notebooks/<name>.ipynb
```

Edit the builder, never the `.ipynb`.

---

## Read them in this order

### 1. [formaldimine_benchmark.ipynb](formaldimine_benchmark.ipynb) — start here
**Does this loop enclose a conical intersection?** The whole method on the one system where the
answer is independently known: FCI in the same basis. Loop transport against a state-averaged gap
scan, on a CI-enclosing loop and two controls, with the nonorthogonal overlap, the gauge fixing and
the two estimators explained from first principles. Establishes the vocabulary every later notebook
uses, and the verdict that CAS(4,4) is the smallest active space where both methods agree.

*If you read only one, read this one.*

### 2. [ethylene_ladder.ipynb](ethylene_ladder.ipynb) — how big must the active space be?
The same two workflows over CAS(2,2) → CAS(12,12) on one system, with **everything else held
fixed** — same basis, same geometry family, above all the same loops. CAS(12,12) is the full valence
space, so it is exact in this basis and the ladder can be *scored*. The answer separates two
questions usually conflated: which active space gets the right answer (CAS(2,2), accidentally) and
which one you could trust without already knowing it (CAS(10,10)).

### 3. [butadiene_ladder.ipynb](butadiene_ladder.ipynb) — the same ladder with no reference
A bigger system, no exact answer available, and the honest consequences. The gap-scan position
scatters over 19° and never settles while loop transport returns the same topology at every rung —
the inversion that gives this project its most interesting result and its least resolved question.
Ends with both methods' *positions* side by side, which is where the disagreement becomes numeric.

### 4. [probing_by_small_loops.ipynb](probing_by_small_loops.ipynb) — is the loop encircling what we think?
**Read this before believing any single π.** A small loop placed *on* an intersection asks, in one
walk, whether loop transport encircles the object the gap scan found. Validated on a control whose
answer was known in advance, then applied across all three systems — where it finds that at ethylene
CAS(2,2) and at every butadiene rung, it does not. Also measures how small an enclosing loop can be
before the state-specific solver gives out (about 1°).

*This is the newest and most consequential notebook; §5–§6 contain a result that cuts against
claims made elsewhere in the repository, and say so.*

---

## Method notebooks — read when the question comes up

### [locating_intersections.ipynb](locating_intersections.ipynb)
**Can a method that returns one bit per loop return a *position*?** Bisection on the loop radius
turns the bit into a distance; several centres triangulate it into a point. The geometry is
validated separately from the physics, the answer comes back as a bracket whose width *is* the
resolution, and the residual from three or more centres refuses to produce a position when no single
degeneracy explains the data. Read it before §3 of the butadiene notebook or §7 of the probing one,
both of which use its output.

### [adaptive_stepping.ipynb](adaptive_stepping.ipynb)
**How should a loop be discretized, and does adapting it pay?** The walk already measures its own
error, so the step size can be steered by it for free. The predicted cost saving did not
materialise — measured at matched quality it is a wash. Its real value is that it walks loops
uniform discretization cannot walk at any affordable N, and reports its own resolution limit instead
of guessing.

### [stepping_comparison.ipynb](stepping_comparison.ipynb)
**One parameter update per point, or optimize each point to convergence?** arXiv:2304.06070 does the
former on a quantum-budget argument; this package does the latter. Both were run on the same loops
so the trade-off is measured rather than asserted.

---

## [summary.ipynb](summary.ipynb) — the cross-system verdict, short
What the three systems say together, and what was built to make them say it. Read it last as a
recap, or first if you want the conclusions before the evidence.

---

## How they tie together

The four numbered notebooks are one argument in sequence:

1. **formaldimine** shows the method works where the truth is known;
2. **ethylene** asks what it costs in active space, against an exact reference;
3. **butadiene** removes the reference and finds the comparator, not the method, is the unstable one;
4. **probing** asks whether the method's answers are about the object everyone assumed — and at two
   of the three systems, they are not.

Each step makes the previous one's claim more precise, and the fourth revises it. The three method
notebooks are the machinery those arguments lean on: **locating** turns bits into positions,
**adaptive stepping** decides where the points go, **stepping comparison** justifies converging each
point at all.

Two cross-cutting conventions worth knowing before reading any of them:

* **cost is counted in CASSCF micro-iterations**, not wall time. Local wall times in this repository
  are indicative only — the laptop sleeps and shares cores; see `docs/compute.md`. Cluster timings
  are dedicated allocations and are labelled where quoted;
* **a refusal is a result.** Every workflow here reports FAILED rather than a Berry phase when a
  continuity check fails, and several conclusions rest on *which* loops could not be walked. A
  refusal excludes nothing about geometry, which is a distinction the notebooks are careful about and
  which one earlier analysis got wrong.

## Where the prose lives instead

The notebooks carry the reasoning and the figures. The documents carry the record:

| for | read |
|---|---|
| consolidated findings across systems | [`docs/findings.md`](../docs/findings.md) |
| the formaldimine numbers in detail | [`docs/results.md`](../docs/results.md) |
| what the workflows do *not* do | [`docs/limitations.md`](../docs/limitations.md) |
| costs, cluster jobs, sizing | [`docs/compute.md`](../docs/compute.md) |
| dated progress log, and how to resume | [`docs/progress.md`](../docs/progress.md) |
| planned work not yet done | [`docs/todo.md`](../docs/todo.md) |
