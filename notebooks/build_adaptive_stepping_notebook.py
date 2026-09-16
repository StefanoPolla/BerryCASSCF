#!/usr/bin/env python
"""Generate notebooks/adaptive_stepping.ipynb -- how step control works and when it pays.

The Jahn-Teller figures are computed live: the model is a 2x2 diagonalization, not electronic
structure, so it costs nothing and lets the reader watch the controller work. Every CASSCF
number is read from results/adaptive/.
"""

import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS: list = []


def md(t): CELLS.append(nbf.v4.new_markdown_cell(t.strip("\n")))
def code(t): CELLS.append(nbf.v4.new_code_cell(t.strip("\n")))


md(r"""
# Adaptive step control: letting the loop choose its own discretization

**The workflow, named in full once:** everything here is **SS-CASSCF Berry-phase loop transport**
— continuation of a *state-specific* CASSCF ground state around a closed nuclear loop, with the Z2
Berry phase read from the sign of the initial–final nonorthogonal overlap. Short form: **loop
transport**. The SA-CASSCF gap scan does not appear in this notebook.

## The question

A loop is walked in $N$ steps. How should $N$ be chosen?

Uniform discretization forces one answer for the whole loop, so $N$ has to be large enough for the
**hardest arc** — and every other arc is then oversampled. On the loops in this project the
difficulty is very uneven: the worst adjacent overlap on butadiene's `B_x` at CAS(8,8) is 0.89
while the median is 0.99. The hard part is a small fraction of the loop.

Two things make step control unusually well suited to this particular problem:

1. **The error indicator is already there and is nearly free.** Every step computes
   $s = |\langle\Psi_{k-1}|\Psi_k\rangle|$ for the continuity check. Evaluating it costs ~0.1 ms
   against a CASSCF solve of 30 ms to 3 minutes. Most adaptive schemes have to pay for their error
   estimate; this one does not.
2. **Unequal spacing costs nothing.** A topological readout needs the chain to be *continuous* and
   *closed*. It does not need the points to be evenly spaced — nothing in the Berry phase refers to
   the spacing at all.

## The controller

Write the per-step **mismatch** as $m = 1 - s$. For a smoothly transported state, $m \sim d^2$ in
the step size $d$, so a step that achieved $m$ with size $d$ should next use

$$ d_\mathrm{new} = d\,\sqrt{m_\mathrm{target}/m} $$

clipped to $[d_\min, d_\max]$ with a cap on growth.

**Steering toward a target is deliberately different from rejecting at a floor.** A hard
accept/reject at the failure threshold oscillates: it accepts until it fails, shrinks, grows back,
fails again. Each rejection costs a *full CASSCF solve* and there is no cheap way to predict one,
so the walk is cheapest when it sits just inside the accept threshold rather than bouncing off it.
The default targets $m = 0.02$ (overlap 0.98) while refusing only below 0.80, so the walk carries
margin instead of skating the limit.

Two details that matter more than they look:

* **The walk lands exactly on $t = 1$.** The final step is clipped so the closing geometry is
  *identical* to the start. Otherwise the endpoint estimator stops being exact, which is the whole
  reason it is worth computing.
* **A rejection restores the previous accepted state.** The trial is discarded entirely. This is
  easy to get subtly wrong and the failure is silent, so it lives in one place.

## Reaching the floor is a result, not a safeguard

If shrinking the step restores continuity, the problem was discretization. If it does not, the
problem is not discretization — and the walk says so rather than returning a sign. That is the
distinction this project previously had to make by hand (`docs/findings.md` §5, the radius-8 versus
radius-6 loops), and it is now discovered automatically.
""")

code(r"""
import json, os, sys
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.getcwd()) if os.path.basename(os.getcwd()) == "notebooks" else os.getcwd()
sys.path.insert(0, ROOT)

from berrycasscf.config import AdaptiveConfig
from berrycasscf.toy import jt_loop_adaptive, jt_loop_berry_phase, jt_gap

plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.3, "figure.autolayout": True})

def load(path):
    p = os.path.join(ROOT, path)
    return json.load(open(p)) if os.path.exists(p) else None

print("Jahn-Teller figures are computed live below (a 2x2 diagonalization, no electronic")
print("structure). Every CASSCF number is read from results/adaptive/.")
""")

md(r"""
## 1. Watching the controller work

Before any quantum chemistry, the same controller is run on the 2×2 linear Jahn–Teller model,
whose Berry phase is known in closed form: $\pi$ for any loop enclosing the origin, $0$ otherwise.
A step here costs a 2×2 diagonalization, so the figures are free — and, more importantly, a failure
seen here can only be a **control** bug, never a CASSCF bug. That separation is the reason this
model exists in the package.

The loop below sits off-centre, so its inner arc passes close to the degeneracy at the origin while
its outer arc is far away. A uniform walk cannot know that. The adaptive walk measures it.
""")

code(r"""
centre, radius = (0.0, 0.85), 1.0
res, walk = jt_loop_adaptive(centre=centre, radius=radius)
ts = np.array(walk.t_values)
pts = np.array([(centre[0] + radius*np.cos(2*np.pi*t), centre[1] + radius*np.sin(2*np.pi*t))
                for t in ts])

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

# --- left: the points the walk chose, on the loop, over the gap it was walking through ---
ax = axes[0]
g = np.linspace(-2.1, 2.1, 400)
X, Y = np.meshgrid(g, g)
ax.contourf(X, Y, 2*np.hypot(X, Y), levels=25, cmap="Blues_r", alpha=0.55)
ax.plot(pts[:, 0], pts[:, 1], "-", color="0.35", lw=0.8, zorder=2)
ax.scatter(pts[:, 0], pts[:, 1], s=26, c="crimson", zorder=3, label=f"{len(ts)-1} points chosen")
ax.plot(0, 0, "*", ms=16, color="gold", mec="k", mew=0.6, zorder=4, label="degeneracy")
ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_title("where the walk spent its solves\n(background: the gap it is walking through)")
ax.legend(loc="lower right", fontsize=7.5)

# --- right: step size around the loop, against the distance to the degeneracy ---
ax = axes[1]
d = walk.step_sizes
t_mid = 0.5*(ts[:-1] + ts[1:])
dist = np.hypot(pts[:, 0], pts[:, 1])
ax.step(ts[1:], d, where="pre", color="crimson", lw=1.6, label="step size $d$")
ax.set_xlabel("loop parameter $t$"); ax.set_ylabel("step size $d$", color="crimson")
ax.tick_params(axis="y", labelcolor="crimson")
ax2 = ax.twinx(); ax2.grid(False)
ax2.plot(ts, dist, color="steelblue", lw=1.4, ls="--", label="distance to degeneracy")
ax2.set_ylabel("distance to degeneracy", color="steelblue")
ax2.tick_params(axis="y", labelcolor="steelblue")
ax.set_title("the controller slows down where the state turns fastest")
fig.suptitle(f"Adaptive walk, loop centred {centre} radius {radius}", y=1.02)
plt.show()

print(f"step sizes span {d.min():.4f} to {d.max():.4f}  (ratio {d.max()/d.min():.1f}x)")
print(f"phase = {'pi' if res.is_nontrivial else '0'}, theory says "
      f"{'pi' if res.expected_nontrivial else '0'}; agrees: {res.agrees_with_theory}")
""")

md(r"""
The two curves are mirror images by construction: the step size is smallest exactly where the loop
passes closest to the degeneracy, and it is not told where the degeneracy is. It infers it from the
overlaps it measures.

## 2. What a rejection looks like

With a deliberately over-large $d_\max$ the walk oversteps and has to back up. The rejected trials
are real CASSCF solves in the chemistry case, so they are charged to the cost — a comparison that
hid them would flatter the method.
""")

code(r"""
cfg_coarse = AdaptiveConfig(d_max=0.5, target_mismatch=0.02, accept_mismatch=0.05)
_, w2 = jt_loop_adaptive(centre=(0.0, 0.0), radius=1.0, cfg=cfg_coarse)

fig, ax = plt.subplots(figsize=(8.4, 3.6))
for e in w2.events:
    colour = "seagreen" if e.accepted else "crimson"
    ax.plot([e.t_from, e.t_to], [e.d, e.d], "-", color=colour, lw=2.2, alpha=0.85)
    ax.plot(e.t_to, e.d, "o" if e.accepted else "x", color=colour, ms=5.5,
            mew=1.6 if not e.accepted else 0.5)
ax.axhline(cfg_coarse.d_min, color="0.4", ls=":", lw=1, label="$d_{min}$ (hard floor)")
ax.set_yscale("log"); ax.set_xlabel("loop parameter $t$"); ax.set_ylabel("step size $d$")
ax.set_title(f"every attempted step: {w2.n_accepted} accepted (green), "
             f"{w2.n_rejected} rejected (red)")
ax.legend(fontsize=8)
plt.show()

print(f"{w2.n_rejected} rejections out of {len(w2.events)} attempts.")
print("A rejected trial is discarded completely: the walk retries from the SAME accepted state,")
print("never from the trial. The previous state is never modified, which is the one correctness")
print("requirement of a rejection and the one whose failure would be silent.")
""")

md(r"""
Notice the shape: the walk oversteps, halves back, and then **settles** rather than oscillating.
That is the controller steering toward $m_\mathrm{target}$ instead of merely rejecting at the
threshold — a pure accept/reject rule would keep bouncing off the limit and pay a full solve each
time it did.

## 3. How close can a loop pass? A predicted limit, and a measurement

There is a hard limit on how near a degeneracy a loop can pass and still be walked. It can be
derived rather than discovered.

A loop of radius $R$ passing at distance $\varepsilon$ sweeps the state's angle at
$\mathrm{d}\theta/\mathrm{d}t \approx 2\pi R/\varepsilon$ near closest approach. Keeping the
per-step angle below $\Delta\theta_\max = 2\arccos(1 - m_\mathrm{accept})$ therefore needs a step of
about $\varepsilon\,\Delta\theta_\max/(2\pi R)$, and the walk gives up when that falls below
$d_\min$:

$$ \boxed{\;\varepsilon_\min \;\approx\; \frac{2\pi R\, d_\min}{\Delta\theta_\max}\;} $$

With the defaults ($d_\min = 10^{-3}$, $m_\mathrm{accept} = 0.1$) this is **0.7% of the loop
radius**. Below that the walk reports the floor instead of a phase.
""")

code(r"""
cfg = AdaptiveConfig()
dtheta_max = 2*np.arccos(1 - cfg.accept_mismatch)
eps_pred = 2*np.pi*1.0*cfg.d_min/dtheta_max

eps_list = [0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.007, 0.005, 0.003, 0.001]
rows = []
for eps in eps_list:
    c = (0.0, 1.0 - eps)
    u = jt_loop_berry_phase(centre=c, radius=1.0, n_points=24)
    r, w = jt_loop_adaptive(centre=c, radius=1.0)
    rows.append(dict(eps=eps,
                     u_ok=bool(u.min_abs_adjacent_overlap >= 0.80),
                     u_minovl=u.min_abs_adjacent_overlap, u_right=u.is_nontrivial,
                     a_floor=w.hit_step_floor, a_closed=w.closed,
                     a_minovl=r.min_abs_adjacent_overlap, a_pts=r.n_points,
                     a_right=r.is_nontrivial))

fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
ax = axes[0]
e = np.array([r["eps"] for r in rows])
ax.semilogx(e, [r["u_minovl"] for r in rows], "o-", color="steelblue",
            label="uniform, $N=24$")
ax.semilogx(e, [r["a_minovl"] if r["a_closed"] else np.nan for r in rows], "s-",
            color="crimson", label="adaptive")
ax.axhline(0.80, color="0.3", ls="--", lw=1, label="continuity threshold")
ax.axvline(eps_pred, color="seagreen", ls=":", lw=1.6,
           label=fr"predicted $\varepsilon_{{min}}$ = {eps_pred:.4f}")
ax.set_xlabel(r"closest approach $\varepsilon$ (loop radius 1)")
ax.set_ylabel("worst adjacent overlap"); ax.set_ylim(0.4, 1.02)
ax.set_title("adaptive holds its margin until the predicted limit")
ax.legend(fontsize=7.5, loc="lower right")

ax = axes[1]
for r in rows:
    y = 1 if r["u_ok"] else 0
    ax.plot(r["eps"], 1.0, "o", ms=8, color="seagreen" if r["u_ok"] else "crimson")
    ok_a = r["a_closed"] and not r["a_floor"]
    ax.plot(r["eps"], 0.0, "s", ms=8, color="seagreen" if ok_a else "crimson")
ax.set_xscale("log"); ax.set_yticks([0, 1])
ax.set_yticklabels(["adaptive", "uniform $N=24$"])
ax.axvline(eps_pred, color="seagreen", ls=":", lw=1.6)
ax.set_ylim(-0.6, 1.6); ax.set_xlabel(r"closest approach $\varepsilon$")
ax.set_title("green = a trustworthy answer, red = refused")
plt.show()

first_floor = min((r["eps"] for r in rows if r["a_floor"]), default=None)
last_ok = min((r["eps"] for r in rows if r["a_closed"] and not r["a_floor"]), default=None)
print(f"predicted floor          : eps = {eps_pred:.4f}")
print(f"closest eps still walked : {last_ok}")
print(f"largest eps that floored : {first_floor}")
print()
print("Uniform N=24 returns the CORRECT sign at every eps here, and fails its own continuity")
print("check from eps = 0.1 inwards. Being right is not the same as being trustworthy:")
print("a user without the answer would have to refuse those runs.")
""")

md(r"""
The measurement lands on the prediction: the walk closes at $\varepsilon = 0.007$ and hits the
floor at $0.005$, against a predicted $0.0070$.

Two things worth taking from the right-hand panel:

* **Adaptive extends the trustworthy range by more than an order of magnitude** in closest
  approach — from $\varepsilon \approx 0.1$ to $\varepsilon \approx 0.007$ — while using 16–23
  points. A uniform walk would need $N \approx 1000$ to hold the same margin at $\varepsilon = 0.007$.
* **Uniform stepping is accidentally right long after it stops being trustworthy.** It returns the
  correct sign at every $\varepsilon$ tested, including ones where its worst overlap has collapsed
  to 0.66. Without knowing the answer in advance there is no way to tell those apart, which is
  exactly the distinction this project keeps insisting on.

This is what makes the bisection of `locating_intersections.ipynb` possible: that search spends its
time on loops that pass deliberately close to a degeneracy.

## 4. The informative failure

A loop through the degeneracy cannot be walked at any step size, and the walk should say so rather
than return a sign that happens to have one.
""")

code(r"""
_, w_through = jt_loop_adaptive(centre=(1.0, 0.0), radius=1.0)   # passes exactly through (0,0)
_, w_near    = jt_loop_adaptive(centre=(0.0, 0.85), radius=1.0)  # passes 0.15 away

fig, ax = plt.subplots(figsize=(8.0, 3.4))
for w, name, col in ((w_near, "passes 0.15 away", "seagreen"),
                     (w_through, "passes THROUGH the degeneracy", "crimson")):
    acc = [(e.t_to, e.d) for e in w.accepted_events]
    ax.step([a[0] for a in acc], [a[1] for a in acc], where="pre", color=col, lw=1.8, label=name)
ax.axhline(AdaptiveConfig().d_min, color="0.3", ls=":", lw=1.2, label="$d_{min}$")
ax.set_yscale("log"); ax.set_xlabel("loop parameter $t$"); ax.set_ylabel("step size $d$")
ax.set_xlim(0, 1); ax.legend(fontsize=8)
ax.set_title("shrinking recovers one loop and cannot recover the other")
plt.show()

print(f"near     : closed={w_near.closed}   floor={w_near.hit_step_floor}")
print(f"through  : closed={w_through.closed}  floor={w_through.hit_step_floor}")
print()
print("The second walk stops partway and reports the floor. It does NOT report a phase.")
print("The message names both possible causes -- the loop passing through a degeneracy, or the")
print("solver switching between nearby stationary solutions -- because the walk cannot tell")
print("them apart and should not pretend to.")
""")

md(r"""
### A failure this caught in the real code

An early version of this walk had a hole in exactly this place. When the step floor stopped the
walk partway round, the run still came back **`status OK`, phase `trivial (0)`** — because the
endpoint checks are skipped when there is no endpoint, every *accepted* step was continuous, and
$|\langle\Psi_\mathrm{last}|\Psi_0\rangle|$ — the single factor carrying the sign of the product
estimator — was never checked at all. The unclosed path was silently treated as a loop.

Two checks now cover it (`loop_closed`, `closing_step_continuous`), and both are regression-tested.
A uniform walk closes by construction and pays nothing for them. It is worth being explicit that
this bug was *not* found by reasoning about the code — it was found because a figure showed a
verdict that could not be right.

## 5. Does it pay on real CASSCF?

The honest answer is **it depends on the loop, and on formaldimine it does not.**

The comparison is a *family against a family*. Adaptive stepping has knobs ($d_\max$,
$m_\mathrm{target}$) and uniform stepping has $N$; comparing one setting of each would prove
nothing, since either can be made to look good by choosing its parameter well. Both are swept.

Cost is counted in **micro-iterations** — CASSCF parameter updates — because wall time carries
per-point overhead that differs between machines. Rejected adaptive trials are charged in full.
""")

code(r"""
form = load("results/adaptive/formaldimine_stepping.json")
buta = load("results/adaptive/butadiene_stepping.json")

def frontier_table(data, title):
    if not data:
        print(f"{title}: not run yet"); return
    print(f"\n{title}")
    print(f"{'loop':>5} {'CAS':>10} {'target ovl':>10} {'uniform':>20} {'adaptive':>22} {'ratio':>7}")
    print("-" * 80)
    for q in data["matched_quality"]:
        u, a = q["uniform_micro"], q["adaptive_micro"]
        ratio = f"{u/a:.2f}x" if (u and a) else "--"
        us = f"{u} ({q['uniform_setting']})" if u else "none passes"
        as_ = f"{a} ({q['adaptive_setting'].split(',')[0]})" if a else "none passes"
        print(f"{q['loop']:>5} {q['cas']:>10} {q['target']:>10.2f} {us:>20} {as_:>22} {ratio:>7}")

frontier_table(form, "FORMALDIMINE / STO-3G  -- cost (micro-iterations) at matched quality")
frontier_table(buta, "BUTADIENE / 6-31G*     -- cost (micro-iterations) at matched quality")
""")

code(r"""
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), squeeze=False)
for ax, data, title in zip(axes[0], (form, buta), ("formaldimine", "butadiene")):
    if not data:
        ax.text(0.5, 0.5, f"{title}\nnot run yet", ha="center", va="center"); ax.axis("off")
        continue
    for method, colour, marker in (("uniform", "steelblue", "o"), ("adaptive", "crimson", "s")):
        rows = [r for r in data["rows"] if r["method"] == method and r["status"] == "OK"
                and r["loop"].endswith("_x")]
        if not rows: continue
        ax.plot([r["micro"] for r in rows], [r["min_overlap"] for r in rows],
                marker, color=colour, ms=6, ls="none", label=method)
    ax.axhline(0.80, color="0.3", ls="--", lw=1, label="continuity threshold")
    ax.set_xscale("log"); ax.set_xlabel("cost (micro-iterations)")
    ax.set_ylabel("worst adjacent overlap (quality)")
    ax.set_title(f"{title}: the CI-enclosing loop\nup and to the left is better")
    ax.legend(fontsize=7.5, loc="lower right")
plt.show()
""")

md(r"""
### Verdict

**On formaldimine adaptive stepping is a wash — between 0.60× and 1.23× — and sometimes worse than
uniform.** That is a negative result and it is reported as one. It is also exactly what should have
been expected: formaldimine's loops have *uniform* difficulty (worst adjacent overlap 0.98 against a
median of 0.99), so there is nothing for step control to exploit, and it pays a small overhead for
carrying margin the loop does not need.

The prediction made before running it — from the spread of adjacent overlaps on saved runs — was
1.3× for formaldimine, 3.1× for butadiene `B_x` at CAS(8,8), and 4.7× for ethylene at CAS(8,8). The
formaldimine number is confirmed; the harder loops are where the case has to be made.

**Where step control does pay is not primarily cost.** It is:

* **loops that uniform stepping cannot walk at any affordable $N$** — §3 above, where the
  trustworthy range extends by more than an order of magnitude in closest approach. This is not a
  speed-up, it is the difference between an answer and a refusal, and it is what makes
  `locating_intersections.ipynb` possible;
* **not having to guess $N$ in advance.** A uniform run needs $N$ chosen before anything is known
  about the loop, and the project's own history contains several runs refused purely for having
  guessed it too low;
* **an automatic diagnosis** of *why* a loop failed — undersampled, or passing through something.
""")

with open(os.path.join(HERE, "adaptive_stepping.ipynb"), "w") as fh:
    nb = nbf.v4.new_notebook(cells=CELLS)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    nbf.write(nb, fh)
print(f"wrote adaptive_stepping.ipynb ({len(CELLS)} cells)")
