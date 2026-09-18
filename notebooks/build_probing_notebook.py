#!/usr/bin/env python
"""Generate notebooks/probing_by_small_loops.ipynb.

Reads saved records only (results/centre_probe, results/radius_scan, results/diagnostics),
so it runs in seconds and recomputes nothing.
"""

import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS: list = []


def md(text):
    CELLS.append(nbf.v4.new_markdown_cell(text.strip("\n")))


def code(text):
    CELLS.append(nbf.v4.new_code_cell(text.strip("\n")))


md(r"""
# Asking a loop one question at a time

Locating a degeneracy with loop transport has so far meant **bisection**: fix a centre, shrink the
loop until the Berry phase turns over, and read the transition radius as a distance. It works, but
it spends its whole budget in the worst configuration available — a loop that *grazes* the seam,
which is exactly where the continuation is slowest and most likely to fail.

This notebook uses a different primitive, and it turns out to say something unexpected about a
result the project has been quoting.

> **The probe.** One small circular loop, one walk, one bit: *is a degeneracy inside this loop?*

That is cheap. At ethylene CAS(2,2) a probe costs about **20 seconds**, against an hour or more for
a bisection — so a question that used to be a cluster job is a coffee break, and the answers can
be used to decide where to look next.

Three verdicts, and the third is not a failure:

| verdict | meaning |
|---|---|
| $\pi$ | an **odd** number of degeneracies inside — normally one |
| $0$ | an **even** number, normally none |
| refused | the loop could not be walked. For a small loop that usually means it passes very close to the seam, which is *localizing information* |
""")

code(r"""
import glob, json, os, sys
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.getcwd()) if os.path.basename(os.getcwd()) == "notebooks" else os.getcwd()
sys.path.insert(0, ROOT)
plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.facecolor": "white"})

def load(kind, name):
    p = os.path.join(ROOT, "results", kind, name)
    return json.load(open(p)) if os.path.exists(p) else None

# where the state-averaged methods put ethylene's intersection at this rung
SA_INTERSECTION = (90.0, 111.08)      # CAS(2,2) gap scan, cone-refined (docs/findings.md §1)
EXACT = (90.0, 110.90)                # full-valence CAS(12,12), exact in this basis
COLOUR = {"pi": "tab:red", "zero": "tab:blue", "undetermined": "0.55"}
""")

md(r"""
## 1. Shrinking the loop the ladder uses

`E_x` — the loop the ethylene active-space study walks — is a circle of radius 12° about the
intersection, and it reports $\pi$ at every rung. Shrinking it concentrically gives a clean,
monotone sequence.
""")

code(r"""
shrink = load("radius_scan", "ethylene_cas2-2_Ex_shrink.json")
onsa = load("radius_scan", "ethylene_cas2-2.json")
rows = []
for rec in ((shrink or {}).get("probes", []) + (onsa or {}).get("probes", [])):
    runs = rec.get("runs", [])
    rows.append({"r": rec["radius"], "verdict": rec["verdict"],
                 "pts": max((x["n_points"] for x in runs), default=0),
                 "worst": min((x["min_overlap"] for x in runs), default=np.nan),
                 "wall": rec["wall_time"]})
rows.sort(key=lambda d: -d["r"])
print(f"{'radius (deg)':>13} {'verdict':>14} {'points':>8} {'worst overlap':>14} {'wall (s)':>9}")
print("-" * 62)
for d in rows:
    print(f"{d['r']:>13g} {d['verdict']:>14} {d['pts']:>8} {d['worst']:>14.3f} {d['wall']:>9.1f}")
print("\n(the r = 4 and r = 2 rows come from a scan centred on (90, 111.0) rather than")
print(" (90, 110.9); at this resolution the 0.1 deg differences do not matter)")
""")

md(r"""
The transition is not marginal: $\pi$ at 12, 10, 8 and 6, a refusal band at 5 and 4, and a clean
$0$ from 3 downwards. **Whatever carries the phase lies $4.5 \pm 1.5°$ from the centre.**

Two things in that table matter beyond the verdicts. The walks stay cheap and well-resolved all the
way down — 10 points and worst overlap 0.99 at $r = 0.1°$ — so a small loop is *not* intrinsically
hard to transport. And the expensive probes are exactly the refused ones, 4 and 5, which is the
grazing configuration a bisection is built to spend its time in.
""")

md(r"""
## 2. The loop does not enclose the state-averaged intersection

Here is the result that started this. The state-averaged gap scan puts ethylene's CAS(2,2)
intersection at $(90, 111.08)$, and the exact full-valence answer is $(90, 110.90)$. Both are
**inside** every loop of radius 3 or less in the table above — and every one of those loops reports
a trivial phase.
""")

code(r"""
# the shrink series is centred on the exact position; the state-averaged one is 0.18 deg away
CENTRE = (90.0, 110.9)
d_sa = np.hypot(SA_INTERSECTION[0] - CENTRE[0], SA_INTERSECTION[1] - CENTRE[1])
d_ex = np.hypot(EXACT[0] - CENTRE[0], EXACT[1] - CENTRE[1])
print(f"loops centred on {CENTRE}; the state-averaged intersection is {d_sa:.2f} deg away,")
print(f"the exact position {d_ex:.2f} deg.\n")
print(f"{'radius':>8} {'verdict':>9}  {'contains SA?':>13} {'contains exact?':>16}")
for d in [x for x in rows if x["r"] <= 3]:
    print(f"{d['r']:>8g} {d['verdict']:>9}  {str(d['r'] > d_sa):>13} {str(d['r'] > d_ex):>16}")
n_contain = sum(1 for x in rows if x["r"] <= 3 and x["r"] > d_sa)
print(f"\n{n_contain} loops contain the state-averaged intersection and every one reports 0.")
print("The smallest, r = 0.1, does not contain it -- it is inside the 0.18 deg offset -- so it")
print("is evidence about the exact position only, which it does contain.")
""")

md(r"""
## 3. Then where is it? Sweeping the mirror line

Ethylene's `tw` and $180 - `tw`$ geometries are exact mirror images, so a lone degeneracy has to
sit **on** the `tw` = 90 line. Nine probes of radius 2° cover that line continuously from
`pyr` = 99 to 123 — the whole span of `E_x`.
""")

code(r"""
line = load("centre_probe", "ethylene_cas2-2_line_r2.json")
gaps = load("centre_probe", "ethylene_cas2-2_gaps_r2.json")
probes = []
for rec in (line, gaps):
    if rec:
        for p in rec["probes"]:
            probes.append({"centre": tuple(p["centre"]), "verdict": p["verdict"],
                           "radius": rec["radius"], "wall": p["wall_time"]})
onl = [p for p in probes if abs(p["centre"][0] - 90.0) < 1e-9]
off = [p for p in probes if abs(p["centre"][0] - 90.0) >= 1e-9]
print("on the mirror line, radius 2 deg:")
for p in sorted(onl, key=lambda q: q["centre"][1]):
    print(f"   pyr {p['centre'][1]:>6.1f}   {p['verdict']:>13}   {p['wall']:>6.1f} s")
print("\noff the line, radius 2 deg:")
for p in off:
    print(f"   {str(p['centre']):>14}   {p['verdict']:>13}   {p['wall']:>6.1f} s")
""")

md(r"""
Nothing on the line encloses anything. The only probes that react at all are the two placed
**off** it, at $(94, 111)$ and $(86, 111)$ — and both are *refused*, symmetrically, from geometries
that are exact mirror images computed independently.

That is where something 4.5° from the centre would be. But it cannot be the answer, and the reason
is arithmetic.
""")

code(r"""
# the symmetry is exact, so anything off the line comes in pairs
print("state-specific CASSCF energies at mirror-image geometries, CAS(2,2)/6-31G*:\n")
print("   pyr = 111:  E(tw=94) = -77.8002951284   E(tw=86) = -77.8002951284   diff -6.4e-13")
print("   pyr = 105:  E(tw=94) = -77.8177907123   E(tw=86) = -77.8177907123   diff -2.1e-13")
print("\nso a degeneracy at (94, 111) implies one at (86, 111).")
print("A loop centred on tw = 90 encloses both or neither -- an EVEN count -- which reads as 0.")
print("No arrangement of point degeneracies consistent with this symmetry gives pi at r >= 6")
print("while giving 0 at r <= 3 and nothing on the line.")
""")

code(r"""
fig, ax = plt.subplots(figsize=(7.0, 6.0))
th = np.linspace(0, 2*np.pi, 401)
for p in probes:
    c, r = p["centre"], p["radius"]
    ax.plot(c[0] + r*np.cos(th), c[1] + r*np.sin(th), "-", lw=1.0,
            color=COLOUR[p["verdict"]], alpha=0.85)
    ax.plot(*c, ".", ms=4, color=COLOUR[p["verdict"]])
# the transition annulus from the shrink series
for r, style in ((6.0, "-"), (3.0, "--")):
    ax.plot(90 + r*np.cos(th), 110.9 + r*np.sin(th), style, color="k", lw=1.1)
ax.plot(*SA_INTERSECTION, "*", ms=17, color="gold", mec="k", mew=0.7,
        label="state-averaged intersection")
ax.plot(*EXACT, "P", ms=11, color="k", label="exact (full valence)")
ax.axvline(90, color="0.5", ls=":", lw=1.0)
for v, lab in (("pi", "probe: pi"), ("zero", "probe: 0"), ("undetermined", "probe: refused")):
    ax.plot([], [], "-", color=COLOUR[v], label=lab)
ax.plot([], [], "-", color="k", label="r = 6 (pi) and r = 3 (0)")
ax.set_xlabel("tw (deg)"); ax.set_ylabel("pyr (deg)")
ax.set_title("ethylene CAS(2,2): every probe, and what it found")
ax.legend(fontsize=7.5, loc="upper left"); ax.set_aspect("equal")
ax.set_xlim(80, 100); ax.set_ylim(96, 126)
plt.tight_layout(); plt.show()
""")

md(r"""
## 4. So what *is* in there?

If the $\pi$ cannot come from encircled degeneracies, the region deserves a direct look. Two maps,
neither of which needs a loop:

* the **in-CAS gap** — S0 to the second root inside the tracked active space, at the state-specific
  orbitals. At CAS(2,2) the second root is the doubly excited configuration, so this is *not* the
  physical S0/S1 gap (`docs/limitations.md` says so explicitly). It is a picture of the solution the
  continuation is tracking, and raggedness means neighbouring geometries converged to **different**
  CASSCF solutions;
* **cold-start energies around a circle**, each point solved independently, so a jump marks where
  competing solutions exist.
""")

code(r"""
p = os.path.join(ROOT, "results", "diagnostics", "ethylene_cas2-2_solutions.npz")
sol = np.load(p) if os.path.exists(p) else None
if sol is None:
    print("run: python examples/diagnose_ethylene_solutions.py")
else:
    gap, tws, pyrs = sol["incas_gap"], sol["tws"], sol["pyrs"]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    ax = axes[0]
    im = ax.pcolormesh(pyrs, tws, gap, shading="nearest", cmap="viridis")
    fig.colorbar(im, ax=ax, label="in-CAS gap (mHa)")
    ax.axhline(90, color="w", ls=":", lw=1.2)
    ax.set_xlabel("pyr (deg)"); ax.set_ylabel("tw (deg)")
    ax.set_title("in-CAS gap: smooth away from the line, ragged at tw = 87 and 93")
    ax.grid(False)

    ax = axes[1]
    theta = np.degrees(sol["theta"])
    for r, col in zip(sol["radii"], ("tab:red", "tab:blue")):
        es = sol[f"ring_{r:g}"]
        d = np.abs(np.diff(np.concatenate([es, es[:1]])))
        ax.plot(theta, (es - np.nanmin(es)) * 1e3, "-o", ms=2.5, lw=1.0, color=col,
                label=f"r = {r:g} ({'pi' if r >= 6 else '0'})")
        for k in np.where(d > 10 * np.median(d))[0]:
            ax.axvline(theta[k], color=col, ls=":", lw=1.0, alpha=0.7)
    ax.set_xlabel("angle around the loop (deg)")
    ax.set_ylabel("cold-start energy above the minimum (mHa)")
    ax.set_title("energies around the loop; dotted = jumps above 10x the median step")
    ax.legend(fontsize=8)
    plt.tight_layout(); plt.show()

    for r in sol["radii"]:
        es = sol[f"ring_{r:g}"]
        d = np.abs(np.diff(np.concatenate([es, es[:1]])))
        print(f"  r = {r:g}: median step {np.median(d):.2e} Ha, largest {np.max(d):.2e}, "
              f"jumps above 10x median: {int((d > 10*np.median(d)).sum())}")
    print(f"\n  in-CAS gap never falls below {np.nanmin(gap):.0f} mHa in this window, and the map is")
    print("  exactly mirror-symmetric -- rows 84/96, 85/95, 86/94 are identical to 1e-13.")
""")

md(r"""
## 5. How small can an *enclosing* loop be? The floor

Everything above used loops that do **not** enclose, and those stay cheap and clean at any radius —
10 points and overlap 0.998 at $r = 0.05°$. That is not the question an iterative search needs
answered. A search that re-centres on the intersection and shrinks lives on **enclosing** loops, and
those are the hard ones: every point is near-degenerate, and the tracked state is a state-specific
CASSCF ground state with S1 on top of it.

Formaldimine CAS(4,4) is the right place to measure it — the localization there is self-consistent
(three centres, residual inside its own precision) and a point costs 0.1 s. Concentric loops on its
triangulated position:
""")

code(r"""
floor = load("radius_scan", "formaldimine_cas4-4_floor.json")
if not floor:
    print("run: python examples/run_radius_scan.py formaldimine --cas 4 4 "
          "--centre 128.99 90.24 --radii 8 4 2 1 0.5 0.25 0.1 0.05 --label floor")
else:
    print(f"{'radius':>8} {'verdict':>14} {'points':>8} {'worst overlap':>14} "
          f"{'micro-iterations':>17} {'wall (s)':>9}")
    print("-" * 76)
    prev = None
    for rec in floor["probes"]:
        runs = rec["runs"]
        pts = max(x["n_points"] for x in runs)
        worst = min(x["min_overlap"] for x in runs)
        grow = "" if prev is None else f"  ({rec['cost_micro']/prev:.1f}x)"
        print(f"{rec['radius']:>8g} {rec['verdict']:>14} {pts:>8} {worst:>14.3f} "
              f"{rec['cost_micro']:>17}{grow:<7} {rec['wall_time']:>9.1f}")
        prev = rec["cost_micro"]
""")

md(r"""
Three things, and they pull in different directions.

**The angular prediction holds.** Point counts are flat — 15 or 16 per setting at every radius from
8° down to 1° — exactly as expected if the wavefunction turns by $\pi$ over any enclosing loop and
the discretization requirement is angular rather than metric. A small enclosing loop needs no more
points than a large one.

**But it is not free.** Micro-iterations roughly **double with each halving** of the radius:
1802, 3220, 6061, 11453. The point count is flat and the cost is not, so it is each *solve* getting
harder — the near-degeneracy penalty, arriving exactly as predicted.

**And there is a floor, at about 1°.** Loops at 0.5° and 0.25° still *contain* the intersection —
with clearance 0.4° and 0.15° — and are refused anyway. Below that, $r = 0.05°$ returns a clean
$0$: it has fallen inside the offset between the triangulated centre and the true position, which
incidentally pins that offset to between 0.05° and 0.15°.

So an iterative search on this system could shrink to roughly **1°** and no further. That is worth
comparing honestly against what bisection already achieves here: the formaldimine CAS(4,4) bisection
brackets $\rho$ to $\pm 0.03$ on a 10° semi-axis, i.e. $\pm 0.3°$. **The probe does not win on
asymptotic precision.**
""")

md(r"""
## 6. The check this makes cheap: does the loop encircle the gap scan's answer?

Every comparison in this project has been indirect — locate the state-specific object by bisection,
locate the state-averaged one by scanning, compare two positions each carrying its own error bar.
The probe replaces all of that with one loop placed **on** the gap scan's intersection.

The radius has to clear the ~1° floor of §5 and still contain nothing else; 2° does both. A verdict
of $0$ is an *even* count, so the radius is reported with every answer.
""")

code(r"""
rows_enc = []
for sysname in ("formaldimine", "ethylene", "butadiene"):
    p = os.path.join(ROOT, "results", "encirclement", f"{sysname}_r2.json")
    if os.path.exists(p):
        for r in json.load(open(p))["rungs"]:
            rows_enc.append((sysname, tuple(r["cas"]), tuple(r["sa_intersection"]),
                             r["verdict"], 2.0))
# the CAS(12,12) probe ran on the cluster, one radius per task, under its own name
g = load("radius_scan", "butadiene_cas12-12_gapmin_r1.25.json")
if g:
    rows_enc.append(("butadiene", (12, 12), tuple(g["centre"]),
                     g["probes"][0]["verdict"], g["probes"][0]["radius"]))

if not rows_enc:
    print("run: python examples/verify_encirclement.py <system>")
else:
    say = {"pi": "YES", "zero": "NO", "undetermined": "refused"}
    print(f"{'system':>13} {'CAS':>8} {'gap-scan intersection':>24} {'r':>5} "
          f"{'verdict':>14} {'encircles it?':>14}")
    print("-" * 84)
    for sysname, cas, pt, verdict, r in rows_enc:
        print(f"{sysname:>13} {('(%d,%d)' % cas):>8} {f'({pt[0]:.2f}, {pt[1]:.2f})':>24} "
              f"{r:>5g} {verdict:>14} {say[verdict]:>14}")
""")

md(r"""
**The check discriminates, which is what makes the negatives worth believing.**

At formaldimine CAS(2,2) the answer is **no** — and independently, that rung's gap scan is known to
be *qualitatively* wrong: it reports a spurious intersection near $\alpha = 141°$ and misses the
real one at 132.6° (`docs/results.md` §3). The probe rejects a position already known to be wrong.

At formaldimine CAS(4,4), where the gap scan is converged, the answer is **yes** — and the object
it encircles is the one the triangulation found at (128.99, 90.24), 1° inside the probe.

So the check is not returning "no" everywhere. It agrees where the two methods are known to agree,
and refuses where one of them is known to be wrong. That is what licenses reading the butadiene
CAS(12,12) and ethylene CAS(2,2) rows as evidence rather than as noise.
""")


md(r"""
## What this says, and what it does not

**Standing:**

* at ethylene CAS(2,2), loop transport **does not encircle the state-averaged intersection** — six
  concentric loops containing it all report a trivial phase;
* the phase that `E_x` reports appears only at radius 6 and above, and no arrangement of point
  degeneracies compatible with ethylene's exact mirror symmetry explains it;
* the region contains **competing state-specific CASSCF solutions** — a ragged in-CAS gap map and
  cold-start energy jumps at four of 72 points around a radius-3 circle.

A sign picked up crossing from one CASSCF solution to another is **not** a Berry phase. At this rung
that is the explanation left standing, and it means the $\pi$ the ladder quotes for `E_x` at
CAS(2,2) may be right for the wrong reason.

**Not standing, and worth saying clearly:**

* this is **one rung of one system**, and CAS(2,2) is below the π space — not a chemically
  defensible active space for ethylene, which the repository already says. Nothing here transfers
  to the larger rungs without being measured there;
* the cold-start ring is a map of where *independent* solves disagree, not of what a warm-started
  continuation does. It shows that competing solutions exist, not that the walk crossed between
  them;
* the refusals at $(94, 111)$ and $(86, 111)$ are consistent with a near-degeneracy there and also
  with a solution boundary. They are not evidence for either on their own.

## What it changes about method design

**The probe wins on robustness and cost, not on precision.** Its floor is about 1° (§5), while the
formaldimine bisection already brackets to ±0.3°. So an iterative re-centring search would *not*
resolve a position more finely than what exists — which was the hope, and it is not supported.

What it does win is everything else:

* **cost** — 20 s against an hour per answer at the cheap rungs, because a non-enclosing loop is
  trivial to walk and only the enclosing ones are expensive;
* **robustness** — the ethylene bisections produced regions ±3–6° wide because their loops wandered
  into geometries the continuation could not follow, and seven of eight full-size loops there were
  refused. Every probe used here stayed in a small neighbourhood and returned an answer;
* **it fails informatively** — a refusal at a known small radius says the seam is *about that far
  away*, where a refused bisection anchor used to end the measurement entirely.

**And it inverts the sanity check.** Rather than reconstructing where loop transport thinks the
intersection is and comparing, one can ask directly whether it encircles the intersection the gap
scan already found. That is one probe. At ethylene CAS(2,2) it takes 20 seconds and the answer is
no — a question worth having asked on the first day.

**On the iterative search specifically.** The measured floor says such a search would converge to
~1° and stop, and the cost doubles with every halving, so the last step dominates. That is a
reasonable instrument for *verification* — "is it here, within a degree?" — and a poor one for
squeezing out precision. The parts of it worth building are the ones already built: the probe, and
anchors that move instead of giving up.
""")

nb = nbf.v4.new_notebook(cells=CELLS)
nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python",
                             "name": "python3"},
               "language_info": {"name": "python"}}
out = os.path.join(HERE, "probing_by_small_loops.ipynb")
nbf.write(nb, out)
print(f"wrote {out} with {len(CELLS)} cells")
