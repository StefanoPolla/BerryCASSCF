#!/usr/bin/env python
"""Generate notebooks/active_space.ipynb.

Reads saved results only; runs no electronic structure, so it executes in seconds.
Companion to results.ipynb, which covers the formaldimine benchmark.
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
# How large an active space does each method actually need?

The formaldimine benchmark (see `results.ipynb`) settled on **CAS(2,2)** for the Berry phase and
**CAS(4,4)** for the state-averaged comparator. Those are small. This notebook asks whether
that is a property of the methods or just of an easy molecule, by running the same two
workflows over a ladder of active spaces on a harder system.

**System: ethylene**, at its twisted-pyramidalized S<sub>1</sub>/S<sub>0</sub> conical
intersection &mdash; the textbook case in nonadiabatic photochemistry. Two rigid coordinates:

| | |
|---|---|
| `tau` | torsion of one CH<sub>2</sub> about the C=C axis (deg) |
| `phi` | umbrella pyramidalization of that same CH<sub>2</sub> (deg) |

Both are rigid rotations, so every bond length and both HCH angles are preserved exactly.

**Everything except the active space is held fixed** &mdash; same basis, same geometry family,
and above all the *same loops*. Letting the loops move with the active space would mean each
CAS was being asked a different question.

Ladder: CAS(2,2), (4,4), (6,6), (8,8), (10,10), (12,12). The last spans the full valence space
&mdash; only the two carbon 1s orbitals stay frozen &mdash; and serves as the in-basis reference.
""")

code(r"""
import glob, json, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

sys.path.insert(0, os.path.abspath(".."))
from berrycasscf.ethylene import CI_REGION, DEFAULT_BASIS
from berrycasscf.report import load_berry_records, berry_table, format_table
from berrycasscf.scan import ScanResult

ROOT = os.path.abspath("..")
ETH = os.path.join(ROOT, "results", "ethylene")

LADDER = [(2, 2), (4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]
REFERENCE_CAS = (12, 12)

plt.rcParams.update({"figure.dpi": 120, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.facecolor": "white"})

def load_scan(ne, ncas, grid="11x11"):
    p = os.path.join(ETH, f"ethylene_scan_cas{ne}-{ncas}_{grid}.npz")
    if not os.path.exists(p):
        return None
    r = ScanResult.load(p)
    return r if np.isfinite(r.e_states).all() else None

def mirror_asymmetry(r):
    g = r.gap * 1000.0
    return float(np.nanmax(np.abs(g - g[::-1, :])))

scans = {cas: load_scan(*cas) for cas in LADDER}
available = [c for c, r in scans.items() if r is not None]
print("scans available:", ", ".join(f"CAS{c}" for c in available))
""")

md(r"""
## 1. A trap worth naming first: active-space drift

Before any result, a methodological point that nearly corrupted this study.

The gap scan warm-starts each grid point from its neighbour, purely for speed. For small active
spaces that is harmless. For larger ones it is not: the *active space itself* can drift along
the scan path, so the answer depends on the route taken to a geometry rather than on the
geometry.

Ethylene gives a clean way to catch this. Geometries at `tau` and `180 - tau` are exact mirror
images &mdash; verified by comparing their full interatomic-distance spectra &mdash; so a
correct gap map **must** be symmetric about `tau = 90`. Any asymmetry is drift, not physics.

With naive warm starting, CAS(8,8)/6-31G\* returned **19.1 mHa** and **29.4 mHa** at two
geometries that are mirror images of each other. A cold start at either gives **26.9 mHa**.

The fix (`ScanConfig.strategy = "best"`) runs both a warm and a cold start at every point and
keeps whichever reaches the lower state-averaged energy. The scans below are checked against
the mirror symmetry, and the residual asymmetry is reported for each.
""")

code(r"""
ref = scans.get(REFERENCE_CAS)
ref_pt = ref.min_gap_point()[:2] if ref is not None else None

def gap_at(r, tau, phi):
    # gap this active space predicts at a given point, in mHa
    i = int(np.argmin(np.abs(r.alphas - tau)))
    j = int(np.argmin(np.abs(r.phis - phi)))
    return r.gap[i, j] * 1e3

print(f"{'active space':>13} {'own min (mHa)':>13} {'at (tau, phi)':>17} "
      f"{'err vs ref':>11} {'gap at ref CI':>14} {'mirror asym':>13}")
print("-" * 89)
for cas in LADDER:
    r = scans.get(cas)
    if r is None:
        print(f"{'CAS%s' % (cas,):>13} {'(not run)':>13}")
        continue
    t, p, g = r.min_gap_point()
    err = "  -" if ref_pt is None else f"{np.hypot(t-ref_pt[0], p-ref_pt[1]):9.1f} d"
    atref = "  -" if ref_pt is None else f"{gap_at(r, *ref_pt):14.3f}"
    print(f"{'CAS%s' % (cas,):>13} {g*1e3:13.3f} {f'({t:.1f}, {p:.1f})':>17} {err} {atref} "
          f"{mirror_asymmetry(r):13.2e}")
print("\nMirror asymmetry at the 1e-3 mHa level or below means the scan is path-independent.")
print("Grid spacing is 4 deg, so a position error of 4 deg is one grid step.")
print("'gap at ref CI' is the sharper measure: what each active space thinks the gap is at the")
print("position the reference puts the intersection. A large value there means that active")
print("space does not see an intersection where there is one.")
""")

md(r"""
## 2. Where each active space puts the intersection
""")

code(r"""
avail = [c for c in LADDER if scans.get(c) is not None]
ncol = 3
nrow = int(np.ceil(len(avail) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.4 * nrow), squeeze=False)
ref = scans.get(REFERENCE_CAS)
ref_pt = ref.min_gap_point()[:2] if ref is not None else None

for k, cas in enumerate(avail):
    ax = axes[k // ncol][k % ncol]
    r = scans[cas]
    T, P = np.meshgrid(r.alphas, r.phis, indexing="ij")
    pcm = ax.pcolormesh(P, T, r.gap * 1e3, shading="nearest", cmap="viridis_r")
    fig.colorbar(pcm, ax=ax, label="mHa")
    t, p, g = r.min_gap_point()
    ax.plot(p, t, "wo", ms=8, mec="k", zorder=6)
    if ref_pt is not None:
        ax.plot(ref_pt[1], ref_pt[0], "kx", ms=10, mew=2, zorder=7)
    ax.set_title(f"CAS{cas}: {g*1e3:.2f} mHa at ({t:.0f}, {p:.0f})", fontsize=9)
    ax.set_xlabel(r"$\phi$ (deg)"); ax.set_ylabel(r"$\tau$ (deg)")
for k in range(len(avail), nrow * ncol):
    axes[k // ncol][k % ncol].set_visible(False)
fig.suptitle("SA-CASSCF S$_1-$S$_0$ gap for ethylene across the active-space ladder\n"
             "white circle = that CAS's minimum;  black cross = reference CAS%s minimum"
             % (REFERENCE_CAS,), y=1.02)
plt.tight_layout(); plt.show()
""")

md(r"""
## 3. The gap through the intersection

A cut along `tau = 90` shows whether each active space produces a genuine cone (linear in the
coordinate, closing to zero) or merely a shallow avoided crossing.
""")

code(r"""
fig, ax = plt.subplots(figsize=(6.4, 3.8))
cmap = plt.get_cmap("viridis")
for k, cas in enumerate(avail):
    r = scans[cas]
    i = int(np.argmin(np.abs(r.alphas - 90.0)))
    ax.plot(r.phis, r.gap[i] * 1e3, "o-", ms=3.5,
            color=cmap(k / max(len(avail) - 1, 1)),
            lw=2.2 if cas == REFERENCE_CAS else 1.3,
            label=f"CAS{cas}" + (" (reference)" if cas == REFERENCE_CAS else ""))
ax.set_xlabel(r"$\phi$ (deg),  at $\tau = 90^\circ$")
ax.set_ylabel(r"$E_1 - E_0$  (mHa)")
ax.set_title("Gap through the intersection, by active space")
ax.legend(fontsize=8); plt.tight_layout(); plt.show()
""")

md(r"""
## 4. Berry phase across the same ladder

The loops are **fixed** for every active space, centred on the intersection located by the
reference CAS. `E_x` encloses it; `E_1` and `E_2` are displaced along `tau` and enclose nothing.
Each is run at two discretizations, so the stability criterion from the formaldimine study
(same phase at every `N`, every run `OK`) applies unchanged.
""")

code(r"""
records = load_berry_records(ETH)
rows = berry_table(records)
if rows:
    print(format_table(rows, ["loop", "CAS", "N", "product", "endpoint",
                              "min|ovl|", "phase", "status"]))
else:
    print("No Berry records yet. Run: python examples/run_active_space_study.py berry")
""")

code(r"""
if rows:
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    colour = {"E_1": "tab:blue", "E_x": "tab:red", "E_2": "tab:green"}
    order = [f"CAS({n},{m})" for n, m in LADDER]
    xs = np.arange(len(order))
    for loop in ["E_1", "E_x", "E_2"]:
        for n, mk in [(13, "o"), (21, "s")]:
            y, x = [], []
            for k, cl in enumerate(order):
                hit = [r for r in rows if r["loop"] == loop and r["CAS"] == cl and r["N"] == n]
                if hit:
                    y.append(hit[0]["product"]); x.append(k)
            if y:
                ax.plot(x, y, mk + "-", color=colour[loop], ms=5, alpha=0.85,
                        label=f"{loop}, N={n}")
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(xs); ax.set_xticklabels([o.replace("CAS", "") for o in order])
    ax.set_xlabel("active space"); ax.set_ylabel(r"product estimator $\Pi$")
    ax.set_title(r"Berry phase vs active space   ($\Pi<0$ is non-trivial)")
    ax.legend(fontsize=7, ncol=3); plt.tight_layout(); plt.show()
""")

md(r"""
### A free consistency check

The two control loops are mirror images of each other: `E_1` is centred at
`tau = 90 - 30` and `E_2` at `tau = 90 + 30`, and the geometries at `tau` and `180 - tau` are
exactly isometric. So **every** quantity computed on `E_1` must equal its counterpart on `E_2`,
digit for digit. Nothing in the continuation knows this &mdash; the two loops are traversed
independently, through different geometries, with independently converged orbitals and
independently chosen CI signs. Agreement is therefore a genuine end-to-end test of the whole
pipeline on a molecule it was not developed on.
""")

code(r"""
if rows:
    print(f"{'CAS':>10}  {'N':>3}   {'Pi(E_1)':>10}  {'Pi(E_2)':>10}   {'|difference|':>12}")
    print("-" * 56)
    worst = 0.0
    for cas in LADDER:
        cl = f"CAS({cas[0]},{cas[1]})"
        for n in sorted({r["N"] for r in rows}):
            a = [r for r in rows if r["CAS"] == cl and r["loop"] == "E_1" and r["N"] == n]
            b = [r for r in rows if r["CAS"] == cl and r["loop"] == "E_2" and r["N"] == n]
            if not (a and b):
                continue
            d = abs(a[0]["product"] - b[0]["product"])
            worst = max(worst, d)
            print(f"{cl:>10}  {n:>3}   {a[0]['product']:>10.6f}  {b[0]['product']:>10.6f}"
                  f"   {d:>12.2e}")
    print(f"\nlargest mirror discrepancy across the ladder: {worst:.2e}")
""")

md(r"""
## 5. Discretization: bigger active spaces need finer loops

Before the verdict, one effect that would otherwise be mistaken for an active-space failure.

The continuity criterion requires every adjacent overlap around the loop to stay above 0.80.
How fast the wavefunction turns as the loop is traversed is itself active-space dependent: a
larger CAS has more room to rearrange, so the same loop at the same `N` produces smaller
adjacent overlaps. Below the threshold the run is reported as FAILED &mdash; correctly, because
the discretization is too coarse to establish continuity &mdash; and refining `N` fixes it.

This is a statement about the **loop discretization**, not about whether the active space can
describe the physics, and the two must not be conflated.
""")

code(r"""
if rows:
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    order = [f"CAS({n},{m})" for n, m in LADDER]
    for n, mk in [(13, "o"), (21, "s"), (31, "^")]:
        xs, ys = [], []
        for k, cl in enumerate(order):
            hit = [r for r in rows if r["loop"] == "E_x" and r["CAS"] == cl and r["N"] == n]
            if hit:
                xs.append(k); ys.append(hit[0]["min|ovl|"])
        if xs:
            ax.plot(xs, ys, mk + "-", ms=5, label=f"N = {n}")
    ax.axhline(0.80, color="crimson", ls="--", lw=1.2, label="continuity threshold")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([o.replace("CAS", "") for o in order])
    ax.set_xlabel("active space"); ax.set_ylabel(r"min $|\langle\Psi_{k-1}|\Psi_k\rangle|$")
    ax.set_title("Loop E$_x$: adjacent overlap falls as the active space grows")
    ax.legend(fontsize=8); plt.tight_layout(); plt.show()
""")

md(r"""
## 6. Verdict: the smallest active space each method needs

Two separate questions, judged with the criteria used throughout this project.

The controls sit outside the scanned region, so they were checked separately: sampling the
SA-CASSCF(2,2) gap on and inside each (centre plus rings of 8 points at radius 6&deg; and
12&deg;) gives a minimum of **81.8 mHa** for both, against **0.177 mHa** at the intersection.
Nothing is enclosed.

**Berry phase.** For each active space, the smallest `N` at which the run passes every check
(all points converged, continuity held, endpoint returned to the same state, both estimators
agreeing), and whether the phase is then correct on all three loops: non-trivial on `E_x`,
trivial on `E_1` and `E_2`.

**SA comparator.** The active space must (a) put its gap minimum inside `E_x`, (b) place it
close to the reference active space's minimum, (c) drive the gap low enough to read as an
intersection, and (d) produce a mirror-symmetric map &mdash; a map that fails the symmetry test
is path-dependent and not trustworthy whatever number it reports.
""")

code(r"""
WANT = {"E_x": "non-trivial (pi)", "E_1": "trivial (0)", "E_2": "trivial (0)"}

def berry_verdict(cas):
    # Smallest N at which all three loops pass every check with the right phase.
    # None  -> rung not run, or not run completely (a missing record must never be
    #          mistaken for a failure); False -> every complete N failed.
    cl = f"CAS({cas[0]},{cas[1]})"
    sub = [r for r in rows if r["CAS"] == cl]
    complete_n = [n for n in sorted({r["N"] for r in sub})
                  if {r["loop"] for r in sub if r["N"] == n} == set(WANT)]
    if not complete_n:
        return None
    for n in complete_n:
        at_n = [r for r in sub if r["N"] == n]
        if all(r["status"] == "OK" for r in at_n) and all(
            r["phase"] == WANT[r["loop"]] for r in at_n
        ):
            return n
    return False

def scan_verdict(cas, tol_deg=8.0, tol_gap=8.0, tol_asym=1e-2):
    r, ref = scans.get(cas), scans.get(REFERENCE_CAS)
    if r is None or ref is None:
        return None
    t, p, g = r.min_gap_point()
    tr, pr, _ = ref.min_gap_point()
    checks = {
        "inside E_x": bool((t - tr) ** 2 + (p - pr) ** 2 <= LOOP_RADIUS ** 2),
        "near reference": abs(t - tr) <= tol_deg and abs(p - pr) <= tol_deg,
        "gap closes": g * 1e3 < tol_gap,
        "mirror-symmetric": mirror_asymmetry(r) < tol_asym,
    }
    return checks

LOOP_RADIUS = 12.0

print(f"{'active space':>13}  {'Berry: smallest working N':>26}   {'SA comparator':>14}   detail")
print("-" * 96)
for cas in LADDER:
    b = berry_verdict(cas)
    sv = scan_verdict(cas)
    btxt = "n/a" if b is None else (f"N = {b}" if b else "fails at every N tested")
    if sv is None:
        stxt, detail = "n/a", ""
    else:
        stxt = "correct" if all(sv.values()) else "WRONG"
        detail = ", ".join(k for k, v in sv.items() if not v) or "all checks pass"
    print(f"{'CAS%s' % (cas,):>13}  {btxt:>26}   {stxt:>14}   {detail}")

ok_b = [c for c in LADDER if berry_verdict(c)]
ok_s = [c for c in LADDER if (lambda v: v is not None and all(v.values()))(scan_verdict(c))]
print()
print(f"Berry phase correct from:   {min(ok_b) if ok_b else 'none of those tested'}")
print(f"SA comparator correct from: {min(ok_s) if ok_s else 'none of those tested'}")
""")

md(r"""
## 7. Conclusions

*Filled in from the numbers above once the full ladder has run &mdash; see
`docs/active_space.md` for the written-up version.*
""")

nb = nbf.v4.new_notebook(cells=CELLS)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
out = os.path.join(HERE, "active_space.ipynb")
nbf.write(nb, out)
print(f"wrote {out} with {len(CELLS)} cells")
