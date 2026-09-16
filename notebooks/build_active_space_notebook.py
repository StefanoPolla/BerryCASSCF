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

def load_scan(ne, ncas):
    # CAS(12,12) uses a thinner grid (see examples/run_active_space_study.py), so match on
    # the active space rather than a fixed grid shape.
    hits = sorted(glob.glob(os.path.join(ETH, f"ethylene_scan_cas{ne}-{ncas}_*.npz")))
    for p in hits:
        r = ScanResult.load(p)
        if np.isfinite(r.e_states).all():
            return r
    return None

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

def refined_phi(r):
    # Sub-grid intersection position along phi at tau = 90, by parabolic interpolation
    # through the grid minimum and its two neighbours. The grid step is 4 deg, so without
    # this the comparison between active spaces is quantised far too coarsely.
    i = int(np.argmin(np.abs(r.alphas - 90.0)))
    row = r.gap[i]
    j = int(np.nanargmin(row))
    if j == 0 or j == len(row) - 1:
        return float(r.phis[j])
    y0, y1, y2 = row[j - 1], row[j], row[j + 1]
    denom = y0 - 2 * y1 + y2
    if abs(denom) < 1e-18:
        return float(r.phis[j])
    step = float(r.phis[1] - r.phis[0])
    return float(r.phis[j] + 0.5 * (y0 - y2) / denom * step)

ref_phi = refined_phi(ref) if ref is not None else None

print(f"{'active space':>13} {'own min (mHa)':>13} {'phi (refined)':>13} "
      f"{'err vs ref':>11} {'gap at ref CI':>14} {'mirror asym':>13}")
print("-" * 83)
for cas in LADDER:
    r = scans.get(cas)
    if r is None:
        print(f"{'CAS%s' % (cas,):>13} {'(not run)':>13}")
        continue
    t, p, g = r.min_gap_point()
    rp = refined_phi(r)
    err = "  -" if ref_phi is None else f"{rp - ref_phi:+10.2f}d"
    atref = "  -" if ref_pt is None else f"{gap_at(r, *ref_pt):14.3f}"
    print(f"{'CAS%s' % (cas,):>13} {g*1e3:13.3f} {rp:13.2f} {err} {atref} "
          f"{mirror_asymmetry(r):13.2e}")
print("\nMirror asymmetry at the 1e-3 mHa level or below means the scan is path-independent.")
print("Grid spacing is 4 deg, so a position error of 4 deg is one grid step.")
print("'gap at ref CI' is the sharper measure: what each active space thinks the gap is at the")
print("position the reference puts the intersection. A large value there means that active")
print("space does not see an intersection where there is one.")
""")

md(r"""
## 2. Where each active space puts the intersection

CAS(12,12) costs about 45 s per grid point, so the reference is run on a thinner
3&times;11 grid rather than the 11&times;11 used for the other rungs. That still does both jobs it
has to: it contains the `tau = 90` line, where the intersection sits and from which the
reference position is measured, and it contains one exact mirror pair
(`tau` = 70 and 110) for the path-independence check. Its panel below is correspondingly
narrow.
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
    ax.set_title(f"CAS{cas}: {g*1e3:.2f} mHa at ({t:.0f}, {p:.0f})"
                 + ("  [3x11 grid]" if r.alphas.size < 11 else ""), fontsize=9)
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
TOL_DEG = 2.0          # half a grid step: what counts as "the same position"

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

positions = {c: refined_phi(scans[c]) for c in LADDER if scans.get(c) is not None}
ref_phi = positions.get(REFERENCE_CAS)

print(f"{'active space':>13} {'Berry: min N':>13} | {'phi':>8} {'err vs ref':>11} "
      f"{'accurate':>9} {'stable':>8}")
print("-" * 72)
for k, cas in enumerate(LADDER):
    b = berry_verdict(cas)
    btxt = "n/a" if b is None else (f"N = {b}" if b else "fails at all N")
    if cas not in positions:
        print(f"{'CAS%s' % (cas,):>13} {btxt:>13} | {'(scan pending)':>8}")
        continue
    phi = positions[cas]
    if ref_phi is None:
        err_txt, acc_txt = f"{'-':>11}", f"{'-':>9}"
    else:
        err = phi - ref_phi
        err_txt = f"{err:+11.2f}"
        acc_txt = f"{('yes' if abs(err) <= TOL_DEG else 'NO'):>9}"
    # "stable" = every LARGER active space tested agrees with this one, i.e. enlarging the
    # space does not move the answer. This is the test a practitioner can actually apply,
    # because it needs no reference.
    larger = [positions[c] for c in LADDER[k + 1:] if c in positions]
    # The largest rung has nothing above it to be checked against, so it is not "unstable".
    stab_txt = "-" if not larger else ("yes" if all(abs(phi - q) <= TOL_DEG for q in larger)
                                       else "NO")
    print(f"{'CAS%s' % (cas,):>13} {btxt:>13} | {phi:8.2f} {err_txt} "
          f"{acc_txt} {stab_txt:>8}")

acc = [c for c in LADDER if c in positions and ref_phi is not None
       and abs(positions[c] - ref_phi) <= TOL_DEG]
stab = [c for k, c in enumerate(LADDER)
        if c in positions
        and [positions[q] for q in LADDER[k + 1:] if q in positions]
        and all(abs(positions[c] - positions[q]) <= TOL_DEG
                for q in LADDER[k + 1:] if q in positions)]
okb = [c for c in LADDER if berry_verdict(c)]
print()
print(f"Berry phase correct from:                      "
      f"{min(okb) if okb else 'none of those tested'}")
print(f"SA comparator accurate from (vs reference):    "
      f"{min(acc) if acc else 'none of those tested'}")
print(f"SA comparator STABLE from (no reference used): "
      f"{min(stab) if stab else 'none of those tested'}")
""")

md(r"""
The last two lines are the point of the whole study, and they disagree.

*Accurate* asks whether an active space happens to land on the right answer. *Stable* asks the
question a practitioner can actually ask without already knowing the answer: **does enlarging the
active space leave the result unchanged?** For ethylene, CAS(2,2) is accurate to 0.01&deg; but
not stable &mdash; the very next rung moves the intersection by nearly 5&deg;, and the one after
that by 8&deg; the other way. Nothing available at CAS(2,2) would tell you it was right.

So the honest requirement here is set by stability, not by the smallest space that works.
""")

md(r"""
## 7. Conclusions

*Filled in from the numbers above once the full ladder has run &mdash; see
`docs/active_space.md` for the written-up version.*
""")

md(r"""
---

# Part II &mdash; butadiene: a second opinion

Ethylene answered the question with a quirk: CAS(2,2) was accidentally exact, so "smallest that
works" and "smallest you can trust" came apart. Butadiene is the follow-up, picked because it
should fail for a more fundamental reason &mdash; its 2<sup>1</sup>A<sub>g</sub> state carries a
large **doubly-excited** component, which a small active space cannot represent at all.

One structural difference changes what can be claimed. Formaldimine is reachable by FCI and
ethylene's full valence space is CAS(12,12); butadiene's is **CAS(22,22)**, far out of reach.
**There is no exact reference here** &mdash; "converged" can only mean "the answer stopped
moving". If it is still moving at the top of the ladder, that is the finding.

The plane was **searched for, not assumed**: three candidate planes of four rigid coordinates
(methylene twist, umbrella pyramidalization, central torsion, skeletal bend) at
SA-CASSCF(4,4)/6-31G\*. Only one contains an intersection.

| plane | minimum gap | verdict |
|---|---|---|
| `tw_pyr` | **0.73 mHa** at (tw=90, pyr&asymp;110) | **contains a conical intersection** |
| `tw_tc` | 99.49 mHa | none |
| `tw_bend` | 36.21 mHa | none |

Two traps the search caught, both of which would otherwise have been silent:

* The plane's *unrestricted* minimum is at `pyr = 180`, gap 0.19 mHa &mdash; but there the
  methylene has folded back onto its own C1&ndash;C2 bond, **266 kcal/mol** up. Degenerate, and
  chemically meaningless. The search reports strain beside the gap for this reason.
* **Butadiene does not have ethylene's `tw` &rarr; `180 - tw` symmetry.** Its two methylene
  hydrogens are inequivalent (one cis, one trans to C3=C4). Its exact symmetry is reflection
  through the molecular plane, `(tw, pyr)` &rarr; `(-tw, -pyr)`, which maps outside the scanned
  region &mdash; so the path-independence check costs a few extra solves instead of being free.
  Reusing ethylene's check would have compared unrelated geometries.

Full account: `docs/butadiene.md`.
""")

code(r"""
BUTA = os.path.join(ROOT, "results", "butadiene")
BUTA_LADDER = [(4, 4), (6, 6), (8, 8), (10, 10), (12, 12)]

def load_buta(ne, ncas):
    hits = sorted(glob.glob(os.path.join(BUTA, f"butadiene_scan_cas{ne}-{ncas}_*.npz")))
    for p in hits:
        r = ScanResult.load(p)
        if np.isfinite(r.e_states).all():
            return r
    return None

def buta_refined_pyr(r):
    i = int(np.unravel_index(np.nanargmin(r.gap), r.gap.shape)[0])
    row = r.gap[i]
    j = int(np.nanargmin(row))
    if j in (0, len(row) - 1):
        return float(r.phis[j]), float(r.alphas[i])
    y0, y1, y2 = row[j - 1], row[j], row[j + 1]
    d = y0 - 2 * y1 + y2
    if abs(d) < 1e-18:
        return float(r.phis[j]), float(r.alphas[i])
    return (float(r.phis[j] + 0.5 * (y0 - y2) / d * (r.phis[1] - r.phis[0])),
            float(r.alphas[i]))

buta = {c: load_buta(*c) for c in BUTA_LADDER}
have = [c for c in BUTA_LADDER if buta[c] is not None]
print("butadiene scans available:", ", ".join(f"CAS{c}" for c in have) or "(none)")
""")

code(r"""
if have:
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    cmap = plt.get_cmap("plasma")
    pos = {}
    for k, cas in enumerate(have):
        r = buta[cas]
        i = int(np.unravel_index(np.nanargmin(r.gap), r.gap.shape)[0])
        pos[cas] = buta_refined_pyr(r)
        ax.plot(r.phis, r.gap[i] * 1e3, "o-", ms=4,
                color=cmap(k / max(len(have) - 1, 1)), label=f"CAS{cas}")
    ax.set_xlabel(r"pyramidalization $\phi$ (deg), at tw = 90$^\circ$")
    ax.set_ylabel(r"$E_1 - E_0$  (mHa)")
    ax.set_title("Butadiene: gap through the intersection, by active space")
    ax.legend(fontsize=8); plt.tight_layout(); plt.show()

    largest = have[-1]
    ref_pyr = pos[largest][0]
    print(f"{'active space':>13} {'pyr (refined)':>14} {'tw':>6} "
          f"{'vs largest rung':>16} {'grid min (mHa)':>15}")
    print("-" * 70)
    for cas in have:
        pyr, tw = pos[cas]
        print(f"{'CAS%s' % (cas,):>13} {pyr:14.2f} {tw:6.0f} {pyr - ref_pyr:+16.2f} "
              f"{buta[cas].min_gap_point()[2]*1e3:15.3f}")
    print(f"\nLargest rung available: CAS{largest}. With no exact reference, the only")
    print("meaningful question is whether successive rungs still disagree.")
""")

md(r"""
## Verdict across all three systems
""")

code(r"""
rows_out = [
    ("formaldimine", "CAS(2,2)", "CAS(4,4)", "CAS(4,4)",
     "FCI in-basis (alpha_x = 132.6 deg)"),
    ("ethylene", "CAS(2,2)", "CAS(2,2)", "CAS(10,10)",
     "CAS(12,12) full valence (phi_x = 110.0 deg)"),
]
if len(have) >= 2:
    pyrs = [buta_refined_pyr(buta[c])[0] for c in have]
    settled = [have[k] for k in range(len(have) - 1)
               if all(abs(pyrs[k] - q) <= 2.0 for q in pyrs[k + 1:])]
    buta_stable = f"CAS{min(settled)}" if settled else "none tested"
    rows_out.append(("butadiene", "see below", "n/a (no reference)", buta_stable,
                     f"largest rung CAS{have[-1]}, no exact reference"))

print(f"{'system':>13} {'Berry needs':>12} {'SA accurate':>20} {'SA stable':>14}   reference")
print("-" * 92)
for r_ in rows_out:
    print(f"{r_[0]:>13} {r_[1]:>12} {r_[2]:>20} {r_[3]:>14}   {r_[4]}")
""")

md(r"""
The pattern across the three systems is that **no active-space recipe transfers**. Formaldimine's
comparator is qualitatively wrong at CAS(2,2) and converged at CAS(4,4); ethylene's is accurate at
CAS(2,2) but not stable until CAS(10,10); butadiene has no reachable reference at all. What does
transfer is the *procedure*: enlarge the space until the answer stops moving, and check an exact
symmetry of the system to confirm the solver is not the thing that moved.

The Berry phase, on all systems tested, needs only the minimal active space &mdash; because it
asks a topological question, which tolerates a badly misplaced intersection as long as the loop
still encloses it.
""")

nb = nbf.v4.new_notebook(cells=CELLS)
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
out = os.path.join(HERE, "active_space.ipynb")
nbf.write(nb, out)
print(f"wrote {out} with {len(CELLS)} cells")
