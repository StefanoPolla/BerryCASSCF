"""State-averaged CASSCF comparator: a two-dimensional gap scan.

This answers the same question as the Berry-phase workflow by a completely different route:
resolve the two lowest singlets with equal-weight SA-CASSCF on a grid covering the region
enclosed by the loop, and look for a point where ``E1 - E0`` goes to zero.

Deliberately *not* implemented: minimum-energy CI optimization. The brief scopes this to a
scan, and a scan is what supports the "is the CI inside the loop?" question without the
extra machinery and failure modes of a constrained optimizer.

Warm starting runs along grid rows in a serpentine order so that consecutive solves are at
neighbouring geometries. This is for speed and SCF robustness only -- unlike the Berry-phase
workflow, nothing here depends on maintaining a continuous gauge.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, asdict

import numpy as np
from pyscf import fci, mcscf

from .casscf import build_mol, run_rhf, transfer_mo, core_count, apply_singlet_constraint
from .config import CasConfig, ScanConfig
from .geometry import Loop, formaldimine_geom


@dataclass
class ScanResult:
    """A rectangular (alpha, phi) grid of state-averaged CASSCF energies."""

    alphas: np.ndarray
    phis: np.ndarray
    e_states: np.ndarray        # (n_alpha, n_phi, nroots)
    converged: np.ndarray       # (n_alpha, n_phi) bool
    cas_label: str
    basis: str
    loop_name: str
    weights: tuple[float, ...]
    wall_time: float = np.nan

    @property
    def gap(self) -> np.ndarray:
        """E1 - E0 in Hartree."""
        return self.e_states[:, :, 1] - self.e_states[:, :, 0]

    def min_gap_point(self) -> tuple[float, float, float]:
        """(alpha, phi, gap) at the grid minimum of the gap."""
        g = np.where(self.converged, self.gap, np.inf)
        i, j = np.unravel_index(np.argmin(g), g.shape)
        return float(self.alphas[i]), float(self.phis[j]), float(g[i, j])

    def resolution(self) -> tuple[float, float]:
        da = float(self.alphas[1] - self.alphas[0]) if self.alphas.size > 1 else np.nan
        dp = float(self.phis[1] - self.phis[0]) if self.phis.size > 1 else np.nan
        return da, dp

    def save(self, path: str) -> None:
        np.savez_compressed(
            path,
            alphas=self.alphas,
            phis=self.phis,
            e_states=self.e_states,
            converged=self.converged,
            cas_label=self.cas_label,
            basis=self.basis,
            loop_name=self.loop_name,
            weights=np.asarray(self.weights),
            wall_time=self.wall_time,
        )

    @classmethod
    def load(cls, path: str) -> "ScanResult":
        d = np.load(path, allow_pickle=False)
        return cls(
            alphas=d["alphas"],
            phis=d["phis"],
            e_states=d["e_states"],
            converged=d["converged"],
            cas_label=str(d["cas_label"]),
            basis=str(d["basis"]),
            loop_name=str(d["loop_name"]),
            weights=tuple(d["weights"].tolist()),
            wall_time=float(d["wall_time"]),
        )

    def summary(self) -> str:
        a, p, g = self.min_gap_point()
        da, dp = self.resolution()
        nconv = int(self.converged.sum())
        return (
            f"SA-CASSCF scan  {self.cas_label}/{self.basis}  loop {self.loop_name}\n"
            f"  grid: {self.alphas.size} x {self.phis.size} "
            f"(d_alpha={da:.3f} deg, d_phi={dp:.3f} deg), "
            f"alpha in [{self.alphas.min():.1f}, {self.alphas.max():.1f}], "
            f"phi in [{self.phis.min():.1f}, {self.phis.max():.1f}]\n"
            f"  converged: {nconv}/{self.converged.size}\n"
            f"  minimum gap {g:.6f} Ha at (alpha, phi) = ({a:.3f}, {p:.3f})"
        )


def run_sa_casscf(
    mol,
    ncas: int,
    nelecas,
    weights,
    mo_guess=None,
    conv_tol: float = 1e-9,
    max_cycle_macro: int = 200,
    mf=None,
):
    """Equal-weight (or specified-weight) state-averaged CASSCF. Returns (e_states, mo, converged)."""
    if mf is None:
        mf = run_rhf(mol)
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    # Spin-adapt *before* wrapping in the state average: the default solver would otherwise
    # return the Ms=0 triplet as a root and the "gap" would be S0/T1, not S0/S1.
    apply_singlet_constraint(mc)
    mc = mc.state_average_(list(weights))
    mc.conv_tol = conv_tol
    mc.max_cycle_macro = max_cycle_macro
    mc.verbose = 0
    mo = mf.mo_coeff if mo_guess is None else np.asarray(mo_guess)
    mc.kernel(mo)
    return np.asarray(mc.e_states, dtype=float), np.asarray(mc.mo_coeff), bool(mc.converged)


def scan_gap(
    loop: Loop,
    cas: CasConfig | None = None,
    scan: ScanConfig | None = None,
    geom_fn=formaldimine_geom,
    progress=None,
    checkpoint: str | None = None,
) -> ScanResult:
    """Two-dimensional SA-CASSCF gap scan over the loop's bounding box.

    If ``checkpoint`` names an existing ``.npz`` written by a previous run with the same grid,
    completed points are reused and only the missing ones are computed.
    """
    cas = cas or CasConfig()
    scan = scan or ScanConfig()
    say = progress or (lambda _m: None)
    t0 = time.time()

    a_lo, a_hi, p_lo, p_hi = loop.bounding_box(margin=scan.margin)
    alphas = np.linspace(a_lo, a_hi, scan.n_alpha)
    phis = np.linspace(p_lo, p_hi, scan.n_phi)

    e_states = np.full((scan.n_alpha, scan.n_phi, scan.nroots), np.nan)
    converged = np.zeros((scan.n_alpha, scan.n_phi), dtype=bool)
    done = np.zeros((scan.n_alpha, scan.n_phi), dtype=bool)

    if checkpoint:
        try:
            prev = ScanResult.load(checkpoint)
            if (
                np.allclose(prev.alphas, alphas)
                and np.allclose(prev.phis, phis)
                and prev.e_states.shape == e_states.shape
            ):
                mask = np.isfinite(prev.e_states).all(axis=2)
                e_states[mask] = prev.e_states[mask]
                converged[mask] = prev.converged[mask]
                done[mask] = True
                say(f"  restored {int(done.sum())} completed grid points from {checkpoint}")
        except (FileNotFoundError, OSError):
            pass

    # For strategy="anchor": solve once, cold, at the centre of the region. Every grid point
    # then gets the same two path-independent guesses (cold, and this anchor transferred in),
    # so no result can depend on the order in which the grid was walked.
    anchor_mol = anchor_mo = None
    if getattr(scan, "strategy", None) == "anchor":
        a_mid = float(alphas[len(alphas) // 2])
        p_mid = float(phis[len(phis) // 2])
        anchor_mol = build_mol(geom_fn(a_mid, p_mid), cas.basis,
                               charge=cas.charge, spin=cas.spin)
        try:
            _, anchor_mo, _ = run_sa_casscf(
                anchor_mol, cas.ncas, cas.nelecas, scan.weights,
                conv_tol=scan.conv_tol, max_cycle_macro=scan.max_cycle_macro,
            )
            say(f"  anchor solved cold at ({a_mid:.3f}, {p_mid:.3f})")
        except Exception as exc:                          # noqa: BLE001
            say(f"  anchor solve failed ({exc}); falling back to cold-only")
            anchor_mol = anchor_mo = None

    prev_mol = None
    prev_mo = None
    for i, alpha in enumerate(alphas):
        # serpentine: reverse every other row so consecutive solves stay neighbours
        col_order = range(scan.n_phi) if i % 2 == 0 else reversed(range(scan.n_phi))
        for j in col_order:
            if done[i, j]:
                continue
            phi = phis[j]
            mol = build_mol(geom_fn(alpha, phi), cas.basis, charge=cas.charge, spin=cas.spin)
            strategy = getattr(scan, "strategy", None) or (
                "warm" if scan.warm_start else "cold"
            )
            attempts = []
            if strategy in ("warm", "best") and prev_mo is not None:
                attempts.append(transfer_mo(prev_mol, prev_mo, mol))
            if strategy == "anchor" and anchor_mo is not None:
                attempts.append(transfer_mo(anchor_mol, anchor_mo, mol))
            if strategy in ("cold", "best", "anchor") or not attempts:
                attempts.append(None)

            best = None
            for guess in attempts:
                try:
                    e, mo, conv = run_sa_casscf(
                        mol, cas.ncas, cas.nelecas, scan.weights,
                        mo_guess=guess,
                        conv_tol=scan.conv_tol,
                        max_cycle_macro=scan.max_cycle_macro,
                    )
                except Exception as exc:                  # noqa: BLE001 - record and continue
                    say(f"    grid point ({alpha:.2f},{phi:.2f}) failed: {exc}")
                    continue
                # Prefer the lower state-averaged energy: that is the better stationary point.
                sa = float(np.dot(scan.weights, e[: len(scan.weights)]))
                if best is None or sa < best[0] - 1e-10:
                    best = (sa, e, mo, conv)

            if best is None:
                converged[i, j] = False
            else:
                _, e, mo, conv = best
                e_states[i, j, :] = e[: scan.nroots]
                converged[i, j] = conv
                prev_mol, prev_mo = mol, mo
            done[i, j] = True
        say(
            f"  row {i + 1}/{scan.n_alpha}  alpha={alpha:7.3f}  "
            f"min gap in row = {np.nanmin(e_states[i, :, 1] - e_states[i, :, 0]):.6f} Ha"
        )
        if checkpoint:
            ScanResult(
                alphas, phis, e_states, converged, cas.cas_label, cas.basis,
                loop.name, tuple(scan.weights), time.time() - t0,
            ).save(checkpoint)

    return ScanResult(
        alphas=alphas,
        phis=phis,
        e_states=e_states,
        converged=converged,
        cas_label=cas.cas_label,
        basis=cas.basis,
        loop_name=loop.name,
        weights=tuple(scan.weights),
        wall_time=time.time() - t0,
    )


def run_fci(mol, nroots: int = 2, mf=None):
    """Spin-adapted FCI in the full orbital space. Returns (e_states, converged).

    FCI is invariant to orbital rotations, so there is nothing to optimize and no active-space
    choice to justify: within a given basis this is the exact answer. Used as the reference
    against which the truncated active spaces are calibrated, and it is the same quantity the
    gap map of arXiv:2304.06070 Fig. 1a reports.
    """
    if mf is None:
        mf = run_rhf(mol)
    solver = fci.FCI(mf)
    solver.nroots = nroots
    solver.spin = 0
    fci.addons.fix_spin_(solver, ss=0)
    e, _ = solver.kernel()
    e = np.atleast_1d(e)
    return np.asarray(e[:nroots], dtype=float), True


def scan_gap_fci(
    loop: Loop,
    cas: CasConfig | None = None,
    scan: ScanConfig | None = None,
    geom_fn=formaldimine_geom,
    progress=None,
    checkpoint: str | None = None,
) -> ScanResult:
    """FCI gap scan over the loop's bounding box.

    Same interface and storage format as :func:`scan_gap`, so the two are directly comparable.
    Only ``cas.basis`` is read from ``cas``; the active space is the whole orbital space.
    Expensive: budget roughly 30 s per point for formaldimine/STO-3G.
    """
    cas = cas or CasConfig()
    scan = scan or ScanConfig()
    say = progress or (lambda _m: None)
    t0 = time.time()

    a_lo, a_hi, p_lo, p_hi = loop.bounding_box(margin=scan.margin)
    alphas = np.linspace(a_lo, a_hi, scan.n_alpha)
    phis = np.linspace(p_lo, p_hi, scan.n_phi)

    e_states = np.full((scan.n_alpha, scan.n_phi, scan.nroots), np.nan)
    converged = np.zeros((scan.n_alpha, scan.n_phi), dtype=bool)
    done = np.zeros((scan.n_alpha, scan.n_phi), dtype=bool)

    if checkpoint:
        try:
            prev = ScanResult.load(checkpoint)
            if (np.allclose(prev.alphas, alphas) and np.allclose(prev.phis, phis)
                    and prev.e_states.shape == e_states.shape):
                mask = np.isfinite(prev.e_states).all(axis=2)
                e_states[mask] = prev.e_states[mask]
                converged[mask] = prev.converged[mask]
                done[mask] = True
                say(f"  restored {int(done.sum())} completed FCI points from {checkpoint}")
        except (FileNotFoundError, OSError):
            pass

    for i, alpha in enumerate(alphas):
        for j, phi in enumerate(phis):
            if done[i, j]:
                continue
            mol = build_mol(geom_fn(alpha, phi), cas.basis, charge=cas.charge, spin=cas.spin)
            try:
                e, conv = run_fci(mol, nroots=scan.nroots)
                e_states[i, j, :] = e
                converged[i, j] = conv
            except Exception as exc:                       # noqa: BLE001
                say(f"    FCI at ({alpha:.2f},{phi:.2f}) failed: {exc}")
            done[i, j] = True
            say(f"    FCI ({alpha:7.3f},{phi:7.3f})  gap = "
                f"{e_states[i, j, 1] - e_states[i, j, 0]:.6f} Ha")
        if checkpoint:
            ScanResult(alphas, phis, e_states, converged, "FCI", cas.basis,
                       loop.name, tuple(scan.weights), time.time() - t0).save(checkpoint)

    return ScanResult(
        alphas=alphas, phis=phis, e_states=e_states, converged=converged,
        cas_label="FCI", basis=cas.basis, loop_name=loop.name,
        weights=tuple(scan.weights), wall_time=time.time() - t0,
    )
