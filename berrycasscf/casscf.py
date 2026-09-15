"""Thin PySCF wrappers: build molecules, solve CASSCF, move orbitals between geometries.

Nothing here is clever; the intent is that every scientific choice is visible.
"""

from __future__ import annotations

import numpy as np
from pyscf import gto, scf, mcscf, lib

from .overlap import CasWavefunction

# Linear-dependence threshold for the S^{+-1/2} construction.
_OAO_EIGENVALUE_FLOOR = 1e-10


def build_mol(geom: str, basis: str, charge: int = 0, spin: int = 0, verbose: int = 0,
              max_memory: int = 4000) -> gto.Mole:
    """Build a Mole in C1 (no point-group symmetry), as the loop breaks symmetry."""
    mol = gto.Mole()
    mol.atom = geom
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.symmetry = False          # C1: the loop breaks symmetry
    mol.unit = "Angstrom"
    mol.verbose = verbose
    mol.max_memory = max_memory
    mol.build()
    return mol


def core_count(mol: gto.Mole, nelecas: int | tuple[int, int]) -> int:
    """Number of doubly occupied core orbitals for a CAS(nelecas, ncas) calculation."""
    n_act = nelecas if isinstance(nelecas, int) else sum(nelecas)
    n_core_elec = mol.nelectron - n_act
    if n_core_elec < 0 or n_core_elec % 2:
        raise ValueError(
            f"Cannot form a closed-shell core: {mol.nelectron} electrons, {n_act} active."
        )
    return n_core_elec // 2


def as_nelecas_pair(nelecas: int | tuple[int, int]) -> tuple[int, int]:
    if isinstance(nelecas, (tuple, list)):
        return (int(nelecas[0]), int(nelecas[1]))
    n = int(nelecas)
    return (n - n // 2, n // 2)


# --- orthonormal-atomic-orbital (OAO) frame ---------------------------------------
# The continuity trick taken from the auto_oo tutorial: MO coefficients are carried
# between geometries in the OAO frame, in which they stay orthonormal.
#
#   C_OAO = S^{1/2} C_AO        (at the source geometry)
#   C_AO  = S^{-1/2} C_OAO      (at the target geometry)
#
# Re-using C_AO directly would give non-orthonormal MOs at the new geometry.

def _overlap_powers(mol: gto.Mole) -> tuple[np.ndarray, np.ndarray]:
    """Return (S^{1/2}, S^{-1/2}) for the AO overlap of ``mol``."""
    s = mol.intor_symmetric("int1e_ovlp")
    w, v = np.linalg.eigh(s)
    if w.min() < _OAO_EIGENVALUE_FLOOR:
        raise np.linalg.LinAlgError(
            f"AO overlap is near-singular (min eigenvalue {w.min():.3e}); "
            "the basis is linearly dependent at this geometry."
        )
    s_half = (v * np.sqrt(w)) @ v.T
    s_inv_half = (v / np.sqrt(w)) @ v.T
    return s_half, s_inv_half


def mo_to_oao(mol: gto.Mole, mo_coeff: np.ndarray) -> np.ndarray:
    """AO-basis MO coefficients -> OAO-basis MO coefficients."""
    s_half, _ = _overlap_powers(mol)
    return s_half @ mo_coeff


def oao_to_mo(mol: gto.Mole, oao_mo_coeff: np.ndarray) -> np.ndarray:
    """OAO-basis MO coefficients -> AO-basis MO coefficients at ``mol``'s geometry."""
    _, s_inv_half = _overlap_powers(mol)
    return s_inv_half @ oao_mo_coeff


def transfer_mo(mol_from: gto.Mole, mo_from: np.ndarray, mol_to: gto.Mole) -> np.ndarray:
    """Carry MOs from one geometry to another through the OAO frame.

    The result is orthonormal at ``mol_to`` by construction, and is the identity map
    when the two geometries coincide.
    """
    return oao_to_mo(mol_to, mo_to_oao(mol_from, mo_from))


def orthonormality_error(mol: gto.Mole, mo_coeff: np.ndarray) -> float:
    """max |C^T S C - I|, a check that a transferred guess is a valid MO set."""
    s = mol.intor_symmetric("int1e_ovlp")
    g = mo_coeff.T @ s @ mo_coeff
    return float(np.abs(g - np.eye(g.shape[0])).max())


# --- SCF and CASSCF ---------------------------------------------------------------

def apply_singlet_constraint(mc) -> None:
    """Restrict a CAS solver to singlets.

    Without this, PySCF's default ``direct_spin1`` solver returns the Ms=0 *triplet* as a
    root. For formaldimine that triplet is below the first excited singlet over much of the
    (alpha, phi) plane, so an unconstrained state average silently resolves the wrong pair of
    states. Every CAS solve in this package is spin-adapted for that reason.
    """
    mc.fcisolver.spin = 0
    mc.fix_spin_(ss=0)


def run_rhf(mol: gto.Mole, dm0: np.ndarray | None = None, conv_tol: float = 1e-10,
            max_cycle: int = 100) -> scf.hf.RHF:
    """RHF, optionally warm-started from a density matrix.

    Its role is mainly to supply integrals to CASSCF; the CASSCF orbital guess comes
    from the continuation, not from here.
    """
    mf = scf.RHF(mol)
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    mf.kernel(dm0=dm0)
    return mf


def run_casscf(
    mol: gto.Mole,
    ncas: int,
    nelecas: int | tuple[int, int],
    mo_guess: np.ndarray | None = None,
    ci0: np.ndarray | None = None,
    conv_tol: float = 1e-10,
    conv_tol_grad: float | None = None,
    max_cycle_macro: int = 100,
    fix_spin: bool = True,
    label: str = "",
    mf: scf.hf.RHF | None = None,
) -> CasWavefunction:
    """State-specific CASSCF for the lowest root, warm-startable in both MOs and CI.

    ``mo_guess`` should already be orthonormal at ``mol`` (use :func:`transfer_mo`).
    """
    if mf is None:
        mf = run_rhf(mol)
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.conv_tol = conv_tol
    if conv_tol_grad is not None:
        mc.conv_tol_grad = conv_tol_grad
    mc.max_cycle_macro = max_cycle_macro
    mc.verbose = 0
    if fix_spin:
        # Keep the solver on the singlet surface; root flipping to a triplet would be a
        # silent state change, which is exactly what the continuation must avoid.
        apply_singlet_constraint(mc)

    mo = mf.mo_coeff if mo_guess is None else np.asarray(mo_guess)
    mc.kernel(mo, ci0=ci0)

    ncore = core_count(mol, nelecas)
    return CasWavefunction(
        mol=mol,
        mo_coeff=np.asarray(mc.mo_coeff),
        ci=np.asarray(mc.ci),
        ncas=ncas,
        nelecas=as_nelecas_pair(mc.nelecas),
        ncore=ncore,
        energy=float(mc.e_tot),
        converged=bool(mc.converged),
        label=label,
        meta={"e_cas": float(mc.e_cas)},
    )


def casci_roots(wfn: CasWavefunction, nroots: int = 2) -> np.ndarray:
    """Spin-adapted CASCI energies at a wavefunction's own converged orbitals.

    This is a **root-flipping risk indicator for the continuation**, not an estimate of the
    true vertical S0/S1 gap. The orbitals are optimized for the ground state alone and the
    excited root is confined to the tracked active space, so for a small CAS it can be
    wildly too large: in CAS(2,2) the second root is the doubly-excited configuration,
    roughly 1 Ha up, while the physical S1 is a single excitation that the active space does
    not even contain. Its job here is to warn when the *tracked* solution is close to
    another root of its own CI problem, which is when continuation can silently jump branch.

    The physically meaningful S0/S1 gap comes from the state-averaged scan in
    :mod:`berrycasscf.scan`.
    """
    mf = scf.RHF(wfn.mol)
    mf.mo_coeff = wfn.mo_coeff
    mf.mo_occ = np.zeros(wfn.mo_coeff.shape[1])
    mf.mo_occ[: wfn.ncore] = 2
    mc = mcscf.CASCI(mf, wfn.ncas, wfn.nelecas)
    apply_singlet_constraint(mc)
    mc.fcisolver.nroots = nroots
    mc.verbose = 0
    mc.kernel(wfn.mo_coeff)
    e = np.atleast_1d(mc.e_tot)
    return np.asarray(e, dtype=float)
