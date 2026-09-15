"""Exact overlap between two CASSCF wavefunctions defined at *different* geometries.

Why this module exists
----------------------
At two neighbouring points of a nuclear loop the wavefunctions are expanded in different
AO bases (the nuclei moved) and different MO bases (the orbitals were re-optimized).
Their overlap is therefore a genuine *nonorthogonal* CI overlap.

PySCF's ``fci.addons.transform_ci_for_orbital_rotation`` / ``fci.addons.overlap`` assume the
two orbital sets span the same space and are related by a unitary. That is false here, so
they are **not** used anywhere in this package.

The method
----------
Let ``|A>`` and ``|B>`` be CASSCF wavefunctions with a common ``ncore`` (doubly occupied) and
``ncas`` (active) orbital count. Write the cross-geometry MO overlap over the occupied-capable
space as blocks ``S_cc``, ``S_cA``, ``S_Ac``, ``S_AA``.

For a determinant pair the occupied-orbital overlap matrix is

    D_IJ = [[S_cc,      S_cA[:, J]],
            [S_Ac[I, :], S_AA[I, J]]]

whose determinant, by the Schur complement identity, is

    det(D_IJ) = det(S_cc) * det( T[I, J] ),    T = S_AA - S_Ac @ inv(S_cc) @ S_cA

so only *minors of the single ncas x ncas matrix* ``T`` depend on the determinant pair.

PySCF orders a determinant as (alpha orbitals ascending)(beta orbitals ascending), so the two
spin blocks factorize with no extra sign, and with ``M_s[I, J] = det(T[occ_I, occ_J])``

    <A|B> = det(S_cc)**2 * sum_{Ia,Ib} c_A[Ia,Ib] * (M_alpha @ c_B @ M_beta.T)[Ia,Ib]

which is a Frobenius inner product. Core orbitals are doubly occupied, hence ``det(S_cc)`` squared.

Everything here is exact: no orthogonality between the two orbital sets is assumed at any point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from pyscf import gto
from pyscf.fci import cistring


@dataclass
class CasWavefunction:
    """A CASSCF solution at one geometry, with everything needed to form overlaps."""

    mol: gto.Mole
    mo_coeff: np.ndarray           # (nao, nmo), AO basis
    ci: np.ndarray                 # (n_str_alpha, n_str_beta)
    ncas: int
    nelecas: tuple[int, int]
    ncore: int
    energy: float = np.nan
    converged: bool = True
    label: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def nocc(self) -> int:
        """Number of orbitals that can be occupied: core + active."""
        return self.ncore + self.ncas

    def occupied_block(self) -> np.ndarray:
        """MO coefficients of the core+active orbitals."""
        return np.asarray(self.mo_coeff)[:, : self.nocc]

    def scaled(self, factor: float) -> "CasWavefunction":
        """Copy with the CI vector multiplied by ``factor`` (used for sign/gauge fixing)."""
        import copy

        new = copy.copy(self)
        new.ci = np.asarray(self.ci) * factor
        return new


@lru_cache(maxsize=None)
def _occ_lists(ncas: int, nelec: int) -> np.ndarray:
    """Ascending occupied-orbital indices for every CI string address. Shape (n_str, nelec)."""
    return np.asarray(cistring.gen_occslst(range(ncas), nelec), dtype=np.int64)


def _minor_determinants(T: np.ndarray, occ: np.ndarray) -> np.ndarray:
    """``M[I, J] = det(T[occ[I], occ[J]])`` for every string pair, vectorized.

    ``occ`` has shape (n_str, k); the result has shape (n_str, n_str).
    For k == 0 (no electrons of this spin) every minor is the empty determinant, 1.
    """
    n_str, k = occ.shape
    if k == 0:
        return np.ones((n_str, n_str))
    # sub[I, J, x, y] = T[occ[I, x], occ[J, y]]
    sub = T[occ[:, None, :, None], occ[None, :, None, :]]
    return np.linalg.det(sub)


def cross_mo_overlap(wfn_a: CasWavefunction, wfn_b: CasWavefunction) -> np.ndarray:
    """MO-basis overlap ``C_a^T S_AO(Ra, Rb) C_b`` over the core+active block.

    ``S_AO`` is the cross-geometry AO overlap and is *not* symmetric.
    """
    s_ao = gto.intor_cross("int1e_ovlp", wfn_a.mol, wfn_b.mol)
    return wfn_a.occupied_block().T @ s_ao @ wfn_b.occupied_block()


def cas_overlap(
    wfn_a: CasWavefunction,
    wfn_b: CasWavefunction,
    return_diagnostics: bool = False,
):
    """Exact nonorthogonal overlap ``<A|B>`` between two CASSCF wavefunctions.

    Both wavefunctions must share ``ncas``, ``nelecas`` and ``ncore``; the geometries,
    AO bases and MO coefficients may differ arbitrarily.

    Returns the scalar overlap, or ``(overlap, diagnostics)`` if ``return_diagnostics``.
    """
    if (wfn_a.ncas, wfn_a.nelecas, wfn_a.ncore) != (wfn_b.ncas, wfn_b.nelecas, wfn_b.ncore):
        raise ValueError(
            "CAS definitions differ: "
            f"{(wfn_a.ncas, wfn_a.nelecas, wfn_a.ncore)} vs "
            f"{(wfn_b.ncas, wfn_b.nelecas, wfn_b.ncore)}"
        )

    ncore, ncas = wfn_a.ncore, wfn_a.ncas
    s_mo = cross_mo_overlap(wfn_a, wfn_b)

    # --- core elimination via the Schur complement --------------------------------
    if ncore > 0:
        s_cc = s_mo[:ncore, :ncore]
        s_ca = s_mo[:ncore, ncore:]
        s_ac = s_mo[ncore:, :ncore]
        s_aa = s_mo[ncore:, ncore:]
        sign_core, logdet_core = np.linalg.slogdet(s_cc)
        if sign_core == 0.0:
            raise np.linalg.LinAlgError(
                "Core-block overlap is singular: the two core spaces are orthogonal "
                "somewhere along the path. The two points are not connected by continuation."
            )
        T = s_aa - s_ac @ np.linalg.solve(s_cc, s_ca)
        core_cond = float(np.linalg.cond(s_cc))
    else:
        sign_core, logdet_core = 1.0, 0.0
        T = s_mo
        core_cond = 1.0

    # Core is doubly occupied -> det(S_cc) enters once per spin.
    log_prefactor = 2.0 * logdet_core
    prefactor_sign = float(sign_core) ** 2  # always +1, kept explicit for clarity

    # --- active-space string minors ------------------------------------------------
    na, nb = wfn_a.nelecas
    occ_a = _occ_lists(ncas, na)
    occ_b = _occ_lists(ncas, nb)
    m_alpha = _minor_determinants(T, occ_a)
    m_beta = _minor_determinants(T, occ_b)

    ci_a = np.asarray(wfn_a.ci).reshape(len(occ_a), len(occ_b))
    ci_b = np.asarray(wfn_b.ci).reshape(len(occ_a), len(occ_b))

    active_part = float(np.sum(ci_a * (m_alpha @ ci_b @ m_beta.T)))
    value = prefactor_sign * np.exp(log_prefactor) * active_part

    if not return_diagnostics:
        return value

    diagnostics = {
        "core_logdet": float(logdet_core),
        "core_cond": core_cond,
        "active_part": active_part,
        "T_min_singular_value": float(np.linalg.svd(T, compute_uv=False).min()),
    }
    return value, diagnostics


def normalized_cas_overlap(wfn_a: CasWavefunction, wfn_b: CasWavefunction) -> float:
    """``<A|B> / sqrt(<A|A><B|B>)``.

    The self-overlaps are 1 for a properly normalized CASSCF solution, so this is a
    cheap guard against an unnormalized CI vector rather than a routine need.
    """
    nab = cas_overlap(wfn_a, wfn_b)
    naa = cas_overlap(wfn_a, wfn_a)
    nbb = cas_overlap(wfn_b, wfn_b)
    return float(nab / np.sqrt(naa * nbb))


def brute_force_cas_overlap(wfn_a: CasWavefunction, wfn_b: CasWavefunction) -> float:
    """Reference implementation: explicit sum over determinant pairs.

    Builds each determinant's occupied-orbital overlap matrix in full (core included)
    instead of eliminating the core analytically. Used only to validate
    :func:`cas_overlap`; cost is O(n_det^2) with no factorization.
    """
    ncore, ncas = wfn_a.ncore, wfn_a.ncas
    s_mo = cross_mo_overlap(wfn_a, wfn_b)
    core = np.arange(ncore)
    na, nb = wfn_a.nelecas
    occ_a = _occ_lists(ncas, na)
    occ_b = _occ_lists(ncas, nb)
    ci_a = np.asarray(wfn_a.ci).reshape(len(occ_a), len(occ_b))
    ci_b = np.asarray(wfn_b.ci).reshape(len(occ_a), len(occ_b))

    def det_block(rows_act, cols_act):
        rows = np.concatenate([core, ncore + rows_act]).astype(int)
        cols = np.concatenate([core, ncore + cols_act]).astype(int)
        return np.linalg.det(s_mo[np.ix_(rows, cols)])

    total = 0.0
    for ia, occ_ia in enumerate(occ_a):
        for ja, occ_ja in enumerate(occ_a):
            d_alpha = det_block(occ_ia, occ_ja)
            if d_alpha == 0.0:
                continue
            for ib, occ_ib in enumerate(occ_b):
                for jb, occ_jb in enumerate(occ_b):
                    if ci_a[ia, ib] == 0.0 or ci_b[ja, jb] == 0.0:
                        continue
                    d_beta = det_block(occ_ib, occ_jb)
                    total += ci_a[ia, ib] * ci_b[ja, jb] * d_alpha * d_beta
    return float(total)
