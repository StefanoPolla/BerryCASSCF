"""Direct search for the smallest S1-S0 gap, to check the cone fit rather than trust it.

The fitted intersection positions elsewhere in this package come from a model: a parabola in
``gap**2``, exact for an ideal cone. That model can be wrong in ways the residual does not always
catch -- a solution discontinuity looks like a steep branch, and a cut that is not conical over
the fitted window biases the apex. This module answers the same question without a model, by
minimizing the gap directly.

Nelder-Mead is used deliberately. Near a conical intersection the gap behaves like
``|displacement|``: continuous but not differentiable at the minimum, so gradient-based methods
have nothing useful to work with. A simplex method only ever compares function values.

Every evaluation is an independent **cold** SA-CASSCF solve, matching the scan default, so the
search is path-independent and reproducible: the same starting point always gives the same
trajectory. Results are cached, since the simplex revisits points.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .casscf import apply_singlet_constraint, build_mol, run_rhf
from .config import CasConfig


@dataclass
class SearchResult:
    """Outcome of a direct gap minimization."""

    x: float
    y: float
    gap: float                       # mHa
    n_evaluations: int
    converged: bool
    message: str
    start: tuple[float, float]
    start_gap: float
    trace: list[tuple[float, float, float]] = field(default_factory=list)

    @property
    def moved(self) -> float:
        """How far the search travelled from its starting point, in the plane's units."""
        return float(np.hypot(self.x - self.start[0], self.y - self.start[1]))

    def summary(self) -> str:
        return (f"start ({self.start[0]:.3f}, {self.start[1]:.3f}) gap {self.start_gap:.4f} mHa"
                f"  ->  ({self.x:.3f}, {self.y:.3f}) gap {self.gap:.4f} mHa"
                f"   [moved {self.moved:.3f}, {self.n_evaluations} evaluations]")


def gap_at(x: float, y: float, geom_fn, cas: CasConfig, weights=(0.5, 0.5)) -> float:
    """Equal-weight state-averaged S1-S0 gap at one point, in mHa. Cold start."""
    from pyscf import mcscf

    mol = build_mol(geom_fn(x, y), cas.basis, charge=cas.charge, spin=cas.spin)
    mc = mcscf.CASSCF(run_rhf(mol), cas.ncas, cas.nelecas)
    apply_singlet_constraint(mc)
    mc = mc.state_average_(list(weights))
    mc.conv_tol = cas.conv_tol
    mc.max_cycle_macro = cas.max_cycle_macro
    mc.verbose = 0
    mc.kernel()
    return float(mc.e_states[1] - mc.e_states[0]) * 1000.0


def minimize_gap(
    geom_fn,
    cas: CasConfig,
    start: tuple[float, float],
    step: tuple[float, float] = (1.0, 1.0),
    max_evaluations: int = 40,
    tol: float = 1e-3,
    progress=None,
) -> SearchResult:
    """Minimize the gap from ``start`` with a bounded number of solves.

    ``step`` sets the initial simplex size in each coordinate; it should be comparable to the
    uncertainty in the starting guess, not to the grid spacing. ``max_evaluations`` is a hard
    budget: each evaluation is a full SA-CASSCF, so the search is meant to *refine* a good guess
    over a few tens of solves, not to find an intersection from nothing.
    """
    say = progress or (lambda _m: None)
    cache: dict[tuple[float, float], float] = {}
    trace: list[tuple[float, float, float]] = []

    def f(v):
        key = (round(float(v[0]), 6), round(float(v[1]), 6))
        if key in cache:
            return cache[key]
        try:
            g = gap_at(key[0], key[1], geom_fn, cas)
        except Exception as exc:                       # noqa: BLE001 - keep the search alive
            say(f"    evaluation at ({key[0]:.3f}, {key[1]:.3f}) failed: {exc}")
            g = float("inf")
        cache[key] = g
        trace.append((key[0], key[1], g))
        say(f"({key[0]:8.3f}, {key[1]:8.3f})  gap = {g:9.4f} mHa")
        return g

    start_gap = f(np.array(start, dtype=float))
    simplex = np.array([
        [start[0], start[1]],
        [start[0] + step[0], start[1]],
        [start[0], start[1] + step[1]],
    ], dtype=float)

    res = minimize(
        f, np.array(start, dtype=float), method="Nelder-Mead",
        options={"initial_simplex": simplex, "maxfev": max_evaluations,
                 "xatol": tol, "fatol": tol, "disp": False},
    )
    return SearchResult(
        x=float(res.x[0]), y=float(res.x[1]), gap=float(res.fun),
        n_evaluations=len(cache), converged=bool(res.success), message=str(res.message),
        start=(float(start[0]), float(start[1])), start_gap=float(start_gap), trace=trace,
    )
