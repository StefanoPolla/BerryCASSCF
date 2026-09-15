"""2x2 linear Jahn-Teller (E x e) model: an analytically known Berry phase.

The point of this module is to exercise the gauge-fixing convention and the two Berry-phase
estimators on a system whose answer is known in closed form, before any quantum chemistry
cost can obscure a sign bug.

The model Hamiltonian in a diabatic basis is

    H(x, y) = kappa * [[x, y], [y, -x]]   (+ a multiple of the identity, which is irrelevant)

whose lower adiabatic state is ``(-sin(theta/2), cos(theta/2))`` with ``theta = atan2(y, x)``.
Transporting it once around the origin sends ``theta -> theta + 2*pi`` and hence flips its sign:
the Berry phase is ``pi`` for any loop enclosing the origin and ``0`` for any loop that does not.

Because the diabatic basis is the same at every point, the overlap between neighbouring states is
just a dot product -- there is no nonorthogonal-CI machinery here. What *is* shared with the CASSCF
workflow is the sequential gauge fixing and the two estimators, which is what we want to test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def jt_hamiltonian(x: float, y: float, kappa: float = 1.0) -> np.ndarray:
    return kappa * np.array([[x, y], [y, -x]], dtype=float)


def jt_ground_state(x: float, y: float, kappa: float = 1.0) -> np.ndarray:
    """Lower adiabatic eigenvector from a black-box diagonalizer (arbitrary sign)."""
    _, vecs = np.linalg.eigh(jt_hamiltonian(x, y, kappa))
    return vecs[:, 0]


def jt_gap(x: float, y: float, kappa: float = 1.0) -> float:
    return float(2.0 * kappa * np.hypot(x, y))


@dataclass
class ToyResult:
    centre: tuple[float, float]
    radius: float
    n_points: int
    product_estimator: float
    endpoint_estimator: float
    min_abs_adjacent_overlap: float
    encloses_origin: bool
    expected_nontrivial: bool

    @property
    def is_nontrivial(self) -> bool:
        return self.product_estimator < 0

    @property
    def agrees_with_theory(self) -> bool:
        return (
            self.is_nontrivial == self.expected_nontrivial
            and np.sign(self.product_estimator) == np.sign(self.endpoint_estimator)
        )


def jt_loop_berry_phase(
    centre: tuple[float, float] = (0.0, 0.0),
    radius: float = 1.0,
    n_points: int = 24,
    kappa: float = 1.0,
    randomize_signs: bool = True,
    seed: int = 0,
) -> ToyResult:
    """Transport the lower adiabatic state around a circle and extract the Z2 Berry phase.

    ``randomize_signs`` multiplies each freshly diagonalized eigenvector by a random +-1,
    imitating the arbitrary overall sign a CASSCF solver returns. A correct implementation
    must be completely insensitive to it.
    """
    rng = np.random.default_rng(seed)
    ts = np.arange(n_points) / n_points
    pts = [
        (centre[0] + radius * np.cos(2 * np.pi * t), centre[1] + radius * np.sin(2 * np.pi * t))
        for t in ts
    ]

    states: list[np.ndarray] = []
    adjacent: list[float] = []
    for k, (x, y) in enumerate(pts):
        v = jt_ground_state(x, y, kappa)
        if randomize_signs:
            v = v * rng.choice([-1.0, 1.0])
        if k > 0:
            raw = float(states[-1] @ v)
            if raw < 0:                     # same sequential gauge fixing as the CASSCF walk
                v = -v
            adjacent.append(abs(raw))
        states.append(v)

    closing = float(states[-1] @ states[0])
    product = float(np.prod(adjacent) * closing)

    # Endpoint estimator: one more transport step onto the identical point.
    v_end = jt_ground_state(*pts[0], kappa)
    if randomize_signs:
        v_end = v_end * rng.choice([-1.0, 1.0])
    if states[-1] @ v_end < 0:
        v_end = -v_end
    endpoint = float(states[0] @ v_end)

    encloses = bool(np.hypot(*centre) < radius)
    return ToyResult(
        centre=centre,
        radius=radius,
        n_points=n_points,
        product_estimator=product,
        endpoint_estimator=endpoint,
        min_abs_adjacent_overlap=float(min(adjacent)) if adjacent else float("nan"),
        encloses_origin=encloses,
        expected_nontrivial=encloses,
    )
