"""BerryCASSCF: a variational Berry-phase conical-intersection detector at CASSCF level.

Two independent workflows answer the same question for a closed loop in nuclear
coordinate space:

1. :func:`berrycasscf.berry.run_loop` -- track the state-specific CASSCF ground state
   around the loop by continuation and read the Z2 Berry phase off the sign of the
   initial-final overlap.
2. :func:`berrycasscf.scan.scan_gap` -- resolve the S0/S1 near-degeneracy inside the
   loop with a state-averaged CASSCF gap scan.

Systems: formaldimine (:mod:`berrycasscf.geometry`), ethylene (:mod:`berrycasscf.ethylene`)
and fulvene (:mod:`berrycasscf.fulvene`). Each exposes a ``<name>_geom(x, y)`` function with
the same signature, so the same drivers work on all of them via their ``geom_fn`` argument.
"""

from .geometry import (
    Loop,
    formaldimine_geom,
    LOOP_CI,
    LOOP_CONTROL_LOWER,
    LOOP_CONTROL_UPPER,
    BENCHMARK_LOOPS,
    REFERENCE_CI_ALPHA_PHI,
    OVERVIEW_REGION,
    SCAN_REGIONS,
)
from .config import CasConfig, ContinuationConfig, ScanConfig
from .overlap import CasWavefunction, cas_overlap, brute_force_cas_overlap
from .casscf import build_mol, run_casscf, transfer_mo
from .continuation import traverse_loop, LoopTraversal
from .berry import analyse, run_loop, BerryResult, TRIVIAL, NONTRIVIAL, UNDETERMINED
from .scan import scan_gap, ScanResult
from .toy import jt_loop_berry_phase
from . import ethylene, fulvene

__version__ = "0.1.0"

__all__ = [
    "Loop", "formaldimine_geom", "LOOP_CI", "LOOP_CONTROL_LOWER", "LOOP_CONTROL_UPPER",
    "BENCHMARK_LOOPS", "REFERENCE_CI_ALPHA_PHI", "OVERVIEW_REGION", "SCAN_REGIONS",
    "CasConfig", "ContinuationConfig", "ScanConfig",
    "CasWavefunction", "cas_overlap", "brute_force_cas_overlap",
    "build_mol", "run_casscf", "transfer_mo",
    "traverse_loop", "LoopTraversal",
    "analyse", "run_loop", "BerryResult", "TRIVIAL", "NONTRIVIAL", "UNDETERMINED",
    "scan_gap", "ScanResult",
    "jt_loop_berry_phase",
    "ethylene", "fulvene",
]
