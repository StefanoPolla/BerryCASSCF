"""Saving and reloading run records, so results can be read without recomputing.

Berry-phase runs are stored as JSON (small, diffable, human-readable); scan grids use
``.npz`` via :meth:`berrycasscf.scan.ScanResult.save`.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, is_dataclass
from typing import Any

import numpy as np


def _plain(obj: Any):
    """Recursively convert numpy scalars/arrays and dataclasses to JSON-safe types."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return _plain(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _plain(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def save_json(obj: Any, path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(_plain(obj), fh, indent=2, sort_keys=False)
    return path


def load_json(path: str) -> Any:
    with open(path) as fh:
        return json.load(fh)


def save_berry_run(result, traversal, path: str) -> str:
    """Store a Berry-phase result together with its full per-point diagnostics."""
    return save_json({"result": result, "traversal": traversal}, path)


def berry_record_exists(path: str) -> bool:
    return os.path.exists(path)
