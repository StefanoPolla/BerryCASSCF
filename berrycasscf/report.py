"""Collect saved Berry-phase records into comparison tables.

Used by the results notebook and by ``examples/summarize.py``.
"""

from __future__ import annotations

import glob
import os
import re

from .store import load_json

_NAME_RE = re.compile(r"(?P<loop>.+)_cas(?P<ne>\d+)-(?P<ncas>\d+)_N(?P<n>\d+)\.json$")


def load_berry_records(directory: str) -> list[dict]:
    """Load every Berry-phase JSON record in ``directory``, newest schema assumed."""
    records = []
    for path in sorted(glob.glob(os.path.join(directory, "*.json"))):
        m = _NAME_RE.search(os.path.basename(path))
        if not m:
            continue
        blob = load_json(path)
        rec = dict(blob["result"])
        rec["_path"] = path
        rec["_traversal"] = blob.get("traversal")
        rec["nelec"] = int(m.group("ne"))
        rec["ncas_from_name"] = int(m.group("ncas"))
        records.append(rec)
    return records


def berry_table(records: list[dict]) -> list[dict]:
    """Flatten records into rows suitable for a table."""
    rows = []
    for r in records:
        rows.append(
            {
                "loop": r["loop_name"],
                "CAS": r["cas_label"],
                "N": r["n_points"],
                "product": r["product_estimator"],
                "endpoint": r["endpoint_estimator"],
                "min|ovl|": r["min_abs_adjacent_overlap"],
                "phase": r["berry_phase"],
                "status": r["status"],
                "t(s)": r.get("wall_time"),
                "warnings": len(r.get("messages") or []),
            }
        )
    rows.sort(key=lambda d: (d["loop"], d["CAS"], d["N"]))
    return rows


def format_table(rows: list[dict], columns: list[str] | None = None) -> str:
    """Plain-text table, so results are readable without pandas."""
    if not rows:
        return "(no records)"
    columns = columns or list(rows[0])

    def cell(v):
        if isinstance(v, float):
            return f"{v:+.4f}" if abs(v) < 1000 else f"{v:.3g}"
        return str(v)

    widths = {c: max(len(c), *(len(cell(r.get(c))) for r in rows)) for c in columns}
    head = "  ".join(c.ljust(widths[c]) for c in columns)
    sep = "  ".join("-" * widths[c] for c in columns)
    body = "\n".join("  ".join(cell(r.get(c)).ljust(widths[c]) for c in columns) for r in rows)
    return f"{head}\n{sep}\n{body}"


def stability_verdict(rows: list[dict], loop: str, cas: str) -> dict:
    """Is the conclusion for one (loop, CAS) stable across discretizations?

    Stability criteria, applied consistently throughout this project:

    1. every traversal at that (loop, CAS) has ``status == "OK"``;
    2. at least two different discretizations ``N`` were run;
    3. all of them return the *same* Berry phase;
    4. the product and endpoint estimators agree in sign in every run (already enforced
       by the OK status, re-checked here).
    """
    subset = [r for r in rows if r["loop"] == loop and r["CAS"] == cas]
    phases = {r["phase"] for r in subset}
    ns = sorted({r["N"] for r in subset})
    all_ok = bool(subset) and all(r["status"] == "OK" for r in subset)
    signs_agree = all(
        (r["product"] < 0) == (r["endpoint"] < 0)
        for r in subset
        if r["endpoint"] is not None
    )
    stable = all_ok and len(ns) >= 2 and len(phases) == 1 and signs_agree
    return {
        "loop": loop,
        "CAS": cas,
        "n_runs": len(subset),
        "N_values": ns,
        "phases": sorted(phases),
        "all_ok": all_ok,
        "estimators_agree": signs_agree,
        "stable": stable,
    }
