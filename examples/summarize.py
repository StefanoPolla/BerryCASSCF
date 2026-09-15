#!/usr/bin/env python
"""Print a table of every saved Berry-phase result, plus the per-(loop, CAS) stability verdict."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from berrycasscf.report import load_berry_records, berry_table, format_table, stability_verdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main() -> int:
    rows = berry_table(load_berry_records(os.path.join(ROOT, "results", "berry")))
    if not rows:
        print("No results yet. Run examples/run_formaldimine_berry.py first.")
        return 1
    print(format_table(rows, ["loop", "CAS", "N", "product", "endpoint",
                              "min|ovl|", "phase", "status", "t(s)"]))
    print("\nStability verdicts (same phase across all N, every run OK, >= 2 discretizations):")
    seen = []
    for r in rows:
        key = (r["loop"], r["CAS"])
        if key in seen:
            continue
        seen.append(key)
        v = stability_verdict(rows, *key)
        mark = "STABLE  " if v["stable"] else "UNSTABLE"
        print(f"  {mark} {v['loop']:5s} {v['CAS']:10s} N={v['N_values']}  phases={v['phases']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
