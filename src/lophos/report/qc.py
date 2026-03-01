from __future__ import annotations

from pathlib import Path

import pandas as pd


def _vc(df: pd.DataFrame, col: str, value: str) -> int:
    """Value-count helper that survives missing columns."""
    if col not in df.columns:
        return 0
    return int((df[col].astype(str) == value).sum())


def write_summary(path: Path, peaks: pd.DataFrame, loops: pd.DataFrame) -> None:
    """Write a lightweight QC summary TSV.

    Stable outputs:
      - bias_call breakdown for peaks and loops
      - evidence_tier breakdown for loops (if available)
      - depth medians (if available)

    This file is meant for quick sanity checks, not deep analysis.
    """
    s: list[tuple[str, object]] = []

    # Peaks
    s.append(("peaks_total", int(len(peaks))))
    s.append(("peaks_maternal", _vc(peaks, "bias_call", "Maternal")))
    s.append(("peaks_paternal", _vc(peaks, "bias_call", "Paternal")))
    s.append(("peaks_balanced", _vc(peaks, "bias_call", "Balanced")))
    s.append(("peaks_undetermined", _vc(peaks, "bias_call", "Undetermined")))
    if "total" in peaks.columns:
        s.append(
            ("peaks_total_median", float(pd.to_numeric(peaks["total"], errors="coerce").median()))
        )

    # Loops
    s.append(("loops_total", int(len(loops))))
    s.append(("loops_maternal", _vc(loops, "bias_call", "Maternal")))
    s.append(("loops_paternal", _vc(loops, "bias_call", "Paternal")))
    s.append(("loops_balanced", _vc(loops, "bias_call", "Balanced")))
    s.append(("loops_undetermined", _vc(loops, "bias_call", "Undetermined")))

    # Evidence tiers (v1.0+)
    if "evidence_tier" in loops.columns:
        s.append(("loops_evidence_direct", _vc(loops, "evidence_tier", "direct")))
        s.append(("loops_evidence_inferred", _vc(loops, "evidence_tier", "inferred")))
        s.append(("loops_evidence_insufficient", _vc(loops, "evidence_tier", "insufficient")))

    # Depth medians
    if "total_pairs" in loops.columns:
        s.append(
            (
                "loops_total_pairs_median",
                float(pd.to_numeric(loops["total_pairs"], errors="coerce").median()),
            )
        )
    elif "informative_pairs" in loops.columns:
        s.append(
            (
                "loops_informative_pairs_median",
                float(pd.to_numeric(loops["informative_pairs"], errors="coerce").median()),
            )
        )

    pd.DataFrame(s, columns=["metric", "value"]).to_csv(path, sep="\t", index=False)
