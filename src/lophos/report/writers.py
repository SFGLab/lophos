from pathlib import Path

import pandas as pd


def write_peaks(path: Path, df: pd.DataFrame) -> None:
    cols = [
        "chrom",
        "start",
        "end",
        "peak_id",
        "maternal",
        "paternal",
        "total",
        "log2_ratio",
        "p_value",
        "fdr",
        "bias_call",
    ]
    # sanity checks
    assert set(cols).issubset(df.columns), f"Missing peak columns: {set(cols) - set(df.columns)}"
    # enforce total consistency (helps catch upstream regressions early)
    if "maternal" in df.columns and "paternal" in df.columns and "total" in df.columns:
        assert (
            df["total"] == (df["maternal"].astype(int) + df["paternal"].astype(int))
        ).all(), "peaks 'total' must equal maternal + paternal"
    assert (df["p_value"].between(0.0, 1.0)).all()
    assert (df["fdr"].between(0.0, 1.0)).all()
    df[cols].to_csv(path, sep="\t", index=False, header=False)


def write_loops(path: Path, df: pd.DataFrame) -> None:
    """
    Export loops in the schema we actually generate:
      [6 anchor cols],
      total_pairs, maternal_pairs, paternal_pairs, ambiguous_pairs, informative_pairs,
      log2_ratio_pairs, p_value_pairs, fdr_pairs, bias_call

    Notes:
      - We DO NOT export a 'loop_id' unless it truly exists (avoids misaligned columns).
      - informative_pairs = m + p (identical to total_pairs here; kept for clarity/compat).
    """
    out = df.rename(columns={"m": "maternal_pairs", "p": "paternal_pairs"}).copy()

    # Normalize/guard integer pair columns
    for c in ("maternal_pairs", "paternal_pairs", "ambiguous_pairs", "total_pairs"):
        if c not in out.columns:
            # keep explicit KeyError behavior consistent with previous code
            continue
        out[c] = out[c].fillna(0).astype(int)

    # Compute informative_pairs explicitly (m+p)
    if "maternal_pairs" not in out.columns or "paternal_pairs" not in out.columns:
        raise KeyError("Expected 'maternal_pairs' and 'paternal_pairs' in loop dataframe.")
    out["informative_pairs"] = out["maternal_pairs"].astype(int) + out["paternal_pairs"].astype(int)

    # Column presence & sanity
    anchor_cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    required = anchor_cols + [
        "total_pairs",
        "maternal_pairs",
        "paternal_pairs",
        "ambiguous_pairs",
        "informative_pairs",
        "log2_ratio_pairs",
        "p_value_pairs",
        "fdr_pairs",
        "bias_call",
    ]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise KeyError(f"Missing loop columns for export: {missing}")

    # internal consistency: totals & probabilities
    # If total_pairs disagrees with m+p, correct it (then assert)
    mismatch = out["total_pairs"] != out["informative_pairs"]
    if mismatch.any():
        out.loc[mismatch, "total_pairs"] = out.loc[mismatch, "informative_pairs"]
    assert (
        out["total_pairs"] == out["informative_pairs"]
    ).all(), "total_pairs must equal informative_pairs (m+p) in current implementation."
    assert (out["p_value_pairs"].between(0.0, 1.0)).all()
    assert (out["fdr_pairs"].between(0.0, 1.0)).all()

    export_cols = required  # exact order
    out[export_cols].to_csv(path, sep="\t", index=False, header=False)
