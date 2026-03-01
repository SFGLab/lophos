from __future__ import annotations

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
    assert (
        df["total"] == (df["maternal"].astype(int) + df["paternal"].astype(int))
    ).all(), "peaks 'total' must equal maternal + paternal"
    assert (df["p_value"].between(0.0, 1.0)).all()
    assert (df["fdr"].between(0.0, 1.0)).all()
    df[cols].to_csv(path, sep="\t", index=False, header=False)


def _normalize_informative_counts(df: pd.DataFrame) -> None:
    """Normalize maternal and paternal pair counts in-place."""
    if "m" in df.columns and "p" in df.columns:
        df.rename(columns={"m": "maternal_pairs", "p": "paternal_pairs"}, inplace=True)
    if "maternal_pairs" not in df.columns or "paternal_pairs" not in df.columns:
        raise KeyError("Expected loop columns 'maternal_pairs'/'paternal_pairs' or 'm'/'p'.")

    df["maternal_pairs"] = df["maternal_pairs"].fillna(0).astype(int)
    df["paternal_pairs"] = df["paternal_pairs"].fillna(0).astype(int)


def _normalize_noninformative_counts(df: pd.DataFrame) -> None:
    """Normalize ambiguous and homozygous pair counts in-place."""
    if "ambiguous_pairs_raw" in df.columns:
        df["ambiguous_pairs"] = df["ambiguous_pairs_raw"].fillna(0).astype(int)
    elif "ambiguous_pairs" in df.columns:
        df["ambiguous_pairs"] = df["ambiguous_pairs"].fillna(0).astype(int)
    else:
        df["ambiguous_pairs"] = 0

    if "homozygous_pairs_raw" in df.columns:
        df["homozygous_pairs"] = df["homozygous_pairs_raw"].fillna(0).astype(int)
    elif "homozygous_pairs" in df.columns:
        df["homozygous_pairs"] = df["homozygous_pairs"].fillna(0).astype(int)
    else:
        df["homozygous_pairs"] = 0


def _compute_total_and_noninformative(df: pd.DataFrame) -> None:
    """Compute total and noninformative pair counts in-place."""
    df["informative_pairs"] = df["maternal_pairs"] + df["paternal_pairs"]

    if "total_pairs" in df.columns:
        df["total_pairs"] = df["total_pairs"].fillna(0).astype(int)
    else:
        df["total_pairs"] = df["informative_pairs"].astype(int)

    if "noninformative_pairs" in df.columns:
        df["noninformative_pairs"] = df["noninformative_pairs"].fillna(0).astype(int)
    else:
        df["noninformative_pairs"] = (df["ambiguous_pairs"] + df["homozygous_pairs"]).astype(int)


def _normalize_ambiguous_frac(df: pd.DataFrame) -> None:
    """Normalize ambiguous fraction in-place."""
    if "ambiguous_frac" in df.columns:
        df["ambiguous_frac"] = df["ambiguous_frac"].astype(float)
    else:
        denom = df["informative_pairs"] + df["noninformative_pairs"]
        denom_safe = denom.where(denom > 0, other=1)
        df["ambiguous_frac"] = df["noninformative_pairs"] / denom_safe


def _normalize_evidence_columns(df: pd.DataFrame) -> None:
    """Normalize evidence-related columns in-place."""
    if "bias_call_direct" not in df.columns:
        df["bias_call_direct"] = df["bias_call"]
    if "bias_call_inferred" not in df.columns:
        df["bias_call_inferred"] = "Undetermined"
    if "bias_call_final" not in df.columns:
        df["bias_call_final"] = df["bias_call"]
    if "evidence_tier" not in df.columns:
        df["evidence_tier"] = "unknown"
    if "inferred_reason" not in df.columns:
        df["inferred_reason"] = "na"

    # keep bias_call aligned to final call for stable downstream behavior
    df["bias_call"] = df["bias_call_final"]


def write_loops(path: Path, df: pd.DataFrame) -> None:
    """Write loop calls as BEDPE-like TSV with a frozen, parseable schema.

    This schema supports:
      - informative (maternal/paternal) evidence for direct loop phasing
      - non-informative evidence (ambiguous + homozygous ties)
      - stats computed on informative counts only (log2_ratio_pairs, p_value_pairs, fdr_pairs)
      - explicit evidence provenance via:
          evidence_tier ∈ {direct, inferred, insufficient}
          bias_call_direct / bias_call_inferred / bias_call_final

    Notes
    -----
    - The core phasing algorithm should populate the evidence columns.
      We still fill defaults defensively to avoid breaking I/O on partial DataFrames.
    - `bias_call` is kept as an alias for the final call for backwards compatibility.
    """
    out = df.copy()

    # --- anchors ---
    anchor_cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    missing_anchors = [c for c in anchor_cols if c not in out.columns]
    if missing_anchors:
        raise KeyError(f"Missing loop anchor columns for export: {missing_anchors}")

    # --- normalize informative counts ---
    _normalize_informative_counts(out)

    # --- normalize non-informative counts ---
    _normalize_noninformative_counts(out)

    # --- compute totals ---
    _compute_total_and_noninformative(out)

    # --- ambiguous fraction ---
    _normalize_ambiguous_frac(out)

    # --- required stats columns ---
    required_stats = ["log2_ratio_pairs", "p_value_pairs", "fdr_pairs"]
    missing_stats = [c for c in required_stats if c not in out.columns]
    if missing_stats:
        raise KeyError(f"Missing loop stats columns for export: {missing_stats}")

    if "bias_call" not in out.columns:
        raise KeyError("Missing loop column 'bias_call' for export.")

    # --- evidence columns (v1.0+) ---
    # Fill defensively if upstream didn't populate them
    _normalize_evidence_columns(out)

    assert (out["p_value_pairs"].between(0.0, 1.0)).all()
    assert (out["fdr_pairs"].between(0.0, 1.0)).all()

    export_cols = anchor_cols + [
        "total_pairs",
        "maternal_pairs",
        "paternal_pairs",
        "ambiguous_pairs",
        "homozygous_pairs",
        "informative_pairs",
        "noninformative_pairs",
        "ambiguous_frac",
        "log2_ratio_pairs",
        "p_value_pairs",
        "fdr_pairs",
        # legacy/final label first, then provenance
        "bias_call",
        "evidence_tier",
        "bias_call_direct",
        "bias_call_inferred",
        "bias_call_final",
        "inferred_reason",
    ]

    missing = [c for c in export_cols if c not in out.columns]
    if missing:
        raise KeyError(f"Missing loop columns for export after normalization: {missing}")

    out[export_cols].to_csv(path, sep="\t", index=False, header=False)
