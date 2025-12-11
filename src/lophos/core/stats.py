import math
from dataclasses import dataclass

import pandas as pd
from scipy.stats import binomtest

from ..constants import PSEUDOCOUNT
from ..utils.fdr import bh_fdr


@dataclass(frozen=True)
class PeakStat:
    m: int
    p: int
    total: int
    log2_ratio: float
    p_value: float


def _log2_ratio(m: int, p: int) -> float:
    return math.log2((m + PSEUDOCOUNT) / (p + PSEUDOCOUNT))


def compute_peak_stats(df_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-peak stats:
      - total = m + p
      - log2_ratio = log2((m+pc)/(p+pc))
      - p_value (two-sided binomial at 0.5); when total==0 => p=1
      - fdr via BH over the *entire* p-value vector (rows with p=1 naturally yield q≈1)
    """
    m = df_counts["maternal"].astype(int)
    p = df_counts["paternal"].astype(int)
    total = m + p

    ratio = [_log2_ratio(int(mi), int(pi)) for mi, pi in zip(m, p, strict=False)]
    pvals = [
        (
            binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue
            if (int(mi) + int(pi)) > 0
            else 1.0
        )
        for mi, pi in zip(m, p, strict=False)
    ]
    qvals = bh_fdr(pvals)

    out = df_counts.copy()
    out["total"] = total
    out["log2_ratio"] = ratio
    out["p_value"] = pvals
    out["fdr"] = qvals

    # light sanity (helps catch any future regressions)
    assert (out["p_value"].between(0.0, 1.0)).all()
    assert (out["fdr"].between(0.0, 1.0)).all()
    # if p==1 on a row, q should be very close to 1 (exactly 1.0 in most BH implementations)
    # we allow tiny FP wiggle but never < 0.99
    assert (out.loc[out["p_value"] == 1.0, "fdr"] >= 0.99).all()

    return out


def compute_loop_stats(df_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute statistics for loop counts.

    Input df must contain:
      - maternal_pairs, paternal_pairs, ambiguous_pairs (ambiguous optional -> filled with zeros)

    Output includes:
      - m, p
      - total_pairs = m + p      (informative only)
      - ambiguous_frac = ambiguous / (m+p+ambiguous)      [safe when denom==0]
      - log2_ratio_pairs = log2((m+pc)/(p+pc))
      - p_value_pairs (two-sided binomial at 0.5 on informative pairs; total_pairs==0 -> p=1)
      - fdr_pairs via BH over p_value_pairs
    """
    df = df_counts.rename(columns={"maternal_pairs": "m", "paternal_pairs": "p"}).copy()

    # Ensure integer types
    df["m"] = df["m"].astype(int)
    df["p"] = df["p"].astype(int)

    # ambiguous_pairs column may not exist in some downstream calls; fill with zeros if absent
    if "ambiguous_pairs" not in df.columns:
        df["ambiguous_pairs"] = 0
    df["ambiguous_pairs"] = df["ambiguous_pairs"].astype(int)

    # Total informative pairs (excluding ambiguous) used for statistical testing
    df["total_pairs"] = df["m"] + df["p"]

    # Ambiguous fraction based on total informative + ambiguous counts
    denom = df["m"] + df["p"] + df["ambiguous_pairs"]
    denom_safe = denom.where(denom > 0, other=1)
    df["ambiguous_frac"] = df["ambiguous_pairs"] / denom_safe

    # Log2 ratio using only informative pairs
    df["log2_ratio_pairs"] = [
        _log2_ratio(int(mi), int(pi)) for mi, pi in zip(df["m"], df["p"], strict=False)
    ]

    # Binomial p-value (two-sided, p=0.5) for informative pairs
    pvals = []
    for mi, pi in zip(df["m"], df["p"], strict=False):
        tot = int(mi) + int(pi)
        pvals.append(
            binomtest(int(mi), tot, 0.5, alternative="two-sided").pvalue if tot > 0 else 1.0
        )
    df["p_value_pairs"] = pvals
    df["fdr_pairs"] = bh_fdr(df["p_value_pairs"])

    # light sanity
    assert (df["p_value_pairs"].between(0.0, 1.0)).all()
    assert (df["fdr_pairs"].between(0.0, 1.0)).all()
    return df
