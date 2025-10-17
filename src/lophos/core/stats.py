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
    m = df_counts["maternal"].astype(int)
    p = df_counts["paternal"].astype(int)
    total = m + p
    ratio = [_log2_ratio(int(mi), int(pi)) for mi, pi in zip(m, p, strict=False)]
    pvals = [
        (
            binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue
            if (mi + pi) > 0
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
    return out


def compute_loop_stats(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Compute statistics for loop counts.

    The input ``df_counts`` should contain columns ``maternal_pairs``,
    ``paternal_pairs`` and ``ambiguous_pairs``.  The returned DataFrame
    includes renamed columns ``m`` and ``p`` for convenience, total counts
    (``total_pairs``), log2 ratios, binomial p-values, FDR-corrected q-values
    (``fdr_pairs``) and an ``ambiguous_frac`` column representing the
    fraction of ambiguous pairs.
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
    # Ambiguous fraction based on total pairs plus ambiguous counts
    denom = df["m"] + df["p"] + df["ambiguous_pairs"]
    # Avoid division by zero by replacing zero denominators with one
    denom_safe = denom.where(denom > 0, other=1)
    df["ambiguous_frac"] = df["ambiguous_pairs"] / denom_safe
    # Log2 ratio using only informative pairs; if both m and p are zero, ratio is 0
    df["log2_ratio_pairs"] = [
        _log2_ratio(int(mi), int(pi)) for mi, pi in zip(df["m"], df["p"], strict=False)
    ]
    # Binomial p-value (two-sided, p=0.5) for informative pairs
    pvals = []
    for mi, pi in zip(df["m"], df["p"], strict=False):
        total = int(mi) + int(pi)
        if total > 0:
            pvals.append(binomtest(int(mi), total, 0.5, alternative="two-sided").pvalue)
        else:
            pvals.append(1.0)
    df["p_value_pairs"] = pvals
    df["fdr_pairs"] = bh_fdr(df["p_value_pairs"])
    return df
