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
        binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue
        if (mi + pi) > 0 else 1.0
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
    df = df_counts.rename(columns={"maternal_pairs": "m", "paternal_pairs": "p"}).copy()
    df["total_pairs"] = df["m"].astype(int) + df["p"].astype(int)
    df["log2_ratio_pairs"] = [_log2_ratio(int(mi), int(pi)) for mi, pi in zip(df["m"], df["p"], strict=False)]
    df["p_value_pairs"] = [
        binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue if (mi + pi) > 0 else 1.0
        for mi, pi in zip(df["m"], df["p"], strict=False)
    ]
    df["fdr_pairs"] = bh_fdr(df["p_value_pairs"])
    return df
