from dataclasses import dataclass
import pandas as pd

@dataclass(frozen=True)
class BiasThresholds:
    min_reads: int = 5
    fdr: float = 0.05
    min_fold: float = 1.5  # used as M >= P * min_fold

def _classify(m: int, p: int, q: float, thr: BiasThresholds) -> str:
    total = m + p
    if total < thr.min_reads:
        return "Undetermined"
    if q <= thr.fdr:
        if m >= max(1, int(p * thr.min_fold)):
            return "Maternal"
        if p >= max(1, int(m * thr.min_fold)):
            return "Paternal"
    return "Balanced"

def call_bias_for_peaks(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.DataFrame:
    df = stats_df.copy()
    df["bias_call"] = [
        _classify(int(m), int(p), float(q), thresholds)
        for m, p, q in zip(df["maternal"], df["paternal"], df["fdr"])
    ]
    return df

def call_bias_for_loops(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.DataFrame:
    df = stats_df.copy()
    df["bias_call"] = [
        _classify(int(m), int(p), float(q), thresholds)
        for m, p, q in zip(df["m"], df["p"], df["fdr_pairs"])
    ]
    return df
