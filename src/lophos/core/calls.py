from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class BiasThresholds:
    """Thresholds controlling allele bias calling.

    Parameters
    ----------
    min_reads : int
        Minimum informative read/pair count (m+p) required to attempt a bias call.
        Below this threshold, the call is set to ``"Undetermined"``.
    fdr : float
        Maximum FDR value (Benjamini–Hochberg corrected p-value) required to
        consider a feature significant.  Non-significant features are called
        ``"Balanced"``.
    min_fold : float
        Minimum fold-change (m vs p) to consider a feature biased.  Must be
        >= 1.0.  A value of 1.0 means any deviation from equality is
        sufficient when ``fdr`` passes.
    min_abs_log2 : float
        Minimum absolute log2 ratio |log2((m+1)/(p+1))| required to call a
        feature biased.  Features with smaller effect sizes are called
        ``"Balanced"``.  Default is 0.0 (no effect-size threshold).
    max_ambiguous_frac : float
        Maximum allowable fraction of ambiguous pairs in a loop before it is
        automatically called ``"Undetermined"``.  Applicable only to loops.
        Defaults to 1.0 (no effect).
    """

    min_reads: int = 5
    fdr: float = 0.05
    min_fold: float = 1.5
    min_abs_log2: float = 0.0
    max_ambiguous_frac: float = 1.0


def _classify(
    m: int,
    p: int,
    q: float,
    thr: BiasThresholds,
    log2_ratio: float | None = None,
    ambiguous_frac: float | None = None,
) -> str:
    """Classify a feature into Maternal/Paternal/Balanced/Undetermined.

    The decision is based on a combination of read counts, FDR threshold,
    minimum fold change, effect-size threshold and ambiguous fraction.

    Parameters
    ----------
    m : int
        Count of maternal reads/pairs.
    p : int
        Count of paternal reads/pairs.
    q : float
        Adjusted p-value (FDR) for the feature.
    thr : BiasThresholds
        Threshold parameters controlling calling behaviour.
    log2_ratio : float | None, optional
        Pre-computed log2 ratio (m/p).  If not provided, it will be
        calculated internally from m and p using the pseudocount defined in
        ``lophos.constants.PSEUDOCOUNT``.  This parameter enables the caller
        to avoid recomputing the ratio.
    ambiguous_frac : float | None, optional
        Fraction of ambiguous pairs for loops.  When provided, loops with
        ``ambiguous_frac > thr.max_ambiguous_frac`` are automatically
        called ``"Undetermined"``.
    """
    total = m + p
    # Low coverage -> Undetermined
    if total < thr.min_reads:
        return "Undetermined"
    # Ambiguous fraction guard for loops
    if ambiguous_frac is not None and ambiguous_frac > thr.max_ambiguous_frac:
        return "Undetermined"
    # Non-significant p-value -> Balanced
    if q > thr.fdr:
        return "Balanced"
    # Compute log2 ratio if not provided
    if log2_ratio is None:
        import math

        from ..constants import PSEUDOCOUNT

        log2_ratio = math.log2((m + PSEUDOCOUNT) / (p + PSEUDOCOUNT))
    # Effect-size threshold
    if abs(float(log2_ratio)) < thr.min_abs_log2:
        return "Balanced"
    # Fold-change threshold
    # m and p always >=0 here.  Use max to ensure at least 1 read on the opposing side
    if m >= max(1, int(p * thr.min_fold)):
        return "Maternal"
    if p >= max(1, int(m * thr.min_fold)):
        return "Paternal"
    # If none of the directional thresholds are met, classify as balanced
    return "Balanced"


def call_bias_for_peaks(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.DataFrame:
    """Apply bias classification to peaks.

    This function reads the maternal/paternal counts, FDR values and log2 ratios
    from ``stats_df`` and assigns a bias call for each peak.  See
    ``BiasThresholds`` for controlling parameters.
    """
    df = stats_df.copy()
    calls: list[str] = []
    log2_ratios = df["log2_ratio"] if "log2_ratio" in df.columns else pd.Series([None] * len(df))
    for m, p, q, r in zip(df["maternal"], df["paternal"], df["fdr"], log2_ratios, strict=False):
        calls.append(
            _classify(
                int(m),
                int(p),
                float(q),
                thresholds,
                log2_ratio=float(r) if r is not None else None,
                ambiguous_frac=None,
            )
        )
    df["bias_call"] = calls
    return df


def call_bias_for_loops(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.DataFrame:
    """Apply bias classification to loops.

    This function reads the informative pair counts (m/p), FDR values,
    log2 ratios and ambiguous fractions from ``stats_df`` and assigns a
    bias call.  See ``BiasThresholds`` for controlling parameters.
    """
    df = stats_df.copy()
    calls: list[str] = []
    log2_ratios = (
        df["log2_ratio_pairs"] if "log2_ratio_pairs" in df.columns else pd.Series([None] * len(df))
    )
    ambiguous_fracs = (
        df["ambiguous_frac"] if "ambiguous_frac" in df.columns else pd.Series([None] * len(df))
    )
    for m, p, q, r, amb in zip(
        df["m"],
        df["p"],
        df["fdr_pairs"],
        log2_ratios,
        ambiguous_fracs,
        strict=False,
    ):
        calls.append(
            _classify(
                int(m),
                int(p),
                float(q),
                thresholds,
                log2_ratio=float(r) if r is not None else None,
                ambiguous_frac=float(amb) if amb is not None else None,
            )
        )
    df["bias_call"] = calls
    return df
