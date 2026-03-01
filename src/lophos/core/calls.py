from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class BiasThresholds:
    """Thresholds controlling allele bias calling.

    Parameters
    ----------
    min_reads : int
        Minimum *informative* read/pair count (m+p) required to attempt a bias call.
        Below this threshold, the call is set to ``"Undetermined"``.
    fdr : float
        Maximum FDR value (Benjamini–Hochberg corrected p-value) required to
        consider a feature significant. Non-significant features are called
        ``"Balanced"``.
    min_fold : float
        Minimum fold-change (m vs p) to consider a feature biased. Must be
        >= 1.0. A value of 1.0 means any deviation from equality is sufficient
        when ``fdr`` passes.
    min_abs_log2 : float
        Minimum absolute log2 ratio |log2((m+pc)/(p+pc))| required to call a
        feature biased. Features with smaller effect sizes are called
        ``"Balanced"``. Default is 0.0 (no effect-size threshold).
    max_ambiguous_frac : float
        Maximum allowable fraction of *non-informative* pairs in a loop before it is
        automatically called ``"Undetermined"``. In v1.0.0 this fraction can include
        both ambiguous and homozygous-like ties depending on the policy used in
        `compute_loop_stats`. Applicable only to loops. Defaults to 1.0 (no effect).
    """

    min_reads: int = 5  # informative reads/pairs
    fdr: float = 0.05
    min_fold: float = 1.5  # practical effect guard
    min_abs_log2: float = 0.0  # additional effect guard (optional)
    max_ambiguous_frac: float = 1.0  # loops-only guard


def _classify(
    m: int,
    p: int,
    q: float,
    thr: BiasThresholds,
    log2_ratio: float | None = None,
    ambiguous_frac: float | None = None,
) -> str:
    """Classify a feature into Maternal/Paternal/Balanced/Undetermined."""
    total = m + p

    # Coverage gate
    if total < thr.min_reads:
        return "Undetermined"

    # Non-informative guard (loops)
    if ambiguous_frac is not None and ambiguous_frac > thr.max_ambiguous_frac:
        return "Undetermined"

    # Significance gate
    if q > thr.fdr:
        return "Balanced"

    # Compute/consume effect size
    if log2_ratio is None:
        import math

        from ..constants import PSEUDOCOUNT

        log2_ratio = math.log2((m + PSEUDOCOUNT) / (p + PSEUDOCOUNT))

    # Effect-size (log2) guard
    if abs(float(log2_ratio)) < thr.min_abs_log2:
        return "Balanced"

    # Practical fold guard using the SAME pseudocount as log2
    from ..constants import PSEUDOCOUNT

    fold = (m + PSEUDOCOUNT) / (p + PSEUDOCOUNT)
    if fold >= thr.min_fold:
        return "Maternal"
    if (1.0 / fold) >= thr.min_fold:
        return "Paternal"

    return "Balanced"


def call_bias_for_peaks(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.DataFrame:
    """Apply bias classification to peaks."""
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

    Expects columns:
      - m, p (informative counts post policy)
      - fdr_pairs
      - log2_ratio_pairs (optional)
      - ambiguous_frac (optional; recommended)

    ambiguous_frac may include homozygous-like ties depending on the policy used
    in compute_loop_stats().
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
        df["m"], df["p"], df["fdr_pairs"], log2_ratios, ambiguous_fracs, strict=False
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


def evidence_tier_direct_for_loops(stats_df: pd.DataFrame, thresholds: BiasThresholds) -> pd.Series:
    """Return a direct-evidence tier label for loops.

    This is NOT the inferred-anchor fallback. It only indicates whether a loop has
    sufficient direct (connectivity) evidence under the thresholds.

    - 'direct'       : total_pairs >= min_reads AND ambiguous_frac <= max_ambiguous_frac
    - 'insufficient' : otherwise

    CLI can later upgrade 'insufficient' to 'inferred' if anchor fallback is enabled.
    """
    if "total_pairs" in stats_df.columns:
        total_pairs = stats_df["total_pairs"].astype(int)
    else:
        total_pairs = (stats_df["m"] + stats_df["p"]).astype(int)

    amb_frac = (
        stats_df["ambiguous_frac"].astype(float)
        if "ambiguous_frac" in stats_df.columns
        else pd.Series([0.0] * len(stats_df), index=stats_df.index)
    )

    ok = (total_pairs >= int(thresholds.min_reads)) & (
        amb_frac <= float(thresholds.max_ambiguous_frac)
    )
    return pd.Series(["direct" if v else "insufficient" for v in ok], index=stats_df.index)
