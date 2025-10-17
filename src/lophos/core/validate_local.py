from __future__ import annotations

from typing import Any

import pandas as pd


def run_local_validation(
    _bam: Any,
    _loops_df: pd.DataFrame,
    loop_calls: pd.DataFrame,
    anchor_pad: int,
    mapq: int,
) -> pd.DataFrame:
    """Approximate local validation for loop bias calls.

    The goal of local validation is to place a loop's maternal/paternal
    imbalance into context relative to other loops.  A simple z-score is
    computed from the global distribution of the difference ``(m - p)`` and
    used to derive a two-sided p-value under a standard normal model.  This
    serves as a proxy for local enrichment when a true local background model
    is unavailable.  Future improvements may bin loops by genomic distance
    and sample matched controls.

    Parameters
    ----------
    _bam : Any
        BAM handle (unused in this approximation).
    _loops_df : pd.DataFrame
        Original loops dataframe (unused in this approximation).
    loop_calls : pd.DataFrame
        Dataframe of loop statistics and bias calls.
    anchor_pad : int
        Anchor padding used during counting (unused here).
    mapq : int
        MAPQ threshold used during counting (unused here).

    Returns
    -------
    pd.DataFrame
        A copy of ``loop_calls`` with two additional columns:
        ``local_enrichment_z`` (float) and ``local_enrichment_p`` (float).
    """
    import numpy as np
    from scipy.stats import norm

    # We ignore anchor_pad, mapq, _bam and _loops_df in this approximation
    _ = (anchor_pad, mapq, _bam, _loops_df)

    out = loop_calls.copy()
    # Ensure required columns exist
    if not {"m", "p"}.issubset(out.columns):
        # Nothing to compute
        out["local_enrichment_z"] = 0.0
        out["local_enrichment_p"] = 1.0
        return out
    # Compute the signed difference between maternal and paternal counts
    diff = (out["m"].astype(float) - out["p"].astype(float)).to_numpy()
    # Use population standard deviation; if std is zero, set to 1 to avoid division by zero
    std = float(np.std(diff, ddof=0))
    if std == 0.0:
        std = 1.0
    mean = float(np.mean(diff))
    # Compute z-scores and two-sided p-values under the normal distribution
    z = (diff - mean) / std
    pvals = 2.0 * norm.sf(np.abs(z))
    out["local_enrichment_z"] = z
    out["local_enrichment_p"] = pvals
    return out
