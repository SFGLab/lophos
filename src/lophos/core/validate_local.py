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

    This is lightweight QC only.

    v1.0+ behavior:
      - if `evidence_tier` is present, the background distribution is estimated
        from loops with `evidence_tier == 'direct'` ONLY (to avoid circularity
        when anchor-fallback inference is enabled).
      - z/p are computed for all loops using that reference distribution.

    Accepts either:
      - `m`, `p` (older/internal naming)
      - `maternal_pairs`, `paternal_pairs` (current normalized naming)
    """
    import numpy as np
    from scipy.stats import norm

    # Currently unused here; keep signature stable.
    _ = (anchor_pad, mapq, _bam, _loops_df)

    out = loop_calls.copy()

    # Select maternal/paternal columns
    if {"m", "p"}.issubset(out.columns):
        mcol, pcol = "m", "p"
    elif {"maternal_pairs", "paternal_pairs"}.issubset(out.columns):
        mcol, pcol = "maternal_pairs", "paternal_pairs"
    else:
        out["local_enrichment_z"] = 0.0
        out["local_enrichment_p"] = 1.0
        out["local_enrichment_ref"] = "none"
        return out

    # Choose background subset (prefer direct-evidence loops)
    ref = out
    ref_label = "all"
    if "evidence_tier" in out.columns:
        direct = out[out["evidence_tier"].astype(str) == "direct"]
        if len(direct) >= 10:
            ref = direct
            ref_label = "direct"
        else:
            ref = out
            ref_label = "all_fallback"

    ref_diff = (ref[mcol].astype(float) - ref[pcol].astype(float)).to_numpy()
    std = float(np.std(ref_diff, ddof=0))
    if std == 0.0:
        std = 1.0
    mean = float(np.mean(ref_diff))

    diff = (out[mcol].astype(float) - out[pcol].astype(float)).to_numpy()
    z = (diff - mean) / std
    pvals = 2.0 * norm.sf(np.abs(z))

    out["local_enrichment_z"] = z
    out["local_enrichment_p"] = pvals
    out["local_enrichment_ref"] = ref_label
    return out
