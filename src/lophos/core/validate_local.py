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
    """
    Placeholder validation that passes data through and adds neutral columns.

    Parameters
    ----------
    _bam : Any
        BAM handle (unused in placeholder).
    _loops_df : pd.DataFrame
        Loops dataframe (unused in placeholder).
    loop_calls : pd.DataFrame
        Loop stats/calls to be annotated.
    anchor_pad : int
        Anchor padding used during counting (unused in placeholder).
    mapq : int
        MAPQ threshold used during counting (unused in placeholder).
    """
    _ = (anchor_pad, mapq, _bam, _loops_df)  # silence linters until implemented

    out = loop_calls.copy()
    out["local_enrichment_z"] = 0.0
    out["local_enrichment_p"] = 1.0
    return out
