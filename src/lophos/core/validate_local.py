import pandas as pd

def run_local_validation(
    _bam,
    _loops_df: pd.DataFrame,
    loop_calls: pd.DataFrame,
    anchor_pad: int,
    mapq: int,
) -> pd.DataFrame:
    _ = (anchor_pad, mapq)  # quiet linters until implemented
    out = loop_calls.copy()
    out["local_enrichment_z"] = 0.0
    out["local_enrichment_p"] = 1.0
    return out
