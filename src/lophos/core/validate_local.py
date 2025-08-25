import pandas as pd

def run_local_validation(_bam, _loops_df: pd.DataFrame, loop_calls: pd.DataFrame, _anchor_pad: int, _mapq: int) -> pd.DataFrame:
    """Placeholder: copy-through with neutral columns."""
    out = loop_calls.copy()
    out["local_enrichment_z"] = 0.0
    out["local_enrichment_p"] = 1.0
    return out
