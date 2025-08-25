import pandas as pd

def run_local_validation(bam, loops_df: pd.DataFrame, loop_calls: pd.DataFrame, anchor_pad: int, mapq: int) -> pd.DataFrame:
    # TODO: implement true local background & z-scores
    loop_calls = loop_calls.copy()
    loop_calls["local_enrichment_z"] = 0.0
    loop_calls["local_enrichment_p"] = 1.0
    return loop_calls
