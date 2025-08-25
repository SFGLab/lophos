from pathlib import Path

import pandas as pd

def write_peaks(path: Path, df: pd.DataFrame) -> None:
    cols = ["chrom","start","end","peak_id","maternal","paternal","total","log2_ratio","p_value","fdr","bias_call"]
    df[cols].to_csv(path, sep="\t", index=False, header=False)

def write_loops(path: Path, df: pd.DataFrame) -> None:
    out = df.rename(columns={"m": "maternal_pairs", "p": "paternal_pairs"})
    export_cols = [
        "chrom1","start1","end1","chrom2","start2","end2","loop_id",
        "maternal_pairs","paternal_pairs","ambiguous_pairs","total_pairs",
        "log2_ratio_pairs","p_value_pairs","fdr_pairs","bias_call"
    ]
    out[export_cols].to_csv(path, sep="\t", index=False, header=False)
