from pathlib import Path
import pandas as pd

BEDPE_COLS = [
    "chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2"
]

def read_bedpe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    df.columns = BEDPE_COLS[: len(df.columns)]
    return df
