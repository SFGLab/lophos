from pathlib import Path
import pandas as pd

BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]

def read_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", header=None)
    df.columns = BED_COLS[: len(df.columns)]
    return df

def to_interval_series(df: pd.DataFrame):
    return df[["chrom", "start", "end"]].itertuples(index=False, name=None)
