from pathlib import Path

import pandas as pd

BEDPE_COLS = [
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "name",
    "score",
    "strand1",
    "strand2",
]


def read_bedpe(path: Path) -> pd.DataFrame:
    # Read all columns; skip comments/track lines
    df = pd.read_csv(path, sep="\t", comment="#", header=None)

    # Must have at least the 6 coordinate columns
    if df.shape[1] < 6:
        raise ValueError(f"BEDPE file {path} has fewer than 6 columns")

    base = BEDPE_COLS
    n = df.shape[1]
    if n <= len(base):
        # Fewer or equal to 10 columns: name what we have
        df.columns = base[:n]
    else:
        # More than 10 columns: keep extras and name them extra1, extra2, ...
        extras = [f"extra{i}" for i in range(1, n - len(base) + 1)]
        df.columns = base + extras

    return df
