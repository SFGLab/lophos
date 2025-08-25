from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
import pandas as pd

@dataclass(frozen=True)
class APAResult:
    matrix: np.ndarray          # aggregate matrix (rows x cols)
    score: float                # simple APA-like score (peak - background)
    window_bp: int
    bin_bp: int

def compute_apa_matrix(
    loops_df: pd.DataFrame,
    window_bp: int = 100_000,
    bin_bp: int = 5_000,
) -> APAResult:
    """
    Placeholder APA: returns a small centered matrix of zeros and score 0.0.
    A real implementation would bin contacts around each loop and average.
    """
    size = max(5, (window_bp // bin_bp) | 1)  # ensure odd
    mat = np.zeros((size, size), dtype=float)
    score = 0.0
    return APAResult(matrix=mat, score=score, window_bp=window_bp, bin_bp=bin_bp)
