from dataclasses import dataclass

import numpy as np
import pandas as pd

@dataclass(frozen=True)
class APAResult:
    matrix: np.ndarray  # type: ignore[type-arg]
    score: float
    window_bp: int
    bin_bp: int

def compute_apa_matrix(
    _loops_df: pd.DataFrame,
    window_bp: int = 100_000,
    bin_bp: int = 5_000,
) -> APAResult:
    """Placeholder APA."""
    size = max(5, (window_bp // bin_bp) | 1)
    mat = np.zeros((size, size), dtype=float)
    return APAResult(matrix=mat, score=0.0, window_bp=window_bp, bin_bp=bin_bp)
