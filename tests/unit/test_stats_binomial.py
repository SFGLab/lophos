import pandas as pd
from lophos.core import stats


def test_peak_stats_shapes():
    df = pd.DataFrame({"peak_id": ["p1", "p2"], "maternal": [10, 0], "paternal": [0, 10]})
    out = stats.compute_peak_stats(df)
    assert {"total", "log2_ratio", "p_value", "fdr"}.issubset(out.columns)
    assert list(out["total"]) == [10, 10]
