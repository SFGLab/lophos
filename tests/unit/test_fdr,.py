from lophos.utils.fdr import bh_fdr

def test_bh_fdr_monotone():
    q = bh_fdr([0.001, 0.01, 0.2, 0.5, 0.9])
    assert all(0.0 <= x <= 1.0 for x in q)
