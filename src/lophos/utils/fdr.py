from collections.abc import Iterable


def bh_fdr(pvals: Iterable[float]) -> list[float]:
    """
    Benjamini–Hochberg step-up procedure (two key properties):
      - q-values are in [0, 1]
      - Monotone (non-decreasing) with respect to sorted p-values

    Implementation detail:
      1) Sort p ascending, compute q_raw = p * n / rank
      2) Enforce monotonicity by applying a cumulative minimum from the *end* (reverse)
      3) Map q back to original order
    """
    p = list(pvals)
    n = len(p)
    indexed = sorted(enumerate(p), key=lambda t: t[1])

    # q_raw in sorted order
    q_raw = [0.0] * n
    for rank, (_idx, pv) in enumerate(indexed, start=1):  # rename idx to _idx (B007)
        q_raw[rank - 1] = (pv * n) / rank

    # monotone non-decreasing adjustment from the end
    q_adj_sorted = [0.0] * n
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, q_raw[i])
        q_adj_sorted[i] = prev

    # map back to original order
    fdr = [0.0] * n
    for new_rank, (orig_idx, _pv) in enumerate(indexed):
        fdr[orig_idx] = q_adj_sorted[new_rank]
    return fdr


def fdr_mask(qvals: Iterable[float], alpha: float) -> list[bool]:
    return [q <= alpha for q in qvals]
