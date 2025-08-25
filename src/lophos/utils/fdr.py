from typing import Iterable, List, Tuple

def bh_fdr(pvals: Iterable[float]) -> List[float]:
    p = list(pvals)
    n = len(p)
    indexed = sorted(enumerate(p), key=lambda t: t[1])
    fdr = [0.0] * n
    prev = 1.0
    for rank, (idx, pv) in enumerate(indexed, start=1):
        q = pv * n / rank
        prev = min(prev, q)
        fdr[idx] = prev
    return fdr

def fdr_mask(qvals: Iterable[float], alpha: float) -> List[bool]:
    return [q <= alpha for q in qvals]
