import math
from dataclasses import dataclass
from typing import Literal

import pandas as pd
from scipy.stats import binomtest

from ..constants import PSEUDOCOUNT
from ..utils.fdr import bh_fdr


@dataclass(frozen=True)
class PeakStat:
    m: int
    p: int
    total: int
    log2_ratio: float
    p_value: float


def _log2_ratio(m: int, p: int) -> float:
    return math.log2((m + PSEUDOCOUNT) / (p + PSEUDOCOUNT))


def compute_peak_stats(df_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-peak stats:
      - total = m + p
      - log2_ratio = log2((m+pc)/(p+pc))
      - p_value (two-sided binomial at 0.5); when total==0 => p=1
      - fdr via BH over the *entire* p-value vector (rows with p=1 naturally yield q≈1)
    """
    m = df_counts["maternal"].astype(int)
    p = df_counts["paternal"].astype(int)
    total = m + p

    ratio = [_log2_ratio(int(mi), int(pi)) for mi, pi in zip(m, p, strict=False)]
    pvals = [
        (
            binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue
            if (int(mi) + int(pi)) > 0
            else 1.0
        )
        for mi, pi in zip(m, p, strict=False)
    ]
    qvals = bh_fdr(pvals)

    out = df_counts.copy()
    out["total"] = total
    out["log2_ratio"] = ratio
    out["p_value"] = pvals
    out["fdr"] = qvals

    # light sanity (helps catch any future regressions)
    assert (out["p_value"].between(0.0, 1.0)).all()
    assert (out["fdr"].between(0.0, 1.0)).all()
    # if p==1 on a row, q should be very close to 1 (exactly 1.0 in most BH implementations)
    # we allow tiny FP wiggle but never < 0.99
    assert (out.loc[out["p_value"] == 1.0, "fdr"] >= 0.99).all()

    return out


HomozygousPolicy = Literal["drop", "split50", "ambiguous"]


def compute_loop_stats(
    df_counts: pd.DataFrame, *, homozygous_policy: HomozygousPolicy = "drop"
) -> pd.DataFrame:
    """
    Compute statistics for loop counts.

    Input df must contain:
      - maternal_pairs, paternal_pairs
      - ambiguous_pairs (optional; filled with zeros)
      - homozygous_pairs (optional; filled with zeros)

    Definitions:
      - "informative" pairs are those assigned maternal or paternal.
      - "non-informative" pairs are ambiguous + homozygous-like (tie) pairs.

    homozygous_policy controls how `homozygous_pairs` contribute:

      - "drop" (default): homozygous pairs are tracked but do NOT contribute to m/p.
      - "ambiguous": homozygous pairs are added to ambiguous (non-informative).
      - "split50": homozygous pairs are split evenly into maternal and paternal *informative* counts.
                   If homozygous is odd, the remainder stays non-informative.

    Output includes (among others):
      - m_raw, p_raw, ambiguous_pairs_raw, homozygous_pairs_raw
      - m, p  (post-policy informative counts)
      - total_pairs = m + p
      - noninformative_pairs
      - ambiguous_frac = noninformative_pairs / (m+p+noninformative_pairs) [safe when denom==0]
      - log2_ratio_pairs, p_value_pairs, fdr_pairs computed on informative pairs only
    """
    df = df_counts.copy()

    # Normalize / ensure required columns exist
    if "ambiguous_pairs" not in df.columns:
        df["ambiguous_pairs"] = 0
    if "homozygous_pairs" not in df.columns:
        df["homozygous_pairs"] = 0

    m_raw = df["maternal_pairs"].astype(int)
    p_raw = df["paternal_pairs"].astype(int)
    amb_raw = df["ambiguous_pairs"].astype(int)
    hom_raw = df["homozygous_pairs"].astype(int)

    if homozygous_policy not in ("drop", "split50", "ambiguous"):
        raise ValueError(
            f"Invalid homozygous_policy={homozygous_policy!r}. "
            "Expected one of: 'drop', 'split50', 'ambiguous'."
        )

    # Apply policy to derive informative m/p and non-informative pools
    if homozygous_policy == "drop":
        m = m_raw
        p = p_raw
        noninfo = amb_raw + hom_raw
    elif homozygous_policy == "ambiguous":
        m = m_raw
        p = p_raw
        noninfo = amb_raw + hom_raw
        # We do NOT destroy *_raw columns; we only treat hom as noninformative through noninfo.
    else:  # split50
        add = (hom_raw // 2).astype(int)
        rem = (hom_raw % 2).astype(int)
        m = m_raw + add
        p = p_raw + add
        # remainder stays non-informative (cannot be deterministically assigned)
        noninfo = amb_raw + rem

    out = df.rename(columns={"maternal_pairs": "m_raw", "paternal_pairs": "p_raw"}).copy()
    out["ambiguous_pairs_raw"] = df["ambiguous_pairs"].astype(int)
    out["homozygous_pairs_raw"] = df["homozygous_pairs"].astype(int)

    out["m"] = m.astype(int)
    out["p"] = p.astype(int)

    out["total_pairs"] = out["m"] + out["p"]
    out["noninformative_pairs"] = noninfo.astype(int)

    denom = out["m"] + out["p"] + out["noninformative_pairs"]
    denom_safe = denom.where(denom > 0, other=1)
    out["ambiguous_frac"] = out["noninformative_pairs"] / denom_safe

    out["log2_ratio_pairs"] = [
        _log2_ratio(int(mi), int(pi)) for mi, pi in zip(out["m"], out["p"], strict=False)
    ]

    pvals: list[float] = []
    for mi, pi in zip(out["m"], out["p"], strict=False):
        tot = int(mi) + int(pi)
        pvals.append(
            binomtest(int(mi), tot, 0.5, alternative="two-sided").pvalue if tot > 0 else 1.0
        )

    out["p_value_pairs"] = pvals
    out["fdr_pairs"] = bh_fdr(out["p_value_pairs"])

    # light sanity
    assert (out["p_value_pairs"].between(0.0, 1.0)).all()
    assert (out["fdr_pairs"].between(0.0, 1.0)).all()

    return out
