from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import Literal

import pandas as pd
import pysam
from scipy.stats import binomtest

from ..constants import PSEUDOCOUNT
from ..io.bam import allele_from_rg, read_is_duplicate
from ..utils.fdr import bh_fdr


@dataclass(frozen=True)
class AnchorFallbackParams:
    """Parameters controlling anchor-based loop inference.

    This is a *proxy* inference for loops when direct connectivity evidence
    (mates/SA) is insufficient.

    Strategy:
      1) Count allele support independently at anchor A and anchor B (peak-style).
      2) Compute stats + bias calls per anchor (BH-FDR across all anchors).
      3) Infer a loop allele only when both anchors are confidently consistent.

    Notes
    -----
    - Counts are resolved per QNAME to avoid double-counting MAT_PAT merged BAMs.
    - If a QNAME has a best-score tie between maternal and paternal, it is treated
      as `homozygous/uninformative` and does not contribute to maternal/paternal
      informative counts (tracked as QC).
    """

    mapq: int = 30
    keep_duplicates: bool = False

    # Calling thresholds for each anchor
    min_reads_anchor: int = 5
    fdr: float = 0.05
    min_fold: float = 1.5
    min_abs_log2: float = 0.0


AnchorName = Literal["A", "B"]


def _iter_overlaps(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    mapq: int,
    keep_dups: bool,
) -> Iterator[pysam.AlignedSegment]:
    """Yield alignments overlapping [start,end) with basic MAPQ/dup filters."""
    for aln in bam.fetch(chrom, max(0, start), end):
        if aln.is_unmapped or aln.mapping_quality < mapq:
            continue
        if not keep_dups and read_is_duplicate(aln):
            continue
        aln_end = aln.reference_end
        if aln_end is None:
            continue
        if not (aln.reference_start < end and aln_end > start):
            continue
        yield aln


def _aln_score(aln: pysam.AlignedSegment) -> tuple[int, int, int]:
    """Higher is better: (MAPQ, AS, -NM)."""
    mq = int(aln.mapping_quality or 0)
    try:
        ascore = int(aln.get_tag("AS"))
    except Exception:
        ascore = 0
    try:
        nm = int(aln.get_tag("NM"))
    except Exception:
        nm = 10**9
    return (mq, ascore, -nm)


def _log2_ratio(m: int, p: int) -> float:
    import math

    return math.log2((m + PSEUDOCOUNT) / (p + PSEUDOCOUNT))


def build_loop_anchors(loops: pd.DataFrame, anchor_pad: int) -> pd.DataFrame:
    """Expand loops BEDPE rows into an anchors table (two rows per loop).

    Input
    -----
    loops : DataFrame with columns chrom1,start1,end1,chrom2,start2,end2

    Output
    ------
    DataFrame with columns:
      loop_index, anchor, chrom, start, end
    where loop_index is the original loops row index (preserved).
    """
    req = {"chrom1", "start1", "end1", "chrom2", "start2", "end2"}
    missing = req - set(loops.columns)
    if missing:
        raise KeyError(f"loops missing required columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for idx, r in loops.iterrows():
        c1, s1, e1 = str(r["chrom1"]), int(r["start1"]), int(r["end1"])
        c2, s2, e2 = str(r["chrom2"]), int(r["start2"]), int(r["end2"])

        rows.append(
            {
                "loop_index": idx,
                "anchor": "A",
                "chrom": c1,
                "start": max(0, s1 - anchor_pad),
                "end": e1 + anchor_pad,
            }
        )
        rows.append(
            {
                "loop_index": idx,
                "anchor": "B",
                "chrom": c2,
                "start": max(0, s2 - anchor_pad),
                "end": e2 + anchor_pad,
            }
        )
    return pd.DataFrame(rows)


def count_anchor_alleles(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    *,
    mapq: int,
    keep_dups: bool,
) -> tuple[int, int, int]:
    """Count (maternal, paternal, homozygous) support in an anchor region.

    Per-QNAME resolution is applied:
      - best score wins
      - best-score tie between maternal and paternal => homozygous (uninformative)

    Returns
    -------
    (m, p, hom)
    """
    per_qname: dict[str, tuple[str, tuple[int, int, int]]] = {}

    for aln in _iter_overlaps(bam, chrom, start, end, mapq, keep_dups):
        qn = aln.query_name
        if not qn:
            continue
        allele = allele_from_rg(aln)
        if allele not in ("maternal", "paternal"):
            continue
        score = _aln_score(aln)

        if qn not in per_qname:
            per_qname[qn] = (allele, score)
            continue

        prev_a, prev_s = per_qname[qn]

        # already homozygous; only replace if strictly better
        if prev_a == "homozygous":
            if score > prev_s:
                per_qname[qn] = (allele, score)
            continue

        if score > prev_s:
            per_qname[qn] = (allele, score)
            continue
        if score < prev_s:
            continue

        # score tie: conflict between maternal/paternal -> homozygous
        if {prev_a, allele} == {"maternal", "paternal"}:
            per_qname[qn] = ("homozygous", score)

    m = sum(1 for a, _ in per_qname.values() if a == "maternal")
    p = sum(1 for a, _ in per_qname.values() if a == "paternal")
    hom = sum(1 for a, _ in per_qname.values() if a == "homozygous")
    return m, p, hom


def count_loop_anchors(
    bam: pysam.AlignmentFile,
    anchors: pd.DataFrame,
    *,
    mapq: int,
    keep_dups: bool,
) -> pd.DataFrame:
    """Count allele support for each anchor row in an anchors table.

    anchors must contain: loop_index, anchor, chrom, start, end
    """
    req = {"loop_index", "anchor", "chrom", "start", "end"}
    missing = req - set(anchors.columns)
    if missing:
        raise KeyError(f"anchors missing required columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for _, r in anchors.iterrows():
        m, p, hom = count_anchor_alleles(
            bam,
            str(r["chrom"]),
            int(r["start"]),
            int(r["end"]),
            mapq=mapq,
            keep_dups=keep_dups,
        )
        rows.append(
            {
                "loop_index": r["loop_index"],
                "anchor": r["anchor"],
                "chrom": r["chrom"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "maternal": int(m),
                "paternal": int(p),
                "homozygous": int(hom),
            }
        )
    return pd.DataFrame(rows)


def _call_anchor_bias(
    m: int,
    p: int,
    q: float,
    *,
    min_reads: int,
    fdr: float,
    min_fold: float,
    min_abs_log2: float,
) -> str:
    """Anchor bias calling: Maternal/Paternal/Balanced/Undetermined."""
    total = m + p
    if total < min_reads:
        return "Undetermined"
    if q > fdr:
        return "Balanced"

    r = _log2_ratio(m, p)
    if abs(r) < min_abs_log2:
        return "Balanced"

    fold = (m + PSEUDOCOUNT) / (p + PSEUDOCOUNT)
    if fold >= min_fold:
        return "Maternal"
    if (1.0 / fold) >= min_fold:
        return "Paternal"
    return "Balanced"


def compute_anchor_stats_and_calls(
    anchors_counts: pd.DataFrame, params: AnchorFallbackParams
) -> pd.DataFrame:
    """Compute stats + calls for anchors (peak-style) with BH-FDR across all anchors."""
    df = anchors_counts.copy()
    if not {"maternal", "paternal"}.issubset(df.columns):
        raise KeyError("anchors_counts must include 'maternal' and 'paternal' columns")

    m = df["maternal"].astype(int)
    p = df["paternal"].astype(int)
    total = m + p
    df["total"] = total
    df["log2_ratio"] = [_log2_ratio(int(mi), int(pi)) for mi, pi in zip(m, p, strict=False)]
    pvals = [
        (
            binomtest(int(mi), int(mi) + int(pi), 0.5, alternative="two-sided").pvalue
            if (int(mi) + int(pi)) > 0
            else 1.0
        )
        for mi, pi in zip(m, p, strict=False)
    ]
    df["p_value"] = pvals
    df["fdr"] = bh_fdr(pvals)

    calls = [
        _call_anchor_bias(
            int(mi),
            int(pi),
            float(qi),
            min_reads=int(params.min_reads_anchor),
            fdr=float(params.fdr),
            min_fold=float(params.min_fold),
            min_abs_log2=float(params.min_abs_log2),
        )
        for mi, pi, qi in zip(df["maternal"], df["paternal"], df["fdr"], strict=False)
    ]
    df["bias_call"] = calls
    return df


def infer_loop_from_anchor_calls(call_a: str, call_b: str) -> tuple[str, str]:
    """Infer loop allele from two anchor bias calls.

    Rules (strict, explainable):
      - Maternal + Maternal => Maternal
      - Paternal + Paternal => Paternal
      - Any Undetermined => Undetermined (insufficient anchor evidence)
      - Maternal vs Paternal => Undetermined (discordant)
      - Balanced with anything => Undetermined (not enough support to infer)
    """
    if call_a == "Undetermined" or call_b == "Undetermined":
        return "Undetermined", "anchor_undetermined"
    if call_a == "Balanced" or call_b == "Balanced":
        return "Undetermined", "anchor_balanced"
    if call_a == call_b and call_a in ("Maternal", "Paternal"):
        return call_a, "anchors_consistent"
    return "Undetermined", "anchors_discordant"


def infer_loops_from_anchors(anchors_called: pd.DataFrame) -> pd.DataFrame:
    """Collapse anchor calls into one inferred loop call per loop_index.

    Input anchors_called must include: loop_index, anchor, bias_call

    Output:
      loop_index, bias_call_inferred, inferred_reason
    """
    req = {"loop_index", "anchor", "bias_call"}
    missing = req - set(anchors_called.columns)
    if missing:
        raise KeyError(f"anchors_called missing required columns: {sorted(missing)}")

    piv = anchors_called.pivot_table(
        index="loop_index", columns="anchor", values="bias_call", aggfunc="first"
    )

    if "A" not in piv.columns or "B" not in piv.columns:
        out = pd.DataFrame({"loop_index": piv.index})
        out["bias_call_inferred"] = "Undetermined"
        out["inferred_reason"] = "missing_anchor"
        return out.reset_index(drop=True)

    inferred = []
    reason = []
    for ca, cb in zip(piv["A"], piv["B"], strict=False):
        call, why = infer_loop_from_anchor_calls(str(ca), str(cb))
        inferred.append(call)
        reason.append(why)

    out = pd.DataFrame(
        {
            "loop_index": piv.index,
            "bias_call_inferred": inferred,
            "inferred_reason": reason,
        }
    )
    return out.reset_index(drop=True)
