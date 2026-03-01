from __future__ import annotations

from collections.abc import Iterator

import pandas as pd
import pysam

from ..io.bam import allele_from_rg, read_is_duplicate


def _iter_overlaps(
    bam: pysam.AlignmentFile, chrom: str, start: int, end: int, mapq: int, keep_dups: bool
) -> Iterator[pysam.AlignedSegment]:
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
    """
    Higher is better.
    (MAPQ, AS, -NM) to prefer:
      - higher mapping quality
      - then higher alignment score
      - then fewer edits
    """
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


def _normalize_allele(a: str | None) -> str | None:
    if a in ("maternal", "paternal"):
        return a
    return None  # treat unknown as not informative for peak phasing


def _update_per_qname_peaks(
    per_qname: dict[str, tuple[str, tuple[int, int, int]]],
    qn: str,
    allele: str,
    score: tuple[int, int, int],
) -> None:
    """Update per_qname dict with best-score tie-breaking logic for peaks mode."""
    if qn not in per_qname:
        per_qname[qn] = (allele, score)
        return

    prev_allele, prev_score = per_qname[qn]

    # If already marked homozygous, keep it unless we see strictly better evidence
    if prev_allele == "homozygous":
        if score > prev_score:
            per_qname[qn] = (allele, score)
        return

    if score > prev_score:
        per_qname[qn] = (allele, score)
        return
    if score < prev_score:
        return

    # Tie: if maternal vs paternal conflict -> homozygous/uninformative
    if {prev_allele, allele} == {"maternal", "paternal"}:
        per_qname[qn] = ("homozygous", score)
    # else: same allele tie -> keep as-is


def count_peaks(
    bam: pysam.AlignmentFile, peaks: pd.DataFrame, mapq: int, window_bp: int, keep_dups: bool
) -> pd.DataFrame:
    """
    Count allele-specific peak support within a fixed window around peak center.

    IMPORTANT for MAT_PAT merged BAMs:
      The same molecule (QNAME) may appear under both maternal and paternal RG.
      We resolve per-QNAME to avoid double counting and artificial bias.

    Tie handling:
      If maternal and paternal records for the same QNAME tie at best score,
      the molecule is treated as "homozygous/uninformative" and does not
      contribute to maternal/paternal counts (tracked as QC in `homozygous`).
    """
    rows = []
    for idx, row in peaks.iterrows():
        chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
        center = (start + end) // 2
        wstart, wend = center - window_bp, center + window_bp

        # Per-QNAME resolution: qname -> (allele, score) or ("homozygous", score)
        per_qname: dict[str, tuple[str, tuple[int, int, int]]] = {}

        for aln in _iter_overlaps(bam, chrom, wstart, wend, mapq, keep_dups):
            qn = aln.query_name
            if not qn:
                continue

            allele = _normalize_allele(allele_from_rg(aln))
            if allele is None:
                continue

            score = _aln_score(aln)
            _update_per_qname_peaks(per_qname, qn, allele, score)

        m = sum(1 for a, _ in per_qname.values() if a == "maternal")
        p = sum(1 for a, _ in per_qname.values() if a == "paternal")
        hom = sum(1 for a, _ in per_qname.values() if a == "homozygous")

        rows.append(
            {
                "peak_id": row.get("name", f"peak_{idx}"),
                "chrom": chrom,
                "start": start,
                "end": end,
                "maternal": m,
                "paternal": p,
                # QC only; stats/writer ignore this unless you choose to surface it later
                "homozygous": hom,
            }
        )

    return pd.DataFrame(rows)
