from typing import Iterator

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

def count_peaks(
    bam: pysam.AlignmentFile, peaks: pd.DataFrame, mapq: int, window_bp: int, keep_dups: bool
) -> pd.DataFrame:
    rows = []
    for idx, row in peaks.iterrows():
        chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
        center = (start + end) // 2
        wstart, wend = center - window_bp, center + window_bp
        m = p = 0
        for aln in _iter_overlaps(bam, chrom, wstart, wend, mapq, keep_dups):
            allele = allele_from_rg(aln)
            if allele == "maternal":
                m += 1
            elif allele == "paternal":
                p += 1
        rows.append(
            {
                "peak_id": row.get("name", f"peak_{idx}"),
                "chrom": chrom,
                "start": start,
                "end": end,
                "maternal": m,
                "paternal": p,
            }
        )
    return pd.DataFrame(rows)
