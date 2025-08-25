from __future__ import annotations

from typing import TypedDict

import pandas as pd
import pysam

from ..io.bam import allele_from_rg, read_is_duplicate


class LoopRow(TypedDict):
    loop_id: str
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    maternal_pairs: int
    paternal_pairs: int
    ambiguous_pairs: int


def _in_interval(pos: int, start: int, end: int) -> bool:
    """Return True if pos is within [start, end)."""
    return start <= pos < end


def count_loops(
    bam: pysam.AlignmentFile,
    loops: pd.DataFrame,
    mapq: int,
    anchor_pad: int,
    keep_dups: bool,
) -> pd.DataFrame:
    """
    Count allele-consistent pairs connecting loop anchors.

    Notes
    -----
    - We assume mates share the same RG tag (common in phased pipelines).
      If that assumption breaks, pairs are counted as 'ambiguous'.
    - Expected columns in `loops`: chrom1, start1, end1, chrom2, start2, end2, [name?].
    - Output columns: maternal_pairs, paternal_pairs, ambiguous_pairs.
    """
    rows: list[LoopRow] = []

    for idx, row in loops.iterrows():
        chr1, s1, e1 = str(row["chrom1"]), int(row["start1"]), int(row["end1"])
        chr2, s2, e2 = str(row["chrom2"]), int(row["start2"]), int(row["end2"])

        a1s, a1e = s1 - anchor_pad, e1 + anchor_pad
        a2s, a2e = s2 - anchor_pad, e2 + anchor_pad

        mm = pp = amb = 0

        for aln in bam.fetch(chr1, max(0, a1s), a1e):
            if aln.is_unmapped or aln.mapping_quality < mapq:
                continue
            if not keep_dups and read_is_duplicate(aln):
                continue
            if not aln.is_paired or aln.next_reference_id < 0:
                continue

            mate_chr = bam.get_reference_name(aln.next_reference_id)
            mate_pos = aln.next_reference_start

            if mate_chr == chr2 and _in_interval(mate_pos, a2s, a2e):
                a = allele_from_rg(aln)
                # In most phased BAMs, mates share RG; if not, treat as ambiguous.
                b = allele_from_rg(aln)
                if a == "maternal" and b == "maternal":
                    mm += 1
                elif a == "paternal" and b == "paternal":
                    pp += 1
                else:
                    amb += 1

        loop_id = str(row["name"]) if "name" in loops.columns else f"loop_{idx}"
        d: LoopRow = {
            "loop_id": loop_id,
            "chrom1": chr1,
            "start1": s1,
            "end1": e1,
            "chrom2": chr2,
            "start2": s2,
            "end2": e2,
            "maternal_pairs": mm,
            "paternal_pairs": pp,
            "ambiguous_pairs": amb,
        }
        rows.append(d)

    return pd.DataFrame(rows)
