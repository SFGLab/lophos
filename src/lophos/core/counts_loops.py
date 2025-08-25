import pandas as pd
import pysam
from ..io.bam import allele_from_rg, read_is_duplicate

def _in_interval(pos: int, start: int, end: int) -> bool:
    return start <= pos < end

def count_loops(bam: pysam.AlignmentFile, loops: pd.DataFrame, mapq: int, anchor_pad: int, keep_dups: bool) -> pd.DataFrame:
    rows = []
    for idx, row in loops.iterrows():
        chr1, s1, e1 = row["chrom1"], int(row["start1"]), int(row["end1"])
        chr2, s2, e2 = row["chrom2"], int(row["start2"]), int(row["end2"])
        a1s, a1e = s1 - anchor_pad, e1 + anchor_pad
        a2s, a2e = s2 - anchor_pad, e2 + anchor_pad
        mm = pp = amb = 0

        # Iterate reads from anchor 1
        for aln in bam.fetch(chr1, max(0, a1s), a1e):
            if aln.is_unmapped or aln.mapping_quality < mapq:
                continue
            if not keep_dups and read_is_duplicate(aln):
                continue
            if aln.next_reference_id < 0:
                continue
            mate_chr = bam.get_reference_name(aln.next_reference_id)
            mate_pos = aln.next_reference_start
            if mate_chr == chr2 and _in_interval(chr2, mate_pos, a2s, a2e):
                a = allele_from_rg(aln)
                b = allele_from_rg(aln) if aln.has_tag("RG") else None  # NOTE: mate RG not directly stored; assume same RG tag for pair
                # In phased alignments, mates should share RG. If not, mark ambiguous.
                if a == "maternal" and b == "maternal":
                    mm += 1
                elif a == "paternal" and b == "paternal":
                    pp += 1
                else:
                    amb += 1
        rows.append({
            "loop_id": row.get("name", f"loop_{idx}"),
            "chrom1": chr1, "start1": s1, "end1": e1,
            "chrom2": chr2, "start2": s2, "end2": e2,
            "maternal_pairs": mm, "paternal_pairs": pp, "ambiguous_pairs": amb
        })
    return pd.DataFrame(rows)
