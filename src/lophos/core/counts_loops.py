from __future__ import annotations

from collections.abc import Iterator
from typing import TypedDict

import pandas as pd
import pysam

from ..io.bam import allele_from_rg, read_is_duplicate


class LoopRow(TypedDict):
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


def _endpoint_in_anchor(
    chrom: str, pos_start: int, pos_end: int, a_chrom: str, a_start: int, a_end: int
) -> bool:
    """
    Return True if an endpoint (chrom, [pos_start,pos_end)) overlaps the anchor interval [a_start,a_end) on a_chrom.
    We accept either endpoint coordinate; for mate-based we only have a point, for SA we may have a short interval.
    """
    if chrom != a_chrom:
        return False
    # Treat endpoint as an interval; overlap if any base intersects.
    return not (pos_end <= a_start or pos_start >= a_end)


def _contact_hits_anchors(
    ep1: tuple[str, int, int],
    ep2: tuple[str, int, int],
    a1: tuple[str, int, int],
    a2: tuple[str, int, int],
) -> bool:
    """
    Return True if the two endpoints (ep1, ep2) land in the two anchors (a1, a2) in either order.
    Endpoints are (chrom, start, end); anchors are (chrom, a_start, a_end).
    """
    c1, s1, e1 = ep1
    c2, s2, e2 = ep2
    achr1, as1, ae1 = a1
    achr2, as2, ae2 = a2

    # ep1 in anchor1 AND ep2 in anchor2 (or swapped)
    return (
        _endpoint_in_anchor(c1, s1, e1, achr1, as1, ae1)
        and _endpoint_in_anchor(c2, s2, e2, achr2, as2, ae2)
    ) or (
        _endpoint_in_anchor(c1, s1, e1, achr2, as2, ae2)
        and _endpoint_in_anchor(c2, s2, e2, achr1, as1, ae1)
    )


# ------------------------ common small helpers ------------------------


def _iter_reads_in_region(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    mapq: int,
    keep_dups: bool,
) -> Iterator[pysam.AlignedSegment]:
    """Yield reads from region with basic MAPQ/dup filtering."""
    for aln in bam.fetch(chrom, max(0, start), end):
        if aln.is_unmapped or aln.mapping_quality < mapq:
            continue
        if not keep_dups and read_is_duplicate(aln):
            continue
        yield aln


# ------------------------ mates mode ------------------------


def _counts_for_single_loop_mates(
    bam: pysam.AlignmentFile,
    chr1: str,
    s1: int,
    e1: int,
    chr2: str,
    s2: int,
    e2: int,
    anchor_pad: int,
    mapq: int,
    keep_dups: bool,
) -> tuple[int, int, int]:
    """
    Count mate-pairs that bridge the two anchors (±pad) in either direction.
    Assumes mates share RG; if not, class as ambiguous.
    """
    a1s, a1e = s1 - anchor_pad, e1 + anchor_pad
    a2s, a2e = s2 - anchor_pad, e2 + anchor_pad

    mm = pp = amb = 0
    seen_qnames: set[str] = set()  # avoid double-counting the same template

    # Fetch reads from both anchors
    for fchr, fs, fe in ((chr1, a1s, a1e), (chr2, a2s, a2e)):
        for aln in _iter_reads_in_region(bam, fchr, fs, fe, mapq, keep_dups):
            if not aln.is_paired or aln.next_reference_id < 0:
                continue

            qn = aln.query_name
            if qn is None:
                continue
            if qn in seen_qnames:
                continue

            mate_chr = bam.get_reference_name(aln.next_reference_id)
            mate_pos = aln.next_reference_start

            hit_a1_to_a2 = (aln.reference_name == chr1) and _in_interval(mate_pos, a2s, a2e)
            hit_a2_to_a1 = (aln.reference_name == chr2) and _in_interval(mate_pos, a1s, a1e)
            if (mate_chr == chr2 and hit_a1_to_a2) or (mate_chr == chr1 and hit_a2_to_a1):
                seen_qnames.add(qn)
                a = allele_from_rg(aln)  # assume mate has same RG in phased pipelines
                if a == "maternal":
                    mm += 1
                elif a == "paternal":
                    pp += 1
                else:
                    amb += 1

    return mm, pp, amb


def _count_loops_mates(
    bam: pysam.AlignmentFile,
    loops: pd.DataFrame,
    mapq: int,
    anchor_pad: int,
    keep_dups: bool,
) -> pd.DataFrame:
    """
    Mate-pair based loop counting (paired-end HiChIP/Hi-C style).
    Counts a pair if one read maps in anchor A (±pad) and its *mate* maps in anchor B (±pad), in either direction.
    """
    rows: list[LoopRow] = []
    for _, row in loops.iterrows():
        chr1, s1, e1 = str(row["chrom1"]), int(row["start1"]), int(row["end1"])
        chr2, s2, e2 = str(row["chrom2"]), int(row["start2"]), int(row["end2"])
        mm, pp, amb = _counts_for_single_loop_mates(
            bam, chr1, s1, e1, chr2, s2, e2, anchor_pad, mapq, keep_dups
        )
        rows.append(
            {
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
        )
    return pd.DataFrame(rows)


# ------------------------ SA mode ------------------------


def _counts_for_single_loop_sa(
    bam: pysam.AlignmentFile,
    chr1: str,
    s1: int,
    e1: int,
    chr2: str,
    s2: int,
    e2: int,
    anchor_pad: int,
    mapq: int,
    keep_dups: bool,
    *,
    sa_min_mapq: int,
    sa_min_seg_len: int,
    sa_min_cis_dist: int,
    sa_allow_trans: bool,
    sa_orientation: str,
    sa_dedup_within_read: bool,
) -> tuple[int, int, int]:
    """
    For each read overlapping either padded anchor, reconstruct adjacent-segment contacts within the read
    (using SA:Z), then count contacts that span the two anchors (in either direction). Allele from RG.
    """
    from .sa_pairs import build_contacts

    a1 = (chr1, s1 - anchor_pad, e1 + anchor_pad)
    a2 = (chr2, s2 - anchor_pad, e2 + anchor_pad)

    mm = pp = amb = 0
    seen_reads: set[str] = set()  # avoid processing same read twice for this loop

    for fchr, fs, fe in ((a1[0], a1[1], a1[2]), (a2[0], a2[1], a2[2])):
        for aln in _iter_reads_in_region(bam, fchr, fs, fe, mapq, keep_dups):
            qn = aln.query_name
            if qn is None:
                continue
            if qn in seen_reads:
                continue
            seen_reads.add(qn)

            contacts = build_contacts(
                aln,
                min_mapq=sa_min_mapq,
                min_seg_len=sa_min_seg_len,
                min_cis_dist=sa_min_cis_dist,
                allow_trans=sa_allow_trans,
                orientation=sa_orientation,
                dedup_within_read=sa_dedup_within_read,
            )
            if not contacts:
                continue

            for c in contacts:
                ep1 = (str(c["chrom1"]), int(c["start1"]), int(c["end1"]))
                ep2 = (str(c["chrom2"]), int(c["start2"]), int(c["end2"]))
                if not _contact_hits_anchors(ep1, ep2, a1, a2):
                    continue
                allele = allele_from_rg(aln)
                if allele == "maternal":
                    mm += 1
                elif allele == "paternal":
                    pp += 1
                else:
                    amb += 1

    return mm, pp, amb


def _count_loops_sa(
    bam: pysam.AlignmentFile,
    loops: pd.DataFrame,
    mapq: int,
    anchor_pad: int,
    keep_dups: bool,
    *,
    sa_min_mapq: int,
    sa_min_seg_len: int,
    sa_min_cis_dist: int,
    sa_allow_trans: bool,
    sa_orientation: str,
    sa_dedup_within_read: bool,
) -> pd.DataFrame:
    """
    SA:Z-based loop counting (long-read chimeric mode).
    For each loop, count reconstructed contacts that span the two anchors (in either direction).
    """
    rows: list[LoopRow] = []
    for _, row in loops.iterrows():
        chr1, s1, e1 = str(row["chrom1"]), int(row["start1"]), int(row["end1"])
        chr2, s2, e2 = str(row["chrom2"]), int(row["start2"]), int(row["end2"])
        mm, pp, amb = _counts_for_single_loop_sa(
            bam,
            chr1,
            s1,
            e1,
            chr2,
            s2,
            e2,
            anchor_pad,
            mapq,
            keep_dups,
            sa_min_mapq=sa_min_mapq,
            sa_min_seg_len=sa_min_seg_len,
            sa_min_cis_dist=sa_min_cis_dist,
            sa_allow_trans=sa_allow_trans,
            sa_orientation=sa_orientation,
            sa_dedup_within_read=sa_dedup_within_read,
        )
        rows.append(
            {
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
        )
    return pd.DataFrame(rows)


# ------------------------ public dispatcher ------------------------


def count_loops(
    bam: pysam.AlignmentFile,
    loops: pd.DataFrame,
    mapq: int,
    anchor_pad: int,
    keep_dups: bool,
    *,
    loop_mode: str = "mates",
    sa_min_mapq: int = 30,
    sa_min_seg_len: int = 50,
    sa_min_cis_dist: int = 1000,
    sa_allow_trans: bool = True,
    sa_orientation: str = "any",
    sa_dedup_within_read: bool = True,
) -> pd.DataFrame:
    """
    Dispatch counting by loop mode.

    Parameters
    ----------
    loop_mode : {'mates','sa'}
        'mates' -> paired-end mate logic (default, previous behavior).
        'sa'    -> SA:Z-based split-read reconstruction (long-read mode).

    SA parameters are used only when loop_mode='sa'.
    """
    mode = (loop_mode or "mates").lower()
    if mode == "mates":
        return _count_loops_mates(
            bam=bam,
            loops=loops,
            mapq=mapq,
            anchor_pad=anchor_pad,
            keep_dups=keep_dups,
        )
    if mode == "sa":
        return _count_loops_sa(
            bam=bam,
            loops=loops,
            mapq=mapq,
            anchor_pad=anchor_pad,
            keep_dups=keep_dups,
            sa_min_mapq=sa_min_mapq,
            sa_min_seg_len=sa_min_seg_len,
            sa_min_cis_dist=sa_min_cis_dist,
            sa_allow_trans=sa_allow_trans,
            sa_orientation=sa_orientation,
            sa_dedup_within_read=sa_dedup_within_read,
        )
    raise ValueError(f"Unknown loop_mode '{loop_mode}'. Expected 'mates' or 'sa'.")
