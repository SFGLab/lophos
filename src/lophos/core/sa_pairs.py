# src/lophos/core/sa_pairs.py
from __future__ import annotations

import pysam

from ..io import bam as bam_io

Endpoint = tuple[str, int, int, str]  # chrom, start, end, strand ('+' or '-')


def _order_key(ep: Endpoint) -> tuple[str, int, int, str]:
    # Sort by chrom, then coordinate. This is acceptable as an approximation of read order.
    chrom, s, e, strand = ep
    return (chrom, s, e, strand)


def _is_convergent(strand1: str, strand2: str) -> bool:
    # Convergent = --> <-- (i.e., + then - when left->right)
    # We’ll enforce this only for short-range cis when requested.
    return (strand1 == "+" and strand2 == "-") or (strand1 == "-" and strand2 == "+")


def _filter_segments_for_endpoints(
    segs: list[dict[str, int | str]],
    *,
    min_mapq: int,
    min_seg_len: int,
) -> list[Endpoint]:
    """Filter by per-segment MAPQ and length on reference; return endpoints."""
    filtered: list[Endpoint] = []
    for s in segs:
        if int(s["mapq"]) < min_mapq:
            continue
        if int(s["ref_len"]) < min_seg_len:
            continue
        chrom = str(s["rname"])
        start0 = int(s["pos0"])
        end0 = start0 + int(s["ref_len"])
        strand = str(s["strand"])
        filtered.append((chrom, start0, end0, strand))
    return filtered


def _accept_pair(
    c1: str,
    s1: int,
    e1: int,
    st1: str,
    c2: str,
    s2: int,
    e2: int,
    st2: str,
    *,
    min_cis_dist: int,
    allow_trans: bool,
    orientation: str,
) -> bool:
    """Inter-chrom/ cis distance / optional orientation gating for a candidate pair."""
    # Inter-chrom filtering
    if c1 != c2:
        return allow_trans

    # cis distance based on segment midpoints
    mid1 = (s1 + e1) // 2
    mid2 = (s2 + e2) // 2
    cis_dist = abs(mid2 - mid1)
    if cis_dist < min_cis_dist:
        # self/near ligation artifacts
        return False

    if orientation == "convergent-short-cis" and cis_dist < 100_000:
        # heuristic: only enforce convergent when short-range
        return _is_convergent(st1, st2)

    return True


def _pair_adjacent_endpoints(
    endpoints: list[Endpoint],
    *,
    min_cis_dist: int,
    allow_trans: bool,
    orientation: str,
    dedup_within_read: bool,
) -> list[dict[str, int | str]]:
    """Sort endpoints and pair ADJACENT only (A-B, B-C, ...), applying filters."""
    endpoints.sort(key=_order_key)

    contacts: list[dict[str, int | str]] = []
    seen: set[tuple[str, int, int, str, int, int]] = set()

    for i in range(len(endpoints) - 1):
        c1, s1, e1, st1 = endpoints[i]
        c2, s2, e2, st2 = endpoints[i + 1]

        if not _accept_pair(
            c1,
            s1,
            e1,
            st1,
            c2,
            s2,
            e2,
            st2,
            min_cis_dist=min_cis_dist,
            allow_trans=allow_trans,
            orientation=orientation,
        ):
            continue

        # dedup within a read
        key = (c1, s1, e1, c2, s2, e2)
        if dedup_within_read and key in seen:
            continue
        seen.add(key)

        contacts.append(
            {
                "chrom1": c1,
                "start1": int(s1),
                "end1": int(e1),
                "chrom2": c2,
                "start2": int(s2),
                "end2": int(e2),
            }
        )
    return contacts


def build_contacts(
    aln: pysam.AlignedSegment,
    *,
    min_mapq: int = 30,
    min_seg_len: int = 50,
    min_cis_dist: int = 1000,
    allow_trans: bool = True,
    orientation: str = "any",  # 'any' | 'convergent-short-cis'
    dedup_within_read: bool = True,
) -> list[dict[str, int | str]]:
    """
    Construct "virtual contacts" from a single long read using its primary+SA segments.

    Returns a list of dicts with keys:
      chrom1, start1, end1, chrom2, start2, end2
    (Coordinates are 0-based half-open.)
    """
    # Gather segments: primary + SA
    segs = bam_io.iter_read_segments(aln)  # List[SegmentDict]
    if not segs:
        return []

    # Filter by per-segment MAPQ and length on reference
    endpoints = _filter_segments_for_endpoints(
        segs,
        min_mapq=min_mapq,
        min_seg_len=min_seg_len,
    )
    if len(endpoints) < 2:
        return []

    # Sort endpoints by (chrom, coord); pair ADJACENT only (A-B, B-C, ...)
    return _pair_adjacent_endpoints(
        endpoints,
        min_cis_dist=min_cis_dist,
        allow_trans=allow_trans,
        orientation=orientation,
        dedup_within_read=dedup_within_read,
    )
