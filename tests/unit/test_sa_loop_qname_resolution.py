from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

import pytest
from lophos.core.counts_loops import _counts_for_single_loop_sa


@dataclass
class FakeAln:
    query_name: str
    reference_name: str
    reference_start: int
    mapping_quality: int
    is_unmapped: bool = False
    is_duplicate: bool = False
    _tags: dict[str, Any] = field(default_factory=dict)

    def has_tag(self, tag: str) -> bool:
        return tag in self._tags

    def get_tag(self, tag: str) -> Any:
        return self._tags[tag]


class FakeBam:
    def __init__(self, reads: list[FakeAln]):
        self._reads = reads

    def fetch(self, chrom: str, start: int, end: int) -> Iterator[FakeAln]:
        # Minimal fetch: return reads on chrom whose start is within [start,end)
        for r in self._reads:
            if r.reference_name == chrom and start <= r.reference_start < end:
                yield r


@pytest.fixture
def monkeypatch_build_contacts(monkeypatch: pytest.MonkeyPatch):
    # Patch build_contacts used inside _counts_for_single_loop_sa.
    # It returns two contacts, one that hits the anchors and one that doesn't.
    from lophos.core import sa_pairs

    def _fake_build_contacts(_aln: FakeAln, **_kwargs: Any) -> list[dict[str, Any]]:
        return [
            # Hitting contact between chr1:100-110 and chr2:200-210
            {
                "chrom1": "chr1",
                "start1": 100,
                "end1": 101,
                "chrom2": "chr2",
                "start2": 200,
                "end2": 201,
            },
            # Non-hitting contact
            {
                "chrom1": "chr1",
                "start1": 9990,
                "end1": 9991,
                "chrom2": "chr2",
                "start2": 9990,
                "end2": 9991,
            },
        ]

    monkeypatch.setattr(sa_pairs, "build_contacts", _fake_build_contacts)


def _run_sa(reads: list[FakeAln]) -> tuple[int, int, int]:
    bam = FakeBam(reads)
    return _counts_for_single_loop_sa(
        bam=bam,  # type: ignore[arg-type]
        chr1="chr1",
        s1=100,
        e1=110,
        chr2="chr2",
        s2=200,
        e2=210,
        anchor_pad=0,
        mapq=0,
        keep_dups=True,
        sa_min_mapq=0,
        sa_min_seg_len=0,
        sa_min_cis_dist=0,
        sa_allow_trans=True,
        sa_orientation="any",
        sa_dedup_within_read=True,
    )


def test_sa_qname_best_alignment_wins(_monkeypatch_build_contacts: None) -> None:
    # Same QNAME appears twice (maternal + paternal). Paternal has better score -> paternal wins.
    reads = [
        FakeAln(
            query_name="q1",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=30,
            _tags={"RG": "maternal", "AS": 10, "NM": 2},
        ),
        FakeAln(
            query_name="q1",
            reference_name="chr2",
            reference_start=205,
            mapping_quality=40,
            _tags={"RG": "paternal", "AS": 12, "NM": 1},
        ),
    ]
    mm, pp, amb = _run_sa(reads)
    assert (mm, pp, amb) == (0, 1, 0)


def test_sa_qname_tie_conflict_becomes_ambiguous(_monkeypatch_build_contacts: None) -> None:
    # Same QNAME appears twice with identical score but conflicting allele -> ambiguous.
    reads = [
        FakeAln(
            query_name="q1",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=40,
            _tags={"RG": "maternal", "AS": 10, "NM": 1},
        ),
        FakeAln(
            query_name="q1",
            reference_name="chr2",
            reference_start=205,
            mapping_quality=40,
            _tags={"RG": "paternal", "AS": 10, "NM": 1},
        ),
    ]
    mm, pp, amb = _run_sa(reads)
    assert (mm, pp, amb) == (0, 0, 1)


def test_sa_counts_one_per_qname(_monkeypatch_build_contacts: None) -> None:
    # Multiple alignments or multiple contacts should not inflate beyond one per QNAME.
    reads = [
        FakeAln(
            query_name="q1",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=30,
            _tags={"RG": "maternal", "AS": 10, "NM": 1},
        ),
        FakeAln(
            query_name="q2",
            reference_name="chr2",
            reference_start=205,
            mapping_quality=30,
            _tags={"RG": "paternal", "AS": 10, "NM": 1},
        ),
        # Duplicate record for q2 with worse score should not change count
        FakeAln(
            query_name="q2",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=10,
            _tags={"RG": "paternal", "AS": 1, "NM": 10},
        ),
    ]
    mm, pp, amb = _run_sa(reads)
    assert (mm, pp, amb) == (1, 1, 0)


def test_sa_multiple_qnames_count_once_each(_monkeypatch_build_contacts: None) -> None:
    reads = [
        FakeAln(
            query_name="q1",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=50,
            _tags={"RG": "maternal", "AS": 20, "NM": 0},
        ),
        FakeAln(
            query_name="q2",
            reference_name="chr2",
            reference_start=205,
            mapping_quality=50,
            _tags={"RG": "paternal", "AS": 20, "NM": 0},
        ),
        FakeAln(
            query_name="q3",
            reference_name="chr1",
            reference_start=105,
            mapping_quality=50,
            _tags={"AS": 20, "NM": 0},  # no RG -> ambiguous
        ),
    ]
    mm, pp, amb = _run_sa(reads)
    assert (mm, pp, amb) == (1, 1, 1)
