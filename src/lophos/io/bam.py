# regex-based allele assignment
import re
from os import PathLike

import pysam

from ..constants import MAT_RG_CANDIDATES, PAT_RG_CANDIDATES


# Compile default patterns for read-group identifiers.  These can be overridden
# at runtime by calling ``set_rg_patterns()``.  The default behaviour is to
# treat any RG tag matching one of the values in MAT_RG_CANDIDATES as
# maternal and any RG tag matching PAT_RG_CANDIDATES as paternal.  Matching
# is case-insensitive.
def _compile_default_pattern(names: set[str]) -> re.Pattern[str]:
    # Join candidate tokens into a single alternation pattern.  Use word
    # boundaries so that e.g. "mat" does not match "maternal" twice.  Note
    # that the default sets are small ("maternal", "mat", "M"), so the
    # resulting pattern is simple.
    escaped = [re.escape(n) for n in names]
    pattern = r"^(?:" + "|".join(escaped) + r")$"
    return re.compile(pattern, re.IGNORECASE)


# Global patterns used by ``allele_from_rg``.  They are compiled once at
# module import but may be overwritten by ``set_rg_patterns`` if the user
# specifies custom patterns via the CLI.
MAT_RG_PATTERN: re.Pattern[str] = _compile_default_pattern(MAT_RG_CANDIDATES)
PAT_RG_PATTERN: re.Pattern[str] = _compile_default_pattern(PAT_RG_CANDIDATES)


def set_rg_patterns(maternal_pattern: str | None, paternal_pattern: str | None) -> None:
    """Override the default regex patterns used to map RG tags to maternal and
    paternal alleles.

    Parameters
    ----------
    maternal_pattern : str | None
        A regular expression pattern matching RG tags that should be assigned
        to the maternal allele.  If ``None``, the default pattern based on
        ``MAT_RG_CANDIDATES`` is retained.
    paternal_pattern : str | None
        A regular expression pattern matching RG tags that should be assigned
        to the paternal allele.  If ``None``, the default pattern based on
        ``PAT_RG_CANDIDATES`` is retained.
    """
    global MAT_RG_PATTERN, PAT_RG_PATTERN
    if maternal_pattern:
        MAT_RG_PATTERN = re.compile(maternal_pattern, re.IGNORECASE)
    if paternal_pattern:
        PAT_RG_PATTERN = re.compile(paternal_pattern, re.IGNORECASE)


def open_bam(path: str | PathLike[str]) -> pysam.AlignmentFile:
    return pysam.AlignmentFile(str(path), "rb")


def read_is_duplicate(aln: pysam.AlignedSegment) -> bool:
    return aln.is_duplicate


def allele_from_rg(aln: pysam.AlignedSegment) -> str | None:
    """Return the allele assignment for a read-group tag.

    The function uses the compiled ``MAT_RG_PATTERN`` and ``PAT_RG_PATTERN``
    regular expressions to decide whether the RG tag corresponds to the
    maternal or paternal haplotype.  If no RG tag is present, or if the tag
    does not match either pattern, ``None`` is returned.
    """
    rg = aln.get_tag("RG") if aln.has_tag("RG") else None
    if rg is None:
        return None
    rg_str = str(rg)
    # If the tag matches the maternal pattern and does not match the paternal
    # pattern, assign as maternal.  In cases where a tag could match both
    # patterns, maternal has precedence to preserve legacy behaviour.
    if MAT_RG_PATTERN.search(rg_str) and not PAT_RG_PATTERN.search(rg_str):
        return "maternal"
    if PAT_RG_PATTERN.search(rg_str):
        return "paternal"
    return None


# ---------------------------------------------------------------------------
# SA:Z-based helpers for long-read (ONT/Pore-C) chimeric contact reconstruction
# ---------------------------------------------------------------------------


# A robust CIGAR parser for reference-consumed length (M, D, N, =, X)
def segment_len_from_cigar(cigar: str) -> int:
    """
    Reference-consumed length from CIGAR: sum of M, D, N, =, X.
    (Insertions I, soft/hard clips S/H, pads P do not consume reference.)
    """
    total = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        if not num:
            # malformed piece; skip gracefully
            continue
        n = int(num)
        if ch in ("M", "D", "N", "=", "X"):
            total += n
        # reset accumulator
        num = ""
    return total


# SA:Z tag format per entry: rname,pos,strand,cigar,mapq,nm;
# POS is 1-based in SA; convert to 0-based here.
_SA_ENTRY_RE = re.compile(
    r"(?P<rname>[^,]+),(?P<pos>[0-9]+),(?P<strand>[+-]),(?P<cigar>[^,]+),(?P<mapq>[0-9]+),(?P<nm>[0-9]+)"
)


def parse_sa_tag(sa_str: str) -> list[dict[str, int | str]]:
    """
    Parse an SA:Z string into a list of segment dicts with keys:
      rname, pos0, strand, cigar, mapq, nm, ref_len

    Notes
    -----
    - SA entries are separated by ';' and may end with a trailing ';'.
    - POS in SA is 1-based; we store 0-based ``pos0``.
    - ``ref_len`` is computed as reference-consumed CIGAR length.
    """
    segs: list[dict[str, int | str]] = []
    if not sa_str:
        return segs
    for part in sa_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        m = _SA_ENTRY_RE.fullmatch(part)
        if not m:
            # skip malformed piece silently
            continue
        rname = m.group("rname")
        pos1 = int(m.group("pos"))
        strand = m.group("strand")
        cigar = m.group("cigar")
        mapq = int(m.group("mapq"))
        nm = int(m.group("nm"))
        pos0 = pos1 - 1  # convert to 0-based
        ref_len = segment_len_from_cigar(cigar)
        segs.append(
            {
                "rname": rname,
                "pos0": pos0,
                "strand": strand,
                "cigar": cigar,
                "mapq": mapq,
                "nm": nm,
                "ref_len": ref_len,
            }
        )
    return segs


def iter_read_segments(aln: pysam.AlignedSegment) -> list[dict[str, int | str]]:
    """
    Gather primary + SA segments for a read into a uniform representation.

    Primary alignment:
      rname   = aln.reference_name
      pos0    = aln.reference_start (0-based)
      strand  = '-' if aln.is_reverse else '+'
      cigar   = aln.cigarstring (if None, fall back to "<len>M" for alignment length)
      mapq    = aln.mapping_quality
      nm      = NM tag if present else 0
      ref_len = reference-consumed length computed from CIGAR

    SA segments are parsed from SA:Z using ``parse_sa_tag``.
    """
    if aln.is_unmapped:
        return []

    # Primary segment
    rname = aln.reference_name
    assert rname is not None  # ensured by is_unmapped check above
    pos0 = int(aln.reference_start)
    strand = "-" if aln.is_reverse else "+"
    # If CIGAR is None (rare), approximate with aligned length as M's
    if aln.cigarstring:
        cigar = aln.cigarstring
    else:
        try:
            alen = int(aln.query_alignment_length)
        except Exception:
            alen = 0
        cigar = f"{alen}M"
    mapq = int(aln.mapping_quality)
    try:
        nm = int(aln.get_tag("NM"))
    except Exception:
        nm = 0
    ref_len = segment_len_from_cigar(cigar)

    segs: list[dict[str, int | str]] = [
        {
            "rname": rname,
            "pos0": pos0,
            "strand": strand,
            "cigar": cigar,
            "mapq": mapq,
            "nm": nm,
            "ref_len": ref_len,
        }
    ]

    # SA:Z segments (optional)
    try:
        sa = aln.get_tag("SA")
    except Exception:
        sa = None
    if sa:
        segs.extend(parse_sa_tag(str(sa)))

    return segs
