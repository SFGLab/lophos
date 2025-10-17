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
