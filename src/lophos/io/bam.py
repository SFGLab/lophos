import os
from os import PathLike
from typing import Optional

import pysam

from ..constants import MAT_RG_CANDIDATES, PAT_RG_CANDIDATES

def open_bam(path: str | PathLike[str]) -> pysam.AlignmentFile:
    return pysam.AlignmentFile(path, "rb")

def read_is_duplicate(aln: pysam.AlignedSegment) -> bool:
    return aln.is_duplicate

def allele_from_rg(aln: pysam.AlignedSegment) -> Optional[str]:
    rg = aln.get_tag("RG") if aln.has_tag("RG") else None
    if rg is None:
        return None
    rg_lower = str(rg).lower()
    if rg_lower in MAT_RG_CANDIDATES:
        return "maternal"
    if rg_lower in PAT_RG_CANDIDATES:
        return "paternal"
    return None
