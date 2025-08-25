from dataclasses import dataclass
from typing import Optional, Iterable, Tuple, List
import logging

logger = logging.getLogger("lophos.motif")

@dataclass(frozen=True)
class MotifCheck:
    has_ctcf_motif: bool
    orientation: Optional[str]  # "convergent", "divergent", "same", or None
    note: str

def check_ctcf_motif_for_anchors(
    anchors: Iterable[Tuple[str, int, int]],
    fasta_path: Optional[str] = None,
    _jaspar_pwm_id: str = "MA0139.1",
) -> List[MotifCheck]:
    """
    Placeholder CTCF motif checker:
    - If no FASTA provided or motif extras not installed, returns 'unknown' but valid.
    - Keeps interface stable for later real motif scanning.
    """
    try:
        if fasta_path is None:
            raise RuntimeError("FASTA not provided")
        # Attempt to import optional deps; if missing, fall back
        import pyfaidx  # noqa: F401
        import Bio  # noqa: F401
        # Real scanning would happen here in a later version.
        # For now, return neutral results to avoid breaking pipelines.
        return [MotifCheck(False, None, "motif_scan_placeholder") for _ in anchors]
    except Exception as e:
        logger.info("Motif check skipped: %s", e)
        return [MotifCheck(False, None, "motif_check_disabled") for _ in anchors]
