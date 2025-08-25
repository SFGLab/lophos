from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass


@dataclass(frozen=True)
class MotifCheck:
    has_ctcf_motif: bool
    orientation: str | None  # "convergent", "divergent", "same", or None
    note: str


def check_ctcf_motif_for_anchors(
    anchors: Iterable[tuple[str, int, int]],
    fasta_path: str | None = None,
    _jaspar_pwm_id: str = "MA0139.1",
) -> list[MotifCheck]:
    """
    Placeholder CTCF motif checker.

    - If no FASTA provided, returns a neutral 'disabled' result.
    - With FASTA given, still returns a placeholder until real scanning is implemented.
    - Public signature is stable so we can drop in a real scanner later.
    """
    if fasta_path is None:
        return [MotifCheck(False, None, "motif_check_disabled") for _ in anchors]

    # Future: scan sequence in FASTA with JASPAR PWM and set orientation.
    return [MotifCheck(False, None, "motif_scan_placeholder") for _ in anchors]
