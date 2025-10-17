# src/lophos/report/summary.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rich.console import Console

PEAKS_COLS = [
    "chrom",
    "start",
    "end",
    "peak_id",
    "maternal",
    "paternal",
    "total",
    "log2_ratio",
    "p_value",
    "fdr",
    "bias_call",
]

# Column schema for loops (.loops.bedpe).  The columns correspond to those
# written by ``lophos.report.writers.write_loops``.  The order of these
# columns is frozen to enable robust downstream parsing and testing.
LOOPS_COLS = [
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "loop_id",
    "maternal_pairs",
    "paternal_pairs",
    "ambiguous_pairs",
    "total_pairs",
    "log2_ratio_pairs",
    "p_value_pairs",
    "fdr_pairs",
    "bias_call",
]


@dataclass(frozen=True)
class SummaryParams:
    out: Path
    prefix: str | None
    fdr: float
    min_reads_peak: int
    min_pairs_loop: int
    write_tsv: bool = True


def _detect_prefix(out: Path) -> str:
    """Infer a common prefix from *.peaks.bed and *.loops.bedpe in `out`."""
    peaks = {p.name.replace(".peaks.bed", "") for p in out.glob("*.peaks.bed")}
    loops = {p.name.replace(".loops.bedpe", "") for p in out.glob("*.loops.bedpe")}
    common = peaks & loops
    if len(common) == 1:
        return next(iter(common))
    raise RuntimeError(
        "Could not uniquely infer prefix. "
        f"peaks={sorted(peaks)} loops={sorted(loops)}; "
        "pass --prefix explicitly."
    )


def _paths(out: Path, prefix: str) -> tuple[Path, Path]:
    pk = out / f"{prefix}.peaks.bed"
    lp = out / f"{prefix}.loops.bedpe"
    if not pk.exists():
        raise FileNotFoundError(f"Missing peaks file: {pk}")
    if not lp.exists():
        raise FileNotFoundError(f"Missing loops file: {lp}")
    return pk, lp


def _read_peaks(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: str, 3: str, 10: str})
    if df.shape[1] != len(PEAKS_COLS):
        raise ValueError(
            f"Unexpected peaks columns: got {df.shape[1]} cols, expected {len(PEAKS_COLS)}"
        )
    df.columns = PEAKS_COLS
    return df


def _read_loops(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: str, 3: str, 14: str})
    if df.shape[1] < len(LOOPS_COLS):
        # Be permissive if extra trailing cols ever appear; take first 15
        raise ValueError(
            f"Unexpected loops columns: got {df.shape[1]} cols, expected >= {len(LOOPS_COLS)}"
        )
    df = df.iloc[:, : len(LOOPS_COLS)]
    df.columns = LOOPS_COLS
    return df


def compute_summary(params: SummaryParams) -> pd.DataFrame:
    console = Console()
    prefix = params.prefix or _detect_prefix(params.out)
    peaks_path, loops_path = _paths(params.out, prefix)

    peaks = _read_peaks(peaks_path)
    loops = _read_loops(loops_path)

    # Totals
    peaks_total = int(len(peaks))
    loops_total = int(len(loops))

    # Significant counts under thresholds
    peaks_signif = int(
        (
            (peaks["fdr"].astype(float) <= params.fdr)
            & (peaks["total"].astype(int) >= params.min_reads_peak)
        ).sum()
    )
    loops_signif = int(
        (
            (loops["fdr_pairs"].astype(float) <= params.fdr)
            & (loops["total_pairs"].astype(int) >= params.min_pairs_loop)
        ).sum()
    )

    # Medians (match earlier QC)
    peaks_total_reads_median = float(peaks["total"].astype(int).median())
    loops_total_pairs_median = float(loops["total_pairs"].astype(int).median())

    # Calls distribution
    pk_calls = peaks["bias_call"].value_counts()
    lp_calls = loops["bias_call"].value_counts()

    def _get(d: pd.Series, key: str) -> int:
        try:
            return int(d.get(key, 0))
        except Exception:
            return 0

    summary_rows = [
        ("peaks_total", peaks_total),
        ("loops_total", loops_total),
        ("peaks_signif", peaks_signif),
        ("loops_signif", loops_signif),
        ("peaks_total_reads_median", f"{peaks_total_reads_median:.3f}"),
        ("loops_total_pairs_median", f"{loops_total_pairs_median:.3f}"),
        ("peaks_calls_Maternal", _get(pk_calls, "Maternal")),
        ("peaks_calls_Paternal", _get(pk_calls, "Paternal")),
        ("peaks_calls_Balanced", _get(pk_calls, "Balanced")),
        ("peaks_calls_Undetermined", _get(pk_calls, "Undetermined")),
        ("loops_calls_Maternal", _get(lp_calls, "Maternal")),
        ("loops_calls_Paternal", _get(lp_calls, "Paternal")),
        ("loops_calls_Balanced", _get(lp_calls, "Balanced")),
        ("loops_calls_Undetermined", _get(lp_calls, "Undetermined")),
        ("fdr_threshold", params.fdr),
        ("min_reads_peak", params.min_reads_peak),
        ("min_pairs_loop", params.min_pairs_loop),
    ]

    df_summary = pd.DataFrame(summary_rows, columns=["metric", "value"])

    # Pretty print block
    console.print("[bold]QC SUMMARY[/bold]")
    for m, v in summary_rows:
        console.print(f"{m}\t{v}")

    # Optional TSV
    if params.write_tsv:
        # Write a deterministic quick QC file name to avoid confusion with legacy naming
        tsv_path = params.out / f"{prefix}.quick_qc.tsv"
        df_summary.to_csv(tsv_path, sep="\t", index=False)

    return df_summary
