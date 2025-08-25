from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from .core import calls, counts_loops, counts_peaks, stats
from .core.calls import BiasThresholds
from .io import bam as bam_io
from .io import bed as bed_io
from .io import bedpe as bedpe_io
from .io.config import load_yaml_if_exists
from .report import qc, writers

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()


@app.callback()
def main() -> None:
    """LOPHOS — Allele-specific phasing of CTCF peaks & loops from phased HiChIP BAMs."""


@app.command("phase")
def phase(
    bam: Annotated[Path, typer.Option(exists=True, help="Phased HiChIP BAM with RG tags")],
    peaks: Annotated[Path, typer.Option(exists=True, help="CTCF peaks (BED)")],
    loops: Annotated[Path, typer.Option(exists=True, help="Loops (BEDPE)")],
    out: Annotated[Path, typer.Option(help="Output prefix (directory will be created)")],
    mapq: Annotated[int, typer.Option(help="Minimum MAPQ to count")] = 30,
    peak_window: Annotated[int, typer.Option(help="Peak summit +/- bp window")] = 500,
    anchor_pad: Annotated[
        int, typer.Option(help="Anchor padding (bp) when matching mates")
    ] = 10_000,
    min_reads_peak: Annotated[int, typer.Option(help="Min total reads to call a peak")] = 5,
    min_pairs_loop: Annotated[int, typer.Option(help="Min informative pairs to call a loop")] = 3,
    fdr: Annotated[float, typer.Option(help="BH-FDR threshold")] = 0.05,
    keep_duplicates: Annotated[bool, typer.Option(help="Keep PCR/optical duplicates")] = False,
    validate_loops: Annotated[str, typer.Option(help="{none,local} (advanced later)")] = "local",
    config: Annotated[
        Path | None, typer.Option(help="YAML config to override/record params")
    ] = None,
) -> None:
    """Phase CTCF peaks and loops using a haplotype-tagged BAM."""
    if config:
        cfg = load_yaml_if_exists(config)
        console.log(f"Loaded config overrides: {cfg}")

    out.parent.mkdir(parents=True, exist_ok=True)
    console.rule("[bold]LOPHOS phasing")

    peaks_df = bed_io.read_bed(peaks)
    loops_df = bedpe_io.read_bedpe(loops)
    bam_handle = bam_io.open_bam(bam)

    peak_counts = counts_peaks.count_peaks(
        bam=bam_handle, peaks=peaks_df, mapq=mapq, window_bp=peak_window, keep_dups=keep_duplicates
    )
    peak_stats = stats.compute_peak_stats(peak_counts)
    peak_calls = calls.call_bias_for_peaks(
        peak_stats, thresholds=BiasThresholds(min_reads=min_reads_peak, fdr=fdr)
    )

    loop_counts = counts_loops.count_loops(
        bam=bam_handle, loops=loops_df, mapq=mapq, anchor_pad=anchor_pad, keep_dups=keep_duplicates
    )
    loop_stats = stats.compute_loop_stats(loop_counts)
    loop_calls = calls.call_bias_for_loops(
        loop_stats, thresholds=BiasThresholds(min_reads=min_pairs_loop, fdr=fdr)
    )

    if validate_loops == "local":
        from .core.validate_local import run_local_validation

        loop_calls = run_local_validation(bam_handle, loops_df, loop_calls, anchor_pad, mapq)

    writers.write_peaks(out.with_suffix(".peaks.bed"), peak_calls)
    writers.write_loops(out.with_suffix(".loops.bedpe"), loop_calls)
    qc.write_summary(out.with_suffix(".summary.tsv"), peak_calls, loop_calls)
    console.log(f"Done. Outputs: {out}.peaks.bed, {out}.loops.bedpe, {out}.summary.tsv")
