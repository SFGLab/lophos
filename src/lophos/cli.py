from pathlib import Path
import typer
from rich.console import Console

from .core.calls import BiasThresholds
from .io.config import load_yaml_if_exists

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()

@app.callback()
def main():
    """LOPHOS — Allele-specific phasing of CTCF peaks & loops from phased HiChIP BAMs."""

@app.command("phase")
def phase(
    bam: Path = typer.Option(..., exists=True, help="Phased HiChIP BAM with RG tags"),
    peaks: Path = typer.Option(..., exists=True, help="CTCF peaks (BED)"),
    loops: Path = typer.Option(..., exists=True, help="Loops (BEDPE)"),
    out: Path = typer.Option(..., help="Output prefix (directory will be created)"),
    mapq: int = typer.Option(30, help="Minimum MAPQ to count"),
    peak_window: int = typer.Option(500, help="Peak summit +/- bp window"),
    anchor_pad: int = typer.Option(10_000, help="Anchor padding (bp) when matching mates"),
    min_reads_peak: int = typer.Option(5, help="Min total reads to call a peak"),
    min_pairs_loop: int = typer.Option(3, help="Min informative pairs to call a loop"),
    fdr: float = typer.Option(0.05, help="BH-FDR threshold"),
    keep_duplicates: bool = typer.Option(False, help="Keep PCR/optical duplicates"),
    validate_loops: str = typer.Option("local", help="{none,local} (advanced later)"),
    config: Path = typer.Option(None, help="YAML config to override/record params"),
):
    """
    Phase CTCF peaks and loops using a haplotype-tagged BAM.
    """
    from .io import bam as bam_io, bed as bed_io, bedpe as bedpe_io
    from .core import counts_peaks, counts_loops, calls, stats
    from .report import writers, qc

    if config:
        cfg = load_yaml_if_exists(config)
        console.log(f"Loaded config overrides: {cfg}")

    out.parent.mkdir(parents=True, exist_ok=True)
    console.rule("[bold]LOPHOS phasing")

    # Load inputs
    peaks_df = bed_io.read_bed(peaks)
    loops_df = bedpe_io.read_bedpe(loops)
    bam_handle = bam_io.open_bam(bam)

    # Count & stats: peaks
    peak_counts = counts_peaks.count_peaks(
        bam=bam_handle, peaks=peaks_df, mapq=mapq, window_bp=peak_window, keep_dups=keep_duplicates
    )
    peak_stats = stats.compute_peak_stats(peak_counts)
    peak_calls = calls.call_bias_for_peaks(
        peak_stats, thresholds=BiasThresholds(min_reads=min_reads_peak, fdr=fdr)
    )

    # Count & stats: loops
    loop_counts = counts_loops.count_loops(
        bam=bam_handle, loops=loops_df, mapq=mapq, anchor_pad=anchor_pad, keep_dups=keep_duplicates
    )
    loop_stats = stats.compute_loop_stats(loop_counts)
    loop_calls = calls.call_bias_for_loops(
        loop_stats, thresholds=BiasThresholds(min_reads=min_pairs_loop, fdr=fdr)
    )

    # Validation (local) — placeholder (works on filtered BAMs)
    if validate_loops == "local":
        from .core.validate_local import run_local_validation
        loop_calls = run_local_validation(bam_handle, loops_df, loop_calls, anchor_pad=anchor_pad, mapq=mapq)

    # Write outputs
    writers.write_peaks(out.with_suffix(".peaks.bed"), peak_calls)
    writers.write_loops(out.with_suffix(".loops.bedpe"), loop_calls)
    qc.write_summary(out.with_suffix(".summary.tsv"), peak_calls, loop_calls)
    console.log(f"Done. Outputs: {out}.peaks.bed, {out}.loops.bedpe, {out}.summary.tsv")
