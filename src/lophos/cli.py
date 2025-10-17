import math
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated, Any, cast

import pandas as pd
import typer
from rich.console import Console

from .core import calls, counts_loops, counts_peaks, stats
from .core.calls import BiasThresholds
from .io import bam as bam_io
from .io import bed as bed_io
from .io import bedpe as bedpe_io
from .io.config import load_yaml_if_exists
from .report import qc, writers
from .report.summary import SummaryParams, compute_summary

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()

# ----------------------------------------------------------------------
# Concurrency helpers
#
# To support multi-threaded counting of peaks and loops, we provide
# helper functions that partition a DataFrame into chunks and process
# each chunk in a separate thread.  Each worker opens its own BAM
# handle to avoid thread-safety issues in pysam.  Results are
# assembled in the original order of the chunks to preserve
# deterministic output.


def _count_worker_generic(
    bam_path: str,
    subdf: pd.DataFrame,
    func: Callable[..., pd.DataFrame],
    param_name: str,
    kwargs: dict[str, Any],
    idx: int,
) -> tuple[int, pd.DataFrame]:
    """Internal helper to execute counting on a sub-dataframe.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.  A new handle will be opened inside
        this function.
    subdf : pandas.DataFrame
        Subset of peaks or loops to process.
    func : callable
        Counting function (e.g., counts_peaks.count_peaks or
        counts_loops.count_loops).
    param_name : str
        Name of the argument expected by ``func`` for the data
        dataframe ("peaks" or "loops").
    kwargs : dict
        Additional keyword arguments passed to ``func`` (e.g., mapq,
        window_bp, anchor_pad, keep_dups).
    idx : int
        Chunk index used to preserve ordering when results are
        reassembled.

    Returns
    -------
    tuple[int, pandas.DataFrame]
        Index and resulting counts DataFrame.
    """
    # Open a fresh BAM handle for this worker
    handle = bam_io.open_bam(bam_path)
    try:
        # Build kwargs including the sub-dataframe under the correct
        # parameter name (peaks or loops)
        call_kwargs = {param_name: subdf}
        call_kwargs.update(kwargs)
        result_df = func(bam=handle, **call_kwargs)
        return idx, result_df
    finally:
        handle.close()


def _count_parallel(
    bam_path: str,
    df: pd.DataFrame,
    func: Callable[..., pd.DataFrame],
    param_name: str,
    threads: int,
    **kwargs: Any,
) -> pd.DataFrame:
    """Partition the input DataFrame and perform counting in parallel.

    If ``threads`` is 1 or the dataframe is empty, this function
    performs the counting sequentially.  Otherwise it splits the
    dataframe into roughly equal chunks and processes each chunk in a
    separate thread.  Each thread opens its own BAM handle via
    ``bam_io.open_bam`` and calls the provided counting function.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file used for counting.
    df : pandas.DataFrame
        Dataframe of peaks or loops to process.
    func : callable
        Counting function (counts_peaks.count_peaks or
        counts_loops.count_loops).
    param_name : str
        Name of the argument in ``func`` that should receive the
        dataframe ("peaks" or "loops").
    threads : int
        Number of worker threads.  If <= 1, sequential execution is
        used.
    **kwargs
        Additional keyword arguments to pass to ``func`` (e.g., mapq,
        window_bp, anchor_pad, keep_dups).

    Returns
    -------
    pandas.DataFrame
        Concatenated results from all chunks, assembled in the
        original chunk order.
    """
    # Sequential path or no work
    if threads <= 1 or len(df) == 0:
        handle = bam_io.open_bam(bam_path)
        try:
            call_kwargs = {param_name: df}
            call_kwargs.update(kwargs)
            result_df = func(bam=handle, **call_kwargs)
            return result_df
        finally:
            handle.close()

    # Determine chunk size.  Use math.ceil to ensure non-empty chunks
    chunk_size = math.ceil(len(df) / threads)
    results: dict[int, Any] = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i in range(threads):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, len(df))
            if start >= end:
                continue
            # Preserve original index to maintain correct IDs (peak_id / loop_id)
            subdf = df.iloc[start:end].copy()
            future = executor.submit(
                _count_worker_generic,
                bam_path,
                subdf,
                func,
                param_name,
                kwargs,
                i,
            )
            futures.append(future)
        for fut in futures:
            idx, part_df = fut.result()
            results[idx] = part_df
    # Concatenate results in order of chunk index to preserve order
    # Import pandas locally to avoid a global dependency for CLI consumers
    return cast(
        pd.DataFrame, pd.concat([results[i] for i in sorted(results.keys())], ignore_index=True)
    )


@app.callback()
def main() -> None:
    """LOPHOS — Allele-specific phasing of CTCF peaks & loops from phased HiChIP BAMs."""


@app.command("phase")
def phase(  # noqa: C901
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
    # New options
    maternal_rgid: Annotated[
        str, typer.Option(help="Regex for maternal RG identifiers")
    ] = "maternal|mat|M",
    paternal_rgid: Annotated[
        str, typer.Option(help="Regex for paternal RG identifiers")
    ] = "paternal|pat|P",
    pseudocount: Annotated[
        float, typer.Option(help="Pseudocount added to M/P in log2 ratio")
    ] = 1.0,
    min_abs_log2: Annotated[
        float, typer.Option(help="Minimum absolute log2 ratio for bias calling")
    ] = 0.0,
    max_ambiguous_frac: Annotated[
        float, typer.Option(help="Max ambiguous fraction to call loops (default 0.5)")
    ] = 0.5,
    primary_only: Annotated[
        bool, typer.Option(help="Filter out ALT/decoy/Un contigs and output .primary.* files")
    ] = False,
    summary: Annotated[bool, typer.Option(help="Run QC summary after phasing")] = False,
    threads: Annotated[
        int, typer.Option(min=1, help="Number of threads for counting (future use)")
    ] = 1,
    log_level: Annotated[
        str, typer.Option(help="Logging level (info, debug, warning, error)")
    ] = "info",
    config: Annotated[
        Path | None, typer.Option(help="YAML config to override/record params")
    ] = None,
) -> None:
    """Phase CTCF peaks and loops using a haplotype-tagged BAM.

    This command counts maternal and paternal reads/pairs supporting peaks and loops,
    performs binomial tests with FDR correction, applies bias calling thresholds
    and writes out bed/bedpe files summarising the results.  Optional features
    include custom RG tag regexes, effect-size thresholds, ambiguous pair
    handling, primary-chromosome filtering and integrated QC summarisation.
    """

    # -----------------------------------------------------------------------
    # 1. Load configuration and resolve parameters
    # -----------------------------------------------------------------------
    defaults: dict[str, Any] = {
        "mapq": 30,
        "peak_window": 500,
        "anchor_pad": 10_000,
        "min_reads_peak": 5,
        "min_pairs_loop": 3,
        "fdr": 0.05,
        "keep_duplicates": False,
        "validate_loops": "local",
        "maternal_rgid": "maternal|mat|M",
        "paternal_rgid": "paternal|pat|P",
        "pseudocount": 1.0,
        "min_abs_log2": 0.0,
        "max_ambiguous_frac": 0.5,
        "primary_only": False,
        "summary": False,
        "threads": 1,
        "log_level": "info",
    }

    def load_cfg(path: Path | None) -> dict[str, Any]:
        if not path:
            return {}
        cfg_data = load_yaml_if_exists(path)
        if cfg_data:
            console.log(f"Loaded config overrides: {cfg_data}")
        return cfg_data or {}

    def resolve_params(cli_vals: dict[str, Any], cfg_vals: dict[str, Any]) -> dict[str, Any]:
        """Return final parameter values. CLI overrides config except
        when the CLI value matches its default and config supplies an override."""
        resolved: dict[str, Any] = {}
        for name, default_val in defaults.items():
            cli_val = cli_vals.get(name, default_val)
            if name in cfg_vals and cli_val == default_val:
                resolved[name] = cfg_vals[name]
            else:
                resolved[name] = cli_val
        return resolved

    # Build dictionaries of CLI values and config values
    cli_dict: dict[str, Any] = {
        "mapq": mapq,
        "peak_window": peak_window,
        "anchor_pad": anchor_pad,
        "min_reads_peak": min_reads_peak,
        "min_pairs_loop": min_pairs_loop,
        "fdr": fdr,
        "keep_duplicates": keep_duplicates,
        "validate_loops": validate_loops,
        "maternal_rgid": maternal_rgid,
        "paternal_rgid": paternal_rgid,
        "pseudocount": pseudocount,
        "min_abs_log2": min_abs_log2,
        "max_ambiguous_frac": max_ambiguous_frac,
        "primary_only": primary_only,
        "summary": summary,
        "threads": threads,
        "log_level": log_level,
    }
    cfg_dict = load_cfg(config)
    params = resolve_params(cli_dict, cfg_dict)

    # -----------------------------------------------------------------------
    # 2. Update module-level constants and RG patterns
    # -----------------------------------------------------------------------
    from . import constants

    constants.PSEUDOCOUNT = float(params["pseudocount"])
    bam_io.set_rg_patterns(
        params["maternal_rgid"] if params["maternal_rgid"] else None,
        params["paternal_rgid"] if params["paternal_rgid"] else None,
    )

    # -----------------------------------------------------------------------
    # 3. Create output directory and announce start
    # -----------------------------------------------------------------------
    out.parent.mkdir(parents=True, exist_ok=True)
    console.rule("[bold]LOPHOS phasing")

    # -----------------------------------------------------------------------
    # 4. Read inputs
    # -----------------------------------------------------------------------
    peaks_df_full = bed_io.read_bed(peaks)
    loops_df_full = bedpe_io.read_bedpe(loops)

    # -----------------------------------------------------------------------
    # 5. Perform counts, stats, calls and optional validation in a helper
    # -----------------------------------------------------------------------
    def perform_phasing() -> tuple[Any, Any]:
        """Count reads/pairs, compute statistics, call biases and run validation."""
        threads_param = int(params["threads"])
        bam_path_str = str(bam)
        # If threads <= 1, fall back to sequential counting using a single BAM handle
        if threads_param <= 1:
            bam_handle = bam_io.open_bam(bam)
            try:
                peak_counts = counts_peaks.count_peaks(
                    bam=bam_handle,
                    peaks=peaks_df_full,
                    mapq=int(params["mapq"]),
                    window_bp=int(params["peak_window"]),
                    keep_dups=bool(params["keep_duplicates"]),
                )
                peak_stats = stats.compute_peak_stats(peak_counts)
                peak_calls = calls.call_bias_for_peaks(
                    peak_stats,
                    thresholds=BiasThresholds(
                        min_reads=int(params["min_reads_peak"]),
                        fdr=float(params["fdr"]),
                        min_fold=1.5,
                        min_abs_log2=float(params["min_abs_log2"]),
                    ),
                )
                loop_counts = counts_loops.count_loops(
                    bam=bam_handle,
                    loops=loops_df_full,
                    mapq=int(params["mapq"]),
                    anchor_pad=int(params["anchor_pad"]),
                    keep_dups=bool(params["keep_duplicates"]),
                )
                loop_stats = stats.compute_loop_stats(loop_counts)
                loop_calls = calls.call_bias_for_loops(
                    loop_stats,
                    thresholds=BiasThresholds(
                        min_reads=int(params["min_pairs_loop"]),
                        fdr=float(params["fdr"]),
                        min_fold=1.5,
                        min_abs_log2=float(params["min_abs_log2"]),
                        max_ambiguous_frac=float(params["max_ambiguous_frac"]),
                    ),
                )
            finally:
                bam_handle.close()
        else:
            # Multi-threaded counting: partition peaks and loops and process in parallel
            # Peaks counts
            peak_counts = _count_parallel(
                bam_path_str,
                peaks_df_full,
                counts_peaks.count_peaks,
                "peaks",
                threads_param,
                mapq=int(params["mapq"]),
                window_bp=int(params["peak_window"]),
                keep_dups=bool(params["keep_duplicates"]),
            )
            peak_stats = stats.compute_peak_stats(peak_counts)
            peak_calls = calls.call_bias_for_peaks(
                peak_stats,
                thresholds=BiasThresholds(
                    min_reads=int(params["min_reads_peak"]),
                    fdr=float(params["fdr"]),
                    min_fold=1.5,
                    min_abs_log2=float(params["min_abs_log2"]),
                ),
            )
            # Loops counts
            loop_counts = _count_parallel(
                bam_path_str,
                loops_df_full,
                counts_loops.count_loops,
                "loops",
                threads_param,
                mapq=int(params["mapq"]),
                anchor_pad=int(params["anchor_pad"]),
                keep_dups=bool(params["keep_duplicates"]),
            )
            loop_stats = stats.compute_loop_stats(loop_counts)
            loop_calls = calls.call_bias_for_loops(
                loop_stats,
                thresholds=BiasThresholds(
                    min_reads=int(params["min_pairs_loop"]),
                    fdr=float(params["fdr"]),
                    min_fold=1.5,
                    min_abs_log2=float(params["min_abs_log2"]),
                    max_ambiguous_frac=float(params["max_ambiguous_frac"]),
                ),
            )
            # Local validation if requested; open a new BAM handle for validation
            if str(params["validate_loops"]) == "local":
                from .core.validate_local import run_local_validation

                # open BAM handle for validation (use one handle; validation is sequential)
                bam_handle_val = bam_io.open_bam(bam)
                try:
                    loop_calls = run_local_validation(
                        bam_handle_val,
                        loops_df_full,
                        loop_calls,
                        int(params["anchor_pad"]),
                        int(params["mapq"]),
                    )
                finally:
                    bam_handle_val.close()
        return peak_calls, loop_calls

    peak_calls_full, loop_calls_full = perform_phasing()

    # -----------------------------------------------------------------------
    # 6. Write outputs for full data
    # -----------------------------------------------------------------------
    writers.write_peaks(out.with_suffix(".peaks.bed"), peak_calls_full)
    writers.write_loops(out.with_suffix(".loops.bedpe"), loop_calls_full)
    qc.write_summary(out.with_suffix(".summary.tsv"), peak_calls_full, loop_calls_full)

    # -----------------------------------------------------------------------
    # 7. Optionally write primary-only outputs
    # -----------------------------------------------------------------------
    if bool(params["primary_only"]):
        import re

        pattern = re.compile(r"^(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|x|y|m)$", re.IGNORECASE)
        exclude_substrings = ["_", "alt", "decoy", "random", "un"]

        def _is_primary(chrom: str) -> bool:
            chrom_l = chrom.lower()
            if any(sub in chrom_l for sub in exclude_substrings):
                return False
            return bool(pattern.match(chrom_l))

        def _is_primary_peak(row: Any) -> bool:
            return _is_primary(str(row["chrom"]))

        def _is_primary_loop(row: Any) -> bool:
            return _is_primary(str(row["chrom1"])) and _is_primary(str(row["chrom2"]))

        peak_calls_primary = peak_calls_full[
            peak_calls_full.apply(_is_primary_peak, axis=1)
        ].reset_index(drop=True)
        loop_calls_primary = loop_calls_full[
            loop_calls_full.apply(_is_primary_loop, axis=1)
        ].reset_index(drop=True)
        writers.write_peaks(out.with_suffix(".primary.peaks.bed"), peak_calls_primary)
        writers.write_loops(out.with_suffix(".primary.loops.bedpe"), loop_calls_primary)

    # -----------------------------------------------------------------------
    # 8. Persist resolved configuration
    # -----------------------------------------------------------------------
    try:
        import json

        run_json_path = out.with_suffix(".run.json")
        with run_json_path.open("w") as fh:
            json.dump(params, fh, indent=2)
        console.log(f"Saved run configuration to {run_json_path}")
    except Exception as exc:
        console.log(f"[red]WARNING:[/red] Failed to write run configuration: {exc}")

    # -----------------------------------------------------------------------
    # 9. Integrated summary (quick QC)
    # -----------------------------------------------------------------------
    if bool(params["summary"]):
        try:
            params_summary = SummaryParams(
                out=out.parent,
                prefix=out.name,
                fdr=float(params["fdr"]),
                min_reads_peak=int(params["min_reads_peak"]),
                min_pairs_loop=int(params["min_pairs_loop"]),
                write_tsv=True,
            )
            compute_summary(params_summary)
        except Exception as e:  # noqa: BLE001
            console.print(f"[red]ERROR computing summary:[/red] {e}")

    console.log(f"Done. Outputs: {out}.peaks.bed, {out}.loops.bedpe, {out}.summary.tsv")


@app.command("summary")
def summary(
    out: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            help="LOPHOS run output directory",
        ),
    ],
    prefix: Annotated[
        str | None,
        typer.Option(
            help="Prefix (e.g., SAMPLE if files are SAMPLE.peaks.bed & SAMPLE.loops.bedpe)"
        ),
    ] = None,
    fdr: Annotated[
        float, typer.Option(min=0.0, help="FDR threshold for 'significant' counts")
    ] = 0.05,
    min_reads_peak: Annotated[
        int, typer.Option(min=0, help="Min (M+P) reads for peaks to count as measurable")
    ] = 5,
    min_pairs_loop: Annotated[
        int, typer.Option(min=0, help="Min (M+P) pairs for loops to count as measurable")
    ] = 3,
    no_tsv: Annotated[
        bool, typer.Option("--no-tsv", help="Do not write <prefix>.qc_summary.tsv")
    ] = False,
) -> None:
    """
    Summarize an existing LOPHOS run (totals, significant features, medians, call breakdown).

    Examples:
      lophos summary --out results/SAMPLE
      lophos summary --out results/SAMPLE --prefix SAMPLE --fdr 0.05
    """
    params = SummaryParams(
        out=out,
        prefix=prefix,
        fdr=fdr,
        min_reads_peak=min_reads_peak,
        min_pairs_loop=min_pairs_loop,
        write_tsv=(not no_tsv),
    )
    try:
        compute_summary(params)
    except Exception as e:  # noqa: BLE001
        console.print(f"[red]ERROR:[/red] {e}")
        # Chain the original exception to satisfy ruff B904
        raise typer.Exit(code=1) from e
