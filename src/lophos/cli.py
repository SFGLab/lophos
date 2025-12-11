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
# ----------------------------------------------------------------------


def _count_worker_generic(
    bam_path: str,
    subdf: pd.DataFrame,
    func: Callable[..., pd.DataFrame],
    param_name: str,
    kwargs: dict[str, Any],
    idx: int,
) -> tuple[int, pd.DataFrame]:
    """Internal helper to execute counting on a sub-dataframe."""
    handle = bam_io.open_bam(bam_path)
    try:
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
    """Partition the input DataFrame and perform counting in parallel."""
    if threads <= 1 or len(df) == 0:
        handle = bam_io.open_bam(bam_path)
        try:
            call_kwargs = {param_name: df}
            call_kwargs.update(kwargs)
            result_df = func(bam=handle, **call_kwargs)
            return result_df
        finally:
            handle.close()

    chunk_size = math.ceil(len(df) / threads)
    results: dict[int, Any] = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i in range(threads):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, len(df))
            if start >= end:
                continue
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
    return cast(
        pd.DataFrame, pd.concat([results[i] for i in sorted(results.keys())], ignore_index=True)
    )


@app.callback()
def main() -> None:
    """LOPHOS — Allele-specific phasing of CTCF peaks & loops from haplotype-tagged BAMs."""


@app.command("phase")
def phase(  # noqa: C901
    bam: Annotated[Path, typer.Option(exists=True, help="Haplotype-tagged BAM (RG labels)")],
    peaks: Annotated[Path, typer.Option(exists=True, help="CTCF peaks (BED)")],
    loops: Annotated[Path, typer.Option(exists=True, help="Loops (BEDPE)")],
    out: Annotated[Path, typer.Option(help="Output prefix (directory will be created)")],
    mapq: Annotated[int, typer.Option(help="Minimum MAPQ to count")] = 30,
    peak_window: Annotated[int, typer.Option(help="Peak summit +/- bp window")] = 500,
    anchor_pad: Annotated[int, typer.Option(help="Anchor padding (bp)")] = 10_000,
    min_reads_peak: Annotated[int, typer.Option(help="Min total reads to call a peak")] = 5,
    min_pairs_loop: Annotated[int, typer.Option(help="Min informative pairs to call a loop")] = 3,
    fdr: Annotated[float, typer.Option(help="BH-FDR threshold")] = 0.05,
    keep_duplicates: Annotated[bool, typer.Option(help="Keep duplicates")] = False,
    validate_loops: Annotated[str, typer.Option(help="{none,local}")] = "local",
    # RG mapping
    maternal_rgid: Annotated[
        str, typer.Option(help="Regex for maternal RG identifiers")
    ] = "maternal|mat|M",
    paternal_rgid: Annotated[
        str, typer.Option(help="Regex for paternal RG identifiers")
    ] = "paternal|pat|P",
    # Effect-size controls
    pseudocount: Annotated[float, typer.Option(help="Pseudocount for log2 ratio")] = 1.0,
    min_abs_log2: Annotated[float, typer.Option(help="Min |log2| for bias calling")] = 0.0,
    max_ambiguous_frac: Annotated[float, typer.Option(help="Max ambiguous fraction (loops)")] = 0.5,
    # Loop mode (new)
    loop_mode: Annotated[
        str,
        typer.Option(
            help="Loop pairing mode: 'mates' (paired-end) or 'sa' (SA:Z chimeric long-read)"
        ),
    ] = "mates",
    # SA:Z mode knobs (used when loop_mode='sa')
    sa_min_mapq: Annotated[int, typer.Option(help="[sa] Min MAPQ per segment")] = 30,
    sa_min_seg_len: Annotated[int, typer.Option(help="[sa] Min aligned segment length (bp)")] = 50,
    sa_min_cis_dist: Annotated[int, typer.Option(help="[sa] Min cis distance to keep (bp)")] = 1000,
    sa_allow_trans: Annotated[bool, typer.Option(help="[sa] Allow trans-chrom contacts")] = True,
    sa_orientation: Annotated[
        str,
        typer.Option(help="[sa] Orientation policy: 'any' or 'convergent-short-cis'"),
    ] = "any",
    sa_dedup_within_read: Annotated[
        bool, typer.Option(help="[sa] Dedup same contact within read")
    ] = True,
    # Misc
    primary_only: Annotated[
        bool, typer.Option(help="Write .primary.* (primary chroms only)")
    ] = False,
    summary: Annotated[bool, typer.Option(help="Run QC summary after phasing")] = False,
    threads: Annotated[int, typer.Option(min=1, help="Threads for counting")] = 1,
    log_level: Annotated[str, typer.Option(help="Logging: info, debug, warning, error")] = "info",
    config: Annotated[
        Path | None, typer.Option(help="YAML config to override/record params")
    ] = None,
) -> None:
    """Phase CTCF peaks and loops using a haplotype-tagged BAM."""
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
        "loop_mode": "mates",
        "sa_min_mapq": 30,
        "sa_min_seg_len": 50,
        "sa_min_cis_dist": 1000,
        "sa_allow_trans": True,
        "sa_orientation": "any",
        "sa_dedup_within_read": True,
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
        """CLI overrides config except when CLI equals default and config provides a value."""
        resolved: dict[str, Any] = {}
        for name, default_val in defaults.items():
            cli_val = cli_vals.get(name, default_val)
            if name in cfg_vals and cli_val == default_val:
                resolved[name] = cfg_vals[name]
            else:
                resolved[name] = cli_val
        return resolved

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
        "loop_mode": loop_mode,
        "sa_min_mapq": sa_min_mapq,
        "sa_min_seg_len": sa_min_seg_len,
        "sa_min_cis_dist": sa_min_cis_dist,
        "sa_allow_trans": sa_allow_trans,
        "sa_orientation": sa_orientation,
        "sa_dedup_within_read": sa_dedup_within_read,
        "primary_only": primary_only,
        "summary": summary,
        "threads": threads,
        "log_level": log_level,
    }
    cfg_dict = load_cfg(config)
    params = resolve_params(cli_dict, cfg_dict)

    # -----------------------------------------------------------------------
    # 2. Update constants and RG patterns
    # -----------------------------------------------------------------------
    from . import constants

    constants.PSEUDOCOUNT = float(params["pseudocount"])
    bam_io.set_rg_patterns(
        params["maternal_rgid"] if params["maternal_rgid"] else None,
        params["paternal_rgid"] if params["paternal_rgid"] else None,
    )

    # -----------------------------------------------------------------------
    # 3. Output dir & banner
    # -----------------------------------------------------------------------
    out.parent.mkdir(parents=True, exist_ok=True)
    console.rule("[bold]LOPHOS phasing")

    # -----------------------------------------------------------------------
    # 4. Read inputs
    # -----------------------------------------------------------------------
    peaks_df_full = bed_io.read_bed(peaks)
    loops_df_full = bedpe_io.read_bedpe(loops)

    # -----------------------------------------------------------------------
    # 5. Perform counts, stats, calls (validation optional)
    # -----------------------------------------------------------------------
    def perform_phasing() -> tuple[Any, Any]:
        threads_param = int(params["threads"])
        bam_path_str = str(bam)

        loop_kwargs_common: dict[str, int | str | bool] = {
            "mapq": int(params["mapq"]),
            "anchor_pad": int(params["anchor_pad"]),
            "keep_dups": bool(params["keep_duplicates"]),
            # new dispatch + SA params (counts_loops should accept these kwargs)
            "loop_mode": str(params["loop_mode"]),
            "sa_min_mapq": int(params["sa_min_mapq"]),
            "sa_min_seg_len": int(params["sa_min_seg_len"]),
            "sa_min_cis_dist": int(params["sa_min_cis_dist"]),
            "sa_allow_trans": bool(params["sa_allow_trans"]),
            "sa_orientation": str(params["sa_orientation"]),
            "sa_dedup_within_read": bool(params["sa_dedup_within_read"]),
        }

        if threads_param <= 1:
            bam_handle = bam_io.open_bam(bam)
            try:
                # Peaks
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
                # Loops
                loop_counts = counts_loops.count_loops(
                    bam=bam_handle,
                    loops=loops_df_full,
                    mapq=int(params["mapq"]),
                    anchor_pad=int(params["anchor_pad"]),
                    keep_dups=bool(params["keep_duplicates"]),
                    loop_mode=str(params["loop_mode"]),
                    sa_min_mapq=int(params["sa_min_mapq"]),
                    sa_min_seg_len=int(params["sa_min_seg_len"]),
                    sa_min_cis_dist=int(params["sa_min_cis_dist"]),
                    sa_allow_trans=bool(params["sa_allow_trans"]),
                    sa_orientation=str(params["sa_orientation"]),
                    sa_dedup_within_read=bool(params["sa_dedup_within_read"]),
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
            # Peaks (parallel)
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
            # Loops (parallel)
            loop_counts = _count_parallel(
                bam_path_str,
                loops_df_full,
                counts_loops.count_loops,
                "loops",
                threads_param,
                **loop_kwargs_common,
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
            # Local validation (optional)
            if str(params["validate_loops"]) == "local":
                from .core.validate_local import run_local_validation

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
    # 6. Write outputs
    # -----------------------------------------------------------------------
    writers.write_peaks(out.with_suffix(".peaks.bed"), peak_calls_full)
    writers.write_loops(out.with_suffix(".loops.bedpe"), loop_calls_full)
    qc.write_summary(out.with_suffix(".summary.tsv"), peak_calls_full, loop_calls_full)

    # -----------------------------------------------------------------------
    # 7. Primary-only outputs
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
    min_reads_peak: Annotated[int, typer.Option(min=0, help="Min (M+P) reads for peaks")] = 5,
    min_pairs_loop: Annotated[int, typer.Option(min=0, help="Min (M+P) pairs for loops")] = 3,
    no_tsv: Annotated[
        bool, typer.Option("--no-tsv", help="Do not write <prefix>.qc_summary.tsv")
    ] = False,
) -> None:
    """Summarize an existing LOPHOS run (totals, significant features, medians, call breakdown)."""
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
        raise typer.Exit(code=1) from e
