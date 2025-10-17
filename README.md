# LOPHOS — LOops & Peaks HaplOtype phasing Suite

> Allele-specific phasing of CTCF **peaks** and **loops** from haplotype-tagged HiChIP BAMs.

---

## Overview (for biologists)

LOPHOS quantifies **maternal vs paternal** support for:

* **CTCF peaks** (from BED) using haplotype-tagged reads around each peak summit.
* **Loops** (from BEDPE) using paired reads connecting the two anchors.

For each feature, LOPHOS:

1. **Counts** maternal (M) and paternal (P) support.
2. **Tests** allele bias with a two-sided **binomial test** (`p = 0.5`).
3. **Adjusts** p-values using **Benjamini–Hochberg (FDR)**.
4. **Calls** each feature as **Maternal / Paternal / Balanced / Undetermined**.

Outputs are plain **BED/BEDPE** tables you can sort, filter, and load in genome browsers. Other formats can be added if needed.

---

## Installation

```bash
# Python ≥ 3.10
git clone https://github.com/<org-or-user>/lophos.git
cd lophos

# Development install (recommended)
pip install -e ".[dev]"
pre-commit install

# Quick checks
python -m lophos --help
pytest -q
```

---

## Quickstart

### Phase peaks & loops

> **Note:** CLI options use **kebab-case** (hyphens), not underscores.

```bash
lophos phase \
  --bam   /path/to/<sample>.bam \
  --peaks /path/to/<peaks>.bed \
  --loops /path/to/<loops>.bedpe \
  --out   results/<sample> \
  --mapq 30 --peak-window 500 --anchor-pad 10000 \
  --min-reads-peak 5 --min-pairs-loop 3 --fdr 0.05 \
  --maternal-rgid "maternal|mat|M" --paternal-rgid "paternal|pat|P" \
  --pseudocount 1.0 --min-abs-log2 0.0 --max-ambiguous-frac 0.5 \
  --validate-loops local

  # Optional flags:
  # --primary-only true  # drop ALT/decoy/unplaced contigs and write .primary.* outputs
  # --summary true       # run QC summary immediately after phasing
```

This produces:

* `results/<sample>.peaks.bed` — allele calls for peaks (see [Schema](docs/SCHEMA.md)).
* `results/<sample>.loops.bedpe` — allele calls for loops (see [Schema](docs/SCHEMA.md)).
* `results/<sample>.summary.tsv` — legacy per-class counts (used internally).
* If `--summary true`, a quick QC table (`results/<sample>.quick_qc.tsv`) is written.
* If `--primary-only true`, additional filtered files `results/<sample>.primary.peaks.bed` and `results/<sample>.primary.loops.bedpe` are created.

### Summarize an existing run (new)

```bash
lophos summary \
  --out results/<sample> \
  --prefix <sample> \
  --fdr 0.05 \
  --min-reads-peak 5 \
  --min-pairs-loop 3
```

This command reads the `*.peaks.bed` and `*.loops.bedpe` files in `--out`, infers the prefix if not provided, and prints a tidy QC report to the console.  The report includes totals, counts of significant features (at the chosen `--fdr` and minimum coverage thresholds), medians and the breakdown of calls.  A quick QC TSV (`<prefix>.quick_qc.tsv`) is written alongside the inputs unless `--no-tsv` is specified.  You can obtain the same quick QC file directly from the phasing step by using `--summary true` on `lophos phase`.

---

## Input requirements

### BAM (haplotype-tagged)

Reads must carry a read group (**RG**) tag indicating parental origin.  By default, any RG tag matching the regular expression `maternal|mat|M` is considered **maternal**, and any tag matching `paternal|pat|P` is considered **paternal** (case‑insensitive).  You can override these patterns at the command line with `--maternal-rgid` and `--paternal-rgid` to accommodate custom pipelines.

### Peaks (BED)

```
chrom  start  end  [name  score  strand]
```

### Loops (BEDPE)

```
chrom1  start1  end1  chrom2  start2  end2  [name  score  strand1  strand2 ...]
```

> **BEDPE with extra columns** is supported. LOPHOS uses the first 6 BEDPE columns and **ignores any extras**.

---

## Output formats

### Peaks (`*.peaks.bed`)

```
chrom  start  end  peak_id  maternal  paternal  total  log2_ratio  p_value  fdr  bias_call
```

### Loops (`*.loops.bedpe`)

```
chrom1  start1  end1  chrom2  start2  end2  loop_id  maternal_pairs  paternal_pairs  ambiguous_pairs  total_pairs  log2_ratio_pairs  p_value_pairs  fdr_pairs  bias_call
```

### Summary (`*.summary.tsv`)

```
metric  value
peaks_total          <int>
peaks_maternal       <int>
peaks_paternal       <int>
peaks_balanced       <int>
peaks_undetermined   <int>
loops_total          <int>
loops_maternal       <int>
loops_paternal       <int>
loops_balanced       <int>
loops_undetermined   <int>
```

Default outputs are **BED/BEDPE** as requested; TSV/CSV/BigBed or track hubs can be added on request.

---

## CLI options (key)

| Option                  | Purpose                                                               | Typical |
| ----------------------- | --------------------------------------------------------------------- | ------- |
| `--mapq`                | Minimum mapping quality to count reads                                | `30`    |
| `--peak-window`         | Peak summit ± bp window for counting reads                            | `500`   |
| `--anchor-pad`          | Padding (bp) added around each loop anchor for mate matching          | `10000` |
| `--min-reads-peak`      | Minimum total (M + P) reads required to call a peak                   | `5`     |
| `--min-pairs-loop`      | Minimum total (M + P) informative pairs required to call a loop       | `3`     |
| `--fdr`                 | BH‑FDR threshold for significance                                      | `0.05`  |
| `--keep-duplicates`     | Count PCR/optical duplicates (set `True` for HiChIP retaining duplicates) | `False` |
| `--validate-loops`      | Additional loop QC: `none` (skip) or `local` (compute z‑score proxy)    | `local` |
| `--maternal-rgid`       | Regex (case‑insensitive) to detect maternal RG tags                    | `"maternal|mat|M"` |
| `--paternal-rgid`       | Regex to detect paternal RG tags                                        | `"paternal|pat|P"` |
| `--pseudocount`         | Pseudocount added to M and P in log2 ratio calculation                  | `1.0`   |
| `--min-abs-log2`        | Minimum absolute log2 ratio (effect size) required to call bias        | `0.0`   |
| `--max-ambiguous-frac`  | Maximum fraction of ambiguous pairs tolerated when calling loops        | `0.5`   |
| `--primary-only`        | Filter out ALT/decoy/unplaced contigs and write `.primary.*` outputs    | `False` |
| `--summary`             | If `True`, run the QC summary at the end of `phase` and write `.quick_qc.tsv` | `False` |
| `--threads`             | Number of threads for counting (reserved for future parallelism)        | `1`     |
| `--log-level`           | Logging verbosity (`info`, `debug`, `warning`, `error`)                | `"info"` |
| `--config`              | YAML configuration file; values override defaults unless specified on CLI | `None` |

---

## Recommended workflow

1. **Sanity‑check RG tags and set patterns**

   ```bash
   samtools view -H <sample>.bam | grep '^@RG'
   ```

   Confirm that your BAM contains read groups identifying parental origin (e.g. `mat`/`pat`).  If your pipeline uses different identifiers, supply case‑insensitive regular expressions via `--maternal-rgid` and `--paternal-rgid` to match your tags instead of editing the source code.

2. **Try one sample** with defaults (above). If many features are **Undetermined**, relax thresholds:

   * `--min-reads-peak 3` and/or `--min-pairs-loop 2`
   * `--mapq 20`
   * `--fdr 0.10`

3. **Sweep parameters** to find sensible settings for your dataset:

   ```bash
   for MQ in 20 30; do
     for PK in 3 5; do
       for LP in 2 3; do
         OUT=results/sweep/<sample>_mq${MQ}_pk${PK}_lp${LP}
         lophos phase --bam <bam> --peaks <bed> --loops <bedpe> \
           --out "$OUT" --mapq $MQ --peak-window 500 --anchor-pad 10000 \
           --min-reads-peak $PK --min-pairs-loop $LP --fdr 0.05 --validate-loops local
       done
     done
   done
   ```

4. **Batch run** once parameters look good (adjust naming as needed):

   ```bash
   BAM_DIR=/path/to/merged_bams
   PEAK_DIR=/path/to/personalized_peaks
   LOOP_DIR=/path/to/personalized_loops/loops
   OUT_DIR=results/batch; mkdir -p "$OUT_DIR"

   find "$BAM_DIR" -name '*.bam' | while read -r BAM; do
     S=$(basename "$BAM" .bam)
     P="$PEAK_DIR/$S.bed"
     L="$LOOP_DIR/$S.bedpe"
     [[ -f "$P" && -f "$L" ]] || { echo "[skip] $S"; continue; }
     lophos phase --bam "$BAM" --peaks "$P" --loops "$L" \
       --out "$OUT_DIR/$S" --mapq 30 --peak-window 500 --anchor-pad 10000 \
       --min-reads-peak 5 --min-pairs-loop 3 --fdr 0.05 --validate-loops local
   done
   ```

5. **Summarize** any completed run:

   ```bash
   lophos summary \
     --out results/<sample> \
     --prefix <sample> \
     --fdr 0.05 \
     --min-reads-peak 5 \
     --min-pairs-loop 3
   ```

   This command reads the `*.peaks.bed` and `*.loops.bedpe` files in `--out`, infers the prefix if not provided, and prints a tidy QC report to the console.  The report includes totals, counts of significant features (at the chosen `--fdr` and minimum coverage thresholds), medians and the breakdown of calls.  A quick QC TSV (`<prefix>.quick_qc.tsv`) is written alongside the inputs unless `--no-tsv` is specified.

   You can also run `lophos phase ... --summary true` to produce the same QC table immediately after phasing—no separate `summary` invocation needed.

---

## Methods (concise)

* **Counts:** per feature, compute maternal (M) and paternal (P) support.
  Peaks: reads within ±`peak-window` around the summit.
  Loops: paired reads bridging padded anchors (`anchor-pad`).
* **Test:** `M ~ Binomial(M+P, 0.5)` (two-sided).
* **Multiple testing:** BH-FDR over features.
* **Calling:**

  * If `M+P` < minimum → **Undetermined**.
  * For loops, if the fraction of ambiguous pairs (`ambiguous_pairs/(M+P+ambiguous_pairs)`) exceeds `--max-ambiguous-frac`, call **Undetermined** regardless of other criteria.
  * Else if `FDR ≤ α` and both an effect‑size threshold and a fold‑change threshold are satisfied:
    * Absolute log2 ratio `|log2((M+PSEUDO)/(P+PSEUDO))| ≥ --min-abs-log2`.
    * And `M ≥ P × min_fold` → **Maternal**, or `P ≥ M × min_fold` → **Paternal** (default `min_fold = 1.5`).
  * Else → **Balanced** (non‑significant or small effect size).

---

## Assumptions & current limitations

* Reads are **haplotype-tagged via RG**. If mates do **not** share RG, loop counting is conservative; explicit mate lookup will be added.
* `--validate-loops local` triggers an approximate local validation: a global z‑score is computed from the distribution of `(M‑P)` across all loops as a proxy for local enrichment.  Future versions will implement true distance‑matched backgrounds and refine the local p‑value.
* Motif orientation checks for CTCF are stubbed; planned as an optional enhancement.

---

## Project structure

```
src/lophos/
  cli.py                   # Typer CLI (phase & summary)
  core/                    # counts, stats, calls, validation (placeholder), APA, motif
  io/                      # BAM / BED / BEDPE / YAML loaders
  report/                  # writers, qc, and summary helpers
tests/                     # unit + integration tests
docs/                      # user docs
examples/                  # example configs/notebooks
```

---

## Roadmap

* Local loop enrichment (background & Z-scores)
* Explicit mate RG lookup in loop counting
* APA-lite matrices
* Config-driven runs & manifests
* Optional CTCF motif orientation checks
* Performance profiling and speedups where needed

---

## Contributing & Conduct

Please see:

* **[CONTRIBUTING](CONTRIBUTING.md)** — development workflow, style, testing
* **[CODE\_OF\_CONDUCT](CODE_OF_CONDUCT.md)** — expected behavior

PRs and issues are welcome.

---

## Citation / Acknowledgment

If this tool supports your work, please cite the repository and acknowledge:
**“Developed by Pranjul Mishra, under the guidance of Dr. Joanna Borkowska and Prof. Dariusz Plewczyński (Structural and Functional Genomics Laboratory).”**

---

## License

MIT — see `LICENSE`.
