# LOPHOS — LOops & Peaks HaplOtype phasing Suite

> Allele-specific phasing of CTCF **peaks** and **loops** from haplotype-tagged HiChIP BAMs.

**Status:** v0.1.1 (pre-alpha)
**Authorship:** Solely developed by **Pranjul Mishra**, under the guidance of **Dr. Joanna Borkowska** and **Prof. Dariusz Plewczyński** (Structural and Functional Genomics Laboratory).

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

```bash
lophos phase \
  --bam   /path/to/<sample>.bam \
  --peaks /path/to/<peaks>.bed \
  --loops /path/to/<loops>.bedpe \
  --out   results/<sample> \
  --mapq 30 --peak_window 500 --anchor_pad 10000 \
  --min_reads_peak 5 --min_pairs_loop 3 --fdr 0.05 \
  --validate-loops local
```

This produces:

* `results/<sample>.peaks.bed`
* `results/<sample>.loops.bedpe`
* `results/<sample>.summary.tsv`

---

## Input requirements

### BAM (haplotype-tagged)

Reads must carry a read group (**RG**) tag indicating parental origin. Recognized tokens (case-insensitive):

* Maternal: `maternal`, `mat`, `M`
* Paternal: `paternal`, `pat`, `P`

If your pipeline uses different tokens, extend the mapping in `src/lophos/constants.py`.

### Peaks (BED)

```
chrom  start  end  [name  score  strand]
```

### Loops (BEDPE)

```
chrom1  start1  end1  chrom2  start2  end2  [name  score  strand1  strand2]
```

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

| Option              | Purpose                              | Typical |
| ------------------- | ------------------------------------ | ------- |
| `--mapq`            | Minimum mapping quality to count     | `30`    |
| `--peak_window`     | Peak summit ±bp window for counting  | `500`   |
| `--anchor_pad`      | Padding around each loop anchor (bp) | `10000` |
| `--min_reads_peak`  | Minimum (M+P) reads to call a peak   | `5`     |
| `--min_pairs_loop`  | Minimum (M+P) pairs to call a loop   | `3`     |
| `--fdr`             | FDR threshold (BH)                   | `0.05`  |
| `--keep_duplicates` | Keep PCR/optical duplicates          | `False` |
| `--validate-loops`  | Extra QC (`none` or `local`)         | `local` |

---

## Recommended workflow

1. **Sanity-check RG tags**

   ```bash
   samtools view -H <sample>.bam | grep '^@RG'
   ```

   Confirm maternal/paternal tokens match the recognized mapping.

2. **Try one sample** with defaults (above).
   If many features are **Undetermined**, relax thresholds:

   * `--min_reads_peak 3` and/or `--min_pairs_loop 2`
   * `--mapq 20`
   * `--fdr 0.10`

3. **Sweep parameters** to find sensible settings for your dataset:

   ```bash
   for MQ in 20 30; do
     for PK in 3 5; do
       for LP in 2 3; do
         OUT=results/sweep/<sample>_mq${MQ}_pk${PK}_lp${LP}
         lophos phase --bam <bam> --peaks <bed> --loops <bedpe> \
           --out "$OUT" --mapq $MQ --peak_window 500 --anchor_pad 10000 \
           --min_reads_peak $PK --min_pairs_loop $LP --fdr 0.05 --validate-loops local
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
       --out "$OUT_DIR/$S" --mapq 30 --peak_window 500 --anchor_pad 10000 \
       --min_reads_peak 5 --min_pairs_loop 3 --fdr 0.05 --validate-loops local
   done
   ```

---

## Methods (concise)

* **Counts:** per feature, compute maternal (M) and paternal (P) support.
  Peaks: reads within ±`peak_window` around the summit.
  Loops: paired reads bridging padded anchors (`anchor_pad`).
* **Test:** `M ~ Binomial(M+P, 0.5)` (two-sided).
* **Multiple testing:** BH-FDR over features.
* **Calling:**

  * If `M+P` < minimum → **Undetermined**
  * Else if `FDR ≤ α` and `M ≥ P × min_fold` → **Maternal**
  * Else if `FDR ≤ α` and `P ≥ M × min_fold` → **Paternal**
  * Else → **Balanced**
    Default `min_fold = 1.5`.

---

## Assumptions & current limitations

* Reads are **haplotype-tagged via RG**. If mates do **not** share RG, loop counting is conservative; explicit mate lookup will be added.
* `--validate-loops local` is a placeholder; full local background/Z-score is planned.
* Motif orientation checks for CTCF are stubbed; planned as an optional enhancement.

---

## Project structure

```
src/lophos/
  cli.py                   # Typer CLI
  core/                    # counts, stats, calls, validation (placeholder), APA, motif
  io/                      # BAM / BED / BEDPE / YAML loaders
  report/                  # writers and summary
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
**“Solely developed by Pranjul Mishra, under the guidance of Dr. Joanna Borkowska and Prof. Dariusz Plewczyński (Structural and Functional Genomics Laboratory).”**

---

## License

MIT — see `LICENSE`.
