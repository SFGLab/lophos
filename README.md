# LOPHOS — LOops & Peaks HaplOtype phasing Suite

> Allele-specific phasing of CTCF **peaks** and **loops** from haplotype-tagged BAMs (paired-end or long-read SA:Z chimeras).

---

## Overview (for biologists)

LOPHOS quantifies **maternal vs paternal** support for:

* **CTCF peaks** (BED) using haplotype-tagged reads around each peak summit.
* **Loops** (BEDPE) using either:

  * **paired mates** bridging the two anchors (**mates** mode), or
  * **split (SA:Z) segments** within a single long read that span the anchors (**sa** mode; ONT/Pore-C-style).

For each feature, LOPHOS:

1. **Counts** maternal (M) and paternal (P) support (ambiguous tracked separately for loops).
2. **Tests** allele bias with a two-sided **binomial test** (`p = 0.5`).
3. **Adjusts** p-values using **Benjamini–Hochberg (FDR)**.
4. **Calls** each feature as **Maternal / Paternal / Balanced / Undetermined**.

Outputs are plain **BED/BEDPE** tables you can sort, filter, and load in genome browsers.

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

### Phase peaks & loops (paired-end “mates” mode, default)

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
  --validate-loops local \
  --loop-mode mates
  # Optional:
  # --primary-only true  # write .primary.* filtered outputs
  # --summary true       # also write <prefix>.quick_qc.tsv
```

### Phase loops from long-read SA:Z chimeras (ONT/Pore-C “sa” mode)

```bash
lophos phase \
  --bam   /path/to/<sample>.bam \
  --peaks /path/to/<peaks>.bed \
  --loops /path/to/<loops>.bedpe \
  --out   results/<sample> \
  --mapq 30 --peak-window 500 --anchor-pad 10000 \
  --min-reads-peak 5 --min-pairs-loop 3 --fdr 0.05 \
  --loop-mode sa \
  --sa-min-mapq 30 \
  --sa-min-seg-len 50 \
  --sa-min-cis-dist 1000 \
  --sa-allow-trans true \
  --sa-orientation any \
  --sa-dedup-within-read true
```

This produces:

* `results/<sample>.peaks.bed` — allele calls for peaks (see schema).
* `results/<sample>.loops.bedpe` — allele calls for loops (see schema).
* `results/<sample>.summary.tsv` — compact run summary.
* If `--summary true`, also `results/<sample>.quick_qc.tsv`.
* If `--primary-only true`, additional: `.primary.peaks.bed`, `.primary.loops.bedpe`.

### Summarize an existing run

```bash
lophos summary \
  --out results/<sample> \
  --prefix <sample> \
  --fdr 0.05 \
  --min-reads-peak 5 \
  --min-pairs-loop 3
```

Prints a tidy report and (by default) writes `<prefix>.quick_qc.tsv`.

---

## Input requirements

### BAM (haplotype-tagged)

Reads must carry a read group (**RG**) tag indicating parental origin.
Defaults (case-insensitive):

* **Maternal:** `maternal|mat|M`
* **Paternal:** `paternal|pat|P`

Override with `--maternal-rgid` / `--paternal-rgid` as needed.

> SA mode expects long-read alignments where split mappings are encoded in **SA:Z** tags. Mates mode expects paired-end flags.

### Peaks (BED)

```
chrom  start  end  [name  score  strand]
```

### Loops (BEDPE)

```
chrom1  start1  end1  chrom2  start2  end2  [name ...]
```

Extra columns are tolerated; only the first 6 are required.

---

## Output formats

### Peaks (`*.peaks.bed`)

```
chrom  start  end  peak_id  maternal  paternal  total  log2_ratio  p_value  fdr  bias_call
```

### Loops (`*.loops.bedpe`)

```
chrom1  start1  end1  chrom2  start2  end2
total_pairs  maternal_pairs  paternal_pairs  ambiguous_pairs  informative_pairs
log2_ratio_pairs  p_value_pairs  fdr_pairs  bias_call
```

> `informative_pairs = maternal_pairs + paternal_pairs` and equals `total_pairs` in the current implementation (ambiguous excluded from stats).

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

---

## Loop modes

LOPHOS supports two loop counting modes:

* **`mates`** (default): Paired-end HiChIP/Hi-C style. A contact is counted when one read maps in anchor A (±`anchor-pad`) and its **mate** maps in anchor B (±`anchor-pad`), in either direction. Allele is taken from RG. Use for standard paired-end data.

* **`sa`** (long-read chimeric; ONT/Pore-C): Contacts are reconstructed from **SA:Z** split alignments within a single read. We pair **adjacent segments** (A–B, B–C, …), filter by per-segment **MAPQ** and **reference length**, drop very short cis links (`< --sa-min-cis-dist`) to avoid self-ligation artifacts, optionally restrict **orientation** for short-range cis, and optionally disallow **trans**. Each accepted contact is tested against loop anchors (±`anchor-pad`). Allele comes from RG; ambiguous haplotypes are tracked.

### SA parameters

* `--sa-min-mapq` (default **30**): per-segment MAPQ cutoff.
* `--sa-min-seg-len` (default **50**): minimum reference-consumed length from CIGAR.
* `--sa-min-cis-dist` (default **1000**): drop cis contacts closer than this (artifact guard).
* `--sa-allow-trans / --no-sa-allow-trans` (default **allow**).
* `--sa-orientation {any, convergent-short-cis}` (default **any**).
* `--sa-dedup-within-read / --no-sa-dedup-within-read` (default **on**).

---

## CLI options (key)

| Option                 | Purpose                                                   | Typical    |     |     |
| ---------------------- | --------------------------------------------------------- | ---------- | --- | --- |
| `--mapq`               | Minimum mapping quality to count reads                    | `30`       |     |     |
| `--peak-window`        | Peak summit ± bp window for counting reads                | `500`      |     |     |
| `--anchor-pad`         | Padding (bp) added around each loop anchor                | `10000`    |     |     |
| `--loop-mode`          | `mates` (paired-end) or `sa` (SA:Z long-read)             | `mates`    |     |     |
| `--sa-*`               | SA-specific knobs (see above)                             | —          |     |     |
| `--min-reads-peak`     | Min (M+P) reads required to call a peak                   | `5`        |     |     |
| `--min-pairs-loop`     | Min informative pairs (M+P) required to call a loop       | `3`        |     |     |
| `--fdr`                | BH-FDR threshold                                          | `0.05`     |     |     |
| `--keep-duplicates`    | Count duplicates                                          | `False`    |     |     |
| `--validate-loops`     | `none` or `local` (z-score proxy)                         | `local`    |     |     |
| `--maternal-rgid`      | Regex (case-insensitive) for maternal RG                  | `"maternal | mat | M"` |
| `--paternal-rgid`      | Regex for paternal RG                                     | `"paternal | pat | P"` |
| `--pseudocount`        | Pseudocount in log2 ratio                                 | `1.0`      |     |     |
| `--min-abs-log2`       | Min absolute log2 ratio for calling bias                  | `0.0`      |     |     |
| `--max-ambiguous-frac` | Max fraction ambiguous pairs tolerated when calling loops | `0.5`      |     |     |
| `--primary-only`       | Write `.primary.*` (primary chroms only)                  | `False`    |     |     |
| `--summary`            | Also write `<prefix>.quick_qc.tsv`                        | `False`    |     |     |
| `--threads`            | Threads for counting                                      | `1`        |     |     |
| `--log-level`          | `info`, `debug`, `warning`, `error`                       | `info`     |     |     |
| `--config`             | YAML overrides (persisted to `<out>.run.json`)            | `None`     |     |     |

---

## Recommended workflow

1. **Sanity-check RG tags**

   ```bash
   samtools view -H <sample>.bam | grep '^@RG'
   ```

   Confirm maternal/paternal identifiers. If different, provide custom regex via `--maternal-rgid` / `--paternal-rgid`.

2. **Choose loop mode**

   * Paired-end data → `--loop-mode mates`.
   * Long-read (ONT/Pore-C) with SA:Z chimeras → `--loop-mode sa` and tune `--sa-*` knobs.

3. **Trial run**, then adjust thresholds if many features are **Undetermined**:

   * Lower `--min-reads-peak` / `--min-pairs-loop`
   * `--mapq 20`
   * `--fdr 0.10`

4. **Batch** across samples (as in your current pipeline).

5. **Summarize** with `lophos summary` or use `--summary true` during `phase`.

---

## Methods (concise)

* **Counts**

  * **Peaks:** reads within ±`peak-window` of the summit (MAPQ/dup filters).
  * **Loops (mates):** one read in anchor A (±`anchor-pad`) and its **mate** in anchor B (±`anchor-pad`).
  * **Loops (sa):** **adjacent SA segments** within a long read form a contact; per-segment MAPQ/length filters; short-cis filter; optional orientation/trans filters; dedup within read; accept if the two endpoints land in the two anchors (either order). Ambiguous haplotypes are tallied separately.

* **Statistics:** two-sided binomial @ 0.5; BH-FDR.

* **Calling:**

  * Below coverage floor → **Undetermined**.
  * Loops with `ambiguous_frac > --max-ambiguous-frac` → **Undetermined**.
  * Otherwise, if `FDR ≤ α` and both effect size (|log2| ≥ `--min-abs-log2`) and a fold guard (≥ 1.5×) pass → **Maternal** or **Paternal**; else **Balanced**.

---

## Assumptions & current limitations

* Reads are **haplotype-tagged via RG**. SA mode propagates RG from the read to each contact.
* `--validate-loops local` uses a global `(M−P)` proxy for z/p; distance-matched backgrounds are planned.
* Motif/orientation checks for CTCF are optional (basic short-cis convergent policy available in SA mode).

---

## Project structure

```
src/lophos/
  cli.py
  core/
    counts_peaks.py
    counts_loops.py         # mates + sa dispatcher
    sa_pairs.py             # SA:Z → contacts (new)
    stats.py
    calls.py
    validate_local.py
  io/
    bam.py                  # RG mapping + SA helpers (new)
    bed.py
    bedpe.py
    config.py
  report/
    writers.py
    qc.py
    summary.py
tests/
docs/
```

---

## Roadmap

* Distance-matched loop validation & APA-lite
* Refinements for SA contact orientation rules
* Config-driven batch/manifest execution
* Golden tiny BAM for CI (peaks + loops, mates + sa)

---

## Contributing & Conduct

See **CONTRIBUTING.md** and **CODE_OF_CONDUCT.md**. PRs and issues welcome.

---

## Citation / Acknowledgment

If this tool supports your work, please cite the repository and acknowledge:
**“Developed by Pranjul Mishra, under the guidance of Dr. Joanna Borkowska and Prof. Dariusz Plewczyński (Structural and Functional Genomics Laboratory).”**

---

## License

MIT — see `LICENSE`.
