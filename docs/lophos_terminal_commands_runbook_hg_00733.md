you got it. here’s a clean, copy-paste cheat sheet of the exact terminal commands we used (and a few safer, corrected variants where you had syntax hiccups). i grouped them by step.

---

## 0) environment & repo (once per machine)

```bash
# choose python (if you use pyenv)
pyenv shell 3.10.14

# from your cloned repo root
pip install -e ".[dev]"
pre-commit install

# quick hygiene (optional but nice)
ruff check . --fix
black .
mypy src/lophos
pytest -q
pre-commit run --all-files
```

---

## 1) define inputs & outputs (HG00733)

```bash
# inputs
BAM="/mnt/evafs/groups/sfglab/jborkowska/ONT_trios/diploid_BAM/NM_filtered_BAM/merged_bams/HG00733_merged_phased.bam"
PEAKS_MAT="/mnt/evafs/groups/sfglab/abhishek/data/ONT_trios/HG00733_MAT_SG/final_output/HG00733_MAT_peaks.narrowPeak"
PEAKS_PAT="/mnt/evafs/groups/sfglab/abhishek/data/ONT_trios/HG00733_PAT_SG/final_output/HG00733_PAT_peaks.narrowPeak"
LOOPS_MAT="/mnt/evafs/groups/sfglab/abhishek/data/ONT_trios/HG00733_MAT_SG/final_output/HG00733_MAT.bedpe"
LOOPS_PAT="/mnt/evafs/groups/sfglab/abhishek/data/ONT_trios/HG00733_PAT_SG/final_output/HG00733_PAT.bedpe"

# paths for merged inputs
MERGE_DIR="/home2/faculty/pmishra/merge_script_lophos"
MERGED_PEAKS="$MERGE_DIR/HG00733_merged_peaks.bed"
MERGED_LOOPS="$MERGE_DIR/HG00733_merged_loops.bedpe"

# output prefix for lophos run
OUT="/home2/faculty/pmishra/lophos_runs/HG00733"
mkdir -p "$OUT"
```

---

## 2) merge peaks & loops

```bash
# (A) merge peaks (take first 3 cols; drop duplicates)
cd "$MERGE_DIR"

python merge_peaks_bed.py \
  --mat "$PEAKS_MAT" \
  --pat "$PEAKS_PAT" \
  --out "$MERGED_PEAKS" \
  --sort-output

# (B) merge loops (undirected; canonicalize anchor order; drop duplicates)
python merge_bedpe_undirected.py \
  --mat "$LOOPS_MAT" \
  --pat "$LOOPS_PAT" \
  --out "$MERGED_LOOPS" \
  --sort-output
```

> both scripts accept headered files; they skip comment/headers and only use the BED/BEDPE coordinate columns.

---

## 3) ensure BAM is indexed

```bash
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  samtools index -@ 8 "$BAM"
fi
```

---

## 4) run LOPHOS

```bash
lophos phase \
  --bam   "$BAM" \
  --peaks "$MERGED_PEAKS" \
  --loops "$MERGED_LOOPS" \
  --out   "$OUT" \
  --mapq 30 \
  --peak-window 500 \
  --anchor-pad 10000 \
  --min-reads-peak 5 \
  --min-pairs-loop 3 \
  --fdr 0.05 \
  --validate-loops local
```

Outputs created:
```bash
"$OUT/HG00733.peaks.bed"
"$OUT/HG00733.loops.bedpe"
"$OUT/HG00733.summary.tsv"
```

---

## 5) quick sanity checks

```bash
# total rows
wc -l "$OUT"/HG00733.*

# first few lines
head -n 3 "$OUT/HG00733.peaks.bed"
head -n 3 "$OUT/HG00733.loops.bedpe"
```

---

## 6) call distributions

```bash
# peaks calls
awk '{print $NF}' "$OUT/HG00733.peaks.bed" | sort | uniq -c

# loops calls
awk '{print $NF}' "$OUT/HG00733.loops.bedpe" | sort | uniq -c
```

---

## 7) significance counts (match CLI thresholds)

```bash
# peaks: FDR in column 10, (M+P) in column 7
awk '$10 <= 0.05 && $7 >= 5 {c++} END{print "peaks_signif:", c+0}' "$OUT/HG00733.peaks.bed"

# loops: FDR = (NF-1), total_pairs = (NF-4)
awk '{
  n=NF; fdr=$(n-1); tot=$(n-4);
  if (fdr <= 0.05 && tot >= 3) c++
} END{print "loops_signif:", c+0}' "$OUT/HG00733.loops.bedpe"
```

---

## 8) direction-consistency checks
(Maternal ↔ positive log2 ratio; Paternal ↔ negative.)

```bash
# peaks: log2 ratio is column 8
awk '($8>0 && $NF=="Maternal") || ($8<0 && $NF=="Paternal"){ok++} END{print "peaks_direction_consistent:", ok+0}' "$OUT/HG00733.peaks.bed"

# loops: log2 ratio is (NF-3)
awk '{n=NF; lr=$(n-3); if ((lr>0 && $NF=="Maternal") || (lr<0 && $NF=="Paternal")) ok++} END{print "loops_direction_consistent:", ok+0}' "$OUT/HG00733.loops.bedpe"
```

---

## 9) top strongly biased features (FDR ≤ 0.01)

```bash
# LOOPS — strongest absolute log2 ratio (any call)
FILE="$OUT/HG00733.loops.bedpe"
awk -v OFS='\t' '{
  n=NF; fdr=$(n-1); lr=$(n-3);
  if (fdr<=0.01) { abs=(lr<0?-lr:lr); print abs, $0 }
}' "$FILE" | sort -k1,1gr | head | cut -f2-

# LOOPS — strongest Maternal
awk -v OFS='\t' '{
  n=NF; fdr=$(n-1); lr=$(n-3); call=$n;
  if (fdr<=0.01 && call=="Maternal") { abs=(lr<0?-lr:lr); print abs, $0 }
}' "$FILE" | sort -k1,1gr | head | cut -f2-

# LOOPS — strongest Paternal
awk -v OFS='\t' '{
  n=NF; fdr=$(n-1); lr=$(n-3); call=$n;
  if (fdr<=0.01 && call=="Paternal") { abs=(lr<0?-lr:lr); print abs, $0 }
}' "$FILE" | sort -k1,1gr | head | cut -f2-

# PEAKS — strongest absolute log2 ratio (FDR ≤ 0.01)
PEAKS="$OUT/HG00733.peaks.bed"
awk -v OFS='\t' '($10<=0.01){ lr=$8; abs=(lr<0?-lr:lr); print abs, $0 }' "$PEAKS" \
  | sort -k1,1gr | head | cut -f2-

# PEAKS — strongest Maternal (FDR ≤ 0.01)
awk -v OFS='\t' '($10<=0.01 && $11=="Maternal"){ lr=$8; abs=(lr<0?-lr:lr); print abs, $0 }' "$PEAKS" \
  | sort -k1,1gr | head | cut -f2-

# PEAKS — strongest Paternal (FDR ≤ 0.01)
awk -v OFS='\t' '($10<=0.01 && $11=="Paternal"){ lr=$8; abs=(lr<0?-lr:lr); print abs, $0 }' "$PEAKS" \
  | sort -k1,1gr | head | cut -f2-
```

> note: earlier you tried `sort -k( NF-1 ),( NF-1 )gr` which is invalid syntax. the snippets above add the absolute value as **first** column and sort on `-k1,1gr`.

---

## 10) quick QC summary (python helper we used)

If you saved the helper as `summary.py` (the small Typer script we drafted), run:

```bash
python summary.py --out "$OUT" --prefix HG00733
```

This prints a QC block and (optionally) writes a small TSV.

---

If you want, I can throw all of this into a `RUNBOOK.md` in the repo so the whole team can reuse it without scrolling through chat.
