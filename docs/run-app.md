# inputs
BAM="/mnt/evafs/groups/sfglab/jborkowska/ONT_trios/diploid_BAM/NM_filtered_BAM/merged_bams/HG00733_merged_phased.bam"
PEAKS="/home2/faculty/pmishra/merge_script_lophos/HG00733_merged_peaks.bed"     # BED3
LOOPS="/home2/faculty/pmishra/merge_script_lophos/HG00733_merged_loops.bedpe"   # BEDPE (+ extra 'src' col OK)

# output
OUT="/home2/faculty/pmishra/lophos_runs/HG00733"


mkdir -p "$OUT"

# (optional) make sure the BAM is indexed
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  samtools index -@ 8 "$BAM"
fi

# run LOPHOS
lophos phase \
  --bam   "$BAM" \
  --peaks "$PEAKS" \
  --loops "$LOOPS" \
  --out   "$OUT" \
  --mapq 30 \
  --peak-window 500 \
  --anchor-pad 10000 \
  --min-reads-peak 5 \
  --min-pairs-loop 3 \
  --fdr 0.05 \
  --validate-loops local
