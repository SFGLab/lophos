# LOPHOS — LOops & Peaks HaplOtype phasing Suite

Allele-specific phasing of CTCF peaks & loops from phased HiChIP BAMs.

## Quickstart
```bash
pip install -e .[dev]
lophos --help


## Example
lophos phase \
  --bam /mnt/.../merged_bams/<sample>.bam \
  --peaks /mnt/.../personalized_peaks/<peaks>.bed \
  --loops /mnt/.../personalized_loops/loops/<loops>.bedpe \
  --out results/<sample> --validate-loops local
