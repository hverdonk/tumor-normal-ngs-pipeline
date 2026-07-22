#!/usr/bin/env bash
set -euo pipefail
sheet=${1:?samplesheet required}; params=${2:?params yaml required}
[[ $(tail -n +2 "$sheet" | cut -d, -f2 | sort | tr '\n' ' ') == 'normal tumor ' ]] || { echo 'Need exactly tumor and normal roles' >&2; exit 1; }
while IFS=, read -r sample role accession r1 r2 rest; do
  [[ -s $r1 && -s $r2 ]] || { echo "Missing FASTQ for $sample" >&2; exit 1; }
done < <(tail -n +2 "$sheet")
for key in reference targets gnomad common_sites truth_vcf truth_bed vep_cache cancer_genes; do
  value=$(awk -F': *' -v k="$key" '$1==k {sub(/^ +/,"",$2); print $2}' "$params")
  [[ -n $value && $value != null && -e $value ]] || { echo "Missing $key: $value" >&2; exit 1; }
done
echo 'Input paths and sample roles passed structural checks.'
