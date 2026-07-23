#!/usr/bin/env bash
set -euo pipefail

sheet=${1:?samplesheet.csv required}
out=${2:?output directory required}
manifest=${3:?checksum output required}
ena_report_url='https://www.ebi.ac.uk/ena/portal/api/filereport'

mkdir -p "$out" "$(dirname "$manifest")"
command -v curl >/dev/null || { echo 'curl not found' >&2; exit 1; }
command -v md5sum >/dev/null || { echo 'md5sum not found' >&2; exit 1; }
command -v sha256sum >/dev/null || { echo 'sha256sum not found' >&2; exit 1; }

tmp=''
cleanup() {
  [[ -z "$tmp" ]] || rm -rf "$tmp"
}
trap cleanup EXIT

while IFS=, read -r sample role accession fq1 fq2 platform library; do
  [[ "$sample" == 'sample' ]] && continue
  report=$(curl --fail --silent --show-error --location --retry 3 \
    --get "$ena_report_url" \
    --data-urlencode "accession=$accession" \
    --data-urlencode 'result=read_run' \
    --data-urlencode 'fields=run_accession,fastq_ftp,fastq_md5' \
    --data-urlencode 'format=tsv')
  row=$(printf '%s\n' "$report" | awk 'NR == 2')
  [[ -n "$row" ]] || { echo "No ENA FASTQ files found for $accession" >&2; exit 1; }

  IFS=$'\t' read -r returned_accession urls md5s <<< "$row"
  [[ "$returned_accession" == "$accession" ]] || {
    echo "ENA returned $returned_accession while requesting $accession" >&2
    exit 1
  }
  IFS=';' read -r -a url_list <<< "$urls"
  IFS=';' read -r -a md5_list <<< "$md5s"
  [[ ${#url_list[@]} -eq 2 && ${#md5_list[@]} -eq 2 ]] || {
    echo "Expected two paired FASTQs and MD5s for $accession" >&2
    exit 1
  }

  tmp=$(mktemp -d "$out/.${sample}.ena.XXXXXX")
  for mate in 1 2; do
    file="${accession}_${mate}.fastq.gz"
    curl --fail --silent --show-error --location --retry 3 \
      "https://${url_list[$((mate - 1))]}" --output "$tmp/$file"
    printf '%s  %s\n' "${md5_list[$((mate - 1))]}" "$tmp/$file" | md5sum --check --status || {
      echo "MD5 verification failed for $file" >&2
      exit 1
    }
  done

  mv "$tmp/${accession}_1.fastq.gz" "$out/${sample}_R1.fastq.gz"
  mv "$tmp/${accession}_2.fastq.gz" "$out/${sample}_R2.fastq.gz"
  rmdir "$tmp"
  tmp=''
done < "$sheet"

sha256sum "$out"/*.fastq.gz > "$manifest"
