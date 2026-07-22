#!/usr/bin/env bash
set -euo pipefail
sheet=${1:?samplesheet.csv required}; out=${2:?output directory required}; manifest=${3:?checksum output required}
mkdir -p "$out" "$(dirname "$manifest")"
command -v fasterq-dump >/dev/null || { echo 'fasterq-dump not found' >&2; exit 1; }
command -v pigz >/dev/null || { echo 'pigz not found' >&2; exit 1; }
tail -n +2 "$sheet" | while IFS=, read -r sample role accession fq1 fq2 platform library; do
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/hcc1395.XXXXXX")
  trap 'rm -rf "$tmp"' EXIT
  prefetch "$accession" --output-directory "$tmp"
  fasterq-dump --split-files --threads 6 --temp "$tmp" --outdir "$tmp" "$tmp/$accession/$accession.sra"
  pigz -p 6 "$tmp/${accession}_1.fastq" "$tmp/${accession}_2.fastq"
  mv "$tmp/${accession}_1.fastq.gz" "$out/${sample}_R1.fastq.gz"
  mv "$tmp/${accession}_2.fastq.gz" "$out/${sample}_R2.fastq.gz"
  rm -rf "$tmp"; trap - EXIT
done
shasum -a 256 "$out"/*.fastq.gz > "$manifest"

