#!/usr/bin/env bash
set -euo pipefail

sheet=${1:?samplesheet.csv required}
out=${2:?output directory required}
manifest=${3:?checksum output required}
ena_report_url='https://www.ebi.ac.uk/ena/portal/api/filereport'

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

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
  log "[$sample/$role] Looking up ENA FASTQ files for $accession"
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
    destination="$out/${sample}_R${mate}.fastq.gz"
    expected_md5=${md5_list[$((mate - 1))]}

    if [[ -f "$destination" ]] && \
      printf '%s  %s\n' "$expected_md5" "$destination" | md5sum --check --status; then
      log "[$sample/$role] Mate $mate already exists and is verified; skipping download"
      continue
    fi

    [[ -f "$destination" ]] && log "[$sample/$role] Existing mate $mate failed MD5; redownloading"
    log "[$sample/$role] Downloading mate $mate ($file)"
    curl --fail --silent --show-error --location --retry 3 \
      --progress-bar "https://${url_list[$((mate - 1))]}" --output "$tmp/$file"
    log "[$sample/$role] Verifying MD5 for mate $mate"
    printf '%s  %s\n' "$expected_md5" "$tmp/$file" | md5sum --check --status || {
      echo "MD5 verification failed for $file" >&2
      exit 1
    }
    log "[$sample/$role] Mate $mate verified"
    mv "$tmp/$file" "$destination"
  done

  log "[$sample/$role] FASTQs ready in $out"
  rmdir "$tmp"
  tmp=''
done < "$sheet"

log "Writing SHA-256 manifest to $manifest"
sha256sum "$out"/*.fastq.gz > "$manifest"
log 'Download complete'
