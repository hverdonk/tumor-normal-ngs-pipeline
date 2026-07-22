#!/usr/bin/env bash
set -euo pipefail
out=${1:-test-data}; mkdir -p "$out/vep"
printf '>chr1\n' > "$out/ref.fa"
awk 'BEGIN{for(i=0;i<1000;i++)printf "ACGT";printf "\n"}' >> "$out/ref.fa"
printf 'chr1\t0\t4000\n' > "$out/targets.bed"
cp "$out/targets.bed" "$out/confident.bed"
for sample in HCC1395 HCC1395BL; do
  for mate in 1 2; do
    { printf '@%s_%s\n' "$sample" "$mate"; printf 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+%s\n' '+'; printf 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n'; } | gzip -c > "$out/${sample}_R${mate}.fq.gz"
  done
done
printf 'sample,role,accession,fastq_1,fastq_2,platform,library\nHCC1395,tumor,TEST_T,%s/HCC1395_R1.fq.gz,%s/HCC1395_R2.fq.gz,ILLUMINA,test\nHCC1395BL,normal,TEST_N,%s/HCC1395BL_R1.fq.gz,%s/HCC1395BL_R2.fq.gz,ILLUMINA,test\n' "$out" "$out" "$out" "$out" > "$out/samplesheet.csv"
printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$out/resource.vcf"
cp "$out/resource.vcf" "$out/truth.vcf"
bgzip -f "$out/resource.vcf"; tabix -f -p vcf "$out/resource.vcf.gz"
bgzip -f "$out/truth.vcf"; tabix -f -p vcf "$out/truth.vcf.gz"
samtools faidx "$out/ref.fa"
gatk CreateSequenceDictionary -R "$out/ref.fa"
echo "Synthetic fixtures created in $out. VEP cache is intentionally empty; use the test only through preprocessing/alignment or provide a cache fixture."

