process FASTQC_PRE {
  tag "$sample pre"
  label 'small'
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  publishDir "${params.outdir}/qc/fastqc/pre", mode:'copy'
  input: tuple val(sample), val(role), val(platform), val(library), path(r1), path(r2)
  output: path '*_fastqc.{html,zip}', emit: reports
  script: """
  fastqc --threads ${task.cpus} ${r1} ${r2}
  """
}

process FASTP {
  tag "$sample"
  label 'medium'
  container 'quay.io/biocontainers/fastp:1.3.6--h43da1c4_0'
  publishDir "${params.outdir}/clean", mode:'copy'
  input: tuple val(sample), val(role), val(platform), val(library), path(r1), path(r2)
  output:
    tuple val(sample), val(role), val(platform), val(library), path("${sample}_R1.clean.fq.gz"), path("${sample}_R2.clean.fq.gz"), emit: reads
    path "${sample}.fastp.{html,json}", emit: reports
  script: """
  fastp -i ${r1} -I ${r2} -o ${sample}_R1.clean.fq.gz -O ${sample}_R2.clean.fq.gz \
    --thread ${task.cpus} --detect_adapter_for_pe --trim_poly_g \
    --disable_quality_filtering --length_required 30 \
    --html ${sample}.fastp.html --json ${sample}.fastp.json
  """
}

process FASTQC_POST {
  tag "$sample post"
  label 'small'
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  publishDir "${params.outdir}/qc/fastqc/post", mode:'copy'
  input: tuple val(sample), val(role), val(platform), val(library), path(r1), path(r2)
  output: path '*_fastqc.{html,zip}', emit: reports
  script: """fastqc --threads ${task.cpus} ${r1} ${r2}"""
}

workflow QC_TRIM {
  take: reads
  main:
    FASTQC_PRE(reads)
    FASTP(reads)
    FASTQC_POST(FASTP.out.reads)
  emit:
    cleaned = FASTP.out.reads
    qc = FASTQC_PRE.out.reports.mix(FASTP.out.reports, FASTQC_POST.out.reports)
}

process ALIGN {
  tag "$sample"
  label 'large'
  container 'hcc1395-wes:1.0.0'
  input:
    tuple val(sample), val(role), val(platform), val(library), path(r1), path(r2)
    path ref_bundle
  output: tuple val(sample), val(role), path("${sample}.sorted.bam"), emit: bam
  script:
  def ref = ref_bundle.find { it.name ==~ /.*\.(fa|fasta)$/ }
  def rg = "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${library}\\tPL:${platform}\\tPU:${sample}.1"
  """
  bwa-mem3 index ${ref}
  bwa-mem3 mem -t ${task.cpus} -R '${rg}' ${ref} ${r1} ${r2} | \
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
  """
}

process MARKDUP_QC {
  tag "$sample"
  label 'medium'
  container 'hcc1395-wes:1.0.0'
  publishDir "${params.outdir}/bam", mode:'copy', pattern:'*.{bam,bai}'
  publishDir "${params.outdir}/qc/alignment", mode:'copy', pattern:'*.{txt,tsv,pdf}'
  input:
    tuple val(sample), val(role), path(bam)
    path ref_bundle
  output:
    tuple val(sample), val(role), path("${sample}.md.bam"), path("${sample}.md.bai"), emit: bam
    path "${sample}.*.{txt,tsv,pdf}", emit: metrics
  script: """
  picard MarkDuplicates I=${bam} O=${sample}.md.bam M=${sample}.duplicates.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  samtools flagstat -@ ${task.cpus} ${sample}.md.bam > ${sample}.flagstat.txt
  samtools stats -@ ${task.cpus} ${sample}.md.bam > ${sample}.samtools_stats.txt
  picard CollectInsertSizeMetrics I=${sample}.md.bam O=${sample}.insert_size.txt H=${sample}.insert_size.pdf
  """
}

process MOSDEPTH {
  tag "$sample"
  label 'medium'
  container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0'
  publishDir "${params.outdir}/qc/coverage", mode:'copy'
  input:
    tuple val(sample), val(role), path(bam), path(bai)
    path confident
  output: path "${sample}.mosdepth*", emit: metrics
  script: """
  mosdepth --threads ${task.cpus} --by ${confident} --thresholds 20,50,100 ${sample}.mosdepth ${bam}
  """
}

workflow ALIGN_QC {
  take: reads; ref; confident
  main:
    ALIGN(reads, ref)
    MARKDUP_QC(ALIGN.out.bam, ref)
    MOSDEPTH(MARKDUP_QC.out.bam, confident)
  emit:
    bams = MARKDUP_QC.out.bam
    metrics = MARKDUP_QC.out.metrics.mix(MOSDEPTH.out.metrics)
}
