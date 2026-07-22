#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_TRIM; ALIGN_QC } from './modules/preprocess'
include { MUTECT2_FILTER } from './modules/calling'
include { ANNOTATE_BENCHMARK } from './modules/reporting'

process MULTIQC_REPORT {
  label 'small'
  container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
  publishDir "${params.outdir}/qc", mode:'copy'
  input: path inputs
  output: path 'multiqc_report.html'; path 'multiqc_data'
  script: """multiqc --force --outdir . ."""
}

workflow {
  if (!params.reference || !params.targets || !params.gnomad || !params.common_sites ||
      !params.truth_vcf || !params.truth_bed || !params.vep_cache)
    error 'Missing required resource path(s); see config/local.example.yaml'
  if (!params.skip_pon && !params.pon) error 'pon is required unless --skip_pon true is explicit'

  rows = Channel.fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header:true)
    .map { r ->
      assert r.role in ['tumor','normal'] : "Invalid role: ${r.role}"
      tuple(r.sample, r.role, r.platform, r.library,
            file(r.fastq_1, checkIfExists:true), file(r.fastq_2, checkIfExists:true))
    }
  rows.map{ it[1] }.collect().map { roles ->
    assert roles.count('tumor') == 1 && roles.count('normal') == 1 : 'Exactly one tumor and normal are required'
  }

  refPath = params.reference.toString()
  dictPath = refPath.replaceFirst(/\.(fa|fasta)$/, '.dict')
  ref = Channel.value([file(refPath, checkIfExists:true), file(refPath+'.fai', checkIfExists:true), file(dictPath, checkIfExists:true)])
  targets = file(params.targets, checkIfExists:true)
  gnomad = Channel.value([file(params.gnomad, checkIfExists:true), file(params.gnomad+'.tbi', checkIfExists:true)])
  common = Channel.value([file(params.common_sites, checkIfExists:true), file(params.common_sites+'.tbi', checkIfExists:true)])
  truth = Channel.value([file(params.truth_vcf, checkIfExists:true), file(params.truth_vcf+'.tbi', checkIfExists:true)])
  confident = file(params.truth_bed, checkIfExists:true)
  vep_cache = file(params.vep_cache, checkIfExists:true)
  cancer_genes = file(params.cancer_genes, checkIfExists:true)
  pon = params.skip_pon ? Channel.value([]) : Channel.value([file(params.pon, checkIfExists:true), file(params.pon+'.tbi', checkIfExists:true)])

  QC_TRIM(rows)
  ALIGN_QC(QC_TRIM.out.cleaned, ref, targets)
  MULTIQC_REPORT(QC_TRIM.out.qc.mix(ALIGN_QC.out.metrics).collect())
  paired = ALIGN_QC.out.bams.collect()
  MUTECT2_FILTER(paired, ref, targets, gnomad, common, pon)
  ANNOTATE_BENCHMARK(MUTECT2_FILTER.out.pass_vcf, ref, truth, confident, vep_cache, cancer_genes)
}
