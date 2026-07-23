#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QC_TRIM; ALIGN_QC } from './modules/preprocess'

workflow {
  if (params.containsKey('targets'))
    error 'targets is not configurable; truth_bed defines both capture and benchmark regions'
  rows = Channel.fromPath(params.samplesheet, checkIfExists:true)
    .splitCsv(header:true)
    .map { r -> tuple(r.sample, r.role, r.platform, r.library,
                      file(r.fastq_1, checkIfExists:true), file(r.fastq_2, checkIfExists:true)) }
  refPath = params.reference.toString()
  dictPath = refPath.replaceFirst(/\.(fa|fasta)$/, '.dict')
  ref = Channel.value([file(refPath,checkIfExists:true), file(refPath+'.fai',checkIfExists:true), file(dictPath,checkIfExists:true)])
  confident = file(params.truth_bed, checkIfExists:true)
  QC_TRIM(rows)
  ALIGN_QC(QC_TRIM.out.cleaned, ref, confident)
}
