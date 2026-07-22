process VEP {
  tag 'functional annotation'
  label 'large'
  container 'ensemblorg/ensembl-vep:release_116.0'
  publishDir "${params.outdir}/variants", mode:'copy'
  input: path vcf; path ref_bundle; path cache
  output: path 'somatic.pass.norm.vep.vcf.gz', emit: vcf
  script:
  def ref = ref_bundle.find { it.name ==~ /.*\.(fa|fasta)$/ }
  """
  vep --offline --cache --dir_cache ${cache} --species ${params.vep_species} \
    --assembly ${params.vep_assembly} --fasta ${ref} --vcf --compress_output bgzip \
    --canonical --hgvs --symbol --numbers --variant_class --sift b --polyphen b \
    --af_gnomade --clin_sig_allele 1 --check_existing --force_overwrite \
    -i ${vcf} -o somatic.pass.norm.vep.vcf.gz
  """
}

process BENCHMARK {
  tag 'SEQC2 v1.2'
  label 'small'
  container 'hcc1395-wes:1.0.0'
  publishDir "${params.outdir}/benchmark", mode:'copy'
  input: path calls; path truth_bundle; path confident; path ref_bundle
  output: path 'benchmark_*.tsv', emit: metrics
  script:
  def truth = truth_bundle.find { it.name.endsWith('.vcf.gz') }
  def ref = ref_bundle.find { it.name ==~ /.*\.(fa|fasta)$/ }
  """
  bcftools norm -f ${ref} -m -any ${truth} -Oz -o truth.norm.vcf.gz
  bcftools index -t truth.norm.vcf.gz
  python ${projectDir}/scripts/benchmark.py ${calls} truth.norm.vcf.gz ${confident} \
    benchmark_metrics.tsv benchmark_strata.tsv
  """
}

process PRIORITIZE {
  tag 'research prioritization'
  label 'small'
  container 'hcc1395-wes:1.0.0'
  publishDir "${params.outdir}/variants", mode:'copy'
  input: path annotated; path genes
  output: path 'prioritized_variants.tsv', emit: table
  script: """
  python ${projectDir}/scripts/prioritize_vep.py ${annotated} ${genes} prioritized_variants.tsv
  """
}

workflow ANNOTATE_BENCHMARK {
  take: calls; ref; truth; confident; cache; genes
  main:
    VEP(calls, ref, cache)
    PRIORITIZE(VEP.out.vcf, genes)
    BENCHMARK(calls, truth, confident, ref)
  emit: annotated = VEP.out.vcf; prioritized = PRIORITIZE.out.table; benchmark = BENCHMARK.out.metrics
}
