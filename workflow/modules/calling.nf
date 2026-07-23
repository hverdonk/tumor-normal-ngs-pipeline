process CALL_AND_FILTER {
  tag 'HCC1395 matched pair'
  label 'large'
  container 'hcc1395-wes:1.0.0'
  publishDir "${params.outdir}/variants", mode:'copy'
  input:
    tuple val(tumor_name), path(tumor_bam), path(tumor_bai),
          val(normal_name), path(normal_bam), path(normal_bai)
    path ref_bundle
    path confident
    path gnomad_bundle
    path common_bundle
    path pon_bundle
  output:
    path 'somatic.pass.norm.vcf.gz', emit: pass
    path 'somatic.pass.norm.vcf.gz.tbi'
    path '*.{vcf.gz,vcf.gz.tbi,tar.gz,table}', emit: intermediates
  script:
  def ref = ref_bundle.find { it.name ==~ /.*\.(fa|fasta)$/ }
  def gnomad = gnomad_bundle.find { it.name.endsWith('.vcf.gz') }
  def common = common_bundle.find { it.name.endsWith('.vcf.gz') }
  def ponVcf = params.skip_pon ? null : pon_bundle.find { it.name.endsWith('.vcf.gz') }
  def ponArg = params.skip_pon ? '' : "--panel-of-normals ${ponVcf}"
  """
  gatk Mutect2 -R ${ref} -I ${tumor_bam} -I ${normal_bam} \
    -tumor ${tumor_name} -normal ${normal_name} -L ${confident} \
    --germline-resource ${gnomad} ${ponArg} \
    --f1r2-tar-gz f1r2.tar.gz -O somatic.unfiltered.vcf.gz
  gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
  gatk GetPileupSummaries -I ${tumor_bam} -V ${common} -L ${confident} -O tumor.pileups.table
  gatk GetPileupSummaries -I ${normal_bam} -V ${common} -L ${confident} -O normal.pileups.table
  gatk CalculateContamination -I tumor.pileups.table -matched normal.pileups.table \
    -O contamination.table --tumor-segmentation segments.table
  gatk FilterMutectCalls -R ${ref} -V somatic.unfiltered.vcf.gz \
    --contamination-table contamination.table --tumor-segmentation segments.table \
    --ob-priors read-orientation-model.tar.gz -O somatic.filtered.vcf.gz
  gatk SelectVariants -R ${ref} -V somatic.filtered.vcf.gz --exclude-filtered \
    --exclude-non-variants -O somatic.pass.vcf.gz
  bcftools norm -f ${ref} -m -any somatic.pass.vcf.gz -Oz -o somatic.pass.norm.vcf.gz
  bcftools index -t somatic.pass.norm.vcf.gz
  """
}

workflow MUTECT2_FILTER {
  take: pair; ref; confident; gnomad; common; pon
  main: CALL_AND_FILTER(pair, ref, confident, gnomad, common, pon)
  emit: pass_vcf = CALL_AND_FILTER.out.pass
}
