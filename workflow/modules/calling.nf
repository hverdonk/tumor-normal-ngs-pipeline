process CALL_AND_FILTER {
  tag 'HCC1395 matched pair'
  label 'large'
  container 'hcc1395-wes:1.0.0'
  publishDir "${params.outdir}/variants", mode:'copy'
  input:
    val pairs
    path ref_bundle
    path targets
    path gnomad_bundle
    path common_bundle
    path pon_bundle
  output:
    path 'somatic.pass.norm.vcf.gz', emit: pass
    path 'somatic.pass.norm.vcf.gz.tbi'
    path '*.{vcf.gz,vcf.gz.tbi,tar.gz,table}', emit: intermediates
  script:
  def tumor = pairs.find{ it[1] == 'tumor' }
  def normal = pairs.find{ it[1] == 'normal' }
  def ref = ref_bundle.find { it.name ==~ /.*\.(fa|fasta)$/ }
  def gnomad = gnomad_bundle.find { it.name.endsWith('.vcf.gz') }
  def common = common_bundle.find { it.name.endsWith('.vcf.gz') }
  if (!tumor || !normal) error 'Could not resolve matched pair'
  def ponVcf = params.skip_pon ? null : pon_bundle.find { it.name.endsWith('.vcf.gz') }
  def ponArg = params.skip_pon ? '' : "--panel-of-normals ${ponVcf}"
  """
  gatk Mutect2 -R ${ref} -I ${tumor[2]} -I ${normal[2]} \
    -tumor ${tumor[0]} -normal ${normal[0]} -L ${targets} \
    --germline-resource ${gnomad} ${ponArg} \
    --f1r2-tar-gz f1r2.tar.gz -O somatic.unfiltered.vcf.gz
  gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
  gatk GetPileupSummaries -I ${tumor[2]} -V ${common} -L ${targets} -O tumor.pileups.table
  gatk GetPileupSummaries -I ${normal[2]} -V ${common} -L ${targets} -O normal.pileups.table
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
  take: pairs; ref; targets; gnomad; common; pon
  main: CALL_AND_FILTER(pairs, ref, targets, gnomad, common, pon)
  emit: pass_vcf = CALL_AND_FILTER.out.pass
}
