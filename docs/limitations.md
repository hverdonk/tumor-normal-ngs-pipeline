# Limitations

- Checked-in summaries are schemas, not results; the public full WES run has not
  been executed here.
- Cell-line drift, library/capture differences, and reference-bundle differences
  can prevent exact reproduction of published performance.
- The SEQC2 confident BED is used as both capture and benchmark territory. This
  keeps evaluation regions identical by construction, but does not model the
  callable territory of a different capture kit.
- Exact normalized matching is less representation-aware than hap.py/vcfeval.
- Query VAF/depth stratification cannot assign missed truth variants to a query
  depth bin. A complete depth analysis should annotate every truth locus from the
  tumor BAM before stratifying FN.
- Mutect2 does not model all copy-number/purity effects, and this scope excludes
  CNVs, structural variants, fusions, MSI, and germline reporting.
- VEP cache content controls ClinVar/frequency availability. SIFT and PolyPhen are
  predictions, not evidence of pathogenicity. CADD/AlphaMissense require separately
  versioned plugin data.
- WES leaves poorly captured and difficult regions incompletely assessed.
