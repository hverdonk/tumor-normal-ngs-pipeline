# Limitations

- Checked-in summaries are schemas, not results; the public full WES run has not
  been executed here.
- Cell-line drift, library/capture differences, and reference-bundle differences
  can prevent exact reproduction of published performance.
- Capture targets limit callable territory; benchmarking only the SEQC2 confident
  BED without intersecting assay-callable regions can penalize uncaptured truth
  loci. A production analysis should benchmark their documented intersection.
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

