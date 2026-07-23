# Interpretation guide

Begin with sample identity and contamination, then coverage parity, then call
metrics. A high apparent precision with poor recall can simply reflect shallow
capture; a high call count can reflect contamination or mapping artifacts.

Review representative TP, FP, and FN loci in IGV with both duplicate-marked BAMs,
the reference, SEQC2 confident BED, and normalized VCFs loaded. Color reads by
strand and inspect mapping/base qualities, clipping, nearby indels, read-end
position, and normal evidence. Useful examples include one high-VAF SNV TP, one
well-supported indel TP, a low-VAF FN, an ambiguous-mapping FP, and a locus with
low normal depth. Record coordinates and screenshots outside the repository if
they contain controlled data.

Likely failure modes include low tumor VAF, insufficient local tumor or normal
coverage, paralogous/repetitive alignment, orientation or strand bias, oxidative
damage, residual germline signal, and non-equivalent indel representation.

The prioritized VEP table is hypothesis generation. ClinVar conflicts, transcript
choice, population ancestry, tumor purity/copy number, and assay validation all
matter. Do not label a variant clinically actionable from this workflow alone.
