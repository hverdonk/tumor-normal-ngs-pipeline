# Methods

## Study design and provenance

HCC1395 (`SRR7890844`) is the tumor and HCC1395BL (`SRR7890845`) its matched
normal. Both are paired-end WES. The workflow requires local checksummed copies
of reads and resources. Reference FASTA, target BED, population VCFs, PoN, VEP
cache, truth VCF, and truth BED must share the same GRCh38 sequence dictionary.

## Read QC and cleanup

FastQC 0.12.1 runs before and after cleanup; MultiQC 1.35 summarizes per-base quality,
adapter content, read length, duplication, and GC distribution. `fastp 1.3.6`
auto-detects paired-end adapters and trims poly-G tails (a known two-color
artifact). Quality filtering is deliberately disabled: hard quality clipping can
remove real alternate-supporting bases and distort low-VAF sensitivity. Reads
shorter than 30 bp after adapter/poly-G cleanup are removed because they map
ambiguously. The fastp JSON contains input/output counts and retained percent.

## Alignment and QC

BWA-MEM3 0.6.0 receives explicit ID, sample, library, platform, and platform-unit
read-group fields. Alignments are coordinate sorted, duplicate-marked (not
removed), and indexed. Samtools flagstat/stats report total, mapped, and properly
paired reads. Picard reports duplicate fraction and insert sizes. Mosdepth over
the capture BED reports mean target coverage and bases at >=20x, >=50x, >=100x.
Tumor and normal rows must be compared directly because normal undercoverage can
impair germline/artifact rejection.

## Somatic calling and filtering

GATK 4.6.2.0 Mutect2 runs jointly with one tumor and one normal, a GRCh38 gnomAD
AF-only resource, and normally a panel of normals. `--skip_pon true` is explicit
and recorded when a compatible PoN is unavailable. The matched normal supplies
site-specific evidence for rare germline variation and sample-specific artifacts.

Tumor and normal pileups at common biallelic sites feed contamination estimation.
F1R2 counts feed the read-orientation prior. FilterMutectCalls consumes both,
after which only PASS variants are retained, split, left-aligned, and normalized.

## Annotation and prioritization

Offline VEP 116.0 annotates gene/transcript, consequence, HGVS coding/protein
change, canonical status, existing IDs, ClinVar significance present in the
cache, gnomAD frequency, SIFT, and PolyPhen. CADD/AlphaMissense may be added only
when licensed/downloaded plugin data and matching versions are documented.
Prioritization uses a separately supplied, versioned cancer-gene symbol list and
should require protein alteration, rarity, cancer-gene membership,
pathogenic ClinVar evidence or damaging predictions, and sufficient tumor depth
and VAF. This is functional prioritization, not clinical interpretation.

## Benchmarking

Calls and SEQC2 v1.2 truth are split and left-normalized against the same FASTA,
then restricted to `High-Confidence_Regions_v1.2.bed`. Exact CHROM/POS/REF/ALT
matches define TP; call-only variants FP; truth-only variants FN. SNVs and indels
are reported separately with precision = TP/(TP+FP), recall = TP/(TP+FN), and
F1 = 2PR/(P+R). VAF and query-depth bands are emitted separately. For a formal
submission, use hap.py or vcfeval as a second, haplotype-aware evaluator.
