# HCC1395 matched tumor-normal WES pipeline

Reproducible discovery, annotation, and SEQC2 benchmarking of somatic SNVs and
small indels from HCC1395 (`SRR7890844`) and HCC1395BL (`SRR7890845`). The
workflow is intentionally limited to research-use functional annotation; it is
not a clinical assay or interpretation.

## What this repository demonstrates

The DSL2 Nextflow workflow performs raw and cleaned FastQC, conservative `fastp`
cleanup, BWA-MEM3 alignment with read groups, coordinate sorting and duplicate
marking, samtools/mosdepth/Picard QC, matched-normal Mutect2, contamination and
read-orientation filtering, normalized PASS extraction, VEP annotation, and
benchmarking against SEQC2 v1.2 inside its high-confidence regions. MultiQC
collects the QC outputs. See [methods](docs/methods.md), [interpretation](docs/interpretation.md),
and [limitations](docs/limitations.md).

No biological metrics are claimed in the checked-in summaries: the full public
WES analysis is compute- and storage-intensive and has not been run in this
repository. Empty, schema-valid summaries are labeled `NOT_RUN`, preventing
example values from being mistaken for observations.

## Inputs and compatibility

Use the exact tumor/normal roles in `config/samplesheet.csv`. The download
script retrieves ENA's gzip-compressed paired FASTQs for the listed public SRA
accessions, verifies ENA's published MD5 checksums, and records SHA-256
checksums for this run:

```bash
scripts/fetch_sra.sh config/samplesheet.csv data/fastq checksums/fastq.sha256
```

Obtain the reference bundle, gnomAD AF-only VCF, common-sites VCF, optional panel
of normals, VEP cache, and SEQC2 v1.2 truth VCF/BED described in
`config/parameters.yaml`. The SEQC2 confident BED is the fixed interval source
for coverage, calling, and benchmarking; a separate capture BED cannot be
configured. These large resources are deliberately not mirrored.
Every FASTA/VCF/BED must be GRCh38 and use an identical sequence dictionary
(`chr1` versus `1` is not interchangeable). Put local paths in a private params
file copied from `config/local.example.yaml`, then generate and verify a manifest:

```bash
cp config/local.example.yaml config/local.yaml
scripts/check_inputs.sh config/samplesheet.csv config/local.yaml
docker build -t hcc1395-wes:1.0.0 containers/
nextflow run workflow/main.nf -profile docker -params-file config/local.yaml -resume
```

The official SEQC2 publication reports the benchmark data at NCBI's
`ReferenceSamples/seqc/Somatic_Mutation_WG` archive. Treat v1.2 VCF and
`High-Confidence_Regions_v1.2.bed` as a matched pair. URLs can move, so the
workflow accepts local, checksum-verified files rather than downloading silently.

## Smoke test

The offline test creates tiny synthetic paired FASTQs, a toy reference, intervals,
and truth VCF. It tests preprocessing/alignment wiring, not biological performance.
Create the pinned Conda environment first (the fixture builder uses bgzip,
tabix, samtools, and GATK):

```bash
scripts/make_test_data.sh test-data
nextflow run workflow/test.nf -profile test,docker
```

Full runs normally need hundreds of GB of scratch space, a populated VEP cache,
and substantial CPU/RAM. Override process resources in a site profile rather
than editing the workflow.

## Outputs

- `results/qc/`: FastQC/MultiQC, fastp, alignment, insert-size, and mosdepth files
- `results/variants/`: unfiltered, filtered, normalized PASS, and annotated VCFs
- `results/variants/prioritized_variants.tsv`: research-only VEP shortlist
- `results/benchmark/`: overall and depth/VAF-stratified SNV/indel metrics
- `results-summary/`: compact hiring-manager-facing tables and figures
- `results/provenance/`: resolved parameters, software/container metadata, and checksums

Run status and exact commands are recorded by Nextflow (`timeline`, `report`,
`trace`, and `dag`). The pipeline stops on missing resources; a panel of normals
may be disabled explicitly, but gnomAD and common-sites resources are required
for a full run.

## Reproducibility notes

The mixed-tool project container and its Conda packages are version-pinned. For regulated or archival
use, mirror images internally and replace tags with immutable OCI digests. Input
content hashes matter more than filenames; retain `checksums/input.sha256` with
each run. The public accessions identify read data but do not guarantee that an
upstream archive will never revise a file.

## References

- Fang et al., *Nature Biotechnology* (2021), SEQC2 reference samples and truth set.
- [GATK Mutect2 documentation](https://gatk.broadinstitute.org/hc/en-us/articles/30332058799003)
- [Ensembl VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)
