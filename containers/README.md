# Container inventory

Build the mixed-tool image before running: `docker build -t hcc1395-wes:1.0.0
containers/`. The workflow also uses these pinned BioContainers/official tags:

The root `environment.yml` contains the Nextflow runner plus the local
utilities needed for `scripts/fetch_sra.sh` and `scripts/make_test_data.sh`.
Tool dependencies are still isolated by process because the latest VEP, fastp,
mosdepth, and HTSlib Bioconda builds currently have mutually exclusive
transitive HTSlib/libdeflate constraints when placed in one Conda environment.

| Purpose | Image |
|---|---|
| Project-image base | `mambaorg/micromamba:2.8.1` |
| FastQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| MultiQC | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1` |
| fastp | `quay.io/biocontainers/fastp:1.3.6--h43da1c4_0` |
| BWA-MEM3 | project image, Conda `bwa-mem3=0.6.0` |
| samtools/bcftools | project image, Conda `samtools=1.24`, `bcftools=1.24` |
| GATK | project image, Conda `gatk4=4.6.2.0` |
| Picard | project image, Conda `picard=3.4.0` |
| mosdepth | `quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_0` |
| VEP | `ensemblorg/ensembl-vep:release_116.0` |
| Mixed alignment/calling/benchmark steps | locally built `hcc1395-wes:1.0.0` |

Tags are pinned for reproducibility; production mirrors should pin resolved image
digests too. Nextflow records the actual image used in its trace.
