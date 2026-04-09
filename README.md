# AmpSeq Indel Analysis Pipeline

A Nextflow (DSL2) pipeline for amplicon sequencing analysis, designed to quantify
CRISPR-induced insertions and deletions (indels) at targeted genomic loci.

## Overview

This pipeline takes paired-end FASTQ files from amplicon sequencing experiments
and produces per-sample indel counts, integration frequencies, and summary plots.

Note: If your reference genome or project directory are not located within your home directory ($HOME), you will need to manually mount them in nextflow.config. See the Configuration section below for details.

## Pipeline Steps

| Step | Tool | Description |
|------|------|-------------|
| 1 | FastQC | Raw read quality control |
| 2 | Trimmomatic | Adapter trimming (PE mode) |
| 3 | FLASH2 | Read merging |
| 4 | BWA-MEM | Alignment to reference genome |
| 5 | Samtools | BAM sorting and indexing |
| 6 | Custom Python | Per-site indel and integration counting |
| 7 | Pandas / Seaborn | Combined Excel report and read depth plot |
| 8 | Custom Python | Per-sample indel visualization |

## Requirements

- Nextflow >= 22.0
- Docker (recommended) or Conda
- Reference genome FASTA with BWA index (e.g. hg38)
- BED file defining target loci

### System Requirements
- Minimum 16GB RAM recommended (BWA-MEM hg38 index requires ~8GB)
- Docker images are built for linux/amd64 and available on Docker Hub

## Docker Images

All pipeline dependencies are containerized and publicly available on Docker Hub:

| Image | Processes |
|-------|-----------|
| `recouto23/ampseq:1.0.0` | FastQC, Trimmomatic, FLASH2, BWA, Samtools, Python QC |
| `recouto23/tsai_indels:1.0.0` | Indel counting (pysam, scikit-bio, pandas) |
| `recouto23/indel_plot:1.0.0` | Indel visualization (matplotlib, seaborn) |

## Usage

```bash
nextflow run main.nf \
  --fastqDir path/to/fastqs/ \
  --outDir Results/ \
  --reference /path/to/hg38.fa \
  --bed targets.bed \
  --tag MyExperiment \
  --threads 4
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastqDir` | `fastq/` | Directory containing input FASTQ files |
| `--outDir` | `Results/` | Output directory |
| `--reference` | *(required)* | Path to reference genome FASTA (must be BWA-indexed) |
| `--bed` | `bed.txt` | BED file of target sites |
| `--tag` | `26BCPXXX` | Experiment tag for output file naming |
| `--threads` | `4` | Threads per process |

## Configuration

Docker is enabled by default in `nextflow.config`. To ensure the reference genome
is accessible inside containers, the user home directory is mounted automatically:

```groovy
docker {
    enabled = true
    runOptions = "-v ${System.getenv('HOME')}:${System.getenv('HOME')}"
}
```

## Outputs

```
Results/
├── 00_fastQC/            # FastQC HTML reports
├── 01_trimmed/           # Trimmed FASTQs
├── 02_merged/            # FLASH2-merged reads
├── 03_mapped/            # Sorted, indexed BAMs
└── 04_indels/
    ├── *_countIndelsReport.tsv
    ├── *_CombinedIndelReport.xlsx
    ├── *_ReadPlot.png
    └── 00_indelPlots/
```

## Expected Input FASTQ Naming

Paired-end FASTQ files should follow Illumina naming conventions:

```
SampleName_S1_L001_R1_001.fastq.gz
SampleName_S1_L001_R2_001.fastq.gz
```
