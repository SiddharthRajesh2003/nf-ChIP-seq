# Nextflow Pipeline for ChIP-seq Data Analysis

A modular Nextflow pipeline for ChIP-seq data analysis, including download, QC, alignment, peak calling, motif analysis, and annotation.

## Table of Contents
- [Overview](#overview)
- Folder Structure
- Getting Started
- Pipeline Modules
- Input Files
- Reference Files
- Results
- Scripts
- Usage
- License

### Overview

This pipeline automates the analysis of ChIP-seq data using Nextflow. It supports SRA download, quality control, read trimming, alignment, peak calling, motif discovery/enrichment, and peak annotation.

### Folder Structure

```
nf-ChIP-seq/
├── main.nf                # Main Nextflow workflow
├── nextflow.config        # Nextflow configuration
├── env.yaml               # Conda environment file
├── modules/               # Nextflow process modules
│   ├── align.nf
│   ├── align_stats.nf
│   ├── build_index.nf
│   ├── filter_bam.nf
│   ├── mkdp.nf
│   ├── motif_analysis.nf
│   ├── motif_discovery.nf
│   ├── motif_enrichment.nf
│   ├── peakcall.nf
│   ├── peak_annotation.nf
│   ├── qc.nf
│   └── trim_reads.nf
├── metadata/              # Sample metadata
│   └── chipseq_metadata.xlsx
├── samples.csv            # Sample sheet for pipeline input
├── generate_sample_csv.ipynb # Notebook for generating sample CSV
├── run_pipeline.sh        # Script to run the pipeline
├── dag-Entry.svg          # Pipeline DAG visualization
├── README.md              # This file
└── .gitignore
```

## Getting Starteed

1. Clone the repository

```bash
git clone https://github.com/SiddharthRajesh2003/nf-ChIP-seq.git
cd nf-ChIP-seq
```

2. Install Nextflow and conda
   
   Follow instructions at Nextflow and Conda

3. Set up the environment

```bash
conda env create -f env.yaml
conda activate chip_env
```

## Pipeline Modules

Each module in modules is a Nextflow process for a specific step:

- [qc.nf](modules/qc.nf): Quality Control (FastQC)
- [trim_reads.nf](modules/trim_reads.nf): Read trimming (Trimgalore)
- [build_index.nf](modules/build_index.nf): Building Bowtie and Samtools index
- [align.nf](modules/align.nf): Alignment using Bowtie2
- [align_stats.nf](modules/align_stats.nf): Alignment statistics with Samtools
- [mkdp.nf](modules/mkdp.nf): PCR Duplicate marking with GATK
- [filter_bam.nf](modules/filter_bam.nf): Filter out blacklisted regions with ngsutilsj
- [peakcall.nf](modules/peakcall.nf): Peak Calling with MACS3.
- [peak_annotation.nf](modules/peak_annotation.nf): Peak Annotation using HOMER
- [motif_analysis.nf](modules/motif_analysis.nf): Motif Analysis using HOMER
- [motif_discovery.nf](modules/motif_discovery.nf): Motif Discovery with MEME if needed
- [motif_enrichment.nf](modules/motif_enrichment.nf): Motif Enrichment analysis using MEME if needed
  

### Input Files

```bash
# samples.csv

sample_id,fastq
SCZ_SRR1234567,path/to/SRR1234567.fastq
```

### Reference files
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P reference/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz -P reference/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -P reference/
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt -P reference/     # Motif DB for Homer/MEME-Chip

# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz
# Download the above file to your local and upload to HPC using scp
# scp your/local/working/directory/hg19-blacklist.v2.bed.gz your.username@gmail.com:/your/HPC/working/directory

gunzip reference/*.gz
```

**Note**: Always use the reference genome and annotation files matching your experimental design.

## Results

Results are organized in the results folder:
- fastqc/ : Quality Control Reports
- bam/ : Aligned Reads
- bam_stats/ : Alignment Statistics
- peaks/ : Results from peak calling
- motif_analysis/ : Motif Analysis outputs
- pipeline_info/ : Pipeline logs and metadata

## Scripts

- [run_pipeline.sh](run_pipeline.sh): Example script to launch the pipeline on HPC.
- [Download samples](#samples-used-for-this-analysis): Samples used for this pipeline.

## Usage

Run the pipeline with:

```bash
nextflow run main.nf \
    --input_csv samples.csv \
    -profile slurm \
    -resume \
    -with-timeline ${OUTPUT_DIR}/pipeline_info/execution_timeline_${SLURM_JOB_ID}.html \
    -with-report ${OUTPUT_DIR}/pipeline_info/execution_report_${SLURM_JOB_ID}.html \
    -with-trace ${OUTPUT_DIR}/pipeline_info/execution_trace_${SLURM_JOB_ID}.txt \
    -with-dag ${OUTPUT_DIR}/pipeline_info/pipeline_dag_${SLURM_JOB_ID}.svg
```

Customize parameters in [nextflow.config](nextflow.config) or via command-line options.

##### Samples Used for this analysis

```bash
base=your/HPC/working/directory
srr="SRR21948648 SRR21948649 SRR21948650 SRR21948651 SRR21948652 SRR21948653 SRR21948654 SRR21948655 SRR21948656 SRR21948657 SRR21948658 SRR21948659 SRR21948660 SRR21948661 SRR21948662 SRR21948663 SRR21948664 SRR21948665 SRR21948666 SRR21948667 SRR21948668 SRR21948669 SRR21948670 SRR21948671 SRR21948672 SRR21948673 SRR21948674 SRR21948675 SRR21948676 SRR21948677"
cd $base

mkdir -p $base/fastq/raw_reads
module load sra-toolkit

for s in $srr; do
        prefetch $s -O sra_files
        fasterq-dump sra_files/$s -O fastq/raw_reads
done
```

## Acknowledgements

- SRA Toolkit
- Bowtie2
- FastQC
- MEME Suite
- HOMER
- JASPAR
  
## License

See LICENSE for details