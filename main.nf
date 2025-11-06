#!usr/bin/env nextflow

nextflow.enable.dsl = 2

params.base = "/N/project/Krolab/Siddharth/Pipelines/CHIP-seq"
params.ref_dir = "${params.base}/reference"
params.input_csv = 'samples.csv'
params.fastq_dir = "${params.base}/fastq"
params.qc_dir_before_trim = "${params.base}/results/fastqc/raw_reads"
params.qc_dir_after_trim = "${params.base}/results/fastqc/trimmed_reads"

include { QualityControl as QC_before_trim } from './modules/qc.nf'
include { QualityControl as QC_after_trim } from './modules/qc.nf'
include { TrimReads } from './modules/trim_reads.nf'

workflow{
    channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map {row -> tuple(row.sample_id, file(row.fastq_path))}
        .set { samples_ch }

    QC_before_trim(samples_ch, params.qc_dir_before_trim)

}