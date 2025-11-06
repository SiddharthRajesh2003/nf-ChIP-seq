#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TrimReads {
    tag "Trimming Reads for ${sample_id}"
    publishDir "${params.fastq_dir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz")

    script:
    """
    trim_galore --output_dir . ${fastq}
    """
}