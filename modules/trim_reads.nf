#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TrimReads {
    tag "Trimming Reads for ${sample_id}"
    publishDir "${params.trimmed_reads}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${fastq.baseName}_trimmed.fq.gz")

    script:
    """
    trim_galore -q 30 --path_to_cutadapt ${params.cutadapt} --gzip ${fastq}
    """
}