#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process EmptyProcess {
    tag "Empty Process"
    publishDir "${params.results}/empty_process", mode: 'copy'

    input: 
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_empty_output.txt"

    script:
    """
    echo "This is an empty process for ${sample_id}" > ${sample_id}_empty_output.txt
    """
}