#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process QualityControl {
    tag "Quality Control on ${sample_id}"
    publishDir "${qc_dir}", mode: "copy"

    input:
    tuple val(sample_id), path(fastq)
    val qc_dir

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${fastq}
    """
}