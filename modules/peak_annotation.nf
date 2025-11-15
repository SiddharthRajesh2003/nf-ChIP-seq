#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AnnotatePeaks {
    tag "Peak Annotation for ${sample_id}"
    publishDir "${params.annotation_dir}", mode: "copy"

    input:
    tuple val(sample_id), path(peaks), path(annotation_db)

    output:
    path "${peaks.baseName}_annotated.bed"

    script:
    """
    
    """
}