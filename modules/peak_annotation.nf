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
    annotatePeaks.pl ${peaks} ${params.ref} \
        -gtf ${params.gtf} \
        -annStats ${sample_id}_annotation_stats.txt \
        -CpG \
        > ${peaks.baseName}_annotated.bed
    """
}