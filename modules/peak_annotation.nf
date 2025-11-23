#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AnnotatePeaks {
    tag "Peak Annotation for ${sample_id}"
    publishDir "${params.annotation_dir}", mode: "copy"

    input:
    tuple val(sample_id), path(peaks)
    path(gtf)

    output:
    tuple path("${sample_id}_annotated.bed"), path("*_stats.txt")

    script:
    """
    annotatePeaks.pl ${peaks} ${params.ref} \
        -gtf ${gtf} \
        -annStats ${sample_id}_annotation_stats.txt \
        -CpG \
        > ${sample_id}_annotated.bed
    """
}