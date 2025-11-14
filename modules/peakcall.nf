#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PeakCalling {
    tag "Peak Calling for ${sample_id}"
    publishDir "${params.peaks_dir}"

    input:
    tuple val(sample_id), path(bam)
    path ref_index

    output:
    path "${sample_id}.broadPeak"

    script:
    """
    mac2 -t ${params.threads} -broad 
    """
}