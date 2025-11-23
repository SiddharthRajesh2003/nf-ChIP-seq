#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PeakCalling {
    tag "Peak Calling for ${sample_id}"
    publishDir "${params.peaks_dir}", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.xls"), path("*.broadPeak"), path("*.gappedPeak")

    script:    
    """
    macs3 callpeak \
        -t ${bam} --broad \
        -f BAM -g ${params.genome_size} \
        -n ${sample_id} \
        --qvalue 0.05 \
        --broad-cutoff 0.1 \
        --nomodel \
        --extsize 200 \
        --keep-dup all
    """
}