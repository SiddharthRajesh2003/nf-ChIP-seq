#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PeakCalling {
    tag "Peak Calling for ${sample_id}"
    publishDir "${params.peaks_dir}"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.xls"), path("*.bed"), path("*.broadPeak"), path("*.gappedPeak")

    script:
    def sample_name = bam.baseName.replace('_sorted.bam', '')
    """
    macs3 -t ${bam} --broad \
        -f BAM -g ${params.genome_size} \
        -n ${sample_name} \
        --qvalue 0.05 \
        --broad-cutoff 0.1 \
        --nomodel \
        --extsize 200 \
        --keep-dup all
    """
}