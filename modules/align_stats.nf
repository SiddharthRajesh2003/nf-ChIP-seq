#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AlignStats {
    tag "Alignment Statistics for ${sample_id}"
    publishDir "${params.stats_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_alignment_stats.txt"

    script:
    """
    samtools flagstat ${bam} > ${bam.baseName}_flagstats.txt

    samtools stats ${bam} > ${bam.baseName}_stats.txt

    samtools idxstats ${bam} > ${bam.baseName}_idxstats.txt
    """
}