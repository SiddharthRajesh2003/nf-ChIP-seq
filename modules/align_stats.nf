#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AlignStats {
    tag "Alignment Statistics for ${sample_id}"
    publishDir "${out_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val out_dir

    output:
    tuple path("*_flagstats.txt"), path("*_stats.txt"), path("*_idxstats.txt")

    script:
    def sample_name = sample_id.split(':')[-1]
    """
    samtools flagstat ${bam} > ${sample_name}_flagstats.txt

    samtools stats ${bam} > ${sample_name}_stats.txt

    samtools idxstats ${bam} > ${sample_name}_idxstats.txt
    """
}