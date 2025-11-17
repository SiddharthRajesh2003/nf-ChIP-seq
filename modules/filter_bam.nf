#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FilterBAM {
    tag "Filtering BAM for ${sample_id}"
    publishDir "${params.filtered}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*_filtered.bam"), path("*_filtered.bai")

    script:
    def sample_name = sample_id.split(':')[-1]  // Get everything after ':'
    """
    ${params.base}/Apps/ngsutilsj bam-filter --mapped --no-qcfail --tag-min MAPQ:30 \
        ${bam} ${sample_name}_filtered0.bam

    ${params.base}/Apps/ngsutilsj bam-filter --bed-exclude ${params.blacklist} \
        ${sample_name}_filtered0.bam ${sample_name}_filtered.bam
    
    samtools index ${sample_name}_filtered.bam -o ${sample_name}_filtered.bai
    """
}