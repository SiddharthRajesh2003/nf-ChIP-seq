#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MarkDuplicates {
    tag "Marking Duplicates for ${sample_id}"
    publishDir "${params.mkdp}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("*_mkdp.bam"), path ("*_marked_duplicates_metrics.txt")

    script:
    def sample_name = sample_id.split(':')[-1]  // Get everything after ':'
    """    
    # Create temp directory for GATK
    mkdir -p ./tmp
    
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_name}_mkdp.bam \\
        -M ${sample_name}_marked_duplicates_metrics.txt \\
        --VALIDATION_STRINGENCY LENIENT \\
        --TMP_DIR ./tmp \\
        --MAX_RECORDS_IN_RAM 5000000 \\
        --REMOVE_DUPLICATES true
    """
}