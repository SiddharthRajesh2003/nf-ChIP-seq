#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BuildIndex {
    tag 'Building Index'
    publishDir "${params.ref_dir}/", mode: 'copy'
    
    container 'quay.io/biocontainers/mulled-v2-258922dd32ede57ed7115a4d4c28f2094cb21eec:4913fcfa10f2f2aff7ece136756af6fe46e6be0d-0'

    input:
    tuple path(fasta), path(gtf)

    output:
    tuple path("*.bt2"), path("*.fai")

    script:
    """
    module load samtools
    
    bowtie2-build ${fasta} ${fasta.baseName}_index

    samtools faidx ${fasta}
    """
}