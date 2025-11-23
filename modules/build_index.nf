#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BuildIndex {
    tag 'Building Index'
    publishDir "${params.ref_dir}/", mode: 'copy'
    
    input:
    tuple path(fasta), path(gtf)

    output:
    tuple path("*.bt2"), path("*.fai")

    script:
    """
    bowtie2-build ${fasta} ${fasta.baseName}_index

    samtools faidx ${fasta}
    """
}