#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process BuildIndex {
    tag 'Building Index'
    publishDir "${params.ref_dir}/", mode: 'copy'
    
    container 'biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1'

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