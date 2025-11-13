#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process BuildIndex {
    tag 'Building Index', mode: 'copy'
    publishDir "${params.ref_dir}/", mode: 'copy'
    
    input:
    tuple path(fasta), path(gtf)

    output:
    path "${fasta.baseName}_index"

    script:
    """
    bowtie2-build ${fasta} ${fasta.baseName}_index
    """
}