#!usr/env/bin nextflow 

nextflow.enable.dsl = 2

process AlignReads {
    tag "Aligning Reads for ${sample_id}"
    publishDir "${params.bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    val ref_index

    output:
    tuple val(sample_id), path()

    script:
    """
    bowtie2 -x ${ref_index} -
    """

}