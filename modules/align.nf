#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AlignReads {
    tag "Aligning Reads for ${sample_id}"
    publishDir "${params.aligned}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path ref_index

    output:
    tuple val(sample_id), path("*.bam")

    script:
    def sample_name = sample_id.split(':')[-1]  // Get everything after ':'
    def index_base = ref_index[0].toString().replaceAll(/\.1\.bt2$/, '')
    """
    bowtie2 -x ${index_base} \
        -p ${params.threads} \
        -U ${fastq} | samtools view -bS -q 25 - | \
        samtools sort -@ ${task.cpus} -o ${sample_name}.bam
    """

}