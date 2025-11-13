#!usr/env/bin nextflow 

nextflow.enable.dsl = 2

process AlignReads {
    tag "Aligning Reads for ${sample_id}"
    publishDir "${params.bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path ref_index

    output:
    tuple val(sample_id), path("*.bam"), path("*.bai")

    script:
    def sample_name = fastq.baseName.replace('.fq.gz', '')
    """
    bowtie2 -x ${ref_index} \
        -p ${params.threads} \
        -U ${fastq} | samtools view -bS -q 25 - | \
        samtools sort -@ ${task.cpus} -o ${sample_name}.bam

    samtools index ${sample_name}.bam
    """

}