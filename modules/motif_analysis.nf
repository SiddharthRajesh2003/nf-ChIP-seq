#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process MotifAnalysis {
    tag "Motif Analysis for ${sample_id}"
    publishDir "${params.motif_dir}/homer", mode: 'copy'

    input:
    tuple val(sample_id), path(peaks)
    path ref

    output:
    tuple val(sample_id), path("*_motifs")

    script:
    """
    cut -f1-3 ${peaks} > peaks.bed

    findMotifsGenome.pl \\
        peaks.bed \
        ${ref} \
        ${sample_id}_homer_motifs \
        -size 200 \
        -mask \
        -p ${task.cpus}
    """
}