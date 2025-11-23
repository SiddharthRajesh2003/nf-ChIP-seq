#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MotifDiscovery {
    tag "Motif Discovery for ${sample_id}"
    publishDir "${params.motif_dir}/meme", mode: 'copy'

    input:
    tuple val(sample_id), path(peaks)
    path genome_fasta
    path genome_fai

    output:
    tuple val(sample_id), path("${sample_id}_motifs")
    
    script:
    """
    # Get peak summits (use top 500 peaks by score)
    sort -k5 -nr ${peaks} | head -${params.top_peaks} > top_peaks.bed

    # Extend peaks to fixed width around summit
    bedtools slop -i top_peaks.bed -g ${genome_fai} -b ${params.motif_window} | \
        bedtools getfasta -fi ${genome_fasta} -bed - -fo peak_sequences.fa

    # Run MEME-ChIP for motif discovery
    meme-chip \\
        -oc ${sample_id}_motifs \\
        -db ${params.motif_database} \\
        -meme-minw ${params.meme_minw} \\
        -meme-maxw ${params.meme_maxw} \\
        -meme-nmotifs ${params.meme_nmotifs} \\
        -meme-p ${task.cpus} \\
        -dreme-e 0.05 \\
        -centrimo-local \\
        -centrimo-score 5.0 \\
        peak_sequences.fa
    """
}