#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AME {
    tag "AME Analysis for ${sample_id}"
    publishDir "${params.motif_dir}/ame", mode: 'copy'
    
    input:
    tuple val(sample_id), path(peaks)
    path genome_fasta
    path genome_fai
    path motif_database
    
    output:
    tuple val(sample_id), path("${sample_id}_ame")
    
    script:
    """
    # Get top peaks
    sort -k5 -nr ${peaks} | head -${params.top_peaks} > top_peaks.bed
    
    # Extend peaks to fixed width
    bedtools slop -i top_peaks.bed -g ${genome_fai} -b ${params.motif_window} | \
        bedtools getfasta -fi ${genome_fasta} -bed - -fo peak_sequences.fa
    
    # Run AME for known motif enrichment
    ame \\
        --oc ${sample_id}_ame \\
        --control --shuffle-- \\
        peak_sequences.fa \\
        ${motif_database}
    """
}