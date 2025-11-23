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
    # Filter out random/unplaced contigs and get top peaks from main chromosomes only
    grep -v "random\\|chrUn\\|_hap\\|_alt" ${peaks} | \\
        sort -k5 -nr | \\
        head -${params.top_peaks} > top_peaks.bed
    
    # Check if we have enough peaks
    num_peaks=\$(wc -l < top_peaks.bed)
    if [ "\$num_peaks" -lt 100 ]; then
        echo "WARNING: Only \$num_peaks peaks on main chromosomes, using all available"
    fi
    
    # Extend peaks to fixed width
    bedtools slop -i top_peaks.bed -g ${genome_fai} -b ${params.motif_window} | \\
        bedtools getfasta -fi ${genome_fasta} -bed  -fo raw_sequences.fa
    
    # Clean up FASTA headers - remove special characters that AME doesn't like
    # Change >chr1:1000-1200 to >seq_1
    awk '/^>/ {print ">seq_" ++i; next} {print}' raw_sequences.fa > peak_sequences.fa
    
    # Validate sequences
    num_seqs=\$(grep -c "^>" peak_sequences.fa)
    echo "Extracted \$num_seqs sequences from main chromosomes"
    
    # Run AME
    ame \\
        peak_sequences.fa \\
        ${motif_database}
    """
}