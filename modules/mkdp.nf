#!usr/env/bin nextflow

nextflow.enable.dsl = 2

process MarkDuplicate {
    tag "Marking Duplicates for ${sample_id}"
    publishDir "${params.bam_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*_marked.bam"), path("*_marked.bai")
    path "*_marked_duplicates_metrics.txt"

    script:
    def sample_name = bam.baseName.replace('.bam', '')
    """
    # Set JAVA memory options (use ~75GB of your 100GB allocation)
    export JAVA_OPTS="-Xmx75g -XX:ParallelGCThreads=4"
    
    # Create temp directory for GATK
    mkdir -p ./tmp
    
    gatk --java-options "-Xmx75g -XX:ParallelGCThreads=4" MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_name}_marked_duplicates.bam \\
        -M ${sample_name}_marked_duplicates_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT \\
        --TMP_DIR ./tmp \\
        --MAX_RECORDS_IN_RAM 5000000
    """
}