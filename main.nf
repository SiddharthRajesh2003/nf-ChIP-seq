#!usr/bin/env nextflow

nextflow.enable.dsl = 2

params.base = "/N/project/Krolab/Siddharth/nf-chip-seq"
params.ref_dir = "${params.base}/reference"
params.ref = "${params.ref_dir}/human_reference.fa"
params.gtf = "${params.ref_dir}/human_reference.gtf"
params.input_csv = 'samples.csv'
params.threads = 8

// Directories for fastq files
params.fastq_dir = "${params.base}/fastq"
params.raw_reads = "${params.fastq_dir}/raw_reads"
params.trimmed_reads = "${params.fastq_dir}/trimmed_reads"
params.results = "${params.base}/results"


// Directories for Quality Control
params.qc_dir_before_trim = "${params.results}/fastqc/raw_reads"
params.qc_dir_after_trim = "${params.results}/fastqc/trimmed_reads"

// Directory for aligned BAM files
params.bam_dir = "${params.results}/bam"

params.aligned = "${params.bam_dir}/original"
params.skip_alignment = false
params.fallback_to_alignment = true

params.mkdp = "${params.bam_dir}/mkdp"
params.skip_mkdp = false
params.fallback_to_mkdp = true


params.filtered = "${params.bam_dir}/filtered"
params.stats_dir = "${params.results}/bam_stats"
params.skip_alignment = false

// Directories for Peak Calling
params.peaks_dir = "${params.results}/peaks"


def helpMessage() {
    log.info"""
    ===============================================
     DNA-seq Analysis Pipeline
    ===============================================
    
    Usage:
      nextflow run main.nf [options]
    
    Required parameters:
      --input           Path to input csv file
      --ref             Path to reference genome FASTA
      --gtf             Path to reference genome GTF
    
    Optional Parameters:
        --skip_alignment   Skip alignment if BAM files are present (default: false)
    """
}

include { QualityControl as QC_before_trim } from './modules/qc.nf'
include { QualityControl as QC_after_trim } from './modules/qc.nf'
include { TrimReads } from './modules/trim_reads.nf'
include { BuildIndex } from './modules/build_index.nf'
include { AlignReads } from './modules/align.nf'
include { AlignStats } from './modules/align_stats.nf'
include { MarkDuplicates } from './modules/mkdp.nf'


def shouldSkipAlignment() {
    if (!params.skip_alignment) {
        return false
    }
    
    def bamDir = file(params.aligned)
    if (!bamDir.exists()) {
        if (params.fallback_to_alignment) {
            log.warn "BAM directory ${params.aligned} does not exist! Will perform alignment instead."
            return false
        } else {
            error "BAM directory ${params.aligned} does not exist and fallback_to_alignment is disabled!"
        }
    }
    
    // Check if there are actually BAM files
    def bamFiles = bamDir.listFiles().findAll { f -> f.name.endsWith('.bam') }
    if (bamFiles.isEmpty()) {
        if (params.fallback_to_alignment) {
            log.warn "No BAM files found in ${params.aligned}! Will perform alignment instead."
            return false
        } else {
            error "No BAM files found in ${params.aligned} and fallback_to_alignment is disabled!"
        }
    }
    
    return true
}

def shouldSkipMKDP() {
    if (!params.skip_mkdp) {
        return false
    }

    def mkdpDir = file(params.mkdp)
    if (!mkdpDir.exists()) {
        if (params.fallback_to_mkdp) {
            log.warn "Mark Duplicates directory ${params.mkdp} does not exist! Will perform Mark Duplicates with Picard"
            return false
        } else {
            error "Mark Duplicates directory ${params.mkdp} does not exist and fallback_to_mkdp is disabled!"
        }
    }

    def mkdpFiles = mkdpDir.listFiles().findAll { f -> f.name.endsWith('_mkdp.bam') }
    if (mkdpFiles.isEmpty()) {
        if (params.fallback_to_mkdp) {
            log.warn "No Mark Duplicates BAM files found in ${params.mkdp}! Will perform Mark Duplicates with Picard"
        } else {
            error "No Mark Duplicates BAM files found in ${params.mkdp} and fallback_to_mkdp is disabled!"
        }
    }
    return true
}

workflow{
    channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map {row -> tuple(row.sample_id, file(row.fastq_path))}
        .set { samples_ch }

    QC_before_trim(samples_ch, params.qc_dir_before_trim)

    sample_id = samples_ch.map { sample_id, _fastq -> sample_id }

    trimmed_ch = TrimReads(samples_ch)

    QC_after_trim(trimmed_ch, params.qc_dir_after_trim)

    ref_index = BuildIndex(tuple(file(params.ref), file(params.gtf)))

    actuallySkipAlignment = shouldSkipAlignment()
    
    if (actuallySkipAlignment) {
        log.info "Skipping alignment as per request and using existing BAM files from ${params.bam_dir}"

        // Create a channel for existing BAM files
        bam_ch = channel
            .fromPath("${params.bam_dir}/*.bam", checkIfExists: true)
            .map { bam -> 
                log.info "Processing existing BAM file: ${bam.name} for sample ${sample_id}"
                tuple(sample_id, bam)
            }
    } else {
        log.info "Performing alignment with Bowtie2"
        
        // Check for existing BAMs that we can reuse
        existing_bam_ch = channel.empty()
        if (file(params.bam_dir).exists()) {
            existing_bam_ch = channel
                .fromPath("${params.bam_dir}/*.bam", checkIfExists: false)
                .map { bam ->
                    log.info "Reusing existing BAM file: ${bam.name} for sample ${sample_id}"
                    tuple(sample_id, bam)
                }
        }
        
        // Perform new alignment
        new_bam_ch = AlignReads(trimmed_ch, ref_index)
        
        // Use cached BAMs if available, otherwise use new BAMs
        bam_ch = existing_bam_ch.mix(new_bam_ch)
    }


    AlignStats(bam_ch)

    actuallySkipMKDP = shouldSkipMKDP()

    if (actuallySkipMKDP) {
        log.info "Skipping Mark Duplicates as per request and using existing Mark Duplicates BAM files from ${params.mkdp}"

        existing_mkdp_ch = channel
            .fromPath("${params.mkdp}/*_mkdp.bam", checkIfExists: true)
            .map { bam ->
                log.info "Processing existing Mark Duplicates BAM file: ${bam.Name}"
                tuple(sample_id, bam)
        }
    } else {
        log.info "Performing Mark Duplicates with Picard"
        // Checking for existing Mark Duplicates BAMs to reuse
        existing_mkdp_ch = channel.empty()
        if (file(params.mkdp).exists()) {
            existing_mkdp_ch = channel
                .fromPath("${params.mkdp}/*_mkdp.bam", checkIfExists: false)
                .map { bam ->
                    log.info "Reusing existing duplicates-marked BAM file: ${bam.baseName} for sample ${sample_id}"
                    tuple(sample_id, bam)
                }
        }

        // Perform Mark Duplicates
        mkdp_ch = MarkDuplicates(bam_ch)

        mkdp_ch = existing_mkdp_ch.mix(mkdp_ch)
    }

    
}