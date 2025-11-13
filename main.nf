#!usr/bin/env nextflow

nextflow.enable.dsl = 2

params.base = "/N/project/Krolab/Siddharth/nf-chip-seq"
params.ref_dir = "${params.base}/reference"
params.ref = "${params.ref_dir}/human_reference.fa"
params.gtf = "${params.ref_dir}/human_reference.gtf"
params.input_csv = 'samples.csv'
params.fastq_dir = "${params.base}/fastq"
params.results = "${params.base}/results"
params.qc_dir_before_trim = "${params.results}/fastqc/raw_reads"
params.qc_dir_after_trim = "${params.results}/fastqc/trimmed_reads"
params.bam_dir = "${params.results}/bam"

params.skip_alignment = false


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

workflow{
    channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map {row -> tuple(row.sample_id, file(row.fastq_path))}
        .set { samples_ch }

    QC_before_trim(samples_ch, params.qc_dir_before_trim)

    trimmed_ch = TrimReads(samples_ch)

    QC_after_trim(trimmed_ch, params.qc_dir_after_trim)

    ref_index = BuildIndex(tuple(file(params.ref), file(params.gtf)))


    actuallySkipAlignment = shouldSkipAlignment()
    
    if (actuallySkipAlignment){
        log.info "Skipping alignment as per request and using existing BAM files from ${params.bam_dir}"

        //Create a channel for existing BAM files
        existing_bam_ch = channel
            .fromPath("${params.bam_dir}/*.bam", checkIfExists: true)
            .map { bam -> 
                log.info "Processing existing BAM file: ${bam.Name}"
                def bai_paths = [
                    file("${bam}.bai"),
                    file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")  // alternative .bai
                ]

                def bai = bai_paths.find { idx -> idx.exists() }
                if (bai) {
                    log.info "Found index for ${bam.name}: ${bai.name}"
                    return tuple(bam, bai)
                } else{
                    log.error "Index file not found for ${bam.name}"
                    log.error "Looked for: ${bai_paths.collect { p -> p.toString() }.join(', ')}"
                    error "Index file not found for ${bam.name}. Please ensure BAM files are properly indexed."
                }
            }
        bam_ch = existing_bam_ch
    } else {
        log.info "Performing alignment with Bowtie2"
        // Check for existing BAMs that we can reuse

        existing_bam_ch = channel.empty()
        if (file(params.bam_dir).exists()) {
            existing_bam_ch = channel
            .fromPath("${params.bam_dir}/*.bam", checkIfExists: true)
            .map { bam ->
                def bai_paths = [
                    file("${bam}.bai"),
                    file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")  // alternative .bai
                ]
                def bai = bai_paths.find { idx -> idx.exists() }
                    if (bai) {
                        log.info "Found existing indexed BAM: ${bam.name}"
                        return tuple(bam, bai)
                    } else {
                        log.warn "Found BAM without index: ${bam.name} - will skip"
                        return null
                    }
            }
            .filter { item -> item!=null}
        }
        // Perform new alignment
        new_bam_ch = AlignReads(trimmed_ch, ref_index)
        
        // Use cached BAMs if available, otherwise use new BAMs
        bam_ch = existing_bam_ch.mix(new_bam_ch)
    }


    AlignStats(bam_ch)


}