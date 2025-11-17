#!usr/bin/env nextflow

nextflow.enable.dsl = 2

params.base = "/N/project/Krolab/Siddharth/nf-chip-seq"
params.ref_dir = "${params.base}/reference"
params.ref = "${params.ref_dir}/hg19.fa"
params.gtf = "${params.ref_dir}/hg19.knownGene.gtf"
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

// Parameters for Trimming
params.cutadapt = "/N/u/sidrajes/Quartz/.conda/envs/cutadapt/bin/cutadapt"

// Directory for aligned BAM files
params.bam_dir = "${params.results}/bam"

params.aligned = "${params.bam_dir}/original"
params.skip_alignment = false
params.fallback_to_alignment = true
params.stats_dir = "${params.results}/bam_stats"
params.stats_out_dir = "${params.stats_dir}/original"


// Parameters for Marking PCR Duplicates in the BAM files
params.mkdp = "${params.bam_dir}/mkdp"
params.skip_mkdp = false
params.fallback_to_mkdp = true
params.mkdp_stats = "${params.stats_dir}/mkdp"


// Parameters for filtering random noise regions from the BAM files
params.filtered = "${params.bam_dir}/filtered"
params.skip_filtering = false
params.blacklist = "${params.ref_dir}/hg19-blacklist.v2.bed"
params.filtered_stats = "${params.stats_dir}/filtered"

// Directories for Peak Calling
params.peaks_dir = "${params.results}/peaks/broad"
params.annotation_dir = "${params.results}/peaks/annotated"
params.genome_size = 'hs'


def helpMessage() {
    log.info"""
    ===============================================
            ChIP-seq Analysis Pipeline
    ===============================================
    
    Usage:
        nextflow run main.nf [options]
    
    Required parameters:
        --input           Path to input csv file
        --ref             Path to reference genome FASTA
        --gtf             Path to reference genome GTF
        --blacklist       Path to blacklist BED file
    
    Optional Parameters:
        --skip_alignment   Skip alignment if BAM files are present (default: false)
        --skip_mkdp        Skip Mark Duplicates if BAM files are present (default: false)
        --skip_filtering   Skip BAM filtering if filtered BAM files are present (default: false)
    """
}

include { QualityControl as QC_before_trim } from './modules/qc.nf'
include { QualityControl as QC_after_trim } from './modules/qc.nf'
include { TrimReads } from './modules/trim_reads.nf'
include { BuildIndex } from './modules/build_index.nf'
include { AlignReads } from './modules/align.nf'
include { AlignStats as OriginalStats } from './modules/align_stats.nf'
include { AlignStats as MKDPStats } from './modules/align_stats.nf'
include { AlignStats as FilteredStats } from './modules/align_stats.nf'
include { MarkDuplicates } from './modules/mkdp.nf'
include { FilterBAM } from './modules/filter_bam.nf'
include { PeakCalling } from './modules/peakcall.nf'
include { AnnotatePeaks } from './modules/peak_annotation.nf'

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

def shouldSkipFiltering() {
    if (!params.skip_filtering) {
        return false
    }

    def filterDir = file(params.filtered)
    if (!filterDir.exists()) {
        if (params.fallback_to_filtering) {
            log.warn "Mark Duplicates directory ${params.mkdp} does not exist! Will perform Mark Duplicates with Picard"
            return false
        } else {
            error "Mark Duplicates directory ${params.mkdp} does not exist and fallback_to_mkdp is disabled!"
        }
    }

    def mkdpFiles = filterDir.listFiles().findAll { f -> f.name.endsWith('_filtered.bam') }
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

    ref_index = BuildIndex(tuple(file(params.ref), file(params.gtf)))
    
    channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map {row -> tuple(row.sample_id, file(row.fastq_path))}
        .set { samples_ch }

    // Create a persistent mapping from short sample name to full sample_id
    // This will be used throughout to maintain the full sample_id
    sample_id_mapping_ch = samples_ch
        .map { sample_id, _fastq -> 
            def short_name = sample_id.split(':')[-1]
            tuple(short_name, sample_id)
        }

    QC_before_trim(samples_ch, params.qc_dir_before_trim)

    trimmed_ch = TrimReads(samples_ch)

    QC_after_trim(trimmed_ch, params.qc_dir_after_trim)

    actuallySkipAlignment = shouldSkipAlignment()
    
    if (actuallySkipAlignment) {
        log.info "Skipping alignment as per request and using existing BAM files from ${params.aligned}"

        // Create a channel for existing BAM files and join with sample_id mapping
        bam_ch = channel
            .fromPath("${params.aligned}/*.bam", checkIfExists: true)
            .map { bam -> 
                def short_name = bam.baseName
                log.info "Found existing BAM file: ${bam.name}"
                tuple(short_name, bam)
            }
            .join(sample_id_mapping_ch)
            .map { _short_name, bam, full_sample_id ->
                log.info "Matched BAM ${bam.name} to sample ${full_sample_id}"
                tuple(full_sample_id, bam)
            }
    } else {
        log.info "Performing alignment with Bowtie2"
        
        // Check for existing BAMs to potentially reuse
        if (file(params.aligned).exists()) {
            existing_bam_ch = channel
                .fromPath("${params.aligned}/*.bam", checkIfExists: false)
                .map { bam ->
                    def short_name = bam.baseName
                    tuple(short_name, bam)
                }
                .join(sample_id_mapping_ch)
                .map { _short_name, bam, full_sample_id ->
                    log.info "Found existing BAM: ${bam.name} for sample ${full_sample_id}"
                    tuple(full_sample_id, bam)
                }
            
            // Find samples that need alignment
            samples_to_align = trimmed_ch
                .join(existing_bam_ch, by: 0, remainder: true)
                .filter { _sample_id, _fastq, existing_bam -> existing_bam == null }
                .map { sample_id, fastq, _existing_bam -> tuple(sample_id, fastq) }
            
            reused_bams = trimmed_ch
                .join(existing_bam_ch, by: 0)
                .map { sample_id, _fastq, bam -> 
                    log.info "Reusing existing BAM for sample ${sample_id}"
                    tuple(sample_id, bam) 
                }
            
            // Perform alignment for samples without existing BAMs
            new_bam_ch = AlignReads(samples_to_align, ref_index)
            
            // Combine reused and newly aligned BAMs
            bam_ch = reused_bams.mix(new_bam_ch)
        } else {
            // No existing BAMs, align everything
            bam_ch = AlignReads(trimmed_ch, ref_index)
        }
    }


    OriginalStats(bam_ch, params.stats_out_dir)

    actuallySkipMKDP = shouldSkipMKDP()

    if (actuallySkipMKDP) {
        log.info "Skipping Mark Duplicates as per request and using existing Mark Duplicates BAM files from ${params.mkdp}"

        mkdp_output_ch = channel
            .fromPath("${params.mkdp}/*_mkdp.bam", checkIfExists: true)
            .map { bam ->
                def short_name = bam.baseName.replaceAll(/_mkdp$/, '')
                tuple(short_name, bam)
            }
            .join(sample_id_mapping_ch)
            .map { short_name, bam, full_sample_id ->
                def metrics = file("${params.mkdp}/${short_name}_marked_duplicates_metrics.txt")
                
                if (!metrics.exists()) {
                    log.warn "Metrics file not found for ${full_sample_id}, creating placeholder"
                    metrics = file("${params.mkdp}/.dummy_metrics.txt")
                }
                
                log.info "Processing existing Mark Duplicates BAM file: ${bam.name} for sample ${full_sample_id}"
                tuple(full_sample_id, bam, metrics)
            }
    } else {
        log.info "Performing Mark Duplicates with Picard"
        
        // Checking for existing Mark Duplicates BAMs to reuse
        existing_mkdp_ch = channel.empty()
        if (file(params.mkdp).exists()) {
            existing_mkdp_ch = channel
                .fromPath("${params.mkdp}/*_mkdp.bam", checkIfExists: false)
                .map { bam ->
                    def short_name = bam.baseName.replaceAll(/_mkdp$/, '')
                    tuple(short_name, bam)
                }
                .join(sample_id_mapping_ch)
                .map { short_name, bam, full_sample_id ->
                    def metrics = file("${params.mkdp}/${short_name}_marked_duplicates_metrics.txt")
                    
                    if (!metrics.exists()) {
                        log.warn "Metrics file not found for ${full_sample_id}, will reprocess"
                        return null
                    }
                    
                    log.info "Reusing existing duplicates-marked BAM file: ${bam.baseName} for sample ${full_sample_id}"
                    tuple(full_sample_id, bam, metrics)
                }
                .filter { f -> f != null }
            
            // Find samples that need MKDP
            samples_needing_mkdp = bam_ch
                .join(existing_mkdp_ch, by: 0, remainder: true)
                .filter { _sample_id, _bam, existing_mkdp, _existing_metrics -> existing_mkdp == null }
                .map { sample_id, bam, _existing_mkdp, _existing_metrics -> tuple(sample_id, bam) }
            
            // Perform Mark Duplicates
            new_mkdp_ch = MarkDuplicates(samples_needing_mkdp)

            mkdp_output_ch = existing_mkdp_ch.mix(new_mkdp_ch)
        } else {
            // No existing MKDP BAMs, process everything
            mkdp_output_ch = MarkDuplicates(bam_ch)
        }
    }

    // Extract just sample_id and bam for stats
    mkdp_ch = mkdp_output_ch.map { sample_id, bam, _metrics -> tuple(sample_id, bam) }
    
    MKDPStats(mkdp_ch, params.mkdp_stats)


    actuallySkipFiltering = shouldSkipFiltering()

    if (actuallySkipFiltering) {
        log.info "Skipping BAM filtering as per request and using existing filtered BAM files from ${params.filtered}"

        filtered_output_ch = channel
            .fromPath("${params.filtered}/*_filtered.bam", checkIfExists: true)
            .map { bam ->
                def short_name = bam.baseName.replaceAll(/_filtered$/, '')
                tuple(short_name, bam)
            }
            .join(sample_id_mapping_ch)
            .map { _short_name, bam, full_sample_id ->
                def bai_paths = [
                    file("${bam}.bai"),
                    file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                ]

                def bai = bai_paths.find { f -> f.exists()}
                if (bai) {
                    log.info "Found index for ${bam.name}: ${bai.name} for sample ${full_sample_id}"
                    return tuple(full_sample_id, bam, bai)
                } else {
                    log.error "Index file not found for ${bam.name}"
                    log.error "Looked for: ${bai_paths.collect{ f -> f.toString()}.join(', ')}"
                    error "Index file not found for ${bam.name}. Please ensure BAM files are properly indexed."
                }
            }
    } else {
        log.info "Performing BAM filtering"
        
        // Checking for existing filtered BAMs to reuse
        existing_filtered_ch = channel.empty()
        if (file(params.filtered).exists()) {
            existing_filtered_ch = channel
                .fromPath("${params.filtered}/*_filtered.bam", checkIfExists: false)
                .map { bam ->
                    def short_name = bam.baseName.replaceAll(/_filtered$/, '')
                    tuple(short_name, bam)
                }
                .join(sample_id_mapping_ch)
                .map { _short_name, bam, full_sample_id ->
                    def bai_paths = [
                        file("${bam}.bai"),
                        file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                    ]
                    def bai = bai_paths.find { f -> f.exists() }
                    if (bai) {
                        log.info "Found existing indexed BAM: ${bam.name} for sample ${full_sample_id}"
                        return tuple(full_sample_id, bam, bai)
                    } else {
                        log.warn "Found BAM without index: ${bam.name} - will reprocess"
                        return null
                    }
                }
                .filter { f -> f != null }
            
            // Find samples that need filtering
            samples_needing_filtering = mkdp_ch
                .join(existing_filtered_ch, by: 0, remainder: true)
                .filter { _sample_id, _mkdp_bam, existing_filtered, _existing_bai -> existing_filtered == null }
                .map { sample_id, mkdp_bam, _existing_filtered, _existing_bai -> tuple(sample_id, mkdp_bam) }
            
            new_filtered_ch = FilterBAM(samples_needing_filtering)

            filtered_output_ch = existing_filtered_ch.mix(new_filtered_ch)
        } else {
            // No existing filtered BAMs, filter everything
            filtered_output_ch = FilterBAM(mkdp_ch)
        }
    }

    // Extract just sample_id and bam for stats
    filtered_ch = filtered_output_ch.map { sample_id, bam, _bai -> tuple(sample_id, bam) }
    
    FilteredStats(filtered_ch, params.filtered_stats)

}