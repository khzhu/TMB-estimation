#!/usr/bin/env nextflow
import groovy.json.JsonSlurper

include { FASTQ_FASTP_FASTQC        } from './subworkflows/fastq_fastp_fastqc/main'
include { FASTQ_ALIGN_MARKDUP_STATS } from './subworkflows/align_markdup_stats/main'
include { INDEX_GENOME              } from '../subworkflows/index_genome/main'
include { BAM_RECALIBRATION         } from './subworkflows/bam_recalibration/main'

// main workflow
workflow {
    log.info """\
    TMB estimation pipeline ${params.release}
    ========================================
    samples file  : ${params.input_json}
    output dir    : ${params.output_dir}
    """

    // Read sample config file
    def jsonSlurper = new JsonSlurper()

    // If params.input_json exists, use that for multi sample processing.
    // Otherwise, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // run samples through the pipeline
    ch_samples = Channel.from(multi_params.collect{ it -> tuple([ id: it.specimen_num, single_end:false ],
                [ file(it.read1, checkIfExists: true), file(it.read2, checkIfExists: true) ]) })
    ch_adapter_fasta = Channel.fromPath(params.adapter_fasta)
    ch_fasta = Channel.value([
        [id:'fasta'],
        file(params.reference_file, checkIfExists: true),
    ])
    ch_target_bed = Channel.fromPath(params.exome_plus_tumor_panel_bed, checkIfExists: true)
    ch_known_indel_sites = Channel.fromPath(params.known_indel_vcf)
    ch_known_indel_sites_tbi = Channel.fromPath(params.known_indel_vcf_tbi)
    ch_known_snp_sites = Channel.fromPath(params.known_snp_vcf)
    ch_known_snp_sites_tbi = Channel.fromPath(params.known_snp_vcf_tbi)
    ch_versions = Channel.empty()

    // With paired-end data
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false
    val_sort_bam      = true

    // index reference genome
    INDEX_GENOME ( ch_fasta )
    ch_versions = ch_versions.mix( INDEX_GENOME.out.versions )

    // Trim raw seqeunce reads
    FASTQ_FASTP_FASTQC ( ch_samples, ch_adapter_fasta, save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
    ch_versions = ch_versions.mix( FASTQ_FASTP_FASTQC.out.versions )

    // Mapping trimmed reads to hg19 genome
    FASTQ_ALIGN_MARKDUP_STATS ( FASTQ_FASTP_FASTQC.out.reads, INDEX_GENOME.out.index,
                                val_sort_bam, ch_fasta)
    ch_versions = ch_versions.mix( FASTQ_ALIGN_MARKDUP_STATS.out.versions )

    BAM_RECALIBRATION ( tuple(FASTQ_ALIGN_MARKDUP_STATS.out.bam, FASTQ_ALIGN_MARKDUP_STATS.out.bam.bai, ch_target_bed),
        ch_fasta, INDEX_GENOME.out.fai, INDEX_GENOME.out.dict, ch_known_snp_sites.join(ch_known_indel_sites),
        ch_known_snp_sites_tbi.join(ch_known_indel_sites_tbi))
    ch_versions = ch_versions.mix( BAM_RECALIBRATION.out.versions )
}