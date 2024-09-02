#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { INDEX_GENOME             } from './subworkflows/index_genome/main'
include { FASTQ_FASTP_FASTQC       } from './subworkflows/fastq_fastp_fastqc/main'
include { ALIGN_MARKDUP_BQSR_STATS } from './subworkflows/align_markdup_bqsr_stats/main'

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
    ch_versions = Channel.empty()

    // index genome reference
    INDEX_GENOME ( [[ id:'g1k_v37', single_end:false ], file(params.reference_file, checkIfExists: true)])

    // Trim raw seqeunce reads with paired-end data
    FASTQ_FASTP_FASTQC ( ch_samples, ch_adapter_fasta, params.save_trimmed_fail,
                        params.save_merged, params.skip_fastp, params.skip_fastqc )
    ch_versions = ch_versions.mix( FASTQ_FASTP_FASTQC.out.versions )

    // Align to the reference genome
    ALIGN_MARKDUP_BQSR_STATS ( FASTQ_FASTP_FASTQC.out.reads,
                INDEX_GENOME.out.index,
                file(params.reference_file, checkIfExists: true),
                params.val_sort_bam,
                file(params.exome_plus_tumor_panel_bed, checkIfExists: true),
                file(params.fai_file, checkIfExists: true),
                file(params.genome_dict, checkIfExists: true),
                file(params.known_snp_vcf, checkIfExists: true),
                file(params.known_snp_vcf_tbi, checkIfExists: true),
                file(params.known_indel_vcf, checkIfExists: true),
                file(params.known_indel_vcf_tbi, checkIfExists: true) )
    ch_versions = ch_versions.mix( ALIGN_MARKDUP_BQSR_STATS.out.versions )
}