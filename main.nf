#!/usr/bin/env nextflow
import groovy.json.JsonSlurper

include { FASTQ_FASTP_FASTQC        } from './subworkflows/fastq_fastp_fastqc/main'
include { BWA_INDEX                 } from '../../modules/bwa/index/main'
include { FASTQ_ALIGN_MARKDUP_STATS } from './subworkflows/align_markdup_stats/main'

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

    // With paired-end data
    save_trimmed_fail = false
    save_merged       = false
    skip_fastqc       = false
    skip_fastp        = false
    val_sort_bam      = true

    FASTQ_FASTP_FASTQC ( ch_samples, ch_adapter_fasta, save_trimmed_fail, save_merged, skip_fastp, skip_fastqc )
    BWA_INDEX ( ch_fasta )
    FASTQ_ALIGN_MARKDUP_STATS ( FASTQ_FASTP_FASTQC.out.reads, BWA_INDEX.out.index, val_sort_bam, ch_fasta )
}