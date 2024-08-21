#!/usr/bin/env nextflow
import groovy.json.JsonSlurper

include {FASTQ_TRIM_FASTQC } from './subworkflows/fastq_trim_fastqc.nf'

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
    skip_fastqc = false
    skip_trim  = false
    FASTQ_TRIM_FASTQC ( ch_samples, ch_adapter_fasta, skip_trim, skip_fastqc )
}