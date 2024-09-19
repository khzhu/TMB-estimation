#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { SNV_MUTECT2              } from './subworkflows/snv_mutect2/main'
include { SNV_STRELKA2             } from './subworkflows/snv_strelka2/main'

// main workflow
workflow {
    log.info """\
    TMB estimation pipeline ${params.release}
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sample sheet  : ${params.input_json}
    output dir    : ${params.output_dir}
    """

    // Read sample config file
    def jsonSlurper = new JsonSlurper()

    // If params.input_json exists, use that for multi sample processing.
    // Otherwise, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // run samples through the pipeline
    samples = Channel.from(multi_params.collect{ it -> tuple([
                id: it.patient_id, tissue: it.tissue, purity: it.purity ],
                [ file(it.bam_tumor, checkIfExists: true), file(it.bam_normal, checkIfExists: true) ],
                [ file(it.bai_tumor, checkIfExists: true), file(it.bai_normal, checkIfExists: true) ] ) })
    ch_versions = Channel.empty()

    bed_files = Channel.fromPath(params.tumor_panel_bed_files, checkIfExists: true)
    ch_input_files = samples.combine(bed_files)
                    .map { meta, input_bams, input_index_files, intervals ->
                        new_meta = meta.clone()
                        new_meta.sid = intervals.baseName != "no_intervals" ? new_meta.id + "_" + intervals.baseName : new_meta.id
                        intervals = intervals.baseName != "no_intervals" ? intervals : []
                        [new_meta, input_bams, input_index_files, intervals]
                    }
    // calling mutect2 somatic variants
    SNV_MUTECT2 (ch_input_files,
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                file(params.gnomad_exome_vcf, checkIfExists: true),
                file(params.gnomad_exome_vcf_tbi, checkIfExists: true),
                file(params.exac_common_vcf, checkIfExists: true),
                file(params.exac_common_vcf_tbi, checkIfExists: true),
                file(params.vep_cache, checkIfExists: true),
                file(params.filter_vcf, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_MUTECT2.out.versions )

    // calling strelka2 somatic variants
    SNV_STRELKA2 (samples,
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                [[ id:'genome'], file(params.tumor_panel_bed_gz, checkIfExists: true),
                                 file(params.tumor_panel_bed_tbi, checkIfExists: true)],
                                 file(params.vep_cache, checkIfExists: true),
                                 file(params.filter_vcf, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_STRELKA2.out.versions )
}