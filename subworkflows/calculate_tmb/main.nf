#!/usr/bin/env nextflow

include { TMB_CALIBRATION  } from '../../modules/tmb/main'

workflow TMB_CALIBER {

    take:
    input_maf        // channel: [mandatory] [ val(meta), path(maf) ]

    main:
    ch_versions         = Channel.empty()

    ch_maf_file = input_maf.collectFile( name: 'result.out', sort: { it.size() }, keepHeader:true, skip:2 )
    TMB_CALIBRATION ( ch_maf_file )
    ch_versions = ch_versions.mix(TMB_CALIBRATION.out.versions)

    emit:
    tmb    = TMB_CALIBRATION.out.tmb         // channel: [ val(meta), path(tmb) ]
    maf    = ch_maf_file                     // channel: [ val(meta), path(maf) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}