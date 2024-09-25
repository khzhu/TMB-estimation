#!/usr/bin/env nextflow

include { TMB_CALIBRATION  } from '../../modules/tmb/main'

workflow TMB_CALIBER {

    take:
    input_maf     // channel: [mandatory] [ val(meta), path(maf) ]

    main:
    ch_versions         = Channel.empty()

    TMB_CALIBRATION ( input_maf )
    ch_versions = ch_versions.mix(TMB_CALIBRATION.out.versions)

    emit:
    tmb    = TMB_CALIBRATION.out.tmb         // channel: [ val(meta), path(tmb) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}