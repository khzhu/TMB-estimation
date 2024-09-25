#!/usr/bin/env nextflow

include { TMB_CALIBRATION  } from '../../modules/tmb/main'

workflow TMB_CALIBER {

    take:
    maf_files     // channel: [mandatory] [ val(meta), path(maf) ]

    main:
    ch_versions         = Channel.empty()

    TMB_CALIBRATION ( maf_files )
    ch_versions = ch_versions.mix(TMB_CALIBRATION.out.versions)

    emit:
    tmb    = TMB_CALIBRATION.out.tmb         // channel: [ val(meta), path(tmb) ]
    maf    = TMB_CALIBRATION.out.maf         // channel: [ val(meta), path(maf) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}