#!/usr/bin/env nextflow

include { TMB_CALIBRATION  } from '../../modules/tmb/main'

workflow TMB_CALIBER {

    take:
    mutect2_maf     // channel: [mandatory] [ val(meta), path(maf) ]
    strelka2_maf    // channel: [mandatory] [ val(meta), path(maf) ]

    main:
    ch_versions         = Channel.empty()

    TMB_CALIBRATION ( mutect2_maf, strelka2_maf )
    ch_versions = ch_versions.mix(TMB_CALIBRATION.out.versions)

    emit:
    tmb    = TMB_CALIBRATION.out.tmb         // channel: [ val(meta), path(tmb) ]
    maf    = TMB_CALIBRATION.out.maf         // channel: [ val(meta), path(maf) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}