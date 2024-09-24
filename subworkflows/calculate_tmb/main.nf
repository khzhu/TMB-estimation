#!/usr/bin/env nextflow

include { TMB_CALIBRATION  } from '../../modules/tmb/main'

workflow TMB_CALIBER {

    take:
    mutect2_maf        // channel: [mandatory] [ val(meta), path(mutect2_maf) ]
    strelka2_maf       // channel: [mandatory] [ val(meta), path(strelka2_maf) ]

    main:
    ch_versions         = Channel.empty()

    mutect2_out = mutect2_maf.collectFile( name: 'mutect2_out.maf', sort: { it.size() }, keepHeader:true, skip:2 )
    strelka2_out = strelka2_maf.collectFile( name: 'strelka2_out.maf', sort: { it.size() }, keepHeader:false, skip:2 )

    TMB_CALIBRATION ( [[id:'ccs'], mutect2_out, strelka2_out ] )
    ch_versions = ch_versions.mix(TMB_CALIBRATION.out.versions)

    emit:
    tmb    = TMB_CALIBRATION.out.tmb         // channel: [ val(meta), path(tmb) ]
    maf    = TMB_CALIBRATION.out.maf         // channel: [ val(meta), path(maf) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}