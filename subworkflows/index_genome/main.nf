#!/usr/bin/env nextflow

//
// FASTA_INDEX: Build aligner specific index for fasta files
//

include { BWA_INDEX           } from '../../modules/bwa/index/main'
include { SAMTOOLS_FAIDX      } from '../../modules/samtools/faidx/main'
include { GATK4_CREATESEQDICT } from '../../modules/gatk4/createseqdict/main'

workflow INDEX_GENOME {

    take:
    ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]

    main:
    ch_versions         = Channel.empty()

    BWA_INDEX ( ch_fasta )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GATK4_CREATESEQDICT ( ch_fasta )
    ch_versions = ch_versions.mix(GATK4_CREATESEQDICT.out.versions)

    emit:
    index    = BWA_INDEX.out.index           // channel: [ val(meta), path(index) ]
    fai      = SAMTOOLS_FAIDX.out.fai        // channel: [ val(meta), path(fai) ]
    dict     = GATK4_CREATESEQDICT.out.dict  // channel: [ val(meta), path(dict) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}