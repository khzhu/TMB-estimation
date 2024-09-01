//
// Base Quality Score Recalibration (BQSR) by GATK4
//

include { GATK4_BASERECALIBRATOR } from '../../modules/gatk4/baserecalibrator/main'

workflow BAM_RECALIBRATION {

    take:
    ch_bam                    // channel: [ path(bam) ]
    ch_bai                    // channel: [ path(bai) ]
    ch_bed                    // channel: [ path(bed) ]
    ch_fasta                  // channel: [ path(fasta) ]
    ch_fai                    // channel: [ path(fai) ]
    ch_dict                   // channel: [ path(dict) ]
    ch_known_sites            // channel: [ path(known_sites) ]
    ch_known_sites_tbi        // channel: [ path(known_sites_tbi) ]

    main:

    ch_versions = Channel.empty()

    GATK4_BASERECALIBRATOR ( [[id: "bqsr" ], ch_bam, ch_bai, ch_bed],
                            ch_fasta, ch_fai, ch_dict,
                            ch_known_sites, ch_known_sites_tbi)
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    emit:
    bam      = GATK4_BASERECALIBRATOR.out.bam    // channel: [ val(meta), path(bam) ]
    bai      = GATK4_BASERECALIBRATOR.out.bai    // channel: [ val(meta), path(bai) ]
    table    = GATK4_BASERECALIBRATOR.out.table  // channel: [ val(meta), path(table) ]
    versions = ch_versions                       // channel: [ versions.yml ]
}