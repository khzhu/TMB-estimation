//
// Base Quality Score Recalibration (BQSR) by GATK4
//

include { GATK4_BASERECALIBRATOR } from '../../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR        } from '../../modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX         } from '../../modules/samtools/index/main'
include { BAM_STATS_SAMTOOLS     } from '../../subworkflows/bam_stats_samtools/main'

workflow BAM_RECALIBRATION {

    take:
    ch_bam_bai_bed      // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    ch_fasta            // channel: [ path(fasta) ]
    ch_fai              // channel: [ path(fai) ]
    ch_dict             // channel: [ path(dict) ]
    ch_known_sites      // channel: [ path(known_sites) ]
    ch_known_sites_tbi  // channel: [ path(known_sites_tbi) ]

    main:

    ch_versions = Channel.empty()

    GATK4_BASERECALIBRATOR ( ch_bam_bai_bed, ch_fasta, ch_fai, ch_dict, ch_known_sites, ch_known_sites_tbi)
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    GATK4_APPLYBQSR ( ch_bam_bai_bed.mix(GATK4_BASERECALIBRATOR.out.table), ch_fasta, ch_fai, ch_dict)

    ch_bqsr = GATK4_APPLYBQSR.out.bam.mix(GATK4_APPLYBQSR.out.cram)

    SAMTOOLS_INDEX ( ch_bqsr )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_aligment_index = ch_bqsr
        .join(SAMTOOLS_INDEX.out.bai,  by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.crai, by: [0], remainder: true)
        .map{meta, bam, cram, bai, crai ->
            if (bai) [ meta, bam, bai ]
            else [ meta, cram, crai ]
        }

    BAM_STATS_SAMTOOLS ( ch_aligment_index, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = GATK4_APPLYBQSR.out.bam           // channel: [ val(meta), path(bam) ]
    cram     = GATK4_APPLYBQSR.out.cram          // channel: [ val(meta), path(cram) ]
    table    = GATK4_BASERECALIBRATOR.out.table  // channel: [ val(meta), path(table) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    crai     = SAMTOOLS_INDEX.out.crai           // channel: [ val(meta), path(crai) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}