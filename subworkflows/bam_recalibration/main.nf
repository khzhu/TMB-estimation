//
// Base Quality Score Recalibration (BQSR) by GATK4
//

include { GATK4_BASERECALIBRATOR } from '../../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR        } from '../../modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX         } from '../../modules/samtools/index/main'
include { SAMTOOLS_CONVERT       } from '../../modules/samtools/convert/main'
include { BAM_STATS_SAMTOOLS     } from '../../subworkflows/bam_stats_samtools/main'

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

    GATK4_BASERECALIBRATOR ( [ch_bam, ch_bai, ch_bed],
                            ch_fasta, ch_fai, ch_dict,
                            ch_known_sites, ch_known_sites_tbi)
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    SAMTOOLS_CONVERT ( Channel.of([ [id: "bqsr", single_end: false],
                        GATK4_BASERECALIBRATOR.out.bam, GATK4_BASERECALIBRATOR.out.bai]),
                        Channel.of([ [ id:'genome' ], ch_fasta]),
                        Channel.of([ [ id:'fai' ], ch_fai ]))
    BAM_STATS_SAMTOOLS ( Channel.of([ [ id:"bamstats", single_end:false ],
                    GATK4_BASERECALIBRATOR.out.bam, GATK4_BASERECALIBRATOR.out.bai]),
                    Channel.of([[ id:'genome' ], ch_fasta]) )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = GATK4_BASERECALIBRATOR.out.bam    // channel: [ val(meta), path(bam) ]
    cram     = SAMTOOLS_CONVERT.out.cram         // channel: [ val(meta), path(cram) ]
    table    = GATK4_BASERECALIBRATOR.out.table  // channel: [ val(meta), path(table) ]
    bai      = GATK4_BASERECALIBRATOR.out.bai    // channel: [ val(meta), path(bai) ]
    crai     = SAMTOOLS_CONVERT.out.crai         // channel: [ val(meta), path(crai) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}