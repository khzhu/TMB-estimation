//
// Alignment to genome with BWA and sort reads by SAMBAMBA
//

include { BWA_MEM                } from '../../modules/bwa/mem/main'
include { SAMTOOLS_SORT          } from '../../modules/samtools/sort/main'
include { SAMBAMBA_MARKDUP       } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT      } from '../../modules/sambamba/flagstat/main'
include { SAMTOOLS_INDEX         } from '../../modules/samtools/index/main'
include { GATK4_BASERECALIBRATOR } from '../../modules/gatk4/baserecalibrator/main'
include { SAMTOOLS_CONVERT       } from '../../modules/samtools/convert/main'
include { SAMTOOLS_STATS         } from '../../modules/samtools/stats/main'

workflow ALIGN_MARKDUP_BQSR_STATS {
    take:
    reads                  // channel (mandatory): [ val(meta), [ path(reads) ] ]
    bwa_index              // channel (mandatory): [ val(meta2), [ path(index) ] ]
    fasta                  // channel (optional) : [ val(meta3), path(fasta) ]
    val_sort_bam           // boolean (mandatory): true or false
    intervals              // channel: [ path(intervals) ]
    fai                    // channel: [ path(fai) ]
    dict                   // channel: [ path(dict) ]
    snp_known_sites        // channel: [ path(known_sites) ]
    snp_known_sites_tbi    // channel: [ path(known_sites_tbi) ]
    indel_known_sites      // channel: [ path(known_sites) ]
    indel_known_sites_tbi  // channel: [ path(known_sites_tbi)

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    BWA_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // Sort bam with samtools
    //
    SAMTOOLS_SORT ( BWA_MEM.out.bam , ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // Run sambamba deduplicate and flagstat
    //
    SAMBAMBA_MARKDUP ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    //
    // Generate sambamba flagstat
    //
    SAMTOOLS_INDEX ( SAMBAMBA_MARKDUP.out.bam )

    //
    // Generate sambamba flagstat
    //
    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )

    //
    // Apply BQSR
    //
    GATK4_BASERECALIBRATOR ( [SAMBAMBA_MARKDUP.out.bam, SAMTOOLS_INDEX.bai, ch_bed],
                                ch_fasta, ch_fai, ch_dict,
                                ch_snp_known_sites, ch_snp_known_sites_tbi
                                ch_indel_known_sites, ch_indel_known_sites_tbi )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    //
    // Convert bam to cram
    //
    SAMTOOLS_CONVERT ( [SAMBAMBA_MARKDUP.out.bam, SAMTOOLS_INDEX.bai], ch_fasta)

    //
    // Generate bamstats
    //
    SAMTOOLS_STATS ( [SAMBAMBA_MARKDUP.out.bam, SAMTOOLS_INDEX.bai], ch_fasta )

    emit:
    bam      = SAMBAMBA_MARKDUP.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai        // channel: [ val(meta), path(bai) ]
    cram     = SAMTOOLS_CONVERT.out.cram     // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_CONVERT.out.crai     // channel: [ val(meta), path(crai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats   // channel: [ val(meta), path(flagstat) ]
    stats    = SAMTOOLS_STATS.out.stats      // channel: [ val(meta), path(stats) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}