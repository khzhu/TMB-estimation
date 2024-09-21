//
// Alignment to genome with BWA and sort reads by SAMBAMBA
//

include { SEQKIT_SPLIT2          } from '../../modules/seqkit/main'
include { BWA_MEM                } from '../../modules/bwa/mem/main'
include { SAMBAMBA_MERGE         } from '../../modules/sambamba/merge/main'
include { SAMTOOLS_SORT          } from '../../modules/samtools/sort/main'
include { SAMBAMBA_MARKDUP       } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT      } from '../../modules/sambamba/flagstat/main'
include { GATK4_BASERECALIBRATOR } from '../../modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR        } from '../../modules/gatk4/applybqsr/main'
include { SAMTOOLS_STATS         } from '../../modules/samtools/stats/main'
include { SAMTOOLS_CONVERT       } from '../../modules/samtools/convert/main'


workflow ALIGN_MARKDUP_BQSR_STATS {
    take:
    reads                  // channel (mandatory): [ val(meta), [ path(reads) ] ]
    bwa_index              // channel (mandatory): [ val(meta2), [ path(index) ] ]
    fasta                  // channel (optional) : [ val(meta2), [ path(fasta) ] ]
    val_sort_bam           // boolean (mandatory): true or false
    intervals              // channel: [ path(intervals) ]
    fai                    // channel: [ val(meta2), path(fai) ]
    dict                   // channel: [ val(meta2), path(dict) ]
    snp_known_sites        // channel: [ path(known_sites) ]
    snp_known_sites_tbi    // channel: [ path(known_sites_tbi) ]
    indel_known_sites      // channel: [ path(known_sites) ]
    indel_known_sites_tbi  // channel: [ path(known_sites_tbi)
    split_reads            // boolean (mandatory): true or false

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    if ( split_reads ) {
        //Split trimmed reads into 2 parts
        SEQKIT_SPLIT2 ( reads )
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
        BWA_MEM ( SEQKIT_SPLIT2.out.reads, bwa_index, [[id:'genome'],fasta], val_sort_bam )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
        SAMBAMBA_MERGE ( BWA_MEM.out.bam )
        // Sort bam with samtools
        SAMTOOLS_SORT ( SAMBAMBA_MERGE.out.bam , [[id:'genome'],fasta] )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    } else {
        BWA_MEM ( reads, bwa_index, [[id:'genome'],fasta], val_sort_bam )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)
        // Sort bam with samtools
        SAMTOOLS_SORT ( BWA_MEM.out.bam , [[id:'genome'],fasta] )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    }

    //
    // Run sambamba deduplicate and flagstat
    //
    SAMBAMBA_MARKDUP ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    //
    // Generate sambamba flagstat
    //
    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )

    //
    // Apply for gatk4 BQSR
    //
    GATK4_BASERECALIBRATOR ( SAMBAMBA_MARKDUP.out.bam,
                            intervals, [[id:'genome'],fasta], fai, dict,
                            snp_known_sites, snp_known_sites_tbi,
                            indel_known_sites, indel_known_sites_tbi)
    GATK4_APPLYBQSR ( SAMBAMBA_MARKDUP.out.bam,
                        GATK4_BASERECALIBRATOR.out.table,
                        intervals, [[id:'genome'],fasta], fai, dict)
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

    //
    // Generate samtool stats
    //
    SAMTOOLS_STATS ( GATK4_APPLYBQSR.out.bam, [[id:'genome'],fasta])

    //
    // Convert to CRAM
    //
    SAMTOOLS_CONVERT ( GATK4_APPLYBQSR.out.bam, [[id:'genome'],fasta])

    emit:
    bam      = GATK4_APPLYBQSR.out.bam        // channel: [ val(meta), path(bam) ]
    bai      = GATK4_APPLYBQSR.out.bai        // channel: [ val(meta), path(bai) ]
    cram     = SAMTOOLS_CONVERT.out.bam       // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_CONVERT.out.crai      // channel: [ val(meta), path(crai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats    // channel: [ val(meta), path(flagstat) ]
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
    versions = ch_versions                    // channel: [ path(versions.yml) ]
}