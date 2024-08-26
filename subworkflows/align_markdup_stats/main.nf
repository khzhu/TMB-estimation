//
// Alignment with BWA
//

include { BWA_MEM           } from '../../modules/bwa/mem/main'
include { SAMBAMBA_MARKDUP  } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT } from '../../modules/sambamba/flagstat/main'

workflow FASTQ_ALIGN_MARKDUP_STATS {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), [ path(index) ] ]
    val_sort_bam    // boolean (mandatory): true or false
    ch_fai          // channel (optional) : [ val(meta3), path(fai) ]

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    BWA_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // Run sambamba deduplicate and flagstat
    //
    SAMBAMBA_MARKDUP ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    //
    // Generate sambamba flagstat
    //
    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_FLAGSTAT.out.versions)

    emit:
    bam      = SAMBAMBA_MARKDUP.out.bam     // channel: [ val(meta), path(bam) ]
    bai      = SAMBAMBA_MARKDUP.out.bai     // channel: [ val(meta), path(bai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats  // channel: [ val(meta), path(flagstat) ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}