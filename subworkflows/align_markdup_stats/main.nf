//
// Alignment with BWA
//
include { BWA_MEM           } from '../../modules/bwa/mem/main'
include { SAMBAMBA_MARKDUP  } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT } from '../../modules/sambamba/flagstat/main'
include { SAMTOOLS_SORT     } from '../../modules/samtools/main'

workflow FASTQ_ALIGN_MARKDUP_STATS {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam    // boolean (mandatory): true or false
    ch_fasta        // channel (optional) : [ val(meta3), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    BWA_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // Run sambamba deduplicate and flagstat
    //

    SAMBAMBA_MARKDUP ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_FLAGSTAT.out.versions)

    SAMTOOLS_SORT ( SAMBAMBA_MARKDUP.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    emit:
    bam      = SAMBAMBA_MARKDUP.out.bam     // channel: [ val(meta), path(bam) ]
    bai      = SAMBAMBA_MARKDUP.out.bai     // channel: [ val(meta), path(bai) ]
    cram     = SAMTOOLS_SORT.out.cram       // channel: [ val(meta), path(cram) ]
    crai     = SAMTOOLS_SORT.out.crai       // channel: [ val(meta), path(crai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats  // channel: [ val(meta), path(flagstat) ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}