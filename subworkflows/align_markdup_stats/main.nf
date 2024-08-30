//
// Alignment to genome with BWA and sort reads by SAMBAMBA
//

include { BWA_MEM           } from '../../modules/bwa/mem/main'
include { SAMBAMBA_SORT     } from '../../modules/sambamba/sort/main'
include { SAMBAMBA_MARKDUP  } from '../../modules/sambamba/markdup/main'
include { SAMBAMBA_FLAGSTAT } from '../../modules/sambamba/flagstat/main'
include { SAMTOOLS_SORT     } from '../../modules/samtools/sort/main'
include { SAMTOOLS_INDEX    } from '../../modules/samtools/index/main'

workflow FASTQ_ALIGN_MARKDUP_STATS {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), [ path(index) ] ]
    ch_fasta         // channel (optional) : [ val(meta3), path(fasta) ]
    val_sort_bam    // boolean (mandatory): true or false

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    BWA_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //
    // Run sambamba sort
    //
    SAMBAMBA_SORT ( BWA_MEM.out.bam )

    //
    // Run sambamba deduplicate and flagstat
    //
    SAMBAMBA_MARKDUP ( SAMBAMBA_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    //
    // Generate sambamba flagstat
    //
    SAMBAMBA_FLAGSTAT ( SAMBAMBA_MARKDUP.out.bam )

    //
    // Sort bam with samtools
    //
    SAMTOOLS_SORT ( SAMBAMBA_MARKDUP.out.bam , ch_fasta)

    //
    // Generate sambamba flagstat
    //
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    emit:
    bam      = SAMTOOLS_SORT.out.bam        // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai       // channel: [ val(meta), path(bai) ]
    flagstat = SAMBAMBA_FLAGSTAT.out.stats  // channel: [ val(meta), path(flagstat) ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}