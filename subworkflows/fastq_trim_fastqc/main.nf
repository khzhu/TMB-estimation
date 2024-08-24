//
// Read QC and trimming
//

include { FASTQC as FASTQC_RAW  } from '../modules/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../modules/fastqc/main'
include { TRIMMOMATIC           } from '../modules/trimmomatic/main'


workflow FASTQ_TRIM_FASTQC {
    take:
    ch_reads              // channel: [ val(meta), path(reads)  ]
    ch_adapter_fasta     // channel: [ path(fasta) ]
    val_skip_trim        // value: boolean
    val_skip_fastqc       // value: boolean

    main:

    ch_versions = Channel.empty()

    ch_fastqc_raw_html = Channel.empty()
    ch_fastqc_raw_zip  = Channel.empty()
    if (!val_skip_fastqc) {
        FASTQC_RAW (
            ch_reads
        )
        ch_fastqc_raw_html = FASTQC_RAW.out.html
        ch_fastqc_raw_zip  = FASTQC_RAW.out.zip
        ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    ch_trim_reads        = ch_reads
    ch_unpaired_reads    = Channel.empty()
    ch_trim_summary      = Channel.empty()
    ch_trim_log          = Channel.empty()
    ch_out_log           = Channel.empty()  
    ch_fastqc_trim_html  = Channel.empty()
    ch_fastqc_trim_zip   = Channel.empty()

    if (!val_skip_trim) {
        TRIMMOMATIC (
            ch_reads
        )
        ch_trim_reads        = TRIMMOMATIC.out.trimmed_reads
        ch_unpaired_reads    = TRIMMOMATIC.out.unpaired_reads
        ch_trim_summary      = TRIMMOMATIC.out.summary
        ch_trim_log          = TRIMMOMATIC.out.trim_log
        ch_out_log           = TRIMMOMATIC.out.out_log
        ch_versions          = ch_versions.mix(TRIMMOMATIC.out.versions.first())

        //
        // Filter empty FastQ files after adapter trimming so FastQC doesn't fail
        //
        if (!val_skip_fastqc) {
            FASTQC_TRIM (
                ch_trim_reads
            )
            ch_fastqc_trim_html = FASTQC_TRIM.out.html
            ch_fastqc_trim_zip  = FASTQC_TRIM.out.zip
            ch_versions      = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        }
    }

    emit:
    reads             = ch_trim_reads         // channel: [ val(meta), path(reads) ]
    unpaired_reads    = ch_unpaired_reads     // channel: [ val(meta), path(unpaired_reads) ]
    trim_log          = ch_trim_log           // channel: [ val(meta), path(trim_log) ]
    out_log           = ch_out_log            // channel: [ val(meta), path(out_log) ]
    trim_summary      = ch_trim_summary       // channel: [ val(meta), path(trim_summary) ]

    fastqc_raw_html  = ch_fastqc_raw_html    // channel: [ val(meta), path(html) ]
    fastqc_raw_zip   = ch_fastqc_raw_zip     // channel: [ val(meta), path(zip) ]
    fastqc_trim_html = ch_fastqc_trim_html   // channel: [ val(meta), path(html) ]
    fastqc_trim_zip  = ch_fastqc_trim_zip    // channel: [ val(meta), path(zip) ]

    versions = ch_versions.ifEmpty(null)     // channel: [ path(versions.yml) ]
}