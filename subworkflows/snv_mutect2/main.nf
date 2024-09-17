//
// GATK4_MUTECT2: Call variants with Mutect2 (2.1 part of GATK 4)
//

include { GATK4_MUTECT2                                         } from '../../modules/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR  } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION                          } from '../../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS                               } from '../../modules/gatk4/filtermutectcalls/main'
include { BCFTOOLS_NORM                                         } from '../../modules/bcftools/norm/main'
include { GATK4_MERGEVCFS                                       } from '../../modules/gatk4/mergevcfs/main'


workflow SNV_MUTECT2 {

    take:
    ch_input_files                    // channel: [mandatory] [ val(meta), path(ch_input_bams), path(ch_input_bams), path(intervals) ]
    ch_fasta                          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_fai                            // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_dict                           // channel: [mandatory] [ val(meta), path(fasta) ]
    germline_resource
    germline_resource_tbi
    pileup_variants
    pileup_variants_tbi

    main:
    ch_versions         = Channel.empty()

    GATK4_MUTECT2 ( ch_input_files,
                    ch_fasta, ch_fai, ch_dict,
                    germline_resource, germline_resource_tbi)
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    GETPILEUPSUMMARIES_TUMOR  ( ch_input_files,
                                ch_fasta, ch_fai, ch_dict,
                                pileup_variants,
                                pileup_variants_tbi,
                                false)

    GETPILEUPSUMMARIES_NORMAL ( ch_input_files,
                                ch_fasta, ch_fai, ch_dict,
                                pileup_variants,
                                pileup_variants_tbi,
                                true)

    GATK4_CALCULATECONTAMINATION ( GETPILEUPSUMMARIES_TUMOR.out.table,
                                    GETPILEUPSUMMARIES_NORMAL.out.table)
    GATK4_FILTERMUTECTCALLS ( GATK4_MUTECT2.out.vcf.combine(GATK4_MUTECT2.out.tbi, by:[0,3]).combine(GATK4_MUTECT2.out.stats, by:[0,3]),
                            GATK4_CALCULATECONTAMINATION.out.contamination,
                            ch_fasta, ch_fai, ch_dict)
    
    BCFTOOLS_NORM ( GATK4_FILTERMUTECTCALLS.out.vcf.combine(GATK4_FILTERMUTECTCALLS.out.tbi, by:[0,3]), ch_fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    GATK4_MERGEVCFS ( BCFTOOLS_NORM.out.vcf, ch_dict )

    emit:
    vcf      = GATK4_MERGEVCFS.out.vcf      // channel: [ val(meta), path("*.vcf.gz") ]
    tbi      = GATK4_MERGEVCFS.out.tbi      // channel: [ val(meta), path("*.tbi") ]
    stats    = GATK4_MUTECT2.out.stats      // channel: [ val(meta), path("*.stats") ]
    f1r2     = GATK4_MUTECT2.out.f1r2
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}