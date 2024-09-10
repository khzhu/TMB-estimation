//
// GATK4_MUTECT2: Call variants with Mutect2 (2.1 part of GATK 4)
//

include { GATK4_MUTECT2                   } from '../../modules/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES_NORMAL } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES_TUMOR  } from '../../modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    } from '../../modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         } from '../../modules/gatk4/filtermutectcalls/main'
include { BCFTOOLS_NORM                   } from '../../modules/bcftools/norm/main'


workflow SNV_MUTECT2 {

    take:
    tuple val(meta), path(bam_tumor)
    tuple val(meta), path(bam_normal)
    tuple val(meta2), path(fasta) // channel: [mandatory] [ val(meta), path(fasta) ]
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)
    path(intervals)
    path(germline_resource)
    path(germline_resource_tbi)
    path(pileup_variants)
    path(pileup_variants_tbi)
    path(panel_of_normals), optional:true

    main:
    ch_versions         = Channel.empty()

    GATK4_MUTECT2 ( bam_tumor,  bam_normal, fasta, fai, dict, intervals, germline_resource, germline_resource_tbi)
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    GATK4_GETPILEUPSUMMARIES_TUMOR  ( bam_tumor, 
                                      intervals, 
                                      fasta, fai, dict, 
                                      pileup_variants, 
                                      pileup_variants_tbi)

    GATK4_GETPILEUPSUMMARIES_NORMAL ( bam_normal, 
                                      intervals, 
                                      fasta, fai, dict,
                                      pileup_variants, 
                                      pileup_variants_tbi)

    GATK4_CALCULATECONTAMINATION ( GATK4_GETPILEUPSUMMARIES_TUMOR.out.table, 
                                GATK4_GETPILEUPSUMMARIES_NORMAL.out.table)
    GATK4_FILTERMUTECTCALLS ( GATK4_MUTECT2.out.vcf, 
                            GATK4_CALCULATECONTAMINATION.out.contamination,
                            fasta, fai, dict)
    
    BCFTOOLS_NORM ( GATK4_FILTERMUTECTCALLS.out.vcf, fasta )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    emit:
    vcf      = BCFTOOLS_NORM.out.vcf        // channel: [ val(meta), path("*.vcf.gz") ]
    tbi      = BCFTOOLS_NORM.out.tbi        // channel: [ val(meta), path("*.tbi") ]
    stats    = GATK4_MUTECT2.out.stats      // channel: [ val(meta), path("*.stats") ]
    f1r2     = GATK4_MUTECT2.out.f1r2
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}
