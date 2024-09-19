//
// GATK4_MUTECT2: Call variants with Mutect2 (2.1 part of GATK 4)
//

include { MANTA_SOMATIC                         } from '../../modules/manta/main'
include { STRELKA_SOMATIC                       } from '../../modules/strelka/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SNV    } from '../../modules/bcftools/norm/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_INDEL  } from '../../modules/bcftools/norm/main'
include { GATK4_MERGEVCFS as STRELKA2_MERGEVCFS } from '../../modules/gatk4/mergevcfs/main'
include { VCF2MAF as STRELKA2_VEP2MAF           } from '../../modules/vcf2maf/main'

workflow SNV_STRELKA2 {

    take:
    ch_input_files  // channel: [mandatory] [ val(meta), path(ch_input_bams), path(ch_index_files) ]
    ch_fasta        // channel: [mandatory] [ val(meta2), path(fasta) ]
    ch_fai          // channel: [mandatory] [ val(meta2), path(fai) ]
    ch_dict         // channel: [mandatory] [ val(meta2), path(dict) ]
    ch_target_bed   // channel: [mandatory] [ val(meta3), path(bed), path(bed_tbi) ]
    vep_cache
    filter_vcf

    main:
    ch_versions         = Channel.empty()

    MANTA_SOMATIC ( ch_input_files,
                    ch_fasta, ch_fai,
                    ch_target_bed)
    ch_versions = ch_versions.mix(MANTA_SOMATIC.out.versions)

    STRELKA_SOMATIC  ( ch_input_files,
                        MANTA_SOMATIC.out.candidate_small_indels_vcf.combine(
                        MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi, by:0),
                        ch_fasta, ch_fai, ch_target_bed)
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    BCFTOOLS_NORM_SNV  ( STRELKA_SOMATIC.out.vcf_snvs.combine(
                        STRELKA_SOMATIC.out.vcf_snvs_tbi, by:0),
                        ch_fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_SNV.out.versions)

    BCFTOOLS_NORM_INDEL ( STRELKA_SOMATIC.out.vcf_indels.combine(
                        STRELKA_SOMATIC.out.vcf_indels_tbi, by:0),
                        ch_fasta)

    ch_input_vcfs = BCFTOOLS_NORM_SNV.out.vcf \
        | combine( BCFTOOLS_NORM_INDEL.out.vcf, by: 0 ) \
        | map { sample, snv_vcf, indel_vcf ->
            tuple( sample, [snv_vcf, indel_vcf] ) }
    STRELKA2_MERGEVCFS ( ch_input_vcfs, ch_dict )
    ch_versions = ch_versions.mix(STRELKA2_MERGEVCFS.out.versions)

    STRELKA2_VEP2MAF ( STRELKA2_MERGEVCFS.out.vcf, ch_fasta, vep_cache, filter_vcf)
    ch_versions = ch_versions.mix(STRELKA2_VEP2MAF.out.versions)

    emit:
    vcf      = STRELKA2_MERGEVCFS.out.vcf    // channel: [ val(meta), path("*.vcf.gz") ]
    tbi      = STRELKA2_MERGEVCFS.out.tbi    // channel: [ val(meta), path("*.tbi") ]
    maf      = STRELKA2_VEP2MAF.out.maf      // channel: [ val(meta), path("*.maf") ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}
