#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { INDEX_GENOME             } from './subworkflows/index_genome/main'
include { FASTQ_FASTP_FASTQC       } from './subworkflows/fastq_fastp_fastqc/main'
include { ALIGN_MARKDUP_BQSR_STATS } from './subworkflows/align_markdup_bqsr_stats/main'
include { SNV_MUTECT2              } from './subworkflows/snv_mutect2/main'
include { SNV_STRELKA2             } from './subworkflows/snv_strelka2/main'
include { TMB_CALIBER              } from './subworkflows/calculate_tmb/main'

// main workflow
workflow {
    log.info """\
    TMB estimation pipeline ${params.release}
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sample sheet  : ${params.input_json}
    output dir    : ${params.output_dir}
    """

    // Read sample config file
    def jsonSlurper = new JsonSlurper()

    // If params.input_json exists, use that for multi sample processing.
    // Otherwise, default to single sample analysis.
    multi_params = params.input_json ? jsonSlurper.parse(new File(params.input_json)).collect{params + it} : [params]
    output_dir = params.output_dir ? params.output_dir : "."
    
    // run samples through the pipeline
    samples = Channel.from(multi_params.collect{ it -> tuple([
                id: it.specimen_id, pid: it.patient_id, single_end:false, tissue: it.tissue, purity: it.purity ],
                [ file(it.read1, checkIfExists: true), file(it.read2, checkIfExists: true) ]) })
    ch_versions = Channel.empty()

    // Index genome reference
    INDEX_GENOME ( [[ id:'genome'], file(params.reference_file, checkIfExists: true)])

    // Trim sequence reads with paired-end data
    FASTQ_FASTP_FASTQC ( samples,
                        file(params.adapter_fasta, checkIfExists:true),
                        params.save_trimmed_fail,
                        params.save_merged,
                        params.skip_fastp,
                        params.skip_fastqc )
    ch_versions = ch_versions.mix( FASTQ_FASTP_FASTQC.out.versions )

    // Align reads to a reference genome
    ALIGN_MARKDUP_BQSR_STATS ( FASTQ_FASTP_FASTQC.out.reads,
                INDEX_GENOME.out.index,
                file(params.reference_file, checkIfExists: true),
                params.val_sort_bam,
                file(params.exome_plus_tumor_panel_bed, checkIfExists: true),
                INDEX_GENOME.out.fai,
                INDEX_GENOME.out.dict,
                file(params.known_snp_vcf, checkIfExists: true),
                file(params.known_snp_vcf_tbi, checkIfExists: true),
                file(params.known_indel_vcf, checkIfExists: true),
                file(params.known_indel_vcf_tbi, checkIfExists: true),
                params.split_reads)
    ch_versions = ch_versions.mix( ALIGN_MARKDUP_BQSR_STATS.out.versions )
    split_sample_bam = ALIGN_MARKDUP_BQSR_STATS.out.bam
                    .branch{
                            tumor   : it[1] =~ /(\_T.sort\.dedup\.recal\.bam)$/
                            normal  : it[1] =~ /(\_N.sort\.dedup\.recal\.bam)$/
                    }
    split_sample_bai = ALIGN_MARKDUP_BQSR_STATS.out.bai
                    .branch{
                            tumor   : it[1] =~ /(\_T.sort\.dedup\.recal\.bai)$/
                            normal  : it[1] =~ /(\_N.sort\.dedup\.recal\.bai)$/
                    }
    ch_paired_bam = split_sample_bam.tumor.combine(split_sample_bam.normal, by: 1)
    ch_paired_bai = split_sample_bai.tumor.combine(split_sample_bai.normal, by: 1)

    bed_files = Channel.fromPath(params.tumor_panel_bed_files, checkIfExists: true)
    ch_input_files = ch_paired_bam.combine(ch_paired_bai, by:1).combine(bed_files)
                    .map { meta, input_bams, input_index_files, intervals ->
                        new_meta = meta.clone()
                        new_meta.sid = intervals.baseName != "no_intervals" ? new_meta.id + "_" + intervals.baseName : new_meta.id
                        intervals = intervals.baseName != "no_intervals" ? intervals : []
                        [new_meta, input_bams, input_index_files, intervals]
                    }
    // calling mutect2 somatic variants
    SNV_MUTECT2 (ch_input_files,
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                [[id:'gnomad'],file(params.gnomad_exome_vcf, checkIfExists: true),
                    file(params.gnomad_exome_vcf_tbi, checkIfExists: true)],
                [[id:'exac'], file(params.exac_common_vcf, checkIfExists: true),
                    file(params.exac_common_vcf_tbi, checkIfExists: true)],
                [[id: 'cosmic'], file(params.cosmic_vcf, checkIfExists: true),
                    file(params.cosmic_vcf_tbi, checkIfExists: true)],
                file(params.vep_cache, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_MUTECT2.out.versions )

    // calling strelka2 somatic variants
    SNV_STRELKA2 (ch_paired_bam.combine(ch_paired_bai, by:1),
                [[ id:'genome'], file(params.reference_file, checkIfExists: true)],
                [[ id:'genome'], file(params.fai_file, checkIfExists: true)],
                [[ id:'genome'], file(params.dict_file, checkIfExists: true)],
                [[ id:'genome'], file(params.tumor_panel_bed_gz, checkIfExists: true),
                                    file(params.tumor_panel_bed_tbi, checkIfExists: true)],
                [[id: 'gnomad'], file(params.gnomad_exome_vcf, checkIfExists: true),
                                    file(params.gnomad_exome_vcf_tbi, checkIfExists: true)],
                [[id: 'cosmic'], file(params.cosmic_vcf, checkIfExists: true),
                                    file(params.cosmic_vcf_tbi, checkIfExists: true)],
                file(params.vep_cache, checkIfExists: true))
    ch_versions = ch_versions.mix( SNV_STRELKA2.out.versions )

    // collecting sample MAF files
    mutect2_maf = SNV_MUTECT2.out.maf
        .map { it -> it[1] }
        .collectFile( name: 'mutect2_merged.maf', keepHeader:true, skip:2, storeDir:params.store_dir )
        .map { [ [ id:'f1'], it ] }

    strelka2_maf = SNV_STRELKA2.out.maf
        .map { it -> it[1] }
        .collectFile( name: 'strelka2_merged.maf', keepHeader:false, skip:2, storeDir:params.store_dir )
        .map { [ [ id:'f1'], it ] }

    TMB_CALIBER ( mutect2_maf, strelka2_maf )
    ch_versions = ch_versions.mix(TMB_CALIBER.out.versions)
}