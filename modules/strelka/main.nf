process STRELKA_SOMATIC {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    input:
    tuple val(meta),  path(input_bams), path(input_index_files)
    tuple val(meta),  path(candidate_small_indels), path(candidate_small_indels_tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(target_bed), path(target_bed_tbi)

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : ""
    def options_manta = candidate_small_indels ? "--indelCandidates ${candidate_small_indels}" : ""
    """

    configureStrelkaSomaticWorkflow.py \\
        $args \\
        --tumor ${input_bams[0]} \\
        --normal ${input_bams[1]} \\
        --referenceFasta $fasta \\
        ${options_target_bed} \\
        ${options_manta} \\
        --runDir strelka

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

    python strelka/runWorkflow.py -m local -j $task.cpus

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.somatic_indels.vcf.gz
    touch ${prefix}.somatic_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.somatic_snvs.vcf.gz
    touch ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
}
