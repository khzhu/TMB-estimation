process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path  fasta
    path  fai
    path  dict
    path  snp_known_sites
    path  snp_known_sites_tbi
    path  indel_known_sites
    path  indel_known_sites_tbi

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"
    def prefix2 = task.ext.prefix2 ?: "${input.baseName}.recal"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def snp_sites_command = snp_known_sites.collect{"--known-sites $it"}.join(' ')
    def indel_sites_command = indel_known_sites.collect{"--known-sites $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibrator  \\
        --input $input \\
        --output ${prefix}.table \\
        --reference $fasta \\
        $interval_command \\
        $snp_sites_command \\
        $indel_sites_command\\
        --tmp-dir . \\
        $args

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input $input \\
        --output ${prefix2}.${input.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file ${prefix}.table \\
        $interval_command \\
        --tmp-dir . \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.table
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
