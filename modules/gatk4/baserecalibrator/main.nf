process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path  fasta
    path  fai
    path  dict
    path  snp_known_sites
    path  indel_known_sites
    path  snp_known_sites_tbi
    path  indel_known_sites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""
    def sites_command = known_sites.collect{"--known-sites $it"}.join(' ')

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
        --known-sites $snp_known_sites \\
        --known-sites $indel_known_sites \\
        --tmp-dir . \\
        $args

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