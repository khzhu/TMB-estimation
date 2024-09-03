process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam_file)
    path intervals
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)
    path  snp_known_sites
    path  snp_known_sites_tbi
    path  indel_known_sites
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
        --input $bam_file \\
        --output ${prefix}.table \\
        --reference $fasta \\
        $interval_command \\
        $snp_sites_command \\
        $indel_sites_command\\
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
