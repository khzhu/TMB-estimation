process SAMBAMBA_MERGE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bams.collect{ "$it"}.join(' ')

    """
    sambamba \\
        merge $args \\
        -t $task.cpus \\
        ${prefix}.bam \\
        $input_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}