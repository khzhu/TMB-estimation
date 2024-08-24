process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(bwa) , emit: index
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        $args \\
        -p bwa/${prefix} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa

    touch bwa/${prefix}.amb
    touch bwa/${prefix}.ann
    touch bwa/${prefix}.bwt
    touch bwa/${prefix}.pac
    touch bwa/${prefix}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}