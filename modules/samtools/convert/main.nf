process SAMTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam")  , emit: bam ,   optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.bai")  , emit: bai ,   optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"

    """
    samtools \\
        faidx \\
        $fasta

    samtools view \\
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        $args \\
        $bam \\
        -o ${prefix}.${output_extension}

    samtools index -@${task.cpus} ${prefix}.${output_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_extension = input.getExtension() == "bam" ? "cram" : "bam"
    def index_extension = output_extension == "bam" ? "bai" : "crai"

    """
    touch ${prefix}.${output_extension}
    touch ${prefix}.${index_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}