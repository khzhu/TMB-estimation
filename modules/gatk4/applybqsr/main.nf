process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals $intervals" : ""

    """
    /gatk/gatk ApplyBQSR \\
        --input $input \\
        --output ${prefix}.${input.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file $bqsr_table \\
        $interval_command \\
        --tmp-dir . \\
        --add-output-sam-program-record \\
        --create-output-bam-md5 \\
        --use-original-qualities \\
        --create-output-bam-index \\
        --QUIET \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}