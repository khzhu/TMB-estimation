process GATK4_CALCULATECONTAMINATION {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(pileup)
    tuple val(meta), path(matched)

    output:
    tuple val(meta), path('*.contamination.table'), emit: contamination
    tuple val(meta), path('*.segmentation.table') , emit: segmentation, optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${pileup.baseName}"
    def matched_command = matched ? "--matched-normal $matched" : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK CalculateContamination] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CalculateContamination \\
        --input $pileup \\
        --output ${prefix}.contamination.table \\
        $matched_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.contamination.table
    touch ${prefix}.segmentation.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
