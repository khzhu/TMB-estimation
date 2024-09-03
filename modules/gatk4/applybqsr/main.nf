process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam_file)
    tuple val(meta), bqsr_table
    path  intervals
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.bam") , emit: bam
    tuple val(meta), path("*.bai") , emit: bai
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.sort.dedup.recal"
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input $bam_file \\
        --output ${prefix}.${bam_file.getExtension()} \\
        --reference $fasta \\
        --bqsr-recal-file $bqsr_table \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
