process VEP {

    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(input_vcf)
    tuple val(meta2), path(fasta)       // Required
    tuple val(meta3), path(gnomad_vcf), path(gnomad_vcf_tbi)
    tuple val(meta4), path(cosmic_vcf), path(cosmic_vcf_tbi)
    path vep_cache                      // Required for VEP running. A default of /.vep is supplied.

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${input_vcf.baseName}.vep"

    """
    gunzip -kd $input_vcf > ${input_vcf.baseName}.vcf

    vep \\
        --fork ${task.cpus} \\
        $args \\
        --custom $gnomad_vcf,gnomAD,vcf,exact,0,AF_POPMAX,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS
        --custom $cosmic_vcf,COSMIC,vcf,exact,0,ID
        --dir $vep_cache \\
        --fasta $fasta \\
        --input_file $input_vcf \\
        --output_file ${prefix}.vcf'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep -version)
    END_VERSIONS
    """
}
