process VCF2MAF {

    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(input_vcf) // Use an uncompressed VCF file!
    path fasta                 // Required
    path vep_cache             // Required for VEP running. A default of /.vep is supplied.

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def vep_cache_cmd = vep_cache       ? "--vep-data $vep_cache" : ""
    def args2 = task.ext.args2          ?: ''
    def tumor_id = tumorID              ?: "${meta.id}_T"
    def normal_id = normalID            ?: "${meta.id}_N"
    // If VEP is present, it will find it and add it to commands.
    // If VEP is not present they will be blank
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if command -v vep &> /dev/null
    then
        VEP_CMD="--vep-path \$(dirname \$(type -p vep))"
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_CMD=""
        VEP_VERSION=""
    fi

    gunzip -k $input_vcf

    vcf2maf.pl \\
        $args \\
        \$VEP_CMD \\
        $args2 \\
        $vep_cache_cmd \\
        --ref-fasta $fasta \\
        --input-vcf $input_vcf \\
        --tumor-id $tumor_id \\
        --normal-id $normal_id \\
        --vcf-tumor-id ${meta.id}_T \\
        --vcf-normal-id ${meta.id}_N \\
        --output-maf ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION\$VEP_VERSION
    END_VERSIONS
    """
}
