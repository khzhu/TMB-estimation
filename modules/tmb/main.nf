//
// TMB_ESTIMATION: Derive TMB estimates from targeted exome capture seqeueuncing data
//

process TMB_CALIBRATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(mutect2_maf), path(strelka2_maf)
	
    output:
    tuple val(meta), path("*.tsv"), emit: tmb
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"             , emit: versions
	

    //Calculate tumor mutation burden per sample
    script:
    def prefix        = task.ext.prefix ?: "tmb"
    def out_maf       = 'all.maf'
    """
    cat $mutect2_maf $strelka2_maf > $out_maf
    Rscript bin/calculate_tmb.R -m $out_maf -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1)
    END_VERSIONS
    """
}