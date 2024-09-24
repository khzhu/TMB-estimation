//
// TMB_ESTIMATION: Derive TMB estimates from targeted exome capture seqeueuncing data
//

process TMB_CALIBRATION {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(maf_file)
	
    output:
    tuple val(meta), path("*.tsv"), emit: tmb
    path "versions.yml"             , emit: versions
	

    //Calculate tumor mutation burden per sample
    script:
    """
	R $projectDir/bin/calculate_tmb.R -m ${maf_file} -o ${prefix}.tmb.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1)
    END_VERSIONS
    """
}