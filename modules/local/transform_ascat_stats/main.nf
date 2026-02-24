process TRANSFORM_ASCAT_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(ascat_purityploidy)

    output:
    tuple val(meta), path("*.brass.summary.stats.txt"), emit: brass_stats
    path "versions.yml"                               , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gender = meta.gender ?: meta.sex ?: "male" 
    
    """
    transform_ascat_purityploidy.py \\
        ${ascat_purityploidy} \\
        ${gender} \\
        ${prefix}.brass.summary.stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gender = meta.gender ?: meta.sex ?: "male" 
    def gender_chr = gender.toLowerCase() == 'female' ? 'X' : 'Y'
    def gender_found = gender.toLowerCase() == 'female' ? 'N' : 'Y'

    """
    cat <<-EOF > ${prefix}.brass.summary.stats.txt
    rho 0.5
    Ploidy 2.0
    GenderChr ${gender_chr}
    GenderChrFound ${gender_found}
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}