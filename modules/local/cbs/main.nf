process CBS {

    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/cbs:latest':
        'mskilab/cbs:latest' }"

    input:
    tuple val(meta), path(tumor_dryclean_cov, stageAs: "tumor_drycleaned_cov.rds"), path(normal_dryclean_cov, stageAs: "normal_drycleaned_cov.rds")
    val(cnsignif)
    val(field)
    val(name)

    output:
    tuple val(meta), path("*.cov.rds")       , emit: cbs_cov_rds
    tuple val(meta), path("*.seg.rds")       , emit: cbs_seg_rds
    tuple val(meta), path("*.nseg.rds")      , emit: cbs_nseg_rds
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def sample_name = name ? name : "${meta.id}"
    def VERSION     = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export RSCRIPT_PATH=\$(echo "\${NEXTFLOW_PROJECT_DIR}/bin/cbsFH_ZC.R")
    Rscript \$RSCRIPT_PATH \\
        -t ${tumor_dryclean_cov} \\
        -n ${normal_dryclean_cov} \\
        --cnsignif ${cnsignif} \\
        -m ${sample_name} \\
        -f ${field}

    mv cov.rds ${prefix}.cov.rds
    mv seg.rds ${prefix}.seg.rds
    mv nseg.rds ${prefix}.nseg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.cov.rds
    touch ${prefix}.seg.rds
    touch ${prefix}.nseg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS
    """
}
