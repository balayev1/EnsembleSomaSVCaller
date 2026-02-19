
process GRIDSS_SOMATIC_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"


    input:
    tuple val(meta), path(gridss_output), path(gridss_output_tbi)
    path(pondir_gridss)     // pon path

    output:
    tuple val(meta), path("*.high_confidence_somatic.vcf")          , emit: somatic_high_vcf,          optional:true
    tuple val(meta), path("*.high_and_low_confidence_somatic.vcf")  , emit: somatic_all_vcf,           optional:true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def output1       = "${prefix}.high_confidence_somatic.vcf"
    def output2       = "${prefix}.high_and_low_confidence_somatic.vcf"
    //def pon           = pondir_gridss ? "cp -s ${pondir_gridss}/* ." : ""
    //def scriptDir     = "dirname \$(which gridss_somatic_filter)".execute().text.trim()
    """

    gridss_somatic_filter \\
    --pondir ${pondir_gridss} \\
    --input ${gridss_output} \\
    --output $output1 \\
    --fulloutput $output2 \\
    -n 1 \\
    -t 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.high_confidence_somatic.vcf
    touch ${prefix}.high_and_low_confidence_somatic.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}