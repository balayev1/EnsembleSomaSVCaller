process VCF_TO_PLAIN {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bcftools=1.20" : null)
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | head -n 1 | awk '{print \$2}'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: vcf.name.replaceFirst(/\.vcf(\.(gz|bgz))?$/, '')
    def output_vcf = "${prefix}.vcf"

    if ("${vcf}" == "${output_vcf}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    if [[ "${vcf}" == *.gz || "${vcf}" == *.bgz ]]; then
        gzip -dc ${vcf} > ${output_vcf}
    else
        cp ${vcf} ${output_vcf}
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "output"
    """
    touch ${prefix}.vcf
    """
}
