process CREATE_DELLY_SAMPLES_TSV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bcftools=1.23.1" : null)
    container "${workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(bcf), path(csi)

    output:
    tuple val(meta), path("${meta.id}_delly_samples.tsv"), emit: tsv

    script:
    """
    bcftools query -l ${bcf} > delly_samples.list

    awk 'NR==1{print \$0"\ttumor"} NR==2{print \$0"\tcontrol"} END{if (NR != 2) {print "Expected exactly 2 DELLY samples, found " NR > "/dev/stderr"; exit 1}}' delly_samples.list > ${meta.id}_delly_samples.tsv
    """
}
