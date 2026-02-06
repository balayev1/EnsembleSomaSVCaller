#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FACETS {
    tag "${meta.id}"
    label 'process_medium'

    container "docker://blcdsdockerregistry/cnv_facets:0.16.0"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path snp_vcf
    path snp_vcf_index
    path targets_bed, optional: true
    path annotation_bed, optional: true

    output:
    tuple val(meta), path("*.csv.gz"),           emit: pileup
    tuple val(meta), path("*.vcf.gz"),           emit: vcf
    tuple val(meta), path("*.cnv.png"),          emit: cnv_plot, optional: true
    tuple val(meta), path("*.cov.pdf"),          emit: cov_plot, optional: true
    tuple val(meta), path("*.spider.pdf"),       emit: spider_plot, optional: true
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Optional parameters with defaults
    def mapq          = args.snp_mapq          ?: 20
    def baq           = args.snp_baq           ?: 20
    def nprocs        = args.snp_nprocs        ?: task.cpus
    def depth_min     = args.depth_min         ?: 25
    def depth_max     = args.depth_max         ?: 4000
    def cval_pre      = args.cval_pre          ?: 25
    def cval_proc     = args.cval_proc         ?: 150
    def nbhd_snp      = args.nbhd_snp          ?: "auto"
    def gbuild        = args.gbuild            ?: "hg38"
    def rnd_seed      = args.rnd_seed          ?: "${prefix}"

    // Flags
    def count_orphans = args.snp_count_orphans ? "-A" : ""
    def unmatched     = args.unmatched         ? "-u" : ""
    def no_cov_plot   = args.no_cov_plot       ? "-np" : ""
    
    // Optional file inputs
    def targets_arg    = targets_bed    ? "-T ${targets_bed}"    : ""
    def annotation_arg = annotation_bed ? "-a ${annotation_bed}" : ""

    """
    cnv_facets.R \\
        -n ${input_normal} \\
        -t ${input_tumor} \\
        -vcf ${snp_vcf} \\
        -o ${prefix} \\
        -mq ${mapq} \\
        -bq ${baq} \\
        -N ${nprocs} \\
        -d ${depth_min} ${depth_max} \\
        -cv ${cval_pre} ${cval_proc} \\
        -snp ${nbhd_snp} \\
        -g ${gbuild} \\
        -s ${rnd_seed} \\
        ${count_orphans} \\
        ${unmatched} \\
        ${no_cov_plot} \\
        ${targets_arg} \\
        ${annotation_arg}

    # Version export
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: \$(cnv_facets.R --version 2>&1 | grep -oP 'Version \\K[0-9.]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv.gz
    touch ${prefix}.vcf.gz
    touch ${prefix}.cnv.png
    touch ${prefix}.cov.pdf
    touch ${prefix}.spider.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: stub
    END_VERSIONS
    """
}