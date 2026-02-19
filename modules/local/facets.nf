#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FACETS {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://blcdsdockerregistry/cnv_facets:0.16.0':'blcdsdockerregistry/cnv_facets:0.16.0'}"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path(snp_vcf)                   // VCF file of SNPs where pileup is to be computed
    path(snp_vcf_index)             // Index file for the SNP VCF
    path(targets_bed)               // BED file of target regions to scan (optional)
    path(annotation_bed)            // BED file of annotation regions with 4th column as feature name (optional)
    val(mapq)                       // Minimum mapping quality for reads to be included in pileup: Default is 20
    val(baq)                        // Minimum base quality for bases to be included in pileup: Default is 20
    val(mindepth)                   // Minimum depth for a SNP in normal sample to be included in pileup: Default is 25
    val(maxdepth)                   // Maximum depth for a SNP in normal sample to be included in pileup: Default is 4000
    val(cval_pre)                   // Lower segmentation size limit: Default is 25 
    val(cval_proc)                  // Upper segmentation size limit: Default is 400
    val(nbhd_snp)                   // If an interval of size nbhd-snp contains more than one SNP, sample a random one: Default is "auto" if paired-end BAM provided, otherwise 250 
    val(gbuild)                     // Reference genome build
    val(count_orphans)              // If TRUE, do not discard anomalous read pairs: Default is FALSE
    val(unmatched)                  // If TRUE, normal sample is unmatched, use tumor reads to call heterozygous SNPs: Default is FALSE
    val(no_cov_plot)                // If TRUE, do not generate coverage plots: Default is FALSE

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
    def mapq_args           = mapq                          ? "-mq ${mapq}" : ""
    def baq_args            = baq                           ? "-bq ${baq}"  : ""
    def depth_args          = (mindepth && maxdepth)        ? "-d ${mindepth} ${maxdepth}" : ""
    def cval_args           = (cval_pre && cval_proc)       ? "-cv ${cval_pre} ${cval_proc}" : ""
    def nbhd_snp_args       = nbhd_snp                      ? "-snp ${nbhd_snp}" : ""
    def gbuild_args         = gbuild                        ? "-g ${gbuild}" : ""
    def rnd_seed_args       = args.rnd_seed                 ? "-s ${prefix}" : ""

    // Flags
    def count_orphans_args  = count_orphans                 ? "-A" : ""
    def unmatched_args      = unmatched                     ? "-u" : ""
    def no_cov_plot_args    = no_cov_plot                   ? "-np" : ""
    
    // Optional file inputs
    def targets_args        = targets_bed                   ? "-T ${targets_bed}"    : ""
    def annotation_args     = annotation_bed                ? "-a ${annotation_bed}" : ""

    """
    cnv_facets.R \\
        -n ${input_normal} \\
        -t ${input_tumor} \\
        -vcf ${snp_vcf} \\
        -o ${prefix} \\
        ${mapq_args} \\
        ${baq_args} \\
        -N ${task.cpus} \\
        ${depth_args} \\
        ${cval_args} \\
        ${nbhd_snp_args} \\
        ${gbuild_args} \\
        ${rnd_seed_args} \\
        ${count_orphans_args} \\
        ${unmatched_args} \\
        ${no_cov_plot_args} \\
        ${targets_args} \\
        ${annotation_args}

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