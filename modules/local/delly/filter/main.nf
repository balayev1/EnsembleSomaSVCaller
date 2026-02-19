process DELLY_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/delly:1.3.3--h4d20210_0' :
        'biocontainers/delly:1.3.3--h4d20210_0' }"
    
    input:
    tuple val(meta), path(bcf), path(bcf_csi), path(sample_tsv)
    val(mode)                       // Filtering mode: "somatic" or "germline"
    val(min_qual)                   // Minimum SV site quality: Default is 300
    val(altaf)                      // Minimum alternative allele frequency: Default is 0.3
    val(min_sv_size)                // Minimum SV size: Default is 0
    val(max_sv_size)                // Maximum SV size: Default is 500000000
    val(ratio_geno)                 // Fraction of genotyped samples: Default is 0.75
    val(pass)                       // If true, only retain SVs that pass all filters
    val(tags)                       // If true, tag filtered sites instead of removing them
    val(cov)                        // Minimum coverage in tumor sample: Default is 10
    val(contamination_in_control)   // Maximum alternative allele frequency in control: Default is 0
    val(geno_qual)                  // Minimum median genotype quality for carriers/non-carriers: Default is 15
    val(rd_del)                     // Maximum read depth fold-change for deletions in carriers vs non-carriers: Default is 0.8
    val(rd_dup)                     // Minimum read depth fold-change for duplications in carriers vs non-carriers: Default is 1.2

    output:
    tuple val(meta), path("*.{bcf}")  , emit: bcf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_delly_filtered"

    def mode_arg = mode ? "--filter ${mode}" : ""
    def min_qual_arg = min_qual ? "--quality ${min_qual}" : ""
    def altaf_arg = altaf ? "--altaf ${altaf}" : ""
    def min_sv_size_arg = min_sv_size ? "--minsize ${min_sv_size}" : ""
    def max_sv_size_arg = max_sv_size ? "--maxsize ${max_sv_size}" : ""
    def ratio_geno_arg = ratio_geno ? "--ratiogeno ${ratio_geno}" : ""
    def pass_arg = pass ? "--pass" : ""
    def tags_arg = tags ? "--tag" : ""

    def sample_tsv_arg = (mode == 'somatic' && sample_tsv) ? "--samples ${sample_tsv}" : ""
    def cov_arg = (mode == 'somatic' && cov) ? "--coverage ${cov}" : ""
    def contamination_in_control_arg = (mode == 'somatic' && contamination_in_control != null) ? "--controlcontamination ${contamination_in_control}" : ""

    def geno_qual_arg = (mode == 'germline' && geno_qual) ? "--gq ${geno_qual}" : ""
    def rd_del_arg = (mode == 'germline' && rd_del) ? "--rddel ${rd_del}" : ""
    def rd_dup_arg = (mode == 'germline' && rd_dup) ? "--rddup ${rd_dup}" : ""

    def bcf_output = "--outfile ${prefix}.bcf"

    """
    delly filter \\
        ${args} \\
        ${bcf_output} \\
        ${mode_arg} \\
        ${min_qual_arg} \\
        ${altaf_arg} \\
        ${min_sv_size_arg} \\
        ${max_sv_size_arg} \\
        ${ratio_geno_arg} \\
        ${pass_arg} \\
        ${tags_arg} \\
        ${sample_tsv_arg} \\
        ${cov_arg} \\
        ${contamination_in_control_arg} \\
        ${geno_qual_arg} \\
        ${rd_del_arg} \\
        ${rd_dup_arg} \\
        ${bcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | head -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    
    def bcf_output = "touch ${prefix}.bcf"

    """
    ${bcf_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$( echo \$(delly --version 2>&1) | sed 's/^.*Delly version: v//; s/ using.*\$//')
    END_VERSIONS
    """
}