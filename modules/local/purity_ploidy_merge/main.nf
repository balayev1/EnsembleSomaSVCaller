process PURITY_PLOIDY_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(ascat_purityploidy), path(facets_vcf), path(sequenza_solutions), val(aceseq_ploidy_purity)

    output:
    tuple val(meta), path("*.purity_ploidy_candidates.tsv"), emit: candidates
    tuple val(meta), path("*.purity_ploidy_consensus.tsv"),  emit: consensus
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aceseqArg = aceseq_ploidy_purity?.toString()?.trim() ? """cmd+=( --aceseq "${aceseq_ploidy_purity}" )""" : ''

    """
    cmd=(
        python3 "\${NEXTFLOW_PROJECT_DIR}/bin/merge_purity_ploidy.py"
        --sample-id "${meta.id}"
        --ascat "${ascat_purityploidy}"
        --facets "${facets_vcf}"
        --sequenza "${sequenza_solutions}"
    )
    ${aceseqArg}
    cmd+=(
        --candidates-out "${prefix}.purity_ploidy_candidates.tsv"
        --consensus-out "${prefix}.purity_ploidy_consensus.tsv"
    )

    "\${cmd[@]}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-EOF > ${prefix}.purity_ploidy_candidates.tsv
    sample_id\tcaller\tbase_caller\toption_rank\tpurity\tploidy\tcluster_label\tin_consensus\tsource_file
    ${meta.id}\tASCAT\tASCAT\t1\t0.35\t2.0\t0\ttrue\tstub
    ${meta.id}\tFACETS\tFACETS\t1\t0.35\t2.0\t0\ttrue\tstub
    EOF

    cat <<-EOF > ${prefix}.purity_ploidy_consensus.tsv
    sample_id\tstatus\tconsensus_cluster_label\tpurity\tploidy\tsupporting_callers\tsupporting_points\tn_supporting_callers\tn_supporting_points\teps\tmin_samples
    ${meta.id}\tconsensus\t0\t0.35\t2.0\tASCAT,FACETS\tASCAT,FACETS\t2\t2\t0.1\t2
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
