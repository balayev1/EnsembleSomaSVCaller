#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SEQUENZAUTILS_GCWIGGLE {
    tag "$meta.id"
    label 'process_low'

    container "docker://drtomc/sequenza-utils:latest"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig.gz"), emit: wig
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    def prefix = task.ext.prefix ?: "${meta.id}"

    def window_arg = args.window ? "${args.window}" : "50"

    """
    sequenza-utils \\
        gc_wiggle \\
        --fasta ${fasta[0]} \\
        -o ${prefix}.wig.gz \\
        -w $window_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}

process SEQUENZA_RUN {
    tag "$meta.id"
    label 'process_high'

    container "docker://sequenza/sequenza"

    input:
    tuple val(meta), path(normalbam), path(normalbai), path(tumourbam), path(tumourbai)
    path fasta
    path wigfile

    output:
    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    prefix = args.prefix ? "${args.prefix}" : "${meta.id}"
    
    def options = []

    // Optional arguments
    if (args.store_seqztmp)  options << "--store_seqztmp"
    if (args.ignore_normal)  options << "--ignore_normal"
    if (args.ratio_priority) options << "--ratio_priority"
    if (args.no_archive)     options << "--no_archive"
    if (meta.gender =~ /(?i)XX/) options << "--x-heterozygous"

    if (args.cellularity)    options << "--cellularity ${args.cellularity}"
    if (args.ploidy)         options << "--ploidy ${args.ploidy}"
    if (args.breaks_file)    options << "--breaks ${args.breaks_file}"

    def cmd_bool_options = options.join(' ')

    def bin_size = args.bin_size ?: "50"
    def cellularity_range = args.cellularity_range ?: "0-1"
    def ploidy_range = args.ploidy_range ?: "1-7"

    """
    # Set up temporary directory
    if [ -z "${args.tmp_dir ?: ''}" ]; then
        export TMPDIR="\$(pwd)/tmp_sequenza_${meta.id}"
    else
        export TMPDIR="${args.tmp_dir}"
    fi

    mkdir -p "\$TMPDIR"

    # Run Sequenza Pipeline
    sequenza-pipeline \\
    --sample-id ${prefix} \\
    --normal-bam $normalbam \\
    --tumor-bam $tumourbam \\
    --reference-gz ${fasta[0]} \\
    --gc_wig $wigfile \\
    --bin ${args.bin_size ?: 50} \\
    --ncpu ${task.cpus} \\
    --mem ${task.memory.toGiga()} \\
    --cellularity-range ${args.cellularity_range ?: "0-1"} \\
    --ploidy-range ${args.ploidy_range ?: "1-7"} \\
    --tmp \$TMPDIR \\
    ${cmd_bool_options}

    # Clean up temporary directory
    if [[ "\$TMPDIR" == *"/tmp_sequenza_${meta.id}"* ]]; then
        rm -rf "\$TMPDIR"
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: \$(Rscript -e "library(sequenza); cat(as.character(packageVersion('sequenza')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: [:]
    prefix = args.prefix ? "${args.prefix}" : "${meta.id}"
    """
    mkdir -p ${prefix}_results
    touch ${prefix}_results/${prefix}_segments.txt
    touch ${prefix}_results/${prefix}_purity_ploidy.txt
    touch ${prefix}_results/${prefix}_alternative_fit.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: \$(Rscript -e "library(sequenza); cat(as.character(packageVersion('sequenza')))")
    END_VERSIONS
    """
}
