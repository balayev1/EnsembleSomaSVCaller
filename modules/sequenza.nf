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
    tuple val(meta), path(normalbam), path(tumourbam), path(normalbai), path(tumourbai)
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
    
    // Arguments
    def x_flag = (meta.gender =~ /(?i)XX/) ? "--x-heterozygous" : ""
    def store_seqz = args.store_seqztmp ? "--store_seqztmp" : ""
    def ignore_normal = args.ignore_normal ? "--ignore_normal" : ""
    def ratio_priority = args.ratio_priority ? "--ratio_priority" : ""
    def no_archive = args.no_archive ? "--no_archive" : ""
    def bin_size = args.bin_size ?: "50"
    def cellularity_range = args.cellularity_range ?: "0-1"
    def ploidy_range = args.ploidy_range ?: "1-7"
    def cellularity_arg = args.cellularity ? "--cellularity ${args.cellularity}" : ""
    def ploidy_arg = args.ploidy ? "--ploidy ${args.ploidy}" : ""
    def breaks_arg = args.breaks_file ? "--breaks ${args.breaks_file}" : ""
    def tmp_arg = args.tmp_dir ? "--tmp ${args.tmp_dir}" : ""

    """
    sequenza-pipeline \\
    --sample-id ${prefix} \\
    --normal-bam $normalbam \\
    --tumor-bam $tumourbam \\
    --normal-bam-index $normalbai \\
    --tumor-bam-index $tumourbai \\
    --reference-gz ${fasta[0]} \\
    --gc_wig $wigfile \\
    --bin $bin_size \\
    --ncpu ${task.cpus} \\
    --mem ${task.memory.toGiga()} \\
    $breaks_arg \\
    $x_flag \\
    $store_seqz \\
    $ignore_normal \\
    $ratio_priority \\
    $cellularity_arg \\
    $ploidy_arg \\
    --cellularity-range $cellularity_range \\
    --ploidy-range $ploidy_range \\
    $no_archive \\
    $tmp_arg

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
