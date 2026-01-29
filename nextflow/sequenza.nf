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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        gc_wiggle \\
        $args \\
        --fasta $fasta \\
        -o ${prefix}.wig.gz

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

process SEQUENZAUTILS_BAM2SEQZ {
    tag "$meta.id"
    label 'process_medium'

    container "docker://drtomc/sequenza-utils:latest"

    input:
    tuple val(meta), path(normalbam), path(tumourbam)
    path fasta
    path wigfile

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        bam2seqz \\
        $args \\
        -n $normalbam \\
        -t $tumourbam \\
        --fasta $fasta \\
        -gc $wigfile \\
        -o ${prefix}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}

process SEQUENZAUTILS_BINNING {
    tag "$meta.id"
    label 'process_medium'

    container "docker://drtomc/sequenza-utils:latest"

    input:
    tuple val(meta), path(seqz)

    output:
    tuple val(meta), path("*.seqz.gz")     , emit: binned_seqz
    tuple val(meta), path("*.seqz.gz.tbi") , emit: index, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:] 
    def prefix = task.ext.prefix ?: "${meta.id}"

    def window_arg = args.window ? "${args.window}" : "50"
    def output_arg = args.output ? "${args.output}" : "${prefix}.bin.seqz.gz"
    def tabix_arg  = args.tabix  ? "${args.tabix}"  : "tabix"

    """
    sequenza-utils \\
        seqz_binning \\
        -s $seqz \\
        -w $window_arg \\
        -o $output_arg \\
        -T $tabix_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bin.seqz.gz
    touch ${prefix}.bin.seqz.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}

process SEQUENZA_RUN {
    tag "${meta.id}_g${gamma}"
    label 'process_medium'

    container "docker://sequenza/sequenza"

    input:
    tuple val(meta), path(seqz), val(gamma)

    output:
    tuple val(meta), path("${prefix}*"), emit: results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    // We define the prefix here so it's available for both the script and the output block
    prefix = args.prefix ? "${args.prefix}_gamma${gamma}" : "${meta.id}_gamma${gamma}"
    
    // Gender logic based on meta
    def female_flag = (meta.gender =~ /(?i)XX/) ? "TRUE" : "FALSE"

    // Map-style argument setup
    def ref_assembly   = args.ref_assembly   ?: "hg38"
    def genome_size    = args.genome_size    ?: 23
    def breaks_method  = args.breaks_method  ?: "het"
    def fit_method     = args.fit_method     ?: "baf"
    def ploidy_file    = args.ploidy_file    ?: "NULL"
    def window         = args.window         ?: 50
    def type           = args.type           ?: "all"
    def maxvar         = args.maxvar         ?: 20
    def min_reads_n    = args.min_reads_normal ?: 10
    def min_reads_b    = args.min_reads_baf    ?: 1
    def height         = args.height         ?: 440
    def width          = args.width          ?: 1440

    """
    Rscript /path/to/your/SEQUENZA_SCRIPT.R \\
        --seqz_file $seqz \\
        --prefix $prefix \\
        --gamma $gamma \\
        --female $female_flag \\
        --ref_assembly $ref_assembly \\
        --genome_size $genome_size \\
        --breaks_method $breaks_method \\
        --fit_method $fit_method \\
        --ploidy_file $ploidy_file \\
        --window $window \\
        --type $type \\
        --maxvar $maxvar \\
        --min.reads.normal $min_reads_n \\
        --min.reads.baf $min_reads_b \\
        --height $height \\
        --width $width

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: \$(Rscript -e "library(sequenza); cat(as.character(packageVersion('sequenza')))")
    END_VERSIONS
    """

    stub:
    prefix = "${meta.id}_gamma${gamma}"
    """
    touch ${prefix}_segments.txt
    touch ${prefix}_purity_ploidy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: 1.0.0
    END_VERSIONS
    """
}