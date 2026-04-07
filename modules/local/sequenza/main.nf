process SEQUENZAUTILS_GCWIGGLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39he88f293_8':
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39he88f293_8' }"

    input:
    tuple val(meta), path(fasta)        // Mandatory: Format should be [meta, fasta]
    val(windowsize)                     // Windowsize for GC-content % calculation: Default is 50

    output:
    tuple val(meta), path("*.wig.gz"), emit: wig
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    def prefix = task.ext.prefix ?: "${meta.id}"

    def window_arg  = windowsize ? "-w ${windowsize}" : ""

    """
    sequenza-utils \\
        gc_wiggle \\
        --fasta ${fasta} \\
        -o ${prefix}.wig.gz \\
        $window_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(sequenza-utils --version 2>&1 | sed 's/sequenza-utils //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.wig.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(sequenza-utils --version 2>&1 | sed 's/sequenza-utils //')
    END_VERSIONS
    """
}

process SEQUENZA_PREP {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39he88f293_8':
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39he88f293_8' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai), val(chr)
    path(fasta)             // Path to reference genome FASTA
    path(gc_wig)            // Path to GC-content wiggle file
    val(hom)                // Threshold to select homozygous positions: Default is 0.9
    val(het)                // Threshold to select heterozygous positions: Default is 0.25
    val(het_f)              // Frequency threshold in forward strand for heterozygous calls: Default is -0.02 (disabled, effective with values >= 0)
    val(qlimit)             // Minimum base quality score: Default is 20
    val(qformat)            // Quality score format: Default is sanger (Phred+33), can also be illumina (Phred+64)
    val(rd_thr)             // Threshold to filter positions by the sum of read depth of the two samples: Default is 20
    val(seqz_bin_size)      // Bin size for binning the seqz file: Default is 50

    output:
    tuple val(meta), path("*.binned.seqz.gz"), emit: seqz_file
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}_${chr}"

    def fasta_arg   = fasta   ? "--fasta ${fasta}" : ""
    def gc_arg      = gc_wig  ? "-gc ${gc_wig}" : ""
    def hom_arg     = hom     ? "--hom ${hom}" : ""
    def het_arg     = het     ? "--het ${het}" : ""
    def het_f_arg   = het_f   ? "--het_f ${het_f}" : ""
    def qlimit_arg  = qlimit  ? "--qlimit ${qlimit}" : ""
    def qformat_arg = qformat ? "--qformat ${qformat}" : ""
    def rd_thr_arg  = rd_thr  ? "-N ${rd_thr}" : ""
    def win_arg     = seqz_bin_size ? "--window ${seqz_bin_size}" : ""

    """
     sequenza-utils bam2seqz \\
        -n ${normal_bam} \\
        -t ${tumor_bam} \\
        ${fasta_arg} \\
        ${gc_arg} \\
        ${hom_arg} \\
        ${het_arg} \\
        ${het_f_arg} \\
        ${qlimit_arg} \\
        ${qformat_arg} \\
        ${rd_thr_arg} \\
        -C ${chr} | \\
    sequenza-utils seqz_binning \\
        ${win_arg} \\
        -s - | gzip > ${prefix}.binned.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(sequenza-utils --version 2>&1 | sed 's/sequenza-utils //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${chr}"
    """
    echo | gzip > ${prefix}.binned.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(sequenza-utils --version 2>&1 | sed 's/sequenza-utils //')
    END_VERSIONS
    """
}

process SEQUENZA_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::sequenza-utils=3.0.0 bioconda::htslib=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39he88f293_8':
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39he88f293_8' }"

    input:
    tuple val(meta), path(seqz_files)

    output:
    tuple val(meta), path("*.merged.seqz.gz"),     emit: merged_seqz
    tuple val(meta), path("*.merged.seqz.gz.tbi"), emit: merged_seqz_tbi
    path "versions.yml"                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\n' ${seqz_files} | sort -V > seqz_files.list

    python3 - <<'PY'
import gzip

paths = []
with open("seqz_files.list") as handle:
    for raw in handle:
        path = raw.strip()
        if path:
            paths.append(path)

with open("merged.seqz", "w") as out:
    wrote_header = False
    for path in paths:
        with gzip.open(path, "rt", errors="replace", newline="") as fh:
            for idx, line in enumerate(fh):
                if idx == 0 and line.startswith("chromosome"):
                    if not wrote_header:
                        out.write(line if line.endswith("\\n") else line + "\\n")
                        wrote_header = True
                    continue
                out.write(line if line.endswith("\\n") else line + "\\n")
PY

    bgzip -c merged.seqz > ${prefix}.merged.seqz.gz

    tabix -s 1 -b 2 -e 2 -S 1 ${prefix}.merged.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        bgzip: \$(echo \$(bgzip -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.merged.seqz.gz
    touch ${prefix}.merged.seqz.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: stub
        bgzip: stub
    END_VERSIONS
    """
}

process SEQUENZA_RUN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "oras://quay.io/balay011/r-sequenza-hg38:3.0.0-copynumber-2e31d59"

    input:
    tuple val(meta), path(seqz_file)
    val(assembly)              // Genome assembly version (hg19/hg38): Default is hg38
    val(rd_thr_tumor)          // Minimum read depth to consider a position in tumor sample: Default is 40
    val(rd_thr_normal)         // Minimum read depth to consider a position in normal sample: Default is 10
    val(purity_range)          // Range of cellularity values to explore: Default is 0.1-1
    val(ploidy_range)          // Range of ploidy values to explore: Default is 1-7

    output:
    tuple val(meta), path("*_alternative_fit.pdf")          , emit: alternative_fit_plot
    tuple val(meta), path("*_alternative_solutions.txt")    , emit: purity_ploidy_est
    tuple val(meta), path("*_confints_CP.txt")              , emit: confints_CP
    tuple val(meta), path("*_mutations.txt")                , emit: mutations
    tuple val(meta), path("*_segments.txt")                 , emit: segments
    tuple val(meta), path("*_chromosome_depths.pdf")        , emit: chromosome_depths_plot
    tuple val(meta), path("*_chromosome_view.pdf")          , emit: chromosome_view_plot
    tuple val(meta), path("*_CN_bars.pdf")                  , emit: CN_bars_plot
    tuple val(meta), path("*_CP_contours.pdf")              , emit: CP_contours_plot
    tuple val(meta), path("*_gc_plots.pdf")                 , emit: gc_plots
    tuple val(meta), path("*_genome_view.pdf")              , emit: genome_view_plot
    tuple val(meta), path("*_model_fit.pdf")                , emit: model_fit_plot
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def gender = meta.gender ?: 'male'
    def assembly_arg  = assembly      ? ", assembly = \"${assembly}\""          : ""
    def tumor_rd_arg  = rd_thr_tumor  ? ", min.reads = ${rd_thr_tumor}"         : ""
    def normal_rd_arg = rd_thr_normal ? ", min.reads.normal = ${rd_thr_normal}" : ""

    """
    #!/usr/bin/env Rscript
    library(sequenza)
    Sys.setenv(VROOM_CONNECTION_SIZE = "2147483647")
    
    purity_tmp <- as.numeric(unlist(strsplit("${purity_range}", "-")))
    cellularity_seq <- seq(from = purity_tmp[1], to = purity_tmp[2], by = 0.01)

    ploidy_tmp <- as.numeric(unlist(strsplit("${ploidy_range}", "-")))
    ploidy_seq <- seq(from = ploidy_tmp[1], to = ploidy_tmp[2], by = 0.1)

    is_female <- FALSE
    if ("${gender}" == "female") {
        is_female <- TRUE
    }

    seqz_input <- "${prefix}.merged.seqz"
    decompress_status <- system2("gzip", c("-dc", "${seqz_file}"), stdout = seqz_input)
    if (decompress_status != 0) {
        stop("Failed to decompress seqz input: ${seqz_file}")
    }

    data <- sequenza.extract(
        file = seqz_input,
        verbose = FALSE,
        parallel = ${task.cpus}
        ${assembly_arg}
        ${tumor_rd_arg}
        ${normal_rd_arg}
    )

    fit <- sequenza.fit(
        data, 
        female = is_female,
        cellularity = cellularity_seq, 
        ploidy = ploidy_seq,
        mc.cores = ${task.cpus}
    )

    sequenza.results(
        sequenza.extract = data, 
        cp.table = fit, 
        sample.id = "${prefix}", 
        out.dir = ".",
        female = is_female
    )

    sequenza_version <- as.character(packageVersion("sequenza"))
    writeLines(
        c(
            paste0('"', "${task.process}", '":'),
            paste0("    sequenza: ", sequenza_version)
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_alternative_fit.pdf
    touch ${prefix}_alternative_solutions.txt
    touch ${prefix}_confints_CP.txt
    touch ${prefix}_mutations.txt
    touch ${prefix}_segments.txt
    touch ${prefix}_chromosome_depths.pdf
    touch ${prefix}_chromosome_view.pdf
    touch ${prefix}_CN_bars.pdf
    touch ${prefix}_CP_contours.pdf
    touch ${prefix}_gc_plots.pdf
    touch ${prefix}_genome_view.pdf
    touch ${prefix}_model_fit.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza: \$(Rscript -e "cat(as.character(packageVersion('sequenza')))")
    END_VERSIONS
    """
}
