process BRASS {
    tag "${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "quay.io/wtsicgp/brass:v6.3.4"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai), path(tumor_bas), \
                     path(normal_bas), path(ascat_summary)  // Mandatory!
    path(fasta)                 // Path to reference genome FASTA
    path(fasta_fai)             // Path to reference genome index FAI
    path(vagrent_dir)           // Path to genome cache directory
    path(depth)                 // Regions of excessive sequencing depth to be ignored
    path(viral)                 // Virus sequences from NCBI
    path(microbe_files)         // Microbe sequence files from NCBI 
    path(gcbins)                // BED file with 5 columns, col 4 number of non-N bases, col 5 GC fraction.
    path(cytoband)              // Path to cytoband file for a species build 
    path(centtel)               // TSV file of usable regions of the chromosome   
    path(brass_np)              // Path to panel of normals file
    path(brass_np_tbi)          // Path to panel of normals index file
    val(species)                // Species name
    val(protocol)               // Sequencing protocol
    val(assembly)               // Assembly version
    val(min_reads)              // Minimum reads to call group: Default is 2
    val(mincn)                  // Minimum CN change for copynumber_flag: Default is 0.3

    output:
    tuple val(meta), path("**/*.bedpe.gz")  , emit: bedpe
    tuple val(meta), path("**/*.vcf.gz")    , emit: vcf
    path("versions.yml")                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"

    def microbe_file   = microbe_files instanceof List ? microbe_files[0] : microbe_files
    def microbe_prefix = microbe_file.name.replaceAll(/\.\d+\.fa\.2bit$/, '')

    def vagr_dir = vagrent_dir instanceof List ? vagrent_dir[0] : vagrent_dir

    """
    # BRASS/HTS expects BAM indexes to be discoverable as <bam>.bai
    if [[ ! -f "${normal_bam}.bai" ]]; then
        ln -sf "${normal_bai}" "${normal_bam}.bai"
    fi
    if [[ ! -f "${tumor_bam}.bai" ]]; then
        ln -sf "${tumor_bai}" "${tumor_bam}.bai"
    fi

    brass.pl \\
        -o ${prefix} \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -g ${fasta} \\
        -d ${depth} \\
        -gc ${vagr_dir}/vagrent.cache.gz \\
        -vi ${viral} \\
        -mi ${microbe_prefix} \\
        -b ${gcbins} \\
        -cb ${cytoband} \\
        -pr ${protocol} \\
        -ct ${centtel} \\
        -ss ${ascat_summary} \\
        -s ${species} \\
        -as ${assembly} \\
        -c ${task.cpus} \\
        -f ${brass_np} \\
        -j ${min_reads} \\
        -cn ${mincn}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        brass: \$(brass.pl -v | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo "stub" | gzip > ${prefix}/${prefix}.bedpe.gz
    echo "stub" | gzip > ${prefix}/${prefix}.vcf.gz
    ln -sf ${prefix}/${prefix}.bedpe.gz .
    ln -sf ${prefix}/${prefix}.vcf.gz .
    touch versions.yml
    """
}
