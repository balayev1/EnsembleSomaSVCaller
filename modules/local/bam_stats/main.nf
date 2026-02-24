process BAM_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "quay.io/wtsicgp/pcap-core:5.8.1"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta_fai)     // Path to reference genome index FAI

    output:
    tuple val(meta), path("${bam}.bas"), emit: bas
    path "versions.yml", emit: versions

    script:
    """
    bam_stats \\
        -i ${bam} \\
        -o ${bam}.bas \\
        -r ${fasta_fai} \\
        --num_threads ${task.cpus}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcap-core: \$(bam_stats -v | awk '/version/ {print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bas
    touch versions.yml
    """
}