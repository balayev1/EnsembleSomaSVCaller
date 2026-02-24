
process CREATE_DELLY_SAMPLES_TSV {
    tag "$meta.id"
    label 'process_low'

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}_delly_samples.tsv"), emit: tsv

    script:
    """
    echo -e "${meta.tumor_id}\ttumor" > ${meta.id}_delly_samples.tsv
    echo -e "${meta.control_id}\tcontrol" >> ${meta.id}_delly_samples.tsv
    """
}