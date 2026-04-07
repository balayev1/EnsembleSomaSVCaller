process VALIDATE_ACESEQ_MANIFEST {
    executor 'local'

    input:
    path samplesheet
    path manifest

    output:
    path 'validated.ok', emit: ready

    script:
    """
    validate_aceseq_manifest.py \
        --samplesheet ${samplesheet} \
        --manifest ${manifest}

    touch validated.ok
    """
}
