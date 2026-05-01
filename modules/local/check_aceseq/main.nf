process VALIDATE_ACESEQ_MANIFEST {
    executor 'local'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
        
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

    stub:
    """
    touch validated.ok
    """
}
