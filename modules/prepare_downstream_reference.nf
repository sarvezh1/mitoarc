process PREPARE_DOWNSTREAM_REFERENCE {

    tag "downstream reference index"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.fai"), emit: indexed_reference

    script:
    """
    set -euo pipefail
    samtools faidx ${reference}
    """
}
