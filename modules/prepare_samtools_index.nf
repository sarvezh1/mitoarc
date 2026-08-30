process PREPARE_SAMTOOLS_INDEX {

    tag "samtools reference index"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    path reference_bundle

    output:
    path 'samtools_reference_bundle.tar'

    script:
    """
    tar -xf ${reference_bundle}
    if [ ! -f reference.fa.fai ]; then
        samtools faidx reference.fa
    fi
    tar -cf samtools_reference_bundle.tar ./*.fa ./*.fa.fai ./reference.*
    """
}
