process PREPARE_REFERENCE {

    tag "reference index"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    path reference

    output:
    path 'reference_bundle.tar'

    script:
    """
    mkdir -p reference
    cp ${reference} reference/reference.fa
    if [ ! -f reference/reference.bwt.2bit.64 ] || [ ! -f reference/reference.0123 ] || [ ! -f reference/reference.ann ] || [ ! -f reference/reference.amb ] || [ ! -f reference/reference.pac ]; then
        bwa-mem2 index -p reference/reference reference/reference.fa
    fi
    tar -cf reference_bundle.tar -C reference .
    """
}
