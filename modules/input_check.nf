process INPUT_CHECK {

    tag "${sample}"

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path(reads1), path(reads2)

    script:
    """
    echo "Validated sample: ${sample}"
    """
}
