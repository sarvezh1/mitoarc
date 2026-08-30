process BWA_MEM2 {

    tag "${sample}"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple val(sample), path(reads1), path(reads2), path(reference_bundle)

    output:
    tuple val(sample), path("${sample}.sam"), emit: alignments

    script:
    """
    tar -xf ${reference_bundle}
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \\
        reference \\
        ${reads1} ${reads2} > ${sample}.sam
    """
}
