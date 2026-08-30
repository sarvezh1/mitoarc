process FLAGSTAT {

    tag "${sample}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/qc/alignment", mode: 'copy'

    input:
    tuple val(sample), path(alignment), path(index)
    val suffix

    output:
    path "${sample}${suffix}.flagstat.txt"

    script:
    """
    samtools flagstat ${alignment} > ${sample}${suffix}.flagstat.txt
    """
}
