process REALIGN_MT_FLAGSTAT {
    tag "${sample} ${representation}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/qc/mtdna_realign", mode: 'copy'
    input:
    tuple val(sample), val(representation), path(alignment), path(index)
    output:
    path "${sample}.${representation}_mt.flagstat.txt"
    script:
    """
    samtools flagstat ${alignment} > ${sample}.${representation}_mt.flagstat.txt
    """
}
