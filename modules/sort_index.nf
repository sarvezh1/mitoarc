process SORT_INDEX {

    tag "${sample}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample), path(alignment)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: sorted_bams

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${alignment}
    samtools index -@ ${task.cpus} ${sample}.sorted.bam ${sample}.sorted.bam.bai
    """
}
