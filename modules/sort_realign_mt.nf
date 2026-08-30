process SORT_REALIGN_MT {
    tag "${sample} ${representation}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/mtDNA/realigned", mode: 'copy', pattern: '*.bam*', saveAs: { filename ->
        filename.contains('.standard_') ? "standard/${filename}" : "shifted/${filename}"
    }
    input:
    tuple val(sample), val(representation), path(alignment)
    output:
    tuple val(sample), val(representation), path("${sample}.${representation}_mt.sorted.bam"), path("${sample}.${representation}_mt.sorted.bam.bai")
    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -o ${sample}.${representation}_mt.sorted.bam ${alignment}
    samtools index -@ ${task.cpus} ${sample}.${representation}_mt.sorted.bam ${sample}.${representation}_mt.sorted.bam.bai
    """
}
