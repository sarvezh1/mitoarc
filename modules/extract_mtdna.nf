process EXTRACT_MTDNA {

    tag "${sample}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/mtDNA", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample), path(alignment), path(index)
    val mt_contig

    output:
    tuple val(sample), path("${sample}.mtDNA.bam"), path("${sample}.mtDNA.bam.bai"), emit: mt_bams

    script:
    """
    if ! samtools idxstats ${alignment} | awk -v contig='${mt_contig}' '\$1 == contig { found=1 } END { exit !found }'; then
        echo "Requested mitochondrial contig '${mt_contig}' is not present in the indexed alignment." >&2
        exit 1
    fi
    samtools view ${alignment} '${mt_contig}' | cut -f1 | sort -u > mt_read_names.txt
    samtools view -bh -N mt_read_names.txt ${alignment} > ${sample}.mtDNA.bam
    samtools index -@ ${task.cpus} ${sample}.mtDNA.bam ${sample}.mtDNA.bam.bai
    """
}
