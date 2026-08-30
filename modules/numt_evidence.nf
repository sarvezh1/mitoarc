process NUMT_EVIDENCE {
    tag "${sample} NUMT evidence"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/numt", mode: 'copy', pattern: '*.numt_*.tsv'

    input:
    tuple val(sample), path(fullref_bam), path(fullref_index), val(std_rep), path(mt_bam), path(mt_index), path(gatk_a), path(gatk_a_index), path(mutserve_a), path(mutserve_a_index), path(mutserve_a_calls), path(gatk_b), path(gatk_b_index), path(mutserve_b), path(mutserve_b_index), path(mutserve_b_calls)
    path truth
    val mt_contig
    val min_mapq

    output:
    path "${sample}.numt_read_evidence.tsv"
    path "${sample}.numt_variant_evidence.tsv"
    path "${sample}.numt_strategy_comparison.tsv"
    path "${sample}.numt_extended_truth_validation.tsv"

    script:
    """
    set -euo pipefail
    numt_evidence \
      --sample '${sample}' \
      --fullref-bam ${fullref_bam} \
      --mt-bam ${mt_bam} \
      --gatk-a ${gatk_a} \
      --mutserve-a ${mutserve_a} \
      --gatk-b ${gatk_b} \
      --mutserve-b ${mutserve_b} \
      --truth ${truth} \
      --mt-contig '${mt_contig}' \
      --min-mapq '${min_mapq}'
    """
}
