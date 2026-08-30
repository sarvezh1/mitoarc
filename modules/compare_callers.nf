process COMPARE_CALLERS {
    tag "${sample} GATK versus mutserve2"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/variants/comparison", mode: 'copy', pattern: '*.caller_*.tsv'

    input:
    tuple val(sample), path(gatk_raw_vcf), path(gatk_raw_index), path(gatk_raw_stats), path(gatk_filtered_vcf), path(gatk_filtered_index), path(mutserve2_vcf), path(mutserve2_index), path(mutserve2_calls)
    path truth

    output:
    path "${sample}.caller_comparison.tsv"
    path "${sample}.caller_summary.tsv"

    script:
    """
    set -euo pipefail
    compare_mtdna_callers \
        '${sample}' \
        ${gatk_raw_vcf} \
        ${gatk_filtered_vcf} \
        ${mutserve2_vcf} \
        ${truth} \
        ${sample}.caller_comparison.tsv \
        ${sample}.caller_summary.tsv
    """
}
