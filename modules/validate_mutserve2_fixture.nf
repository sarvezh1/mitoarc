process VALIDATE_MUTSERVE2_FIXTURE {
    tag "${sample} mutserve2 fixture truth"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/variants/mutserve2", mode: 'copy', pattern: '*.mutserve2.truth_*.tsv'

    input:
    tuple val(sample), path(vcf), path(index), path(calls)
    path truth

    output:
    path "${sample}.mutserve2.truth_comparison.tsv"
    path "${sample}.mutserve2.truth_variants.tsv"

    script:
    """
    set -euo pipefail
    compare_calls_to_truth \
        '${sample}' \
        ${calls} \
        ${truth} \
        ${sample}.mutserve2.truth_comparison.tsv \
        ${sample}.mutserve2.truth_variants.tsv
    """
}
