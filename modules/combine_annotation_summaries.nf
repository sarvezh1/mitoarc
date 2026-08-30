process COMBINE_ANNOTATION_SUMMARIES {
    tag "${sample} annotation summary"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/annotation", mode: 'copy', pattern: '*.annotation_summary.tsv'

    input:
    tuple val(sample), path(gatk_summary), path(mutserve_summary)

    output:
    path "${sample}.annotation_summary.tsv"

    script:
    """
    set -euo pipefail
    head -n 1 ${gatk_summary} > ${sample}.annotation_summary.tsv
    tail -n +2 ${gatk_summary} >> ${sample}.annotation_summary.tsv
    tail -n +2 ${mutserve_summary} >> ${sample}.annotation_summary.tsv
    """
}

