process COMBINE_CONTAMINATION {
    tag "${sample} contamination summary"
    container 'quay.io/biocontainers/haplocheck@sha256:1d4e006a69dab62cc0304509be83cc8c11f932cbf6a88fe4aacbb072e8eebacd'
    publishDir "${params.outdir}/interpretation/contamination", mode: 'copy', pattern: '*.contamination.tsv'

    input:
    tuple val(sample), path(gatk_assessment), path(mutserve_assessment)

    output:
    path "${sample}.contamination.tsv"

    script:
    """
    set -euo pipefail
    head -n 1 ${gatk_assessment} > ${sample}.contamination.tsv
    tail -n +2 ${gatk_assessment} >> ${sample}.contamination.tsv
    tail -n +2 ${mutserve_assessment} >> ${sample}.contamination.tsv
    """
}

