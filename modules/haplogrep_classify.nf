process HAPLOGREP_CLASSIFY {
    tag "${sample} ${caller} HaploGrep"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/haplogroup/native", mode: 'copy', pattern: '*.haplogrep.native.tsv'

    input:
    tuple val(sample), val(caller), path(vcf)
    val tree
    val het_level

    output:
    tuple val(sample), val(caller), path("${sample}.${caller}.haplogroup.tsv"), path("${sample}.${caller}.haplogrep.native.tsv"), emit: assignment

    script:
    """
    set -euo pipefail
    haplogrep_call '${sample}' '${caller}' ${vcf} '${tree}' '${het_level}' \
        ${sample}.${caller}.haplogrep.native.tsv \
        ${sample}.${caller}.haplogroup.tsv
    """
}

