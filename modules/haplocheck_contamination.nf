process HAPLOCHECK_CONTAMINATION {
    tag "${sample} ${caller} Haplocheck"
    container 'quay.io/biocontainers/haplocheck@sha256:1d4e006a69dab62cc0304509be83cc8c11f932cbf6a88fe4aacbb072e8eebacd'
    publishDir "${params.outdir}/interpretation/contamination/native", mode: 'copy', pattern: '*.haplocheck.native.tsv'
    publishDir "${params.outdir}/interpretation/contamination/native", mode: 'copy', pattern: '*.haplocheck.native.raw.txt'

    input:
    tuple val(sample), val(caller), path(vcf)

    output:
    tuple val(sample), val(caller), path("${sample}.${caller}.contamination.tsv"), emit: assessment
    path "${sample}.${caller}.haplocheck.native.tsv"
    path "${sample}.${caller}.haplocheck.native.raw.txt"

    script:
    """
    set -euo pipefail
    haplocheck_call '${sample}' '${caller}' ${vcf} \
        ${sample}.${caller}.haplocheck.native.tsv \
        ${sample}.${caller}.contamination.tsv
    """
}

