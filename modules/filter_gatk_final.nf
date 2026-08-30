process FILTER_GATK_FINAL {
    tag "${sample} canonical FilterMutectCalls"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/gatk", mode: 'copy', pattern: '*.gatk.filtered.vcf.gz*'

    input:
    tuple val(sample), path(raw_vcf), path(raw_index), path(raw_stats)
    path gatk_reference_bundle
    path metrics

    output:
    tuple val(sample), path("${sample}.gatk.filtered.vcf.gz"), path("${sample}.gatk.filtered.vcf.gz.tbi"), emit: filtered

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    gatk FilterMutectCalls \\
        -R standard_mt.fa \\
        -V ${raw_vcf} \\
        --stats ${raw_stats} \\
        --mitochondria-mode \\
        -O ${sample}.gatk.filtered.vcf.gz
    """
}
