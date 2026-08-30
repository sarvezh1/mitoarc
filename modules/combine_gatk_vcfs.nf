process COMBINE_GATK_VCFS {
    tag "${sample} canonical GATK VCF"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/gatk", mode: 'copy', pattern: '*.gatk.raw.vcf.gz*'

    input:
    tuple val(sample), path(standard_vcf), path(standard_index), path(standard_stats), path(lifted_vcf), path(lifted_index), path(lifted_stats)
    path gatk_reference_bundle

    output:
    tuple val(sample), path("${sample}.gatk.raw.vcf.gz"), path("${sample}.gatk.raw.vcf.gz.tbi"), path("${sample}.gatk.raw.vcf.gz.stats"), emit: raw

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    gatk MergeMutectStats \\
        -stats ${standard_stats} \\
        -stats ${lifted_stats} \\
        -output ${sample}.gatk.raw.vcf.gz.stats
    gatk GatherVcfs \\
        -I ${lifted_vcf} \\
        -I ${standard_vcf} \\
        -O ${sample}.gatk.raw.vcf.gz
    gatk IndexFeatureFile -I ${sample}.gatk.raw.vcf.gz
    gatk ValidateVariants -R standard_mt.fa -V ${sample}.gatk.raw.vcf.gz
    """
}
