process LIFTOVER_MT_VCF {
    tag "${sample} shifted VCF liftover"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/gatk/intermediate", mode: 'copy', pattern: '*.lifted.vcf.gz*'

    input:
    tuple val(sample), val(representation), path(shifted_vcf), path(shifted_index)
    path gatk_reference_bundle

    output:
    tuple val(sample), path("${sample}.shifted.lifted.vcf.gz"), path("${sample}.shifted.lifted.vcf.gz.tbi"), emit: lifted

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    gatk LiftoverVcf \\
        -I ${shifted_vcf} \\
        -O ${sample}.shifted.lifted.vcf.gz \\
        -CHAIN shifted_to_standard.chain \\
        -REJECT rejected.vcf \\
        -R standard_mt.fa
    if grep -qv '^#' rejected.vcf; then
        echo 'Liftover rejected shifted variants:' >&2
        cat rejected.vcf >&2
        exit 1
    fi
    gatk IndexFeatureFile -I ${sample}.shifted.lifted.vcf.gz
    gatk ValidateVariants -R standard_mt.fa -V ${sample}.shifted.lifted.vcf.gz
    """
}
