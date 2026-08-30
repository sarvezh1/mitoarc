process GATK_MUTECT2 {
    tag "${sample} ${representation} Mutect2"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/gatk/intermediate", mode: 'copy', pattern: '*.vcf.gz*'

    input:
    tuple val(sample), val(input_representation), path(alignment), path(index)
    path gatk_reference_bundle
    path metrics
    val representation

    output:
    tuple val(sample), val(representation), path("${sample}.${representation}.gatk.raw.vcf.gz"), path("${sample}.${representation}.gatk.raw.vcf.gz.tbi"), path("${sample}.${representation}.gatk.raw.vcf.gz.stats"), emit: raw

    script:
    def ref = representation == 'standard' ? 'standard_mt' : 'shifted_mt'
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    interval=\$(cat ${ref}.interval)
    gatk Mutect2 \\
        -R ${ref}.fa \\
        -I ${alignment} \\
        -L "\${interval}" \\
        --mitochondria-mode \\
        -O ${sample}.${representation}.gatk.raw.vcf.gz
    test -s ${sample}.${representation}.gatk.raw.vcf.gz
    test -s ${sample}.${representation}.gatk.raw.vcf.gz.tbi
    """
}
