process GATK_NUMT_AUDIT {
    tag "${sample} GATK 4.6.2.0 NUMT audit"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/numt", mode: 'copy'

    input:
    tuple val(sample), path(raw_vcf), path(raw_index), path(raw_stats), path(filtered_vcf), path(filtered_index), val(std_rep), path(mt_bam), path(mt_index)

    output:
    path "${sample}.gatk_4.6.2.0_numt_audit.tsv"

    script:
    """
    set -euo pipefail
    gatk_numt_audit \
      --sample '${sample}' \
      --raw ${raw_vcf} \
      --filtered ${filtered_vcf} \
      --mt-bam ${mt_bam}
    """
}
