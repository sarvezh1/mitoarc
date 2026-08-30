process NUMT_FIXTURE_PREFLIGHT {
    tag "${sample} NUMT fixture pre-flight"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants/numt", mode: 'copy'

    input:
    tuple val(sample), path(fullref_bam), path(fullref_index), path(extracted_bam), path(extracted_index), val(std_rep), path(standard_bam), path(standard_index), val(shift_rep), path(shifted_bam), path(shifted_index)
    path design
    val mt_contig

    output:
    path "${sample}.numt_fixture_validation.tsv"

    script:
    """
    set -euo pipefail
    numt_fixture_preflight \
      --sample '${sample}' \
      --design ${design} \
      --fullref ${fullref_bam} \
      --extracted ${extracted_bam} \
      --standard ${standard_bam} \
      --shifted ${shifted_bam} \
      --mt-contig '${mt_contig}'
    """
}
