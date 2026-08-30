process PREPARE_REPORTING_DATA {
    tag "${sample} reporting data contract"
    container 'mitoarc-reporting:1.0.0'
    publishDir "${params.outdir}/report", mode: 'copy', pattern: '*.reporting_summary.tsv'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.reporting_summary.json'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.mitoplot_data.tsv'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.reporting_features.tsv'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.reporting_variants.tsv'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.reporting_caller_comparison.tsv'

    input:
    tuple val(sample), path(depth), path(metrics), path(source_comparison), path(caller_summary),
          path(haplogroups), path(haplogroup_comparison), path(contamination), path(annotation_summary),
          path(gatk_annotation), path(mutserve2_annotation), path(gatk_consensus_validation),
          path(mutserve2_consensus_validation), path(input_validation)
    path canonical_features
    val reference_label
    val mt_contig
    val coverage_mapq
    val homoplasmy_threshold
    val mutserve2_threshold
    val numt_mode
    val pipeline_version
    val nextflow_version

    output:
    tuple val(sample), path("${sample}.reporting_summary.tsv"), path("${sample}.reporting_summary.json"),
          path("${sample}.mitoplot_data.tsv"), path("${sample}.reporting_features.tsv"),
          path("${sample}.reporting_variants.tsv"), path("${sample}.reporting_caller_comparison.tsv"),
          emit: contract

    script:
    """
    set -euo pipefail
    prepare_reporting_data \
      --sample '${sample}' \
      --reference-label '${reference_label}' \
      --mt-contig '${mt_contig}' \
      --coverage-mapq '${coverage_mapq}' \
      --homoplasmy-threshold '${homoplasmy_threshold}' \
      --mutserve2-threshold '${mutserve2_threshold}' \
      --numt-mode '${numt_mode}' \
      --pipeline-version '${pipeline_version}' \
      --nextflow-version '${nextflow_version}' \
      --input-validation ${input_validation} \
      --metrics ${metrics} \
      --depth ${depth} \
      --caller-comparison ${source_comparison} \
      --caller-summary ${caller_summary} \
      --haplogroups ${haplogroups} \
      --haplogroup-comparison ${haplogroup_comparison} \
      --contamination ${contamination} \
      --annotation-summary ${annotation_summary} \
      --annotation ${gatk_annotation} \
      --annotation ${mutserve2_annotation} \
      --consensus-validation ${gatk_consensus_validation} \
      --consensus-validation ${mutserve2_consensus_validation} \
      --features ${canonical_features} \
      --output-summary ${sample}.reporting_summary.tsv \
      --output-json ${sample}.reporting_summary.json \
      --output-plot-data ${sample}.mitoplot_data.tsv \
      --output-features ${sample}.reporting_features.tsv \
      --output-variants ${sample}.reporting_variants.tsv \
      --output-comparison ${sample}.reporting_caller_comparison.tsv
    """
}
