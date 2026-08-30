process VALIDATE_REPORTING {
    tag "${sample} reporting validation"
    container 'mitoarc-reporting:1.0.0'
    publishDir "${params.outdir}/report/validation", mode: 'copy'

    input:
    tuple val(sample), path(summary), path(plot_data), path(features), path(variants), path(comparison),
          path(metrics), path(annotation_summary), path(gatk_annotation), path(mutserve2_annotation),
          path(source_comparison), path(svgs), path(pngs), path(plot_manifest), path(html)

    output:
    path "${sample}.reporting_validation.tsv"

    script:
    """
    set -euo pipefail
    validate_mtdna_reporting \
      --sample '${sample}' \
      --summary ${summary} \
      --metrics ${metrics} \
      --annotation-summary ${annotation_summary} \
      --annotation ${gatk_annotation} \
      --annotation ${mutserve2_annotation} \
      --source-comparison ${source_comparison} \
      --comparison ${comparison} \
      --variants ${variants} \
      --plot-data ${plot_data} \
      --features ${features} \
      --plot-manifest ${plot_manifest} \
      --html ${html} \
      --svg ${sample}.caller_comparison.svg \
      --svg ${sample}.coverage_thresholds.svg \
      --svg ${sample}.mitoplot.svg \
      --svg ${sample}.mtDNA_depth_profile.svg \
      --svg ${sample}.vaf_landscape.svg \
      --svg ${sample}.variant_consequences.svg \
      --png ${sample}.caller_comparison.png \
      --png ${sample}.coverage_thresholds.png \
      --png ${sample}.mitoplot.png \
      --png ${sample}.mtDNA_depth_profile.png \
      --png ${sample}.vaf_landscape.png \
      --png ${sample}.variant_consequences.png \
      --output ${sample}.reporting_validation.tsv
    """
}

