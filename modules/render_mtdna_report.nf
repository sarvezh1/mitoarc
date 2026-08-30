process RENDER_MTDNA_REPORT {
    tag "${sample} integrated report"
    container 'mitoarc-reporting:1.0.0'
    publishDir "${params.outdir}/report", mode: 'copy', pattern: '*.mtDNA_report.html'

    input:
    tuple val(sample), path(summary), path(variants), path(comparison), path(svgs), path(pngs)
    path icon

    output:
    tuple val(sample), path("${sample}.mtDNA_report.html"), emit: report

    script:
    """
    set -euo pipefail
    render_mtdna_report \
      --sample '${sample}' \
      --summary ${summary} \
      --variants ${variants} \
      --comparison ${comparison} \
      --mitoplot ${sample}.mitoplot.png \
      --depth-plot ${sample}.mtDNA_depth_profile.png \
      --threshold-plot ${sample}.coverage_thresholds.png \
      --vaf-plot ${sample}.vaf_landscape.png \
      --consequence-plot ${sample}.variant_consequences.png \
      --caller-plot ${sample}.caller_comparison.png \
      --icon ${icon} \
      --output ${sample}.mtDNA_report.html
    """
}
