process GENERATE_REPORT_PLOTS {
    tag "${sample} reporting plots"
    container 'mitoarc-reporting:1.0.0'
    publishDir "${params.outdir}/plots", mode: 'copy', pattern: '*.svg'
    publishDir "${params.outdir}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.plot_manifest.tsv'

    input:
    tuple val(sample), path(summary), path(summary_json), path(plot_data), path(features),
          path(variants), path(comparison)

    output:
    tuple val(sample), path("${sample}.*.svg"), path("${sample}.*.png"),
          path("${sample}.plot_manifest.tsv"), emit: plots

    script:
    """
    set -euo pipefail
    export MTDNA_REPORTING_RENDERER_VERSION=1.0.1
    render_mtdna_plots \
      --sample '${sample}' \
      --summary ${summary} \
      --plot-data ${plot_data} \
      --features ${features} \
      --variants ${variants} \
      --comparison ${comparison} \
      --output-prefix ${sample}
    """
}
