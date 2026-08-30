process REPORTING_FUNCTIONAL_FIXTURE {
    tag "canonical reporting functional fixture"
    container 'mitoarc-reporting:1.0.0'
    publishDir "${params.outdir}/plots/functional", mode: 'copy', pattern: '*.svg'
    publishDir "${params.outdir}/plots/functional", mode: 'copy', pattern: '*.png'
    publishDir "${params.outdir}/report/functional", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/report/functional", mode: 'copy', pattern: '*.reporting_summary.tsv'
    publishDir "${params.outdir}/report/functional", mode: 'copy', pattern: '*.reporting_validation.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.mitoplot_data.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.reporting_features.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.reporting_variants.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.reporting_caller_comparison.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.plot_manifest.tsv'
    publishDir "${params.outdir}/report/functional/data", mode: 'copy', pattern: '*.reporting_summary.json'

    input:
    path fixture_vcf
    path expected
    path canonical_reference
    path canonical_features
    path no_numt_evidence
    path icon
    val coverage_mapq
    val homoplasmy_threshold
    val mutserve2_threshold
    val pipeline_version
    val nextflow_version

    output:
    path "FUNCTIONAL.mtDNA_report.html"
    path "FUNCTIONAL.reporting_summary.tsv"
    path "FUNCTIONAL.reporting_validation.tsv"
    path "FUNCTIONAL.*.svg"
    path "FUNCTIONAL.*.png"
    path "FUNCTIONAL.reporting_summary.json"
    path "FUNCTIONAL.mitoplot_data.tsv"
    path "FUNCTIONAL.reporting_features.tsv"
    path "FUNCTIONAL.reporting_variants.tsv"
    path "FUNCTIONAL.reporting_caller_comparison.tsv"
    path "FUNCTIONAL.plot_manifest.tsv"

    script:
    """
    set -euo pipefail
    export MTDNA_REPORTING_RENDERER_VERSION=1.0.1
    mtdna_annotate \
      --sample FUNCTIONAL --caller functional_fixture \
      --vcf ${fixture_vcf} --reference ${canonical_reference} \
      --canonical-reference ${canonical_reference} --features ${canonical_features} \
      --numt-evidence ${no_numt_evidence} --homoplasmy-threshold '${homoplasmy_threshold}' \
      --output FUNCTIONAL.functional_fixture.annotated.tsv \
      --summary FUNCTIONAL.functional_fixture.annotation_summary.tsv \
      --metadata FUNCTIONAL.functional_fixture.annotation_metadata.tsv

    mtdna_consensus \
      --sample FUNCTIONAL --caller functional_fixture --vcf ${fixture_vcf} \
      --reference ${canonical_reference} --canonical-reference ${canonical_reference} \
      --homoplasmy-threshold '${homoplasmy_threshold}' --heteroplasmy-mode reference \
      --output FUNCTIONAL.functional_fixture.consensus.fa \
      --validation FUNCTIONAL.functional_fixture.consensus_validation.tsv

    make_reporting_functional_fixture \
      --sample FUNCTIONAL --coverage-mapq '${coverage_mapq}' \
      --depth FUNCTIONAL.reporting_depth.tsv --metrics FUNCTIONAL.mtdna_metrics.tsv

    prepare_reporting_data \
      --sample FUNCTIONAL --reference-label NC_012920.1.fa --mt-contig chrM \
      --coverage-mapq '${coverage_mapq}' --homoplasmy-threshold '${homoplasmy_threshold}' \
      --mutserve2-threshold '${mutserve2_threshold}' --numt-mode none \
      --pipeline-version '${pipeline_version}' --nextflow-version '${nextflow_version}' \
      --metrics FUNCTIONAL.mtdna_metrics.tsv --depth FUNCTIONAL.reporting_depth.tsv \
      --annotation-summary FUNCTIONAL.functional_fixture.annotation_summary.tsv \
      --annotation FUNCTIONAL.functional_fixture.annotated.tsv \
      --consensus-validation FUNCTIONAL.functional_fixture.consensus_validation.tsv \
      --features ${canonical_features} \
      --output-summary FUNCTIONAL.reporting_summary.tsv \
      --output-json FUNCTIONAL.reporting_summary.json \
      --output-plot-data FUNCTIONAL.mitoplot_data.tsv \
      --output-features FUNCTIONAL.reporting_features.tsv \
      --output-variants FUNCTIONAL.reporting_variants.tsv \
      --output-comparison FUNCTIONAL.reporting_caller_comparison.tsv

    render_mtdna_plots \
      --sample FUNCTIONAL --summary FUNCTIONAL.reporting_summary.tsv \
      --plot-data FUNCTIONAL.mitoplot_data.tsv --features FUNCTIONAL.reporting_features.tsv \
      --variants FUNCTIONAL.reporting_variants.tsv \
      --comparison FUNCTIONAL.reporting_caller_comparison.tsv --output-prefix FUNCTIONAL

    render_mtdna_report \
      --sample FUNCTIONAL --summary FUNCTIONAL.reporting_summary.tsv \
      --variants FUNCTIONAL.reporting_variants.tsv \
      --comparison FUNCTIONAL.reporting_caller_comparison.tsv \
      --mitoplot FUNCTIONAL.mitoplot.png \
      --depth-plot FUNCTIONAL.mtDNA_depth_profile.png \
      --threshold-plot FUNCTIONAL.coverage_thresholds.png \
      --vaf-plot FUNCTIONAL.vaf_landscape.png \
      --consequence-plot FUNCTIONAL.variant_consequences.png \
      --caller-plot FUNCTIONAL.caller_comparison.png \
      --icon ${icon} \
      --output FUNCTIONAL.mtDNA_report.html

    validate_mtdna_reporting \
      --sample FUNCTIONAL --summary FUNCTIONAL.reporting_summary.tsv \
      --metrics FUNCTIONAL.mtdna_metrics.tsv \
      --annotation-summary FUNCTIONAL.functional_fixture.annotation_summary.tsv \
      --annotation FUNCTIONAL.functional_fixture.annotated.tsv \
      --comparison FUNCTIONAL.reporting_caller_comparison.tsv \
      --variants FUNCTIONAL.reporting_variants.tsv \
      --plot-data FUNCTIONAL.mitoplot_data.tsv \
      --features FUNCTIONAL.reporting_features.tsv \
      --plot-manifest FUNCTIONAL.plot_manifest.tsv \
      --html FUNCTIONAL.mtDNA_report.html \
      --svg FUNCTIONAL.caller_comparison.svg \
      --svg FUNCTIONAL.coverage_thresholds.svg \
      --svg FUNCTIONAL.mitoplot.svg \
      --svg FUNCTIONAL.mtDNA_depth_profile.svg \
      --svg FUNCTIONAL.vaf_landscape.svg \
      --svg FUNCTIONAL.variant_consequences.svg \
      --png FUNCTIONAL.caller_comparison.png \
      --png FUNCTIONAL.coverage_thresholds.png \
      --png FUNCTIONAL.mitoplot.png \
      --png FUNCTIONAL.mtDNA_depth_profile.png \
      --png FUNCTIONAL.vaf_landscape.png \
      --png FUNCTIONAL.variant_consequences.png \
      --expected ${expected} --output FUNCTIONAL.reporting_validation.tsv
    """
}
