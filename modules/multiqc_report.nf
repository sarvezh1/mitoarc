process MULTIQC_REPORT {
    tag "technical QC aggregation"
    container 'quay.io/biocontainers/multiqc@sha256:e8821df086710ff4fe535814e70eb6390d4685dd6e962fc743bff58b35b0e603'
    publishDir "${params.outdir}/report", mode: 'copy', overwrite: true, pattern: 'multiqc_report.html'
    publishDir "${params.outdir}/report", mode: 'copy', overwrite: true, pattern: 'multiqc_report_data'

    input:
    path qc_files

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    set -euo pipefail
    mkdir -p qc_input
    for source in ${qc_files}; do
      cp -L "\$source" qc_input/
    done
    multiqc --force --no-ansi --no-version-check --no-ai \
      --cl-config 'show_analysis_paths: false' \
      --template simple \
      --title 'MitoArc technical QC' \
      --filename multiqc_report.html --outdir . --data-dir qc_input
    sanitize_multiqc_data --data multiqc_report_data
    """
}
