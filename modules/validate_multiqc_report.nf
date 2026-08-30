process VALIDATE_MULTIQC_REPORT {
    tag "MultiQC validation"
    container 'quay.io/biocontainers/multiqc@sha256:e8821df086710ff4fe535814e70eb6390d4685dd6e962fc743bff58b35b0e603'
    publishDir "${params.outdir}/report/validation", mode: 'copy', overwrite: true

    input:
    path report
    path data

    output:
    path "multiqc_reporting_validation.tsv"

    script:
    """
    set -euo pipefail
    validate_multiqc_report --report ${report} --data ${data} --output multiqc_reporting_validation.tsv
    """
}
