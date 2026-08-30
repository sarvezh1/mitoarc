process FASTQC {

    tag "${sample}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', pattern: '*.zip'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path(reads1), path(reads2), path('*.html'), path('*.zip')

    script:
    """
    fastqc --outdir . ${reads1} ${reads2}
    """
}
