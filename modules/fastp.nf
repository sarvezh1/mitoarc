process FASTP {

    tag "${sample}"
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_3'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.html'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/reads", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample), path(reads1), path(reads2), path(fastqc_html), path(fastqc_zip)

    output:
    tuple val(sample), path("${sample}.cleaned_R1.fastq.gz"), path("${sample}.cleaned_R2.fastq.gz"), emit: cleaned_reads
    path "${sample}.fastp.html"
    path "${sample}.fastp.json"

    script:
    """
    fastp \\
        --in1 ${reads1} \\
        --in2 ${reads2} \\
        --out1 ${sample}.cleaned_R1.fastq.gz \\
        --out2 ${sample}.cleaned_R2.fastq.gz \\
        --html ${sample}.fastp.html \\
        --json ${sample}.fastp.json
    """
}
