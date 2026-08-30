process ANNOTATE_MTDNA {
    tag "${sample} ${caller} mitochondrial annotation"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/annotation", mode: 'copy', pattern: '*.annotated.tsv'
    publishDir "${params.outdir}/interpretation/annotation/native", mode: 'copy', pattern: '*.annotation_metadata.tsv'

    input:
    tuple val(sample), val(caller), path(vcf), path(numt_evidence)
    path gatk_reference_bundle
    path canonical_reference
    path features
    val homoplasmy_threshold

    output:
    tuple val(sample), val(caller), path("${sample}.${caller}.annotated.tsv"), path("${sample}.${caller}.annotation_summary.tsv"), emit: annotation
    path "${sample}.${caller}.annotation_metadata.tsv"

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    mtdna_annotate \
        --sample '${sample}' \
        --caller '${caller}' \
        --vcf ${vcf} \
        --reference standard_mt.fa \
        --canonical-reference ${canonical_reference} \
        --features ${features} \
        --numt-evidence ${numt_evidence} \
        --homoplasmy-threshold '${homoplasmy_threshold}' \
        --output ${sample}.${caller}.annotated.tsv \
        --summary ${sample}.${caller}.annotation_summary.tsv \
        --metadata ${sample}.${caller}.annotation_metadata.tsv
    """
}

