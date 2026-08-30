process NORMALIZE_MUTSERVE2 {
    tag "${sample} canonical mutserve2 VCF"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/variants/mutserve2", mode: 'copy', pattern: '*.mutserve2.raw.vcf.gz*'
    publishDir "${params.outdir}/variants/mutserve2", mode: 'copy', pattern: '*.mutserve2.calls.tsv'

    input:
    tuple val(sample), path(native_vcf), path(native_tsv), path(native_raw_tsv), path(metadata)
    path gatk_reference_bundle
    val mt_contig

    output:
    tuple val(sample), path("${sample}.mutserve2.raw.vcf.gz"), path("${sample}.mutserve2.raw.vcf.gz.tbi"), path("${sample}.mutserve2.calls.tsv"), emit: canonical

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    printf '%s\n' '${sample}' > sample_name.txt
    bcftools reheader --samples sample_name.txt --output reheader.vcf.gz ${native_vcf}
    bcftools norm \
        --fasta-ref standard_mt.fa \
        --check-ref e \
        --multiallelics -any \
        reheader.vcf.gz \
      | bcftools sort --output-type z --output-file ${sample}.mutserve2.raw.vcf.gz
    bcftools index --tbi ${sample}.mutserve2.raw.vcf.gz

    test "\$(bcftools query --list-samples ${sample}.mutserve2.raw.vcf.gz)" = '${sample}'
    test "\$(bcftools query --format '%CHROM\n' ${sample}.mutserve2.raw.vcf.gz | sort -u)" = '${mt_contig}'
    mutserve2_vcf_to_tsv '${sample}' ${sample}.mutserve2.raw.vcf.gz ${sample}.mutserve2.calls.tsv
    """
}
