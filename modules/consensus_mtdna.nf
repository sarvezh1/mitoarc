process CONSENSUS_MTDNA {
    tag "${sample} ${caller} mitochondrial consensus"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/consensus", mode: 'copy', pattern: '*.consensus.fa'
    publishDir "${params.outdir}/interpretation/consensus", mode: 'copy', pattern: '*.consensus_validation.tsv'

    input:
    tuple val(sample), val(caller), path(vcf)
    path gatk_reference_bundle
    path canonical_reference
    val homoplasmy_threshold
    val heteroplasmy_mode

    output:
    tuple val(sample), val(caller), path("${sample}.${caller}.consensus.fa"), path("${sample}.${caller}.consensus_validation.tsv"), emit: consensus

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    mtdna_consensus \
        --sample '${sample}' \
        --caller '${caller}' \
        --vcf ${vcf} \
        --reference standard_mt.fa \
        --canonical-reference ${canonical_reference} \
        --homoplasmy-threshold '${homoplasmy_threshold}' \
        --heteroplasmy-mode '${heteroplasmy_mode}' \
        --output ${sample}.${caller}.consensus.fa \
        --validation ${sample}.${caller}.consensus_validation.tsv
    """
}

