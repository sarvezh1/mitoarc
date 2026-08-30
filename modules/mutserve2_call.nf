process MUTSERVE2_CALL {
    tag "${sample} mutserve2"
    container 'mitoarc-mutserve2:2.0.3'
    publishDir "${params.outdir}/variants/mutserve2/native", mode: 'copy', pattern: '*.mutserve2.*'

    input:
    tuple val(sample), val(representation), path(alignment), path(index)
    path gatk_reference_bundle
    val mt_contig
    val heteroplasmy_threshold

    output:
    tuple val(sample), path("${sample}.mutserve2.native.vcf.gz"), path("${sample}.mutserve2.native.tsv"), path("${sample}.mutserve2.native_raw.tsv"), path("${sample}.mutserve2.metadata.tsv"), emit: native_calls

    script:
    """
    set -euo pipefail
    tar -xf ${gatk_reference_bundle}
    mutserve call \
        --reference standard_mt.fa \
        --output ${sample}.mutserve2.native.vcf.gz \
        --threads ${task.cpus} \
        --level ${heteroplasmy_threshold} \
        --contig-name ${mt_contig} \
        --no-ansi \
        --write-raw \
        ${alignment}

    test -s ${sample}.mutserve2.native.vcf.gz
    test -s ${sample}.txt
    test -s ${sample}_raw.txt
    mv ${sample}.txt ${sample}.mutserve2.native.tsv
    mv ${sample}_raw.txt ${sample}.mutserve2.native_raw.tsv

    printf 'key\tvalue\n' > ${sample}.mutserve2.metadata.tsv
    printf 'sample\t%s\n' '${sample}' >> ${sample}.mutserve2.metadata.tsv
    printf 'mutserve2_version\t2.0.3\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'release_url\thttps://github.com/seppinho/mutserve/releases/tag/v2.0.3\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'release_archive_sha256\t59b85d263906b73e3b8c9bbd21c29557af8e50e5d17d3926850b1877ce19ce2c\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'container\tmitoarc-mutserve2:2.0.3\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'input_representation\tstandard_mt_realign_bam\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'mt_contig\t%s\n' '${mt_contig}' >> ${sample}.mutserve2.metadata.tsv
    printf 'heteroplasmy_threshold\t%s\n' '${heteroplasmy_threshold}' >> ${sample}.mutserve2.metadata.tsv
    printf 'map_quality_default\t20\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'base_quality_default\t20\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'alignment_quality_default\t30\n' >> ${sample}.mutserve2.metadata.tsv
    printf 'baq_default\tfalse\n' >> ${sample}.mutserve2.metadata.tsv
    """
}
