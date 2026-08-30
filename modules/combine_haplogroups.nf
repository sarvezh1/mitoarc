process COMBINE_HAPLOGROUPS {
    tag "${sample} haplogroup comparison"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/haplogroup", mode: 'copy', pattern: '*.haplogroup*.tsv'

    input:
    tuple val(sample), path(gatk_standard), path(mutserve_standard)
    val tree

    output:
    path "${sample}.haplogroups.tsv"
    path "${sample}.haplogroup_comparison.tsv"

    script:
    def treeVersion = tree.toString().tokenize('@')[-1]
    """
    set -euo pipefail
    compare_haplogroups '${sample}' ${gatk_standard} ${mutserve_standard} \
        /opt/haplogrep3/trees/phylotree-rcrs/${treeVersion}/tree.xml \
        ${sample}.haplogroups.tsv ${sample}.haplogroup_comparison.tsv '${tree}'
    """
}
