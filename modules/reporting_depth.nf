process REPORTING_DEPTH {
    tag "${sample} reporting coverage"
    container 'quay.io/biocontainers/mosdepth@sha256:52f0d8cb1e61714a8f4e7ccb0f720f3e9dd2928fc4367eb78a42e1cf7c24658d'
    publishDir "${params.outdir}/report/data", mode: 'copy', pattern: '*.reporting_depth.tsv'

    input:
    tuple val(sample), path(alignment), path(index)
    val mt_contig
    val coverage_mapq

    output:
    tuple val(sample), path("${sample}.reporting_depth.tsv"), emit: depth

    script:
    """
    set -euo pipefail
    mosdepth --threads ${task.cpus} --mapq ${coverage_mapq} --chrom '${mt_contig}' ${sample}.reporting ${alignment}
    printf 'sample\tposition\tcoverage\n' > ${sample}.reporting_depth.tsv
    gzip -cd ${sample}.reporting.per-base.bed.gz | \
      awk -v OFS='\t' -v sample='${sample}' '{ for (position = \$2 + 1; position <= \$3; position++) print sample, position, \$4 }' \
      >> ${sample}.reporting_depth.tsv
    """
}

