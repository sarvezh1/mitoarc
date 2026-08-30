process COMBINE_MTDNA_METRICS {

    tag "${sample}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1'
    publishDir "${params.outdir}/qc/mtdna", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), path(mt_coverage), path(autosomal_coverage)
    val coverage_mapq

    output:
    path "${sample}.mtdna_metrics.tsv"

    script:
    """
    awk -F '\t' -v OFS='\t' -v sample='${sample}' -v coverage_mapq='${coverage_mapq}' '
    FILENAME == ARGV[1] && FNR == 2 {
        mt_contig=\$2; mean_mt=\$3; median_mt=\$4; min_mt=\$5; max_mt=\$6
        pct_1=\$7; pct_10=\$8; pct_100=\$9; pct_500=\$10; pct_1000=\$11
    }
    FILENAME == ARGV[2] && FNR == 2 {
        autosomal_contigs=\$2; mean_auto=\$4; median_auto=\$5; auto_status=\$6
    }
    END {
        print "sample", "mt_contig", "mean_mt_depth", "median_mt_depth", "min_mt_depth", "max_mt_depth", \
              "pct_mt_ge_1x", "pct_mt_ge_10x", "pct_mt_ge_100x", "pct_mt_ge_500x", "pct_mt_ge_1000x", \
              "autosomal_contigs", "mean_autosomal_depth", "median_autosomal_depth", "mtdna_copy_number_proxy", "coverage_mapq", "copy_number_proxy_status"

        proxy="NA"
        status="unavailable_autosomal_depth"
        if (mean_mt != "NA" && mean_auto != "NA" && (mean_auto + 0) > 0) {
            proxy=sprintf("%.6f", 2 * mean_mt / mean_auto)
            status="available_coverage_derived_proxy"
        } else if (mean_auto != "NA" && (mean_auto + 0) == 0) {
            status="unavailable_zero_autosomal_depth"
        } else if (auto_status != "") {
            status="unavailable_" auto_status
        }

        print sample, mt_contig, mean_mt, median_mt, min_mt, max_mt, pct_1, pct_10, pct_100, pct_500, pct_1000, \
              autosomal_contigs, mean_auto, median_auto, proxy, coverage_mapq, status
    }' ${mt_coverage} ${autosomal_coverage} > ${sample}.mtdna_metrics.tsv
    """
}
