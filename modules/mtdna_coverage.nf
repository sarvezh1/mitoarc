process MTDNA_COVERAGE {

    tag "${sample}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1'
    publishDir "${params.outdir}/qc/mtdna", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), path(alignment), path(index)
    val mt_contig
    val coverage_mapq

    output:
    tuple val(sample), path("${sample}.mtDNA_coverage.tsv")

    script:
    """
    mosdepth --threads ${task.cpus} --mapq ${coverage_mapq} --chrom '${mt_contig}' ${sample}.mtDNA ${alignment}
    gzip -cd ${sample}.mtDNA.per-base.bed.gz > per_base.tsv

    awk '
    {
        span = \$3 - \$2
        depth = \$4 + 0
        total += span
        weighted_sum += span * depth
        histogram[depth] += span
        if (minimum == "" || depth < minimum) minimum = depth
        if (maximum == "" || depth > maximum) maximum = depth
        if (depth >= 1) covered_1 += span
        if (depth >= 10) covered_10 += span
        if (depth >= 100) covered_100 += span
        if (depth >= 500) covered_500 += span
        if (depth >= 1000) covered_1000 += span
    }
    END {
        for (depth in histogram) print depth, histogram[depth] > "depth_histogram.unsorted.tsv"
        print total, weighted_sum, minimum, maximum, covered_1, covered_10, covered_100, covered_500, covered_1000 > "depth_totals.tsv"
    }' per_base.tsv

    sort -n -k1,1 depth_histogram.unsorted.tsv > depth_histogram.tsv

    awk -v sample='${sample}' -v mt_contig='${mt_contig}' '
    NR == FNR {
        total=\$1; weighted_sum=\$2; minimum=\$3; maximum=\$4
        covered_1=\$5; covered_10=\$6; covered_100=\$7; covered_500=\$8; covered_1000=\$9
        lower_rank=int((total + 1) / 2)
        upper_rank=int((total + 2) / 2)
        next
    }
    {
        cumulative += \$2
        if (lower_depth == "" && cumulative >= lower_rank) lower_depth=\$1
        if (upper_depth == "" && cumulative >= upper_rank) upper_depth=\$1
    }
    END {
        print "sample\tmt_contig\tmean_mt_depth\tmedian_mt_depth\tmin_mt_depth\tmax_mt_depth\tpct_mt_ge_1x\tpct_mt_ge_10x\tpct_mt_ge_100x\tpct_mt_ge_500x\tpct_mt_ge_1000x"
        if (total > 0) {
            printf "%s\t%s\t%.6f\t%.6f\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\\n", \
                sample, mt_contig, weighted_sum / total, (lower_depth + upper_depth) / 2, minimum, maximum, \
                100 * covered_1 / total, 100 * covered_10 / total, 100 * covered_100 / total, \
                100 * covered_500 / total, 100 * covered_1000 / total
        } else {
            printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\\n", sample, mt_contig
        }
    }' depth_totals.tsv depth_histogram.tsv > ${sample}.mtDNA_coverage.tsv
    """
}
