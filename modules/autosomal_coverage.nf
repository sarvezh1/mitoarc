process AUTOSOMAL_COVERAGE {

    tag "${sample}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h05c3d44_1'
    publishDir "${params.outdir}/qc/mtdna", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(sample), path(alignment), path(index), path(autosomal_regions), path(autosomal_contigs), path(autosomal_mapping_context)
    val coverage_mapq

    output:
    tuple val(sample), path("${sample}.autosomal_coverage.tsv")

    script:
    """
    printf 'sample\tautosomal_contigs\ttotal_autosomal_bases\tmean_autosomal_depth\tmedian_autosomal_depth\tstatus\n' > ${sample}.autosomal_coverage.tsv

    if [ ! -s ${autosomal_regions} ]; then
        printf '${sample}\tNA\t0\tNA\tNA\tno_autosomal_contigs_detected\n' >> ${sample}.autosomal_coverage.tsv
    else
        selected_count=\$(awk -F '\t' 'NR == 2 { print \$1 }' ${autosomal_mapping_context})
        mapped_selected_count=\$(awk -F '\t' 'NR == 2 { print \$2 }' ${autosomal_mapping_context})
        paste -sd, ${autosomal_contigs} > selected_contigs.csv

        # A retained full header can mechanically advertise all autosomes even
        # when a reduced BAM contains mapped data on only a subset. Whole-genome
        # mosdepth is then both memory-wasteful and biologically uninterpretable.
        if [ "\${selected_count}" -ge 22 ] && [ "\${mapped_selected_count}" -lt "\${selected_count}" ]; then
            printf '${sample}\t%s\t0\tNA\tNA\tincomplete_mapped_autosomal_context\n' "\$(cat selected_contigs.csv)" >> ${sample}.autosomal_coverage.tsv
            exit 0
        fi

        mosdepth --threads ${task.cpus} --mapq ${coverage_mapq} --no-per-base --by ${autosomal_regions} ${sample}.autosomal ${alignment}
        gzip -cd ${sample}.autosomal.regions.bed.gz > autosomal_depths.tsv

        sort -k5,5n -k2,2n autosomal_depths.tsv > autosomal_depths.sorted.tsv
        total_autosomal_bases=\$(awk '{s += \$3 - \$2} END {print s + 0}' autosomal_depths.tsv)
        lower_rank=\$(( (total_autosomal_bases + 1) / 2 ))
        upper_rank=\$(( (total_autosomal_bases + 2) / 2 ))
        awk -v sample='${sample}' -v contigs="\$(cat selected_contigs.csv)" -v lower_rank="\${lower_rank}" -v upper_rank="\${upper_rank}" '
        {
            span = \$3 - \$2
            total_bases += span
            weighted_sum += span * \$5
            depth = \$5 + 0
            cumulative += span
            if (lower_depth == "" && cumulative >= lower_rank) lower_depth = depth
            if (upper_depth == "" && cumulative >= upper_rank) upper_depth = depth
        }
        END {
            if (total_bases > 0) {
                printf "%s\t%s\t%d\t%.6f\t%.6f\tavailable\\n", sample, contigs, total_bases, weighted_sum / total_bases, (lower_depth + upper_depth) / 2
            } else {
                printf "%s\t%s\t0\tNA\tNA\tno_autosomal_bases\\n", sample, contigs
            }
        }' autosomal_depths.sorted.tsv >> ${sample}.autosomal_coverage.tsv
    fi
    """
}
