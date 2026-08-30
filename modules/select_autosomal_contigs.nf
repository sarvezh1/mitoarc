process SELECT_AUTOSOMAL_CONTIGS {

    tag "${sample}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(sample), path(alignment), path(index)
    val autosomal_contigs

    output:
    tuple val(sample), path(alignment), path(index), path('autosomal_regions.bed'), path('autosomal_contigs.txt'), path('autosomal_mapping_context.tsv')

    script:
    """
    samtools idxstats ${alignment} > idxstats.tsv

    if [ -n '${autosomal_contigs}' ]; then
        printf '%s\n' '${autosomal_contigs}' | tr ',' '\n' > requested_contigs.txt
        awk 'NR == FNR { present[\$1]=1; next } ! (\$1 in present) { print \$1 }' idxstats.tsv requested_contigs.txt > missing_contigs.txt
        if [ -s missing_contigs.txt ]; then
            echo "Requested autosomal contig(s) not present in alignment: \$(paste -sd, missing_contigs.txt)" >&2
            exit 1
        fi
        cp requested_contigs.txt autosomal_contigs.txt
    else
        awk '\$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ && \$2 > 0 { print \$1 }' idxstats.tsv > autosomal_contigs.txt
    fi

    awk 'NR == FNR { selected[\$1]=1; next } (\$1 in selected) && \$2 > 0 { print \$1 "\t0\t" \$2 "\t" \$1 }' \
        autosomal_contigs.txt idxstats.tsv > autosomal_regions.bed

    selected_count=\$(awk 'NF { n++ } END { print n+0 }' autosomal_contigs.txt)
    mapped_selected_count=\$(awk 'NR == FNR { selected[\$1]=1; next } (\$1 in selected) && \$3 > 0 { n++ } END { print n+0 }' autosomal_contigs.txt idxstats.tsv)
    printf 'selected_contig_count\tmapped_selected_contig_count\n%s\t%s\n' "\${selected_count}" "\${mapped_selected_count}" > autosomal_mapping_context.tsv
    """
}
