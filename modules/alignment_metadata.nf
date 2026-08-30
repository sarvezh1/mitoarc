process ALIGNMENT_METADATA {
    tag "${sample} FASTQ alignment metadata"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/qc/input", mode: 'copy', pattern: '*.input_alignment_validation.tsv'

    input:
    tuple val(sample), path(alignment), path(index)
    val reference
    val mt_contig

    output:
    tuple val(sample), path(alignment), path(index), path("${sample}.input_alignment_validation.tsv"), emit: alignment_records

    script:
    """
    set -euo pipefail
    samtools idxstats ${alignment} > idxstats.tsv
    autosomes=\$(awk '\$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ && \$2 > 0 { n++ } END { print n+0 }' idxstats.tsv)
    nuclear=\$(awk '\$1 != "*" && \$1 !~ /^(chrM|MT|M)\$/ && \$2 > 0 { n++ } END { print n+0 }' idxstats.tsv)
    if [ "\${autosomes}" -ge 22 ]; then context='FULL_REFERENCE'; numt='EVALUABLE_FULL_REFERENCE'; auto='available_full_reference'; elif [ "\${autosomes}" -gt 0 ]; then context='PARTIAL_NUCLEAR_CONTEXT'; numt='NUMT_CONTEXT_INCOMPLETE'; auto='unavailable_incomplete_autosomal_context'; else context='MT_ONLY'; numt='NOT_EVALUABLE_MT_ONLY_INPUT'; auto='unavailable_no_nuclear_context'; fi
    printf 'sample\tinput_type\tsource_format\tsort_order\tindex_status\theader_sample_status\tmt_contig\tmt_contig_present\tcanonical_autosomes_present\tnuclear_contig_count\talignment_context\treference_compatibility\tnumt_evaluation_capability\tautosomal_depth_capability\tstatus\n' > ${sample}.input_alignment_validation.tsv
    printf '%s\tFASTQ\tFASTQ\tcoordinate\tgenerated\tNA\t%s\t%s\t%s\t%s\t%s\tPASS\t%s\t%s\tPASS\n' '${sample}' '${mt_contig}' "\$(awk -v mt='${mt_contig}' '\$1 == mt { print "true"; found=1 } END { if (!found) print "false" }' idxstats.tsv)" "\${autosomes}" "\${nuclear}" "\${context}" "\${numt}" "\${auto}" >> ${sample}.input_alignment_validation.tsv
    """
}
