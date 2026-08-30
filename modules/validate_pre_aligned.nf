process VALIDATE_PREALIGNED {
    tag "${sample} ${input_mode}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/qc/input", mode: 'copy', pattern: '*.input_alignment_validation.tsv'

    input:
    tuple val(sample), val(input_mode), path(source), path(reference), path(input_index)

    output:
    tuple val(sample), val(input_mode), path("${sample}.full_reference.bam"), path("${sample}.full_reference.bam.bai"), val(reference), path("${sample}.input_alignment_validation.tsv"), emit: standard_alignment

    script:
    """
    set -euo pipefail
    source_format='${input_mode}'
    samtools quickcheck -v ${source}
    samtools view -H ${input_mode == 'CRAM' ? '-T ' + reference : ''} ${source} > header.sam
    index_status='provided'
    if [ "${input_index}" = "no_index" ] || [ ! -s "${input_index}" ]; then
        index_status='generated'
        if [ "${input_mode}" = 'CRAM' ]; then
            samtools index -@ ${task.cpus} ${source}
        else
            samtools index -@ ${task.cpus} ${source}
        fi
    fi

    sort_order=\$(awk '\$1 == "@HD" { for (i=1; i<=NF; i++) if (\$i ~ /^SO:/) { sub(/^SO:/, "", \$i); print \$i; exit } }' header.sam)
    if [ "\${sort_order}" != 'coordinate' ]; then
        echo "${input_mode} for sample '${sample}' must be coordinate-sorted (header SO:coordinate); found '\${sort_order:-missing}'. Re-sort the source before running the pipeline." >&2
        exit 1
    fi

    if [ "${input_mode}" = 'CRAM' ]; then
        if ! samtools view -T ${reference} -c ${source} >/dev/null; then
            echo "CRAM '${source}' cannot be decoded with the supplied reference '${reference}'." >&2
            exit 1
        fi
    fi

    samtools idxstats ${source} > idxstats.tsv
    samtools faidx ${reference}
    # Validate contigs that carry mapped records. BAM headers may legitimately
    # retain unused decoy @SQ entries that are absent from an analysis-set FASTA.
    awk 'BEGIN { bad=0 } FNR==NR { ref[\$1]=\$2; next } \$1 != "*" && \$3 > 0 { seen[\$1]=1; if (!(\$1 in ref)) { print "alignment contig absent from supplied reference: " \$1 > "/dev/stderr"; bad=1 } else if (ref[\$1] != \$2) { print "contig length mismatch: " \$1 " alignment=" \$2 " reference=" ref[\$1] > "/dev/stderr"; bad=1 } } END { exit bad }' ${reference}.fai idxstats.tsv

    mt_present=\$(awk -v mt='${params.mt_contig}' '\$1 == mt { print "true"; found=1 } END { if (!found) print "false" }' idxstats.tsv)
    if [ "\${mt_present}" != 'true' ]; then
        echo "Requested mitochondrial contig '${params.mt_contig}' is absent from ${input_mode} '${source}'." >&2
        exit 1
    fi
    canonical_autosomes=\$(awk '\$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])\$/ && \$2 > 0 { n++ } END { print n+0 }' idxstats.tsv)
    nuclear_count=\$(awk '\$1 != "*" && \$1 !~ /^(chrM|MT|M)\$/ && \$2 > 0 { n++ } END { print n+0 }' idxstats.tsv)
    if [ "\${canonical_autosomes}" -eq 0 ]; then context='MT_ONLY'; numt_capability='NOT_EVALUABLE_MT_ONLY_INPUT'; auto_capability='unavailable_no_nuclear_context'
    elif [ "\${canonical_autosomes}" -lt 22 ]; then context='PARTIAL_NUCLEAR_CONTEXT'; numt_capability='NUMT_CONTEXT_INCOMPLETE'; auto_capability='unavailable_incomplete_autosomal_context'
    else context='FULL_REFERENCE'; numt_capability='EVALUABLE_FULL_REFERENCE'; auto_capability='available_full_reference'; fi

    sm_values=\$(awk -F'\t' '\$1 ~ /^@RG/ { for (i=1; i<=NF; i++) if (\$i ~ /^SM:/) { sub(/^SM:/, "", \$i); print \$i } }' header.sam | sort -u)
    sm_count=\$(printf '%s\n' "\${sm_values}" | awk 'NF { n++ } END { print n+0 }')
    if [ "\${sm_count}" -eq 0 ]; then sm_status='WARNING_MISSING_SM'
    elif [ "\${sm_count}" -gt 1 ]; then echo "${input_mode} '${source}' contains multiple biological SM values: \$(printf '%s' "\${sm_values}" | paste -sd, -)." >&2; exit 1
    elif [ "\$(printf '%s\n' "\${sm_values}")" = '${sample}' ]; then sm_status='PASS_EXACT_SINGLE_SM'
    else sm_status='WARNING_SM_DIFFERS_SAMPLESHEET_AUTHORITY'; fi

    if [ "${input_mode}" = 'CRAM' ]; then
        samtools view -T ${reference} -b ${source} > ${sample}.full_reference.bam
    else
        if [ "\$(basename ${source})" != '${sample}.full_reference.bam' ]; then ln -s ${source} ${sample}.full_reference.bam; fi
    fi
    samtools index -@ ${task.cpus} ${sample}.full_reference.bam ${sample}.full_reference.bam.bai
    printf 'sample\tinput_type\tsource_format\tsort_order\tindex_status\theader_sample_status\tmt_contig\tmt_contig_present\tcanonical_autosomes_present\tnuclear_contig_count\talignment_context\treference_compatibility\tnumt_evaluation_capability\tautosomal_depth_capability\tstatus\n' > ${sample}.input_alignment_validation.tsv
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\ttrue\t%s\t%s\t%s\tPASS\t%s\t%s\tPASS\n' '${sample}' '${input_mode}' "\${source_format}" "\${sort_order}" "\${index_status}" "\${sm_status}" '${params.mt_contig}' "\${canonical_autosomes}" "\${nuclear_count}" "\${context}" "\${numt_capability}" "\${auto_capability}" >> ${sample}.input_alignment_validation.tsv
    """
}
