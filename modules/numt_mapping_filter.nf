process NUMT_MAPPING_FILTER {
    tag "${sample} Strategy B MAPQ ${min_mapq}"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    publishDir "${params.outdir}/variants/numt", mode: 'copy', pattern: '*.numt_mapping_audit.tsv'

    input:
    tuple val(sample), path(alignment), path(index)
    val mt_contig
    val min_mapq

    output:
    tuple val(sample), path("${sample}.strategy_b.fullref.bam"), path("${sample}.strategy_b.fullref.bam.bai"), emit: filtered_bam
    path "${sample}.numt_mapping_audit.tsv", emit: audit

    script:
    """
    set -euo pipefail
    # Emit one compact record per alignment, spill the name grouping through
    # external sort, and reduce one read name at a time. This preserves the
    # original read-level decision exactly without retaining every read name in
    # memory for high-depth human inputs.
    samtools view ${alignment} | awk -F '\t' -v OFS='\t' -v mt='${mt_contig}' -v minq='${min_mapq}' '
      {
        secondary=int(\$2/256)%2
        supplementary=int(\$2/2048)%2
        eligible=(!secondary && !supplementary && \$3 == mt && \$5 >= minq)
        print \$1, (eligible ? 1 : 0), (eligible ? \$5 : "NA")
      }
    ' | sort -S 512M -T . -k1,1 | awk -F '\t' -v OFS='\t' '
      function emit_decision() {
        if (read_name == "") return
        decision=keep ? "KEEP" : "FLAG"
        reason=keep ? "PRIMARY_CHRM_MAPQ_GE_THRESHOLD" : "NO_PRIMARY_CHRM_MAPQ_GE_THRESHOLD"
        print read_name, decision, reason, (keep ? maxq : "NA")
      }
      {
        if (\$1 != read_name) {
          emit_decision()
          read_name=\$1
          keep=0
          maxq=0
        }
        if (\$2 == 1) {
          keep=1
          if (\$3 > maxq) maxq=\$3
        }
      }
      END { emit_decision() }
    ' > decisions.tsv
    awk -F '\t' '\$2 == "KEEP" {print \$1}' decisions.tsv > keep_names.txt
    samtools view -bh -N keep_names.txt ${alignment} > ${sample}.strategy_b.fullref.bam
    samtools index -@ ${task.cpus} ${sample}.strategy_b.fullref.bam ${sample}.strategy_b.fullref.bam.bai
    printf 'sample\tread_name\tstrategy_b_decision\treason\tmax_primary_chrM_mapq\tthreshold\n' > ${sample}.numt_mapping_audit.tsv
    awk -F '\t' -v OFS='\t' -v sample='${sample}' -v minq='${min_mapq}' '{print sample,\$1,\$2,\$3,\$4,minq}' decisions.tsv >> ${sample}.numt_mapping_audit.tsv
    """
}
