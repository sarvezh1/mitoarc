process VALIDATE_GATK_FIXTURE {
    tag "${sample} GATK fixture truth"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    publishDir "${params.outdir}/variants/gatk", mode: 'copy', pattern: '*.truth_comparison.tsv'

    input:
    tuple val(sample), path(raw_vcf), path(raw_index), path(raw_stats)
    path truth

    output:
    path "${sample}.gatk.truth_comparison.tsv"

    script:
    """
    set -euo pipefail
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${raw_vcf} | sort -u > observed.tsv
    awk -F '\\t' 'NR > 1 {print \$1,\$2,\$3,\$4}' OFS='\\t' ${truth} | sort -u > expected.tsv
    awk 'FILENAME == ARGV[1] {e[\$0]=1; next} {if (e[\$0]) tp++; else fp++} END {print tp+0, fp+0}' expected.tsv observed.tsv > counts.observed
    awk 'FILENAME == ARGV[1] {o[\$0]=1; next} {if (!o[\$0]) fn++} END {print fn+0}' observed.tsv expected.tsv > counts.missed
    printf 'sample\\ttrue_positives\\tfalse_positives\\tfalse_negatives\\n' > ${sample}.gatk.truth_comparison.tsv
    printf '%s\\t%s\\t%s\\t%s\\n' '${sample}' \$(awk '{print \$1}' counts.observed) \$(awk '{print \$2}' counts.observed) \$(cat counts.missed) >> ${sample}.gatk.truth_comparison.tsv
    """
}
