process REALIGN_MT {
    tag "${sample} ${representation}"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'
    input:
    tuple val(sample), path(reads1), path(reads2), path(reference_bundle)
    val representation
    output:
    tuple val(sample), val(representation), path("${sample}.${representation}_mt.sam")
    script:
    def prefix = representation == 'standard' ? 'standard_mt' : 'shifted_mt'
    """
    set -euo pipefail
    tar -xf ${reference_bundle}
    bwa-mem2 mem -t ${task.cpus} -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' ${prefix} ${reads1} ${reads2} > ${sample}.${representation}_mt.sam
    """
}
