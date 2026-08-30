process PREPARE_MT_READS {
    tag "${sample} mtDNA reads"
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
    input:
    tuple val(sample), path(alignment), path(index)
    output:
    tuple val(sample), path("${sample}.mtDNA.R1.fastq.gz"), path("${sample}.mtDNA.R2.fastq.gz"), emit: prepared_reads
    script:
    """
    set -euo pipefail
    samtools collate -@ ${task.cpus} -Ou ${alignment} \\
      | samtools fastq -@ ${task.cpus} -n \\
          -1 ${sample}.mtDNA.R1.fastq.gz \\
          -2 ${sample}.mtDNA.R2.fastq.gz \\
          -0 /dev/null -s /dev/null -
    """
}
