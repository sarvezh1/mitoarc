process FINALIZE_ALIGNMENT_CONTRACT {
    tag "${sample} standardized alignment"

    input:
    tuple val(sample), val(input_mode), path(alignment), path(index), val(reference), val(context_hint), path(metadata)

    output:
    tuple val(sample), val(input_mode), path("${sample}.contract.bam"), path("${sample}.contract.bam.bai"), val(reference), val(context_hint), path(metadata), emit: contract

    script:
    """
    set -euo pipefail
    if [ "\$(basename ${alignment})" != '${sample}.contract.bam' ]; then ln -s ${alignment} ${sample}.contract.bam; fi
    if [ "\$(basename ${index})" != '${sample}.contract.bam.bai' ]; then ln -s ${index} ${sample}.contract.bam.bai; fi
    """
}
