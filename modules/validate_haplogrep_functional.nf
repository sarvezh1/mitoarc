process VALIDATE_HAPLOGREP_FUNCTIONAL {
    tag "${sample} HaploGrep functional fixture"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/haplogroup/functional", mode: 'copy'

    input:
    tuple val(sample), val(caller), path(standardized), path(native_output)
    path expected

    output:
    path "${sample}.haplogrep_functional_validation.tsv"
    path "${sample}.functional_fixture.haplogroup.tsv"

    script:
    """
    set -euo pipefail
    validate_haplogrep_fixture ${standardized} ${expected} \
        ${sample}.haplogrep_functional_validation.tsv
    """
}
