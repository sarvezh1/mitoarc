process VALIDATE_ANNOTATION_FUNCTIONAL {
    tag "FUNCTIONAL mitochondrial annotation and consensus"
    container 'mitoarc-haplogrep:3.3.2'
    publishDir "${params.outdir}/interpretation/annotation/functional", mode: 'copy', pattern: '*.annotated.tsv'
    publishDir "${params.outdir}/interpretation/annotation/functional", mode: 'copy', pattern: '*.annotation_*.tsv'
    publishDir "${params.outdir}/interpretation/consensus/functional", mode: 'copy', pattern: '*.consensus.fa'
    publishDir "${params.outdir}/interpretation/consensus/functional", mode: 'copy', pattern: '*.consensus_validation.tsv'

    input:
    path fixture_vcf
    path expected
    path canonical_reference
    path features
    path no_numt_evidence
    val homoplasmy_threshold

    output:
    path "FUNCTIONAL.functional_fixture.annotated.tsv"
    path "FUNCTIONAL.functional_fixture.annotation_summary.tsv"
    path "FUNCTIONAL.functional_fixture.annotation_metadata.tsv"
    path "FUNCTIONAL.reference.consensus.fa"
    path "FUNCTIONAL.reference.consensus_validation.tsv"
    path "FUNCTIONAL.iupac.consensus.fa"
    path "FUNCTIONAL.iupac.consensus_validation.tsv"
    path "FUNCTIONAL.annotation_consensus_functional_validation.tsv"

    script:
    """
    set -euo pipefail
    mtdna_annotate \
        --sample FUNCTIONAL \
        --caller functional_fixture \
        --vcf ${fixture_vcf} \
        --reference ${canonical_reference} \
        --canonical-reference ${canonical_reference} \
        --features ${features} \
        --numt-evidence ${no_numt_evidence} \
        --homoplasmy-threshold '${homoplasmy_threshold}' \
        --output FUNCTIONAL.functional_fixture.annotated.tsv \
        --summary FUNCTIONAL.functional_fixture.annotation_summary.tsv \
        --metadata FUNCTIONAL.functional_fixture.annotation_metadata.tsv

    mtdna_consensus \
        --sample FUNCTIONAL --caller functional_fixture \
        --vcf ${fixture_vcf} --reference ${canonical_reference} \
        --canonical-reference ${canonical_reference} \
        --homoplasmy-threshold '${homoplasmy_threshold}' \
        --heteroplasmy-mode reference \
        --output FUNCTIONAL.reference.consensus.fa \
        --validation FUNCTIONAL.reference.consensus_validation.tsv

    mtdna_consensus \
        --sample FUNCTIONAL --caller functional_fixture \
        --vcf ${fixture_vcf} --reference ${canonical_reference} \
        --canonical-reference ${canonical_reference} \
        --homoplasmy-threshold '${homoplasmy_threshold}' \
        --heteroplasmy-mode iupac \
        --output FUNCTIONAL.iupac.consensus.fa \
        --validation FUNCTIONAL.iupac.consensus_validation.tsv

    validate_mtdna_functional \
        --annotated FUNCTIONAL.functional_fixture.annotated.tsv \
        --expected ${expected} \
        --canonical-reference ${canonical_reference} \
        --reference-consensus FUNCTIONAL.reference.consensus.fa \
        --iupac-consensus FUNCTIONAL.iupac.consensus.fa \
        --reference-validation FUNCTIONAL.reference.consensus_validation.tsv \
        --iupac-validation FUNCTIONAL.iupac.consensus_validation.tsv \
        --output FUNCTIONAL.annotation_consensus_functional_validation.tsv
    """
}
