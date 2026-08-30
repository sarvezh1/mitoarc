process PREPARE_GATK_MT_REFERENCES {
    tag "GATK mt reference assets"
    container 'broadinstitute/gatk:4.6.2.0'

    input:
    path mt_reference_bundle

    output:
    path 'gatk_mt_reference_bundle.tar'

    script:
    """
    set -euo pipefail
    tar -xf ${mt_reference_bundle}

    samtools faidx standard_mt.fa
    samtools faidx shifted_mt.fa
    gatk CreateSequenceDictionary -R standard_mt.fa -O standard_mt.dict
    gatk CreateSequenceDictionary -R shifted_mt.fa -O shifted_mt.dict

    mt_length=\$(awk '\$1 == "mt_length" {print \$2}' mt_shift_metadata.tsv)
    shift_offset=\$(awk '\$1 == "shift_offset" {print \$2}' mt_shift_metadata.tsv)
    standard_contig=\$(awk '\$1 == "original_contig" {print \$2}' mt_shift_metadata.tsv)
    shifted_contig=\$(awk '\$1 == "shifted_contig" {print \$2}' mt_shift_metadata.tsv)
    standard_end=\$((mt_length))
    shifted_sensitive_start=\$((mt_length - shift_offset + 1))

    # The standard call covers the non-rotated part (offset+1..length).
    # The shifted call covers the final shifted segment, which maps to
    # canonical coordinates 1..offset across the circular origin.
    printf '%s:%s-%s\n' "\${standard_contig}" "\$((shift_offset + 1))" "\${standard_end}" > standard_mt.interval
    printf '%s:%s-%s\n' "\${shifted_contig}" "\${shifted_sensitive_start}" "\${mt_length}" > shifted_mt.interval

    # Picard LiftOver consumes intervals on the UCSC chain target and emits
    # query coordinates. Target is therefore shifted mtDNA and query is the
    # canonical standard mtDNA. Two blocks describe the rotation explicitly.
    {
        printf 'chain 0 %s %s + 0 %s %s %s + %s %s 1\\n' "\${shifted_contig}" "\${mt_length}" "\$((mt_length - shift_offset))" "\${standard_contig}" "\${mt_length}" "\${shift_offset}" "\${mt_length}"
        printf '%s\\n' "\$((mt_length - shift_offset))"
        printf '\\n'
        printf 'chain 0 %s %s + %s %s %s %s + 0 %s 2\\n' "\${shifted_contig}" "\${mt_length}" "\$((mt_length - shift_offset))" "\${mt_length}" "\${standard_contig}" "\${mt_length}" "\${shift_offset}"
        printf '%s\\n' "\${shift_offset}"
        printf '\\n'
    } > shifted_to_standard.chain

    tar -cf gatk_mt_reference_bundle.tar \
        standard_mt.fa standard_mt.fa.fai standard_mt.dict \
        shifted_mt.fa shifted_mt.fa.fai shifted_mt.dict \
        mt_shift_metadata.tsv standard_mt.interval shifted_mt.interval \
        shifted_to_standard.chain
    """
}
