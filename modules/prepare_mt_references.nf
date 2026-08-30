process PREPARE_MT_REFERENCES {

    tag "MT reference preparation"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    path reference
    val mt_contig
    val shift_offset

    output:
    path 'mt_reference_bundle.tar'

    script:
    """
    set -euo pipefail
    awk -v requested='${mt_contig}' -v standard='${mt_contig}' -v shifted='${mt_contig}_shifted' -v offset='${shift_offset}' '
        BEGIN { in_contig = 0; found = 0; sequence = "" }
        /^>/ {
            name = substr(\$0, 2)
            split(name, name_fields, /[[:space:]]/)
            name = name_fields[1]
            in_contig = (name == requested)
            if (in_contig) { found = 1 }
            next
        }
        in_contig { gsub(/[[:space:]]/, ""); sequence = sequence \$0 }
        END {
            if (!found) {
                print "ERROR: Requested mitochondrial contig is not present in the reference FASTA: " requested > "/dev/stderr"
                exit 2
            }
            n = length(sequence)
            if (n == 0) { print "ERROR: mitochondrial sequence is empty" > "/dev/stderr"; exit 1 }
            if (offset !~ /^[0-9]+\$/ || offset <= 0 || offset >= n) {
                print "ERROR: shift_offset must be > 0 and < mitochondrial contig length (" n "); got " offset > "/dev/stderr"
                exit 1
            }
            print ">" standard > "standard_mt.fa"
            print sequence > "standard_mt.fa"
            print ">" shifted > "shifted_mt.fa"
            print substr(sequence, offset + 1) substr(sequence, 1, offset) > "shifted_mt.fa"
            print "original_contig\t" standard > "mt_shift_metadata.tsv"
            print "shifted_contig\t" shifted >> "mt_shift_metadata.tsv"
            print "mt_length\t" n >> "mt_shift_metadata.tsv"
            print "shift_offset\t" offset >> "mt_shift_metadata.tsv"
        }
    ' '${reference}'

    bwa-mem2 index -p standard_mt standard_mt.fa
    bwa-mem2 index -p shifted_mt shifted_mt.fa
    tar -cf mt_reference_bundle.tar standard_mt.* shifted_mt.* mt_shift_metadata.tsv
    """
}
