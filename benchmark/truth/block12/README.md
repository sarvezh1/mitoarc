# Controlled truth design

The 36 true mitochondrial alleles comprise six nominal VAF targets crossed
with six fixed archetypes. Sites are selected deterministically from the
canonical 16,569 bp `chrM` sequence. Non-boundary sites are more than one read
length apart; boundary sites are within 1 kb of the origin/end and are also
more than one read length apart in circular distance. Candidate sites require
local base diversity and avoid a same-base predecessor for indel anchors.

The technical placement avoids deliberately selecting known clinical alleles;
this is not a pathogenicity panel. Exact positions are frozen in
`block12_truth_panel.nominal.tsv` once generated. Replicates reuse these
conceptually matched positions while independently sampling fragments with
distinct seeds, isolating stochastic sampling from sequence-context effects.

Controlled NUMT challenges use synthetic nuclear contigs containing exact
mitochondrial homolog segments and deterministic unique flanks. NUMT_FALSE
pairs combine an ALT-bearing homolog read with a nuclear-unique mate;
AMBIGUOUS pairs use homolog-compatible mates. True mitochondrial alleles in
NUMT_LIKE and AMBIGUOUS_CONTEXT panel strata remain TRUE_MT records. These are
controlled mapping challenges, not complete biological NUMT models.

The benchmark reference also contains deterministic neutral scaffolds named
`chr1` through `chr22`. They satisfy the frozen FASTQ alignment contract needed
to exercise Strategy B/C; they do not represent complete human autosomes and
must never be described as complete WGS context. Accuracy conclusions are
limited to the explicitly controlled homolog and challenge records.
