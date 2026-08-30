# Mitochondrial annotation and caller-specific consensus

MitoArc provides descriptive functional annotation and consensus generation. It
does not assign clinical significance, remove variants, fuse callers, or select
an authoritative callset.

## Annotation resource and method

The annotation layer is a lightweight pipeline-native annotator (`mtdna_annotate`
version 1.0.1) backed by the NCBI RefSeq record `NC_012920.1`. The pinned assets
are:

- `resources/mitochondrial/NC_012920.1.fa`: the verbatim 16,569 bp RefSeq FASTA;
- `resources/mitochondrial/NC_012920.1.features.tsv`: a compact derivation of the
  CDS, rRNA, tRNA, and circular D-loop intervals in the accession's feature
  table; and
- `resources/mitochondrial/NC_012920.1.provenance.tsv`: retrieval URLs, dates,
  source checksums, derivation notes, and the NCBI genetic-code source.

NCBI's molecular-data policy states that NCBI itself places no restriction on
use or distribution of its molecular data, while noting that NCBI cannot grant
rights asserted by original submitters. The small accession-specific resources
are therefore redistributed with explicit provenance and NCBI attribution. See
<https://www.ncbi.nlm.nih.gov/home/about/policies/> and
<https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1>.

This approach was selected instead of a whole-genome VEP cache because the
required scope is one fixed mitochondrial molecule and four feature classes.
It avoids a multi-gigabyte dependency while keeping the source accession and
derivation auditable. The annotator uses Python's standard library and the
pinned vertebrate mitochondrial genetic code (NCBI translation table 2).

For each ALT allele, the annotator reports the overlapping gene and feature,
region type, descriptive Sequence Ontology-style consequence, and a protein
change for simple coding SNVs when derivable. Coding indels are described as
frameshift or in-frame events, but no amino-acid change is fabricated. Caller
depth and VAF are retained with explicit source columns. The exact process-level
invocation is:

```text
mtdna_annotate \
  --sample SAMPLE --caller CALLER --vcf CALLER.vcf.gz \
  --reference standard_mt.fa \
  --canonical-reference NC_012920.1.fa \
  --features NC_012920.1.features.tsv \
  --numt-evidence SAMPLE.numt_variant_evidence.tsv \
  --homoplasmy-threshold 0.95 \
  --output SAMPLE.CALLER.annotated.tsv \
  --summary SAMPLE.CALLER.annotation_summary.tsv \
  --metadata SAMPLE.CALLER.annotation_metadata.tsv
```

Annotation runs only when `standard_mt.fa` is exactly 16,569 bases and its
sequence is identical to the pinned `NC_012920.1` sequence. Contig naming does
not affect compatibility. A truncated, synthetic, alternate, or otherwise
non-identical sequence retains caller metrics and NUMT evidence but has `NA`
biological fields and `annotation_status=UNAVAILABLE`; the pipeline does not
fail.

## Variant state and NUMT evidence

`--homoplasmy_threshold` defaults to `0.95`. VAF greater than or equal to the
threshold is `NEAR_HOMOPLASMIC`; lower measured VAF is `HETEROPLASMIC`; missing
VAF is `UNKNOWN`. This is descriptive, not a pathogenicity class. The 0.95
default follows the convention used by the gnomAD mitochondrial callset, which
defines 95--100% alternate fraction as homoplasmic. See
<https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/>.

When `--numt_strategy comparison` is active, the variant key
`sample/caller/contig/position/ref/alt` joins the NUMT `strategy_c_status` into
the annotation. No variant is removed. An evaluated call without an exact
evidence row is `NO_MATCHING_EVIDENCE`. With comparison mode disabled, every
row is `NOT_EVALUATED`.

## Consensus algorithm

`mtdna_consensus` version 1.0.0 independently consumes the filtered canonical
GATK VCF and canonical mutserve2 VCF. It starts from the unshifted
`standard_mt.fa`, considers only `PASS` or unfiltered records, checks every REF
allele before applying it, rejects overlapping or conflicting edits, and
applies edits in descending coordinate order so indels have deterministic
coordinates and length changes.

Calls at or above `--homoplasmy_threshold` replace REF. The conservative default
`--consensus_heteroplasmy_mode reference` leaves lower-VAF sites at REF.
`--consensus_heteroplasmy_mode iupac` represents heteroplasmic SNVs with the
IUPAC code for the set containing REF and all lower-VAF single-base ALT alleles:

| Alleles | Code | Alleles | Code |
|---|---:|---|---:|
| A/C | M | A/G | R |
| A/T | W | C/G | S |
| C/T | Y | G/T | K |
| A/C/G | V | A/C/T | H |
| A/G/T | D | C/G/T | B |
| A/C/G/T | N | | |

Heteroplasmic indels and other non-SNV alleles remain REF because no honest
single-character ambiguity representation exists. Near-homoplasmic indels may
change consensus length and the validation TSV records the length delta.
Records whose contig name contains `shifted`, REF mismatches, overlapping edits,
or multiple near-homoplasmic ALT conflicts produce `validation_status=FAILED_QC`.
The FASTA header records sample, caller, unshifted reference identity, threshold,
heteroplasmy mode, and canonical-reference compatibility.

## Functional fixture

`test_data/annotation/FUNCTIONAL.annotation.vcf` contains only synthetic calls
against the pinned rCRS sequence; they are test alleles, not patient findings.
The five sites exercise the control region (73A>G), RNR1 rRNA (1555A>G), TRNL1
tRNA (3243A>G), an ND1 synonymous SNV (3312C>T, `p.Pro2=`), and an ND1 missense
SNV (3314T>C, `p.Met3Thr`). VAFs exercise both states and both consensus modes.
Expected annotations and bases are pinned in `FUNCTIONAL.expected.tsv`.

## Outputs

```text
results/interpretation/
|-- annotation/
|   |-- SAMPLE.gatk.annotated.tsv
|   |-- SAMPLE.mutserve2.annotated.tsv
|   |-- SAMPLE.annotation_summary.tsv
|   |-- native/SAMPLE.CALLER.annotation_metadata.tsv
|   `-- functional/
`-- consensus/
    |-- SAMPLE.gatk.consensus.fa
    |-- SAMPLE.gatk.consensus_validation.tsv
    |-- SAMPLE.mutserve2.consensus.fa
    |-- SAMPLE.mutserve2.consensus_validation.tsv
    `-- functional/
```

The two callers remain separate throughout. Upstream VCFs are read-only inputs.
