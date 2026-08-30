# Pre-aligned BAM/CRAM entry modes

The pipeline accepts exactly one samplesheet schema per run:

```text
sample,fastq_1,fastq_2
sample,bam
sample,cram
```

Pre-aligned inputs are validated before downstream work. Files must be non-empty,
pass `samtools quickcheck`, have a readable header with `SO:coordinate`, contain
the requested mitochondrial contig, and be compatible with the supplied FASTA.
Missing BAI/CRAI files are generated in the Nextflow work directory; the source
directory is never modified. CRAM is decoded once to a standardized BAM in the
work directory.

Reference preparation has separate execution contracts. FASTQ mode prepares the
whole alignment reference with BWA-MEM2 indexes and a samtools FASTA index before
whole-genome alignment. BAM and CRAM modes schedule only
`PREPARE_DOWNSTREAM_REFERENCE`, which creates a samtools FASTA index in the
Nextflow work/cache area. They never schedule whole-reference BWA-MEM2 indexing.
The original FASTA remains the authoritative reference for compatibility checks,
CRAM decoding, mitochondrial reference extraction, annotation, and consensus.

The standardized internal alignment contract is `(sample, input_mode,
full_reference_bam, bam_index, reference, alignment_context_status,
input_validation_metadata)`. Existing scientific modules consume its legacy
three-field projection, so caller, NUMT, annotation, consensus, and plot logic
are not duplicated or changed.

Contig context is classified as `FULL_REFERENCE` (all canonical autosomes),
`PARTIAL_NUCLEAR_CONTEXT` (some canonical autosomes), or `MT_ONLY`. The input
validation TSV records the corresponding NUMT and autosomal-depth capabilities.
MT-only input remains eligible for mtDNA extraction and mtDNA analyses, while
autosomal depth/copy-number are unavailable and NUMT evidence is explicitly
`NOT_EVALUABLE_MT_ONLY_INPUT`.

The early-validation harness covers missing-input and mixed-schema failures.
`tests/prepare_input_fixtures.sh` derives a CRAM and an mtDNA-only BAM from the
existing tiny alignment using samtools; it downloads no data. Corrupt, unsorted,
incompatible-reference, missing-mt-contig, and multiple-SM inputs are rejected
by `VALIDATE_PREALIGNED` before downstream analysis.

Examples:

```bash
nextflow run main.nf -profile docker --input samplesheet_fastq.csv --reference GRCh38.fa
nextflow run main.nf -profile docker --input samplesheet_bam.csv --reference GRCh38.fa
nextflow run main.nf -profile docker --input samplesheet_cram.csv --reference GRCh38.fa
```
