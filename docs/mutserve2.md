# Independent mutserve2 call path

MitoArc runs mutserve2 as an independent mitochondrial caller and provides a
descriptive comparison with the GATK baseline. It does not
perform caller fusion, consensus generation, caller ranking, annotation, or
threshold tuning.

## Version and container provenance

- mutserve2: `2.0.3`, the current stable (non-prerelease) upstream release
  published on 2024-12-23.
- Release: <https://github.com/seppinho/mutserve/releases/tag/v2.0.3>
- Release archive SHA-256:
  `59b85d263906b73e3b8c9bbd21c29557af8e50e5d17d3926850b1877ce19ce2c`
- Container definition: `containers/mutserve2/Dockerfile`
- Local image name: `mitoarc-mutserve2:2.0.3`
- Runtime base:
  `eclipse-temurin:17.0.13_11-jre-jammy@sha256:7b8b90b35214757840f79fc78794bf4b06cb6ac92a8853429de16366e837c302`

The upstream project does not publish a dedicated maintained v2.0.3 caller
container. The pinned mtDNA-Server 2 container maintained by the same group
contains mutserve2, but also bundles unrelated analysis and reporting tools.
The local multi-stage definition therefore installs only the verified official
mutserve2 release into a fixed JRE image.

Build the image once before using the Docker profile:

```bash
docker build --pull --tag mitoarc-mutserve2:2.0.3 containers/mutserve2
```

## Input and method

mutserve2 documents a sorted and indexed BAM/CRAM as its input and calls the
mitochondrial contig directly. The pipeline supplies the canonical
standard-mtDNA realignment BAM together with the matching
`standard_mt.fa`. This BAM avoids the original whole-genome alignment context
and is already coordinate sorted, indexed, read-group tagged, and realigned to
the canonical mitochondrial sequence. A shifted-reference mutserve2 branch is
not part of the documented mutserve2 method and is not added.

The native minimum heteroplasmy level is `0.01`. It is exposed as
`--mutserve2_heteroplasmy_threshold`, defaults to `0.01`, and is written to each
sample metadata file. Other native calling defaults remain unchanged: mapping
quality 20, base quality 20, alignment quality 30, and BAQ disabled.

The exact caller command is:

```bash
mutserve call \
    --reference standard_mt.fa \
    --output <sample>.mutserve2.native.vcf.gz \
    --threads <task.cpus> \
    --level <mutserve2_heteroplasmy_threshold> \
    --contig-name <mt_contig> \
    --no-ansi \
    --write-raw \
    <sample>.standard_mt.sorted.bam
```

mutserve2 v2.0.3 emits calls marked `PASS` and has no separate filtered-call
product for this invocation. The pipeline does not invent an additional
filtering layer. It retains the native VCF, native call table, and native raw
evidence table. A separate bcftools 1.20 process reheaders the VCF sample to the
pipeline sample ID, splits multiallelic records if present, checks every REF
allele against `standard_mt.fa`, sorts, bgzip-compresses, and tabix-indexes the
canonical VCF.

`alt_depth` is `NA` in the normalized calls table because mutserve2's VCF
reports `DP` and `AF`, but not an exact allele-depth (`AD`) field. The pipeline
does not fabricate a caller-equivalent AD by multiplying rounded AF by DP.

## Workflow architecture

```text
standard mtDNA realignment BAM + canonical reference bundle
  -> MUTSERVE2_CALL
       -> native VCF, call TSV, raw-evidence TSV, metadata TSV
  -> NORMALIZE_MUTSERVE2
       -> canonical raw VCF + TBI, normalized calls TSV
  -> VALIDATE_MUTSERVE2_FIXTURE (when --fixture_truth is set)
       -> truth counts and per-truth-variant metrics

frozen canonical GATK raw/filtered VCFs + canonical mutserve2 VCF + truth
  -> COMPARE_CALLERS
       -> descriptive caller comparison and caller summary TSVs
```

GATK call presence in the comparison is based on the frozen canonical raw VCF;
`gatk_status` is read from the frozen filtered VCF. No combined or final call is
created. Caller comparison runs for every sample. When `--fixture_truth` is not
provided, `truth_status` and truth-derived TP/FP/FN fields are `NA`; shared and
caller-only counts remain available.

## Frozen fixture result

The unchanged truth set is `chrM:25 T>C` and `chrM:200 T>C`.

| Caller | Variant | Status | Depth | VAF |
| --- | --- | --- | ---: | ---: |
| GATK | chrM:25 T>C | PASS | 24 | 0.929 |
| mutserve2 | chrM:25 T>C | PASS | 24 | 1.00 |
| GATK | chrM:200 T>C | PASS | 12 | 0.929 |
| mutserve2 | chrM:200 T>C | PASS | 12 | 1.00 |

Both independent baselines have TP=2, FP=0, FN=0. There are two shared calls,
zero GATK-only calls, and zero mutserve2-only calls. The AF values are retained
as caller-specific estimates and are not expected to be numerically identical.

## Outputs

```text
results/variants/
├── mutserve2/
│   ├── <sample>.mutserve2.raw.vcf.gz
│   ├── <sample>.mutserve2.raw.vcf.gz.tbi
│   ├── <sample>.mutserve2.calls.tsv
│   ├── <sample>.mutserve2.truth_comparison.tsv
│   ├── <sample>.mutserve2.truth_variants.tsv
│   └── native/
│       ├── <sample>.mutserve2.native.vcf.gz
│       ├── <sample>.mutserve2.native.tsv
│       ├── <sample>.mutserve2.native_raw.tsv
│       └── <sample>.mutserve2.metadata.tsv
└── comparison/
    ├── <sample>.caller_comparison.tsv
    └── <sample>.caller_summary.tsv
```

This implementation establishes independent caller baselines only. Comparative
performance requires later real-data or simulation benchmarking and is not a
novelty claim for either caller or this comparison.
