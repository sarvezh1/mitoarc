# Validation summary

MitoArc v1.0.0 has two complementary validation layers: an end-to-end run on real HG002 sequencing data and a controlled truth benchmark. These are technical research validations, not clinical validation.

## Real HG002

The real-data run used genuine GIAB HG002 source data. A reduced indexed BAM retained complete chrM and selected real nuclear context. It exercised the BAM path end to end and generated variant, QC, interpretation, MultiQC, plot, and self-contained HTML reporting outputs.

- Mean mtDNA depth: 131,324×
- Haplogroup: H5a7 from both GATK- and mutserve2-derived callsets
- Contamination: unavailable (`INSUFFICIENT_DATA`)
- Copy-number proxy: unavailable because mapped autosomal context was incomplete

This was not complete-WGS nuclear context and not full-human FASTQ validation. No HG002 reads, alignments, variants, or human reference are distributed here.

## Controlled benchmark

The benchmark used three replicates, depths of 100/250/500/1000/2000/5000×, 36 truth alleles per generated input, 648 input truth observations, deterministic seeds, and exact normalized truth matching.

H1 NUMT comparisons use `OBSERVABLE_NUMT_ASSOCIATED_FP`: GATK Strategy A was TP 193 / FP 33 / FN 23 (recall 0.894), while Strategy B was TP 172 / FP 0 / FN 44 (recall 0.796).

H4 caller comparisons use `TOTAL_EXACT_TRUTH_FP`: GATK Strategy A was TP 548 / FP 33 / FN 100, and mutserve2 Strategy A was TP 311 / FP 17 / FN 337. These describe the evaluated frozen configurations and do not establish a universal caller ranking.

GATK recall rose from 0.741 at 100× to 0.917 at 5000×. This does not define an optimal depth. The matched circular test found no incremental boundary-specific recovery under the tested fixture.

## Interpretation constraints

- Strategy A is native/baseline caller behavior.
- Strategy B applies mapping-confidence filtering using full-reference context.
- Strategy C is an evidence classification framework; its fragment-level allele-support implementation is SNV-scoped, and indels may be unresolved.
- mutserve2 beta insertion/deletion options were disabled in the evaluated/default v1.0.0 configuration.
- Caller disagreement without truth is not evidence that either call is erroneous.
- Synthetic NUMT contexts, fixed truth positions, three replicates, and low-depth VAF quantization limit generalization.
- Adaptive downsampling was not implemented or validated, and no optimal depth was established.

## Source data

The `source_data/` directory contains compact frozen tables for the benchmark design, primary scoped results, NUMT observability, depth/resource trends, caller comparison, and HG002 summary. Public claims should preserve the metric-scope column.
