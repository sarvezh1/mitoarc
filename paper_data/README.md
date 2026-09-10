# MitoArc paper data

These tab-separated files contain machine-readable values underlying selected benchmark results reported in the MitoArc manuscript and Supplementary Information. They correspond to MitoArc version 1.0.0.

- `software_versions.tsv`: verified software and version manifest.
- `caller_depth_metrics.tsv`: caller TP, FP, FN, precision, recall, F1, and achieved-VAF error by depth. FP uses `TOTAL_EXACT_TRUTH_FP`.
- `numt_observability.tsv`: designed, generated, mapped, extracted, and observable NUMT challenge counts by depth and replicate.
- `numt_strategy_metrics.tsv`: Strategy A/B TP, observable NUMT-associated FP, FN, and recall for technical replicates and pooled results. FP uses `OBSERVABLE_NUMT_ASSOCIATED_FP`.
- `circular_matched_metrics.tsv`: matched standard-only and standard-plus-shifted circular-coordinate results.
- `resource_metrics.tsv`: elapsed time, CPU hours, peak RSS, work-directory size, and output size for complete physical workflow executions.

The three benchmark replicates are technical stochastic replicates, not biological replicates. VAF error is calculated against achieved VAF among recovered truth alleles. Strategy C is an evidence-classification view rather than a separate caller and is not represented as a separate caller dataset. Resource measurements include caller and NUMT-comparison branches and are not caller-specific.
