# Controlled benchmark

This package contains the compact definitions, schemas, deterministic generator, analysis harness, frozen source summaries, and regression tests used for MitoArc's controlled technical benchmark. Generated reads, run directories, logs, per-run work products, and large results are intentionally excluded.

## Design

- Three deterministic replicates with distinct seeds
- Six depths: 100, 250, 500, 1000, 2000, and 5000×
- 36 truth alleles per generated input and 648 input truth observations
- SNV and indel classes, fixed VAF/context positions, controlled NUMT-like contexts, and boundary/non-boundary positions
- Exact normalized allele matching

Synthetic NUMT challenge labels are controlled technical contexts; they do not reproduce the full diversity of biological NUMTs. The benchmark establishes no optimal depth and is not clinical validation.

## Package layout

- `config/`: experiment-grid and portable orchestration settings
- `manifests/`: manifest schema and compact examples
- `truth/`: truth schema and frozen controlled truth definitions
- `scripts/benchmark_harness.py`: validation, deterministic generation, run orchestration, normalization, exact matching, resource accounting, and summarization
- `scripts/`: frozen panel/analysis utilities retained with stable identifiers where renaming would break source-data provenance
- `tests/`: unit and regression checks
- `summaries/`: compact frozen source summaries
- `../docs/validation/source_data/`: release-facing claim tables with explicit metric scopes

Generated state belongs under `benchmark/generated/`, `benchmark/runs/`, `benchmark/logs/`, `benchmark/metrics/`, and `benchmark/plots/`; these paths are excluded from the public release.

## Harness

```bash
python3 benchmark/scripts/benchmark_harness.py --help
python3 benchmark/scripts/benchmark_harness.py validate-manifest \
  --manifest benchmark/manifests/example_manifest.tsv
python3 -m unittest discover -s benchmark/tests -p 'test_*.py'
```

Manifest paths are resolved relative to the manifest directory. Executable rows require a reference, input, and truth table. Generated BAM/FASTQ files receive provenance TSVs, source inputs are never overwritten, and random sampling uses the manifest's explicit integer seed.

## Strategy mapping and metric scope

- Strategy A maps to native `--numt_strategy none` behavior.
- Strategies B and C require `--numt_strategy comparison`; B is the mapping-confidence-filtered call path and C is the evidence classification table.
- Strategy C fragment support is SNV-scoped; non-SNV evidence can be unresolved.
- GATK `standard_plus_shifted` uses the canonical filtered VCF. `standard_only` is a distinct raw/filter-aware benchmark view.
- mutserve2 consumes the standard mitochondrial realignment in v1.0.0. Its beta insertion/deletion options are not enabled.

NUMT comparisons use `OBSERVABLE_NUMT_ASSOCIATED_FP`. Pooled caller comparisons use `TOTAL_EXACT_TRUTH_FP`. Under these scopes, the authoritative pooled values are GATK Strategy A TP 548 / FP 33 / FN 100 and mutserve2 Strategy A TP 311 / FP 17 / FN 337.

The matched circular test found no incremental boundary-specific recovery under the tested fixture. Caller disagreement without truth is not treated as error.
