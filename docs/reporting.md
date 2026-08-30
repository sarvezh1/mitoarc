# Reporting outputs

The reporting workflow renders existing machine-readable pipeline results as
portable technical summaries and publication-oriented figures. It does not
alter variant calls, merge callers, select a preferred caller, resolve NUMT
strategy choices, or add clinical interpretation.

## Reproducible environments

MultiQC runs from the immutable image
`quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0` pinned by digest. The
portable `simple` template is generated with network version checks and optional
online integrations disabled. FastQC archives, fastp JSON, and samtools
flagstat files are copied into a clean staging directory so analysis paths are
not exposed. Runtime-local prefixes are removed from the parsed MultiQC source
records, and the portable data directory is retained beside the report.

Static plots and the integrated report use the local
`mitoarc-reporting:1.0.0` image. Its base image is pinned by digest and
all Python dependencies are version-pinned: Python 3.11, Matplotlib 3.10.5,
NumPy 2.3.2, Pillow 11.3.0, contourpy 1.3.3, cycler 0.12.1, fonttools 4.63.0,
kiwisolver 1.5.0, packaging 26.3, pyparsing 3.3.2, python-dateutil
2.9.0.post0, and six 1.17.0. Rendering uses the non-interactive Matplotlib Agg
backend, a fixed style, and deterministic SVG hashing.

Build the plotting image before using the Docker profile:

```text
docker build -t mitoarc-reporting:1.0.0 containers/reporting
```

## Data contracts

For each sample, reporting writes:

- `SAMPLE.reporting_summary.tsv` and a JSON equivalent: one auditable row with
  coverage, copy-number proxy, caller counts, caller-specific variant states,
  caller comparison, haplogroups, contamination, NUMT counts, annotation and
  consensus status, thresholds, reference compatibility, and tool versions;
- `SAMPLE.reporting_depth.tsv`: explicit canonical-position depth used by the
  linear and circular coverage tracks;
- `SAMPLE.reporting_variants.tsv`: the annotation rows retained independently
  by caller, with the caller-comparison relationship appended;
- `SAMPLE.reporting_caller_comparison.tsv`: shared and caller-only sites with
  both callers' original status, depth, and VAF fields;
- `SAMPLE.reporting_features.tsv`: the canonical mitochondrial feature track,
  emitted only for a sequence-compatible `NC_012920.1` reference;
- `SAMPLE.mitoplot_data.tsv`: one coverage row per position, expanded only at
  variant positions to keep each caller row inspectable; and
- `SAMPLE.plot_manifest.tsv`: plotter version, reference compatibility, genome
  length, circular-view status, callers, and plotted variant positions.

Per-base reporting depth is computed from the existing extracted mtDNA BAM with
the same `--coverage_mapq` setting used by the coverage metrics. The process
uses mosdepth 0.3.8 pinned by digest and publishes the resulting TSV; plotting
does not parse BAM files.

## Figures

Every sample receives SVG and 180-dpi PNG versions of:

- a linear, unsmoothed mtDNA depth profile;
- percentage of bases at 1×, 10×, 100×, 500×, and 1000×;
- caller-specific VAF by canonical position, with marker shape encoding
  heteroplasmic versus near-homoplasmic state and a separate outline for
  NUMT-risk evidence;
- supported functional-consequence counts by caller;
- shared, GATK-only, and mutserve2-only call counts; and
- a circular mitochondrial genome view.

The circular view is drawn locally in polar coordinates. Protein-coding genes,
rRNAs, tRNAs, and the control region are separate feature classes; strand is
represented by inner and outer feature bands. Coverage is a log-scaled radial
track used only to make the wide dynamic range visible. Variants occupy
caller-specific radii, marker shape preserves variant state, and a red outline
marks NUMT-associated mapping evidence. Position 1 and the standard reference
breakpoint are explicitly marked. No genomic smoothing or caller fusion is
performed.

For a truncated or sequence-incompatible reference, the ordinary depth and
threshold figures remain available as clearly labelled technical figures.
Biological consequence and circular views instead contain an unavailable
status panel. They are never presented as canonical human mitochondrial
interpretation.

## Integrated report

`SAMPLE.mtDNA_report.html` contains embedded CSS, small vanilla JavaScript for
table filtering and sorting, and base64-embedded PNG figures. It has no CDN or
external asset requirement and can be opened directly from disk. Its sections
cover sample overview, mtDNA QC, circular view, coverage, variant landscape,
caller comparison, NUMT evidence, haplogroups, contamination, functional
annotation, caller-preserving variant rows, caller-specific consensus status,
technical QC, and reproducibility metadata.

The HTML is a downstream view. All scientific values remain separately
available in the source or reporting-contract TSV/JSON files. The report never
reads terminal logs. `multiqc_report.html` remains a separate technical-QC
artifact.

## Validation fixture

When `--annotation_functional_fixture` and
`--annotation_functional_expected` are supplied, a deterministic 16,569 bp
reporting fixture is generated under `results/plots/functional/` and
`results/report/functional/`. Its depth curve includes a reproducible coverage
drop and 100 uncovered bases. The existing canonical test VCF places variants
at positions 73, 1555, 3243, 3312, and 3314, covering control-region, rRNA,
tRNA, synonymous, and missense annotations as well as heteroplasmic and
near-homoplasmic states. These are synthetic test alleles, not patient
findings.

Automated validation checks upstream metric equality, exact preservation of
annotation rows, caller disagreement and VAF preservation, NUMT status counts,
continuous canonical coordinates, expected functional positions and features,
SVG XML validity, PNG signatures and dimensions, embedded report assets, and
absence of private paths or external asset dependencies.

## Output layout

```text
results/
|-- plots/
|   |-- SAMPLE.mtDNA_depth_profile.svg / .png
|   |-- SAMPLE.coverage_thresholds.svg / .png
|   |-- SAMPLE.vaf_landscape.svg / .png
|   |-- SAMPLE.variant_consequences.svg / .png
|   |-- SAMPLE.caller_comparison.svg / .png
|   |-- SAMPLE.mitoplot.svg / .png
|   `-- functional/
`-- report/
    |-- SAMPLE.mtDNA_report.html
    |-- SAMPLE.reporting_summary.tsv
    |-- data/
    |-- functional/
    |-- validation/
    |-- multiqc_report.html
    `-- multiqc_report_data/
```
