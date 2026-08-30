#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
harness="$root/benchmark/scripts/benchmark_harness.py"
manifest="$root/benchmark/manifests/mini_benchmark.tsv"
config="$root/benchmark/config/mini_benchmark.nextflow.config"
reference="$root/test_data/fixture-reference.fa"
truth="$root/benchmark/truth/mini_truth.tsv"
normalized_truth="$root/benchmark/generated/mini_truth.normalized.tsv"
prefix="$root/benchmark/generated/MINI"

python3 "$harness" validate-truth --truth "$truth" --reference "$reference" --output "$normalized_truth"

if [[ ! -s "${prefix}_R1.fastq.gz" || ! -s "${prefix}_R2.fastq.gz" || ! -s "${prefix}.provenance.tsv" || ! -s "${prefix}.achieved_truth.tsv" ]]; then
    python3 "$harness" generate-fastq \
        --benchmark-id mini_input_d100_vaf020_vaf080_r1 \
        --sample MINI --reference "$reference" --truth "$normalized_truth" \
        --output-prefix "$prefix" --depth 100 --seed 314159 \
        --read-length 100 --fragment-length 180
fi

python3 "$harness" validate-manifest --manifest "$manifest"

first_id=mini_d100_vaf020_A_gatk_circ_r1
second_id=mini_d100_vaf080_A_mutserve_std_r1

if [[ ! -d "$root/benchmark/runs/$first_id" && ! -d "$root/benchmark/runs/$second_id" ]]; then
    python3 "$harness" run --manifest "$manifest" --config "$config"
else
    echo 'ERROR: mini run directories already exist; preserve them and invoke the harness with --resume explicitly' >&2
    exit 2
fi

python3 "$harness" normalize-calls --manifest "$manifest"
python3 "$harness" match --manifest "$manifest"
python3 "$harness" resources --manifest "$manifest" --attempt 1

python3 "$harness" run --manifest "$manifest" --config "$config" --resume
python3 "$harness" resources --manifest "$manifest" --attempt 2

for benchmark_id in "$first_id" "$second_id"; do
    resource="$root/benchmark/metrics/$benchmark_id/benchmark_resource_metrics.attempt-002.tsv"
    awk -F '\t' 'NR == 2 && $13 == 0 { ok=1 } END { exit !ok }' "$resource"
done

python3 "$harness" summarize --manifest "$manifest"
echo 'MINI_BENCHMARK_PASS'
