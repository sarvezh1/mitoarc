#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"

preview_dir=$(mktemp -d)
trap 'rm -r -- "$preview_dir"' EXIT

preview() {
    local mode=$1
    local samplesheet=$2
    NXF_OFFLINE=true nextflow run main.nf -preview -profile docker \
        --input "$samplesheet" \
        --reference test_data/fixture-reference.fa \
        --shift_offset 250 \
        --outdir "$preview_dir/${mode}-results" \
        -work-dir "$preview_dir/${mode}-work" \
        -with-dag "$preview_dir/${mode}.dag.dot" >/dev/null
}

assert_present() {
    local dag=$1
    local process=$2
    rg -q "label=\"${process}\"" "$dag"
}

assert_absent() {
    local dag=$1
    local process=$2
    if rg -q "label=\"${process}\"" "$dag"; then
        echo "Unexpected process ${process} in ${dag}" >&2
        exit 1
    fi
}

preview fastq test_data/samplesheet.csv
assert_present "$preview_dir/fastq.dag.dot" PREPARE_REFERENCE
assert_present "$preview_dir/fastq.dag.dot" PREPARE_SAMTOOLS_INDEX
assert_present "$preview_dir/fastq.dag.dot" BWA_MEM2
assert_absent "$preview_dir/fastq.dag.dot" PREPARE_DOWNSTREAM_REFERENCE

for mode in bam cram mt_only; do
    case "$mode" in
        bam) samplesheet=test_data/samplesheet_bam.csv ;;
        cram) samplesheet=test_data/samplesheet_cram.csv ;;
        mt_only) samplesheet=test_data/samplesheet_mt_only.csv ;;
    esac
    preview "$mode" "$samplesheet"
    assert_present "$preview_dir/${mode}.dag.dot" PREPARE_DOWNSTREAM_REFERENCE
    assert_absent "$preview_dir/${mode}.dag.dot" PREPARE_REFERENCE
    assert_absent "$preview_dir/${mode}.dag.dot" PREPARE_SAMTOOLS_INDEX
    assert_absent "$preview_dir/${mode}.dag.dot" BWA_MEM2
done

echo 'Input reference-routing checks passed.'
