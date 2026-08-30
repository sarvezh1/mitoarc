#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"

samtools_image='quay.io/biocontainers/samtools:1.20--h50ea8bc_1'
failure_dir=/tmp/mitoarc-input-failure-fixtures
source_bam="$root/test_data/fixtures/TEST.sorted.bam"
reference="$root/test_data/fixture-reference.fa"
mkdir -p "$failure_dir"

samtools() {
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --volume "$root:$root" \
    --volume "$failure_dir:$failure_dir" \
    --workdir "$root" \
    "$samtools_image" samtools "$@"
}

expect_failure() {
  local name=$1 sheet=$2 pattern=$3 supplied_reference=${4:-$reference}
  local command=(nextflow run main.nf -profile docker --input "$sheet"
    --shift_offset 250 --outdir "/tmp/mitoarc-input-${name}-results"
    -work-dir "/tmp/mitoarc-input-${name}-work")
  if [ "$supplied_reference" != 'NONE' ]; then
    command+=(--reference "$supplied_reference")
  fi
  if "${command[@]}" 2>&1 | tee "/tmp/mitoarc-input-${name}.log"; then
    echo "expected failure did not occur: ${name}" >&2
    exit 1
  fi
  rg -q "$pattern" "/tmp/mitoarc-input-${name}.log"
  printf 'PASS\t%s\t%s\n' "$name" "$pattern"
}

printf 'sample,bam\nTEST,%s/missing.bam\n' "$failure_dir" > "$failure_dir/missing_bam.csv"
printf 'sample,cram\nTEST,%s/missing.cram\n' "$failure_dir" > "$failure_dir/missing_cram.csv"
printf 'not a BAM file\n' > "$failure_dir/corrupt.bam"
printf 'sample,bam\nTEST,%s/corrupt.bam\n' "$failure_dir" > "$failure_dir/corrupt_bam.csv"

samtools view -h "$source_bam" > "$failure_dir/unsorted.source.sam"
awk 'BEGIN { OFS="\t" } $1 == "@HD" { for (i=1; i<=NF; i++) if ($i ~ /^SO:/) $i="SO:unsorted" } { print }' \
  "$failure_dir/unsorted.source.sam" > "$failure_dir/unsorted.sam"
samtools view -b -o "$failure_dir/unsorted.bam" "$failure_dir/unsorted.sam"
samtools index "$failure_dir/unsorted.bam"
printf 'sample,bam\nTEST,%s/unsorted.bam\n' "$failure_dir" > "$failure_dir/unsorted_bam.csv"

samtools view -H "$source_bam" \
  | awk '$1 != "@SQ" || $2 == "SN:chr1"' > "$failure_dir/no_mt.header.sam"
samtools view -b -o "$failure_dir/no_mt.bam" "$failure_dir/no_mt.header.sam"
samtools index "$failure_dir/no_mt.bam"
printf 'sample,bam\nTEST,%s/no_mt.bam\n' "$failure_dir" > "$failure_dir/no_mt.csv"

samtools view -H "$source_bam" | awk '$1 !~ /^@RG/' > "$failure_dir/multiple_sm.header.sam"
printf '@RG\tID:RG_A\tSM:ALPHA\n@RG\tID:RG_B\tSM:BETA\n' >> "$failure_dir/multiple_sm.header.sam"
samtools reheader "$failure_dir/multiple_sm.header.sam" "$source_bam" > "$failure_dir/multiple_sm.bam"
samtools index "$failure_dir/multiple_sm.bam"
printf 'sample,bam\nTEST,%s/multiple_sm.bam\n' "$failure_dir" > "$failure_dir/multiple_sm.csv"

awk '/^>/ { print; next } !changed { first=substr($0,1,1); replacement=(first == "A" ? "C" : "A"); print replacement substr($0,2); changed=1; next } { print }' \
  "$reference" > "$failure_dir/incompatible-reference.fa"
printf 'sample,cram\nTEST,%s\n' "$root/test_data/fixtures/TEST.sorted.cram" > "$failure_dir/cram.csv"

expect_failure missing_bam "$failure_dir/missing_bam.csv" 'BAM file does not exist'
expect_failure missing_cram "$failure_dir/missing_cram.csv" 'CRAM file does not exist'
expect_failure corrupt_bam "$failure_dir/corrupt_bam.csv" 'quickcheck|corrupt\.bam|VALIDATE_PREALIGNED'
expect_failure unsorted_bam "$failure_dir/unsorted_bam.csv" 'must be coordinate-sorted.*unsorted'
expect_failure cram_without_reference "$failure_dir/cram.csv" 'Missing required parameter: --reference' NONE
expect_failure incompatible_cram_reference "$failure_dir/cram.csv" 'MD5|cannot be decoded|reference' "$failure_dir/incompatible-reference.fa"
expect_failure missing_mt_contig "$failure_dir/no_mt.csv" "Requested mitochondrial contig 'chrM' is absent"
expect_failure multiple_sm "$failure_dir/multiple_sm.csv" 'multiple biological SM values: ALPHA,BETA'
expect_failure mixed test_data/samplesheet_mixed_invalid.csv 'Invalid samplesheet schema'
expect_failure missing_samplesheet test_data/samplesheet_missing_bam.csv 'Samplesheet does not exist'

echo 'Pre-aligned input validation checks passed.'
