#!/usr/bin/env bash
set -euo pipefail

check_compatibility() {
    local reference_fai=$1
    local idxstats=$2
    awk 'BEGIN { bad=0 } FNR==NR { ref[$1]=$2; next } $1 != "*" && $3 > 0 { if (!($1 in ref)) { bad=1 } else if (ref[$1] != $2) { bad=1 } } END { exit bad }' "$reference_fai" "$idxstats"
}

fixture_dir=$(mktemp -d)
trap 'rm -r -- "$fixture_dir"' EXIT

printf 'chr1\t100\t0\t100\t101\nchrM\t16569\t0\t100\t101\n' > "$fixture_dir/reference.fai"

# An unused header-only decoy must not make an otherwise compatible BAM fail.
printf 'chr1\t100\t3\t0\nchrM\t16569\t5\t0\ndecoy1\t50\t0\t0\n*\t0\t0\t0\n' > "$fixture_dir/unused-decoy.idxstats"
check_compatibility "$fixture_dir/reference.fai" "$fixture_dir/unused-decoy.idxstats"

# A mapped record on a missing contig must still fail.
printf 'chr1\t100\t3\t0\nchrM\t16569\t5\t0\ndecoy1\t50\t1\t0\n*\t0\t0\t0\n' > "$fixture_dir/used-decoy.idxstats"
if check_compatibility "$fixture_dir/reference.fai" "$fixture_dir/used-decoy.idxstats"; then
    echo 'mapped missing contig was incorrectly accepted' >&2
    exit 1
fi

# A mapped contig with a conflicting length must still fail.
printf 'chr1\t101\t3\t0\nchrM\t16569\t5\t0\n*\t0\t0\t0\n' > "$fixture_dir/wrong-length.idxstats"
if check_compatibility "$fixture_dir/reference.fai" "$fixture_dir/wrong-length.idxstats"; then
    echo 'mapped contig length mismatch was incorrectly accepted' >&2
    exit 1
fi

echo 'Pre-aligned mapped-contig reference compatibility checks passed.'
