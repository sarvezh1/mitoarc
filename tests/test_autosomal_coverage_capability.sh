#!/usr/bin/env bash
set -euo pipefail

fixture_dir=$(mktemp -d)
trap 'rm -r -- "$fixture_dir"' EXIT

for n in $(seq 1 22); do
    printf 'chr%s\n' "$n"
done > "$fixture_dir/full-selected.txt"

for n in $(seq 1 22); do
    mapped=0
    if [ "$n" -eq 1 ] || [ "$n" -eq 2 ] || [ "$n" -eq 17 ]; then mapped=10; fi
    printf 'chr%s\t1000\t%s\t0\n' "$n" "$mapped"
done > "$fixture_dir/reduced.idxstats"

selected_count=$(awk 'NF { n++ } END { print n+0 }' "$fixture_dir/full-selected.txt")
mapped_selected_count=$(awk 'NR == FNR { selected[$1]=1; next } ($1 in selected) && $3 > 0 { n++ } END { print n+0 }' \
    "$fixture_dir/full-selected.txt" "$fixture_dir/reduced.idxstats")

if [ "$selected_count" -ne 22 ] || [ "$mapped_selected_count" -ne 3 ]; then
    echo 'reduced full-header fixture was not detected correctly' >&2
    exit 1
fi
if ! { [ "$selected_count" -ge 22 ] && [ "$mapped_selected_count" -lt "$selected_count" ]; }; then
    echo 'reduced full-header context was incorrectly considered coverage-capable' >&2
    exit 1
fi

printf 'chr1\nchr2\nchr17\n' > "$fixture_dir/explicit-subset.txt"
subset_count=$(awk 'NF { n++ } END { print n+0 }' "$fixture_dir/explicit-subset.txt")
subset_mapped=$(awk 'NR == FNR { selected[$1]=1; next } ($1 in selected) && $3 > 0 { n++ } END { print n+0 }' \
    "$fixture_dir/explicit-subset.txt" "$fixture_dir/reduced.idxstats")
if [ "$subset_count" -ne 3 ] || [ "$subset_mapped" -ne 3 ]; then
    echo 'explicit mapped subset fixture was not preserved' >&2
    exit 1
fi
if [ "$subset_count" -ge 22 ] && [ "$subset_mapped" -lt "$subset_count" ]; then
    echo 'explicit smaller-region workflow was incorrectly disabled' >&2
    exit 1
fi

echo 'Autosomal coverage capability guard checks passed.'
