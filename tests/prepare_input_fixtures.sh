#!/usr/bin/env bash
set -euo pipefail

# Uses the pinned samtools container; it never downloads data.
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
mkdir -p test_data/fixtures
source_bam=${1:-test_data/fixtures/TEST.sorted.bam}
reference=${2:-test_data/fixture-reference.fa}
samtools_image='quay.io/biocontainers/samtools:1.20--h50ea8bc_1'

samtools() {
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        --volume "$root:$root" \
        --workdir "$root" \
        "$samtools_image" samtools "$@"
}

samtools quickcheck "$source_bam"
samtools view -bh "$source_bam" chrM > test_data/fixtures/TEST.mt_only.full_header.bam
samtools view -H test_data/fixtures/TEST.mt_only.full_header.bam \
    | awk '$1 != "@SQ" || $2 == "SN:chrM"' \
    > test_data/fixtures/TEST.mt_only.header.sam
samtools reheader test_data/fixtures/TEST.mt_only.header.sam \
    test_data/fixtures/TEST.mt_only.full_header.bam \
    > test_data/fixtures/TEST.mt_only.bam
samtools index test_data/fixtures/TEST.mt_only.bam
samtools view -C -T "$reference" -o test_data/fixtures/TEST.sorted.cram "$source_bam"
samtools index test_data/fixtures/TEST.sorted.cram
rm -f test_data/fixtures/TEST.mt_only.full_header.bam test_data/fixtures/TEST.mt_only.header.sam

cat > test_data/samplesheet_mt_only.csv <<'EOF'
sample,bam
TEST,test_data/fixtures/TEST.mt_only.bam
EOF

echo "Created deterministic CRAM and mtDNA-only BAM fixtures under test_data/fixtures/"
