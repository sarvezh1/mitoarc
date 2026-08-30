#!/usr/bin/env bash
set -euo pipefail

fixture_dir=$(mktemp -d)
trap 'rm -r -- "$fixture_dir"' EXIT

records="$fixture_dir/records.sam"
legacy="$fixture_dir/legacy.tsv"
streaming="$fixture_dir/streaming.tsv"

printf '%s\n' \
    $'readA\t0\tchrM\t1\t10\t1M\t*\t0\t0\tA\tI' \
    $'readB\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI' \
    $'readA\t0\tchrM\t2\t25\t1M\t*\t0\t0\tA\tI' \
    $'readC\t256\tchrM\t3\t60\t1M\t*\t0\t0\tA\tI' \
    $'readD\t2048\tchrM\t4\t60\t1M\t*\t0\t0\tA\tI' \
    $'readE\t0\tchrM\t5\t30\t1M\t*\t0\t0\tA\tI' \
    $'readF\t0\tchrM\t6\t20\t1M\t*\t0\t0\tA\tI' \
    $'readE\t0\tchrM\t7\t40\t1M\t*\t0\t0\tA\tI' \
    $'readG\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI' > "$records"

awk -F '\t' -v OFS='\t' -v mt='chrM' -v minq='20' '
  {
    name[$1]=1
    secondary=int($2/256)%2
    supplementary=int($2/2048)%2
    if (!secondary && !supplementary && $3 == mt && $5 >= minq) {
      keep[$1]=1
      if ($5 > maxq[$1]) maxq[$1]=$5
    }
  }
  END {
    for (n in name) {
      decision=(n in keep) ? "KEEP" : "FLAG"
      reason=(n in keep) ? "PRIMARY_CHRM_MAPQ_GE_THRESHOLD" : "NO_PRIMARY_CHRM_MAPQ_GE_THRESHOLD"
      print n, decision, reason, ((n in keep) ? maxq[n] : "NA")
    }
  }
' "$records" | sort -k1,1 > "$legacy"

awk -F '\t' -v OFS='\t' -v mt='chrM' -v minq='20' '
  {
    secondary=int($2/256)%2
    supplementary=int($2/2048)%2
    eligible=(!secondary && !supplementary && $3 == mt && $5 >= minq)
    print $1, (eligible ? 1 : 0), (eligible ? $5 : "NA")
  }
' "$records" | sort -S 1M -T "$fixture_dir" -k1,1 | awk -F '\t' -v OFS='\t' '
  function emit_decision() {
    if (read_name == "") return
    decision=keep ? "KEEP" : "FLAG"
    reason=keep ? "PRIMARY_CHRM_MAPQ_GE_THRESHOLD" : "NO_PRIMARY_CHRM_MAPQ_GE_THRESHOLD"
    print read_name, decision, reason, (keep ? maxq : "NA")
  }
  {
    if ($1 != read_name) {
      emit_decision()
      read_name=$1
      keep=0
      maxq=0
    }
    if ($2 == 1) {
      keep=1
      if ($3 > maxq) maxq=$3
    }
  }
  END { emit_decision() }
' > "$streaming"

cmp "$legacy" "$streaming"
echo 'NUMT mapping filter streaming equivalence checks passed.'
