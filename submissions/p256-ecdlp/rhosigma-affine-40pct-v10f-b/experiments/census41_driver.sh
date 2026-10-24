#!/bin/bash
# census41_driver.sh <dp_bits> <trials_total> <lam_num> <lam_den> <mult> <out.jsonl> [workers]
set -e
B=$1; T=$2; LN=$3; LD=$4; M=$5; OUT=$6; W=${7:-16}
PER=$(( (T + W - 1) / W ))
: > "$OUT"
for i in $(seq 0 $((W-1))); do
  SEED=$(( 0x41f000 + i * 7919 + B * 104729 + LN * 31 + M * 17 ))
  ( nice -n 10 ./census_v10f41 "$B" "$PER" "$SEED" "$LN" "$LD" "$M" >> "$OUT" ) &
done
wait
wc -l < "$OUT"
