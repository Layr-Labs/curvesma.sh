#!/bin/bash
# sparse-arm census: arm A dense (dp_bits=0) control, arm B sparse (dp_bits=4,6,8)
# each arm = 8 processes x NTR trials, different seeds; JSON lines -> out/
set -u
cd "$(dirname "$0")"
BIN=./census
NTR=${NTR:-500}
NPROC=8
mkdir -p out
for dp in 0 4 6 8; do
  for p in 0 1 2 3 4 5 6 7; do
    seed=$(( 771000 + dp*1000 + p ))
    $BIN $dp $NTR $seed 6 10 > out/arm${dp}_p${p}.json 2> out/arm${dp}_p${p}.err &
  done
done
wait
echo DONE > out/ALLDONE

# 8 processes x 500 trials per arm, dp_bits in {0 (arm A dense control), 4, 6, 8 (arm B sparse)}.
# Aggregation + Clopper-Pearson: experiments/make_cert_v10h.py, written to evidence/census_sparse.json.
