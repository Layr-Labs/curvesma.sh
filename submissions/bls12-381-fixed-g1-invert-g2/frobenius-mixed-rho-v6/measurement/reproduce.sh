#!/usr/bin/env bash
# Reproduce the v11 step measurement: per-step counts, restart, PRF, DP record, REF, dump/replica at DP 2 and 60.
# Usage: bash reproduce.sh <work-dir> [native]   (native: 2^20 rounds per seed). Needs CURVESMASH, CORE, PY.
set -euo pipefail
work=${1:?work dir}
here=$(cd "$(dirname "$0")" && pwd)
core=${CORE:-$HOME/experimental/yukon/curvesmash-core}
py=${PY:-$HOME/experimental/yukon/work-gt/.venv/bin/python}
S_FRAC=0.447265625
measure() { # <dir> <modes>: callgrind counts of the walk (N1 = 256, N2 = 512) and the mode loops
  LIB=blst N1=256 N2=512 MODES="$2" MODE_N1=4096 MODE_N2=8192 "$core/measure/measure.sh" "$here/step.c" blst "$1/inputs.bin" "$1/out"
}
analyze() { # <dir> <s> <label>
  python3 -I "$here/analyze.py" "$1/out" --n1 256 --n2 512 --s "$2" --label "$3"
}
for seed in 1 2; do
  dir="$work/seed$seed"
  mkdir -p "$dir"
  "$py" -I "$here/gen_inputs.py" --seed "$seed" --s "$S_FRAC" --dp-bits 60 --out "$dir/prod"
  measure "$dir/prod" "mul sqr restartbench storebench prfbench"
  analyze "$dir/prod" "$S_FRAC" "prod-seed$seed"
  for kind in add sqr; do
    s=0
    [ "$kind" = sqr ] && s=1
    "$py" -I "$here/gen_inputs.py" --seed "$seed" --s "$s" --dp-bits 60 --out "$dir/forced-$kind"
    measure "$dir/forced-$kind" "mul sqr"
    analyze "$dir/forced-$kind" "$s" "forced-$kind-seed$seed"
  done
  LIB=blst N1=256 N2=512 REF_FUNCS=ref_mul MODES="" "$core/measure/measure.sh" "$here/step.c" blst "$dir/prod/inputs.bin" "$dir/ref"
  for dp in 2 60; do
    "$py" -I "$here/gen_inputs.py" --seed "$seed" --s "$S_FRAC" --dp-bits "$dp" --out "$dir/dp$dp"
    CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$dir/dp$dp/inputs.bin" 512 "$dir/dp$dp/nat" dump
    cp "$dir/dp$dp/nat/native-512.stdout" "$dir/dp$dp/dump-512.txt"
    cp "$dir/dp$dp/nat/native-512.stderr" "$dir/dp$dp/native-512.stderr"
    "$py" -I "$here/replica.py" "$dir/dp$dp" 512
  done
  if [ "${2:-}" = "native" ]; then
    CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$dir/prod/inputs.bin" 1048576 "$dir/native"
  fi
done
CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$work/seed1/prod/inputs.bin" 0 "$work/vectors" vectors
"$py" -I "$here/check_vectors.py" "$work/vectors/native-0.stdout"
