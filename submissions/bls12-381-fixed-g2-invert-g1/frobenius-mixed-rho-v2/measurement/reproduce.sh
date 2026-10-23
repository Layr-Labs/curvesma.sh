#!/usr/bin/env bash
# Reproduce the blst step measurement (organizer image via measure.sh) and the replica check.
# Usage: bash reproduce.sh <work-dir>
# Needs: Docker (linux/arm64), py_ecc (work-gt/.venv), and the organizer clone at $CURVESMASH
#        (default ../curvesmash next to curvesmash-core). Run from the package directory.
set -euo pipefail
work=${1:?work dir}
here=$(cd "$(dirname "$0")" && pwd)
core=${CORE:-$HOME/experimental/yukon/curvesmash-core}
py=${PY:-$HOME/experimental/yukon/work-gt/.venv/bin/python}
for label in "s46-seed1 0.46 1" "s46-seed2 0.46 2" "s00-seed1 0.0 1"; do
  set -- $label
  dir="$work/$1"
  mkdir -p "$dir"
  "$py" "$here/gen_inputs.py" --seed "$3" --s "$2" --out "$dir"
  LIB=blst N1=64 N2=128 MODES="mul sqr" MODE_N1=4096 MODE_N2=8192 \
    "$core/measure/measure.sh" "$core/examples/blst-gt-rho/step.c" blst "$dir/inputs.bin" "$dir/out"
  python3 -I "$here/analyze.py" "$dir/out" --s "$2" --label "$1"
  cp "$dir/out/dump-128.txt" "$dir/out/native-128.stderr" "$dir/"
  "$py" "$here/replica.py" "$dir" 128
done
