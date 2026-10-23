#!/usr/bin/env bash
# Reproduce the blst step measurement (v10: each trail ends at its first distinguished point; the walk restarts
# with the restart element Y_j = g^(f_j) of the restart table, j = mix((w + 1) GAMMA + t) mod TR, t the walk's restart
# counter; a restarted state that is itself a distinguished point is stored and restarted again, chained restarts):
# production S = 458 at the header DP_BITS = 60 (walk-spec section 6), the forced per-step runs (all-adder and
# all-square at DP_BITS = 60), the restart and distinguished-point record costs (modes restartbench and storebench),
# the replica checks (production DP_BITS = 60; restart check at DP_BITS = 4) and the sortnet table.
# Usage: bash reproduce.sh <work-dir> [native]
#   Docker (linux/arm64), the organizer image (CURVESMASH, default ../curvesmash next to the core),
#   the core's measure/measure.sh (CORE, default ../curvesmash-core) and py_ecc (PY).
#   With the optional argument "native", also runs the uninstrumented kernel for 2^20 rounds at
#   DP_BITS = 60 to recount cycle events and operations (about 10 minutes per seed).
set -euo pipefail
work=${1:?work dir}
here=$(cd "$(dirname "$0")" && pwd)
core=${CORE:-$HOME/experimental/yukon/curvesmash-core}
py=${PY:-$HOME/experimental/yukon/work-gt/.venv/bin/python}
S_FRAC=0.447265625   # 458 / 1024: the production square count
measure() { # <dir>: walk runs at N1/N2 and the mode loops on <dir>/inputs.bin
  LIB=blst N1=256 N2=512 MODES="mul sqr restartbench storebench" MODE_N1=4096 MODE_N2=8192 \
    "$core/measure/measure.sh" "$here/step.c" blst "$1/inputs.bin" "$1/out"
}
analyze() { # <dir> <s> <label>
  python3 -I "$here/analyze.py" "$1/out" --n1 256 --n2 512 --s "$2" --label "$3"
}
for seed in 1 2; do
  dir="$work/seed$seed"
  mkdir -p "$dir"
  # production: header DP_BITS = 60 (no distinguished points in measurement runs)
  "$py" "$here/gen_inputs.py" --seed "$seed" --s "$S_FRAC" --dp-bits 60 --out "$dir"
  measure "$dir"
  analyze "$dir" "$S_FRAC" "seed$seed"
  cp "$dir/out/dump-512.txt" "$dir/out/native-512.stderr" "$dir/"
  "$py" "$here/replica.py" "$dir" 512
  # forced per-step runs at DP_BITS = 60: all-adder (S = 0) and all-square (S = 1024)
  for kind in add sqr; do
    case $kind in
      add) sfrac=0 ;;
      sqr) sfrac=1 ;;
    esac
    fdir="$dir/forced-$kind"
    mkdir -p "$fdir"
    "$py" "$here/gen_inputs.py" --seed "$seed" --s "$sfrac" --dp-bits 60 --out "$fdir"
    measure "$fdir"
    analyze "$fdir" "$sfrac" "$kind-seed$seed"
  done
  # restart check at DP_BITS = 4 (restarts and chained restarts exercised): kernel dump and replica
  rdir="$dir/dp4"
  mkdir -p "$rdir"
  "$py" "$here/gen_inputs.py" --seed "$seed" --s "$S_FRAC" --dp-bits 4 --out "$rdir"
  CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$rdir/inputs.bin" 512 "$rdir/nat" dump
  cp "$rdir/nat/native-512.stdout" "$rdir/dump-512.txt"
  cp "$rdir/nat/native-512.stderr" "$rdir/native-512.stderr"
  "$py" "$here/replica.py" "$rdir" 512
  # sortnet table (walk-spec section 9.2): the standalone bitonic bench, not step.c
  LIB=blst N1=64 N2=128 MODES="sortnet sortscan" MODE_N1=4 MODE_N2=8 \
    "$core/measure/measure.sh" "$here/sortnet/sortnet.c" blst "$dir/inputs.bin" "$dir/sortnet"
  if [ "${2:-}" = "native" ]; then
    CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$dir/inputs.bin" 1048576 "$dir/native"
  fi
done
python3 -I "$here/sortnet/check_bitonic.py" > "$here/sortnet/check_bitonic.out"
