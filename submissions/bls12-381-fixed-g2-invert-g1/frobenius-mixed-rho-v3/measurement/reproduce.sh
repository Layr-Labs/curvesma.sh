#!/usr/bin/env bash
# Reproduce the blst step measurement at the production square count (S = 458), the interval
# analysis and the replica check. Usage: bash reproduce.sh <work-dir> [native]
#   Docker (linux/arm64), the organizer image (CURVESMASH, default ../curvesmash next to the core),
#   the core's measure/measure.sh (CORE, default ../curvesmash-core) and py_ecc (PY).
#   With the optional argument "native", also runs the uninstrumented kernel for 2^20 rounds to
#   recount cycle events and operations (about 10 minutes).
set -euo pipefail
work=${1:?work dir}
here=$(cd "$(dirname "$0")" && pwd)
core=${CORE:-$HOME/experimental/yukon/curvesmash-core}
py=${PY:-$HOME/experimental/yukon/work-gt/.venv/bin/python}
S_FRAC=0.447265625   # 458 / 1024: the production square count
for seed in 1 2; do
  dir="$work/seed$seed"
  mkdir -p "$dir"
  "$py" "$here/gen_inputs.py" --seed "$seed" --s "$S_FRAC" --out "$dir"
  LIB=blst N1=256 N2=512 MODES="mul sqr" MODE_N1=4096 MODE_N2=8192 \
    "$core/measure/measure.sh" "$here/step.c" blst "$dir/inputs.bin" "$dir/out"
  python3 -I "$here/analyze.py" "$dir/out" --n1 256 --n2 512 --s "$S_FRAC" --label "seed$seed"
  cp "$dir/out/dump-512.txt" "$dir/out/native-512.stderr" "$dir/"
  "$py" "$here/replica.py" "$dir" 512
  # c2 excludes distinguished-point work: the same inputs with the DP density set to zero (header DP_BITS = 60)
  mkdir -p "$dir/nodp"
  "$py" -c "
import struct
b = bytearray(open('$dir/inputs.bin', 'rb').read())
assert struct.unpack_from('<Q', b, 32)[0] == 4
struct.pack_into('<Q', b, 32, 60)
open('$dir/nodp/inputs.bin', 'wb').write(bytes(b))"
  LIB=blst N1=256 N2=512 MODES="mul sqr" MODE_N1=4096 MODE_N2=8192 \
    "$core/measure/measure.sh" "$here/step.c" blst "$dir/nodp/inputs.bin" "$dir/nodp/out"
  python3 -I "$here/analyze.py" "$dir/nodp/out" --n1 256 --n2 512 --s "$S_FRAC" --label "seed$seed-nodp"
  if [ "${2:-}" = "native" ]; then
    CORE="$core" bash "$here/native_run.sh" "$here/step.c" "$dir/inputs.bin" 1048576 "$dir/native"
  fi
done
