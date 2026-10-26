#!/usr/bin/env bash
set -euo pipefail
seed=${1:-1}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
out=$(mktemp -d "${TMPDIR:-/tmp}/rho-sigma2.XXXXXX")
python3 -I "$here/walk_ref.py" gen "$out/inputs.txt" "$seed"
docker build --platform linux/arm64 -t curvesma-reference-count "$repo/scripts/reference_costs"
docker run --rm --platform linux/arm64 --network none \
  -v "$here:/harness:ro" -v "$out:/work" \
  curvesma-reference-count sh /harness/run.sh
python3 -I "$here/walk_ref.py" check "$out/inputs.txt" 8 "$out/dump-8.txt" "$out/native-8.stderr"
python3 -I "$here/summarize.py" "$out" --steps-per-n 65536 --C 1000 --n1 4 --n2 8 \
  --mode single:1 --mode-n1 1024 --mode-n2 2048
echo "Raw profiles retained at $out"
