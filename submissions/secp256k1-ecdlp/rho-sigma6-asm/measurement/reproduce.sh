#!/usr/bin/env bash
# Reproduce the per-step instruction counts from the repository root:
#   bash submissions/secp256k1-ecdlp/rho-sigma6-asm/measurement/reproduce.sh [seed]
# Needs Docker (linux/arm64). Uses the organizer's reference-count image. No library
# checkout is needed: all field arithmetic is in secp_field.h and secp_sub.h.
set -euo pipefail
seed=${1:-1}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
out=$(mktemp -d "${TMPDIR:-/tmp}/rho-sigma6-asm.XXXXXX")
python3 -I "$here/walk_ref.py" gen "$out/inputs.txt" "$seed"   # a few minutes
python3 -I "$here/frac_check.py" "$seed"
docker build --platform linux/arm64 -t curvesma-reference-count "$repo/scripts/reference_costs"
docker run --rm --platform linux/arm64 --network none \
  -v "$here:/harness:ro" -v "$out:/work" \
  curvesma-reference-count sh /harness/run.sh
python3 -I "$here/walk_ref.py" check "$out/inputs.txt" 8 "$out/dump-8.txt"
python3 -I "$here/summarize.py" "$out" --steps-per-n 65536 --C 1000 --n1 4 --n2 8 \
  --mode single:1 --mode-n1 1024 --mode-n2 2048
echo "Raw profiles retained at $out"
