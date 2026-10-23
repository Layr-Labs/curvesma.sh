#!/usr/bin/env bash
# Reproduce the per-step instruction counts from the repository root:
#   bash submissions/bls12-381-g1-ecdlp-plain/rho-sigma6-asm/measurement/reproduce.sh [seed]
# Needs Docker (linux/arm64) and git. Uses the organizer's reference-count image.
set -euo pipefail
seed=${1:-1}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
out=$(mktemp -d "${TMPDIR:-/tmp}/rho-sigma6-asm.XXXXXX")
git clone --quiet https://github.com/supranational/blst.git "$out/blst"
git -C "$out/blst" checkout --quiet f7a88067ed98385ebc6cf41ccecf21ead706a5b7
mkdir "$out/work"
python3 -I "$here/walk_ref.py" gen "$out/work/inputs.txt" "$seed"   # about 2 minutes
python3 -I "$here/frac_check.py" "$seed"
docker build --platform linux/arm64 -t curvesma-reference-count "$repo/scripts/reference_costs"
docker run --rm --platform linux/arm64 --network none \
  -v "$out/blst:/blst:ro" -v "$here:/harness:ro" -v "$out/work:/work" \
  curvesma-reference-count sh /harness/run.sh
python3 -I "$here/walk_ref.py" check "$out/work/inputs.txt" 16 "$out/work/dump-16.txt" "$out/work/native-16.stderr"
python3 -I "$here/summarize.py" "$out/work" --steps-per-n 16384 --C 2000 --n1 8 --n2 16 \
  --mode single:1 --mode-n1 2048 --mode-n2 4096
echo "Raw profiles retained at $out/work"
