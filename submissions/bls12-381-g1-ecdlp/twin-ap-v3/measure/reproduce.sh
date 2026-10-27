#!/usr/bin/env bash
# Host-side driver, modelled on scripts/reference_costs/reproduce.sh. Needs Docker with linux/arm64 support
# (native AArch64, or binfmt QEMU on x86-64; the counts are Callgrind's guest-instruction counts either way).
# Usage: bash reproduce.sh <challenge repository root>
set -euo pipefail
root=$(cd "$1" && pwd); here=$(cd "$(dirname "$0")" && pwd)
out=$(mktemp -d "${TMPDIR:-/tmp}/twin-ap-v3-measure.XXXXXX")
git clone --quiet https://github.com/supranational/blst.git "$out/blst"
git -C "$out/blst" checkout --quiet f7a88067ed98385ebc6cf41ccecf21ead706a5b7
mkdir "$out/work"
docker build --platform linux/arm64 -t curvesma-reference-count "$root/scripts/reference_costs"   # the organizers' image and tag
docker run --rm --platform linux/arm64 --network none -v "$out/blst:/blst:ro" \
  -v "$root/scripts/reference_costs:/harness:ro" -v "$here:/m:ro" -v "$out/work:/work" \
  curvesma-reference-count sh /m/run.sh | tee "$out/work/output.txt"
# the organizers' published bls_field totals (scripts/reference_costs/measurement.json) must reproduce exactly
grep -q "^/work/bls_field : Ir(4096)=1568788 Ir(8192)=3137556 per-iteration=383 remainder=0$" "$out/work/output.txt"
grep -q "^/work/side_output lanes : .* per-iteration=1329 remainder=0$" "$out/work/output.txt"
grep -q "^/work/side_fused_test : .* per-iteration=1093 remainder=0$" "$out/work/output.txt"
test "$(grep -c '^selfcheck ok' "$out/work/output.txt")" = 3
echo "reproduced: organizers' bls_field totals, the 1329-instruction blst-layer and the 1093-instruction fused side output; self-checks passed"
