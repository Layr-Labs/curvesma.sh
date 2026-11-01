#!/usr/bin/env bash
# Re-measure the lane-loop side-output kernel side_lanes with the organizers' own pipeline
# (scripts/reference_costs: pinned image, GCC 14.2.0, Valgrind Callgrind 3.24.0,
# blst f7a88067 with its AArch64 assembly, flags -O3 -fno-tree-vectorize
# -mgeneral-regs-only). Usage, from the repository root, Docker with linux/arm64:
#   bash submissions/bls12-381-g1-ecdlp/<candidate>/measure/reproduce.sh .
# Prints the Collected totals; the price is (Ir(8192) - Ir(4096)) / 4096.
set -euo pipefail
root=$(cd "${1:-.}" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
out=$(mktemp -d "${TMPDIR:-/tmp}/side-kernel.XXXXXX")
git clone --quiet https://github.com/supranational/blst.git "$out/blst"
git -C "$out/blst" checkout --quiet f7a88067ed98385ebc6cf41ccecf21ead706a5b7
python3 "$here/gen_side_lanes.py" "$out/blst" | cmp - "$here/side_lanes.S"   # the shipped .S is the generated one
mkdir "$out/results"
docker build --platform linux/arm64 -t curvesma-reference-count "$root/scripts/reference_costs"
docker run --rm --platform linux/arm64 --network none \
  -v "$out/blst:/blst:ro" -v "$root/scripts/reference_costs:/harness:ro" \
  -v "$here:/m:ro" -v "$out/results:/work" \
  curvesma-reference-count sh /m/run_measure_lanes.sh | tee "$out/results/measure.log"
python3 - "$out/results" <<'PY'
import re, sys, pathlib
d = pathlib.Path(sys.argv[1]); c = {}
for k in ("bls_field", "side_lanes"):
    c[k] = [int(re.search(r"Collected : (\d+)", (d / f"{k}-{n}.stderr").read_text()).group(1)) for n in (4096, 8192)]
    print(k, c[k], "per iteration:", (c[k][1] - c[k][0]) / 4096)
assert (c["bls_field"][1] - c["bls_field"][0]) / 4096 == 383
assert all("selfcheck ok" in (d / f"side_lanes-{n}.stdout").read_text() for n in (4096, 8192))
print("OK: results retained at", d)
PY
