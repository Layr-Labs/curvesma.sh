#!/usr/bin/env bash
set -euo pipefail
harness=$(cd "$(dirname "$0")" && pwd)
output=$(mktemp -d "${TMPDIR:-/tmp}/curvesma-reference.XXXXXX")
git clone --quiet https://github.com/bitcoin-core/secp256k1.git "$output/secp"
git -C "$output/secp" checkout --quiet e74a0187aac1cac734d5ad23ca5e40ceb6ddd4a8
git clone --quiet https://github.com/supranational/blst.git "$output/blst"
git -C "$output/blst" checkout --quiet f7a88067ed98385ebc6cf41ccecf21ead706a5b7
mkdir "$output/results"
docker build --platform linux/arm64 -t curvesma-reference-count "$harness"
docker run --rm --platform linux/arm64 --network none \
  -v "$output/secp:/secp:ro" -v "$output/blst:/blst:ro" \
  -v "$harness:/harness:ro" -v "$output/results:/work" \
  curvesma-reference-count sh /harness/run.sh
python3 "$harness/../reference_operation_costs.py" --results "$output/results"
echo "Raw profiles and pinned checkouts retained at $output"
