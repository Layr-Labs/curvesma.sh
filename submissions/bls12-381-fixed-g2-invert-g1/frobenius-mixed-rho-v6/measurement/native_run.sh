#!/usr/bin/env bash
# Uninstrumented run of step.c in the organizer image: native_run.sh <kernel.c> <inputs.bin> <rounds> <out-dir> [mode]
# EXTRA_CFLAGS adds compile flags (test builds only, e.g. -DCAP_STEPS_OVERRIDE=1024). Needs <core>/.cache/blst
# (created by measure.sh) and the organizer Dockerfile at $CURVESMASH. Outside the core tree set CORE to curvesmash-core.
set -euo pipefail
kernel=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
inputs=$2
rounds=$3
out=$4
mode=${5:-}
extra=${EXTRA_CFLAGS:-}
core=${CORE:-$(cd "$(dirname "$0")/../.." && pwd)}
repo=${CURVESMASH:-$core/../curvesmash}
mkdir -p "$out"
cp "$inputs" "$out/inputs.bin"
docker build -q --platform linux/arm64 -t curvesma-reference-count "$repo/scripts/reference_costs" > /dev/null
docker run --rm --platform linux/arm64 --network none \
  -e KERNEL="$(basename "$kernel")" -e ROUNDS="$rounds" -e MODE="$mode" -e EXTRA="$extra" \
  -v "$core/.cache/blst:/blst:ro" -v "$(dirname "$kernel"):/harness:ro" \
  -v "$(cd "$out" && pwd):/work" \
  curvesma-reference-count sh -c '
    set -eu
    cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
    cc -c /blst/build/assembly.S -o /work/assembly.o
    cc -O3 -fno-tree-vectorize -mgeneral-regs-only $EXTRA -DREF_ONLY -o /work/kernel_native "/harness/$KERNEL" /work/blst.o /work/assembly.o
    /work/kernel_native /work/inputs.bin "$ROUNDS" $MODE > /work/native-$ROUNDS.stdout 2> /work/native-$ROUNDS.stderr
    tail -1 /work/native-$ROUNDS.stderr'
