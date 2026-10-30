#!/bin/sh
# Runs inside the organizers' reference-cost image (see reproduce_lanes.sh).
# Mounts: /blst (pinned blst), /harness (scripts/reference_costs), /m (this dir), /work (results).
set -eu
cc --version | head -1
valgrind --version
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/bls_field.c /work/blst.o /work/assembly.o -o /work/bls_field
sha256sum /work/bls_field /m/side_lanes.S /m/side_lanes_kernel.c /m/gen_side_lanes.py
cc -c /m/side_lanes.S -o /work/side_lanes.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /m/side_lanes_kernel.c /work/side_lanes.o /work/blst.o /work/assembly.o -o /work/side_lanes
objdump -d /work/side_lanes > /work/side_lanes.dis
for kind in bls_field side_lanes; do
  for n in 4096 8192; do
    valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$kind-$n.callgrind /work/$kind "$n" > /work/$kind-$n.stdout 2> /work/$kind-$n.stderr
    echo "== $kind $n"; cat /work/$kind-$n.stdout; grep -E "Collected|refs" /work/$kind-$n.stderr
  done
done
