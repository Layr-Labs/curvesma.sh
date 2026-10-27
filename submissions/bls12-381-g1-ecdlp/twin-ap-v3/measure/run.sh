#!/bin/sh
# twin-ap-v2 measurement, run inside the organizers' image (scripts/reference_costs/Dockerfile) with the same flags as
# scripts/reference_costs/run.sh. Mounts: /blst = blst at f7a88067ed98385ebc6cf41ccecf21ead706a5b7 (read-only),
# /harness = scripts/reference_costs (read-only), /m = this directory (read-only), /work = output.
set -eu
cc --version | head -1
valgrind --version
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/bls_field.c /work/blst.o /work/assembly.o -o /work/bls_field
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /m/fp_ops.c /work/blst.o /work/assembly.o -o /work/fp_ops
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /m/side_output.c -o /work/side_output.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /work/side_output.o /work/assembly.o -o /work/side_output
sha256sum /work/bls_field /work/fp_ops /work/side_output
ir() { valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/cg.out "$@" 2>/work/cg.err >/work/cg.stdout; sed -n 's/.*I *refs: *//p' /work/cg.err | tr -d ,; }
row() { a=$(ir "$@" 4096); b=$(ir "$@" 8192); echo "$* : Ir(4096)=$a Ir(8192)=$b per-iteration=$(( (b - a) / 4096 )) remainder=$(( (b - a) % 4096 ))"; }
row /work/bls_field
for op in mul sqr add sub; do row /work/fp_ops $op; done
for mode in lanes call exported; do row /work/side_output $mode; done
# outputs of the side-output program (self-check line), and the expected organizer totals
/work/side_output lanes 8192 | head -1
/work/side_output exported 8192 | head -1
