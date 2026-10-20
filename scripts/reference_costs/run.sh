#!/bin/sh
set -eu
cc --version | head -1
valgrind --version
cc -O3 -fno-tree-vectorize -mgeneral-regs-only -o /work/field /harness/field.c
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/fp12.c /work/blst.o /work/assembly.o -o /work/fp12
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/bls_field.c /work/blst.o /work/assembly.o -o /work/bls_field
for kind in field fp12 bls_field; do
  for n in 4096 8192; do
    valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$kind-$n.callgrind /work/$kind "$n" > /work/$kind-$n.stdout 2> /work/$kind-$n.stderr
    tail -8 /work/$kind-$n.stderr
  done
done
