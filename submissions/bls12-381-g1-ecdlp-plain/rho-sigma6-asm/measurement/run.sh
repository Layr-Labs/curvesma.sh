#!/bin/sh
# Runs inside the organizer's reference-count image (scripts/reference_costs/Dockerfile)
# with the pinned blst checkout at /blst, this directory at /harness and outputs in /work.
# blst (used only for the batch inversion and input conversion) is built exactly as in
# the organizer's harness; the kernel uses the organizer flags. Every instruction counts.
set -eu
cc --version | head -1
valgrind --version
FLAGS="-O3 -fno-tree-vectorize -mgeneral-regs-only"
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc $FLAGS -o /work/kernel /harness/step.c /work/blst.o /work/assembly.o
for n in 8 16; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/step-$n.callgrind \
    /work/kernel /work/inputs.txt "$n" > /work/step-$n.stdout 2> /work/step-$n.stderr
  tail -1 /work/step-$n.stderr
  grep -h "^summary:" /work/step-$n.callgrind
done
for n in 2048 4096; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/mode-single-$n.callgrind \
    /work/kernel /work/inputs.txt "$n" single > /dev/null 2> /work/mode-single-$n.stderr
  grep -h "^summary:" /work/mode-single-$n.callgrind
done
/work/kernel /work/inputs.txt 16 dump > /work/dump-16.txt 2> /work/native-16.stderr
cat /work/native-16.stderr
