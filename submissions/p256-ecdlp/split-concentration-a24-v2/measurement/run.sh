#!/bin/sh
set -eu
cc --version | head -1
valgrind --version
cc -O3 -fno-tree-vectorize -mgeneral-regs-only -o /work/kernel /harness/step.c
for n in 4 8; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/step-$n.callgrind \
    /work/kernel /work/inputs.txt "$n" > /work/step-$n.stdout 2> /work/step-$n.stderr
  tail -1 /work/step-$n.stderr
  grep -h "^summary:" /work/step-$n.callgrind
done
for n in 1024 2048; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/mode-single-$n.callgrind \
    /work/kernel /work/inputs.txt "$n" single > /dev/null 2> /work/mode-single-$n.stderr
  grep -h "^summary:" /work/mode-single-$n.callgrind
done
/work/kernel /work/inputs.txt 8 dump > /work/dump-8.txt 2> /work/native-8.stderr
cat /work/native-8.stderr
