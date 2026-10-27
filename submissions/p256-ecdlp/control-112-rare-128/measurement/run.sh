#!/bin/sh
# Runs inside the organizer's reference-count image (scripts/reference_costs/Dockerfile)
# with this directory at /harness and outputs in /work. Organizer compiler flags.
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
# Anchor round: the same source with ANCHOR_PERIOD=2 (one AND immediate differs), so Ir(2)-Ir(1) is an anchor round.
cc -O3 -fno-tree-vectorize -mgeneral-regs-only -DANCHOR_PERIOD=2 -o /work/kernel-a2 /harness/step.c
for n in 1 2; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/anchor-$n.callgrind \
    /work/kernel-a2 /work/inputs.txt "$n" > /dev/null 2> /work/anchor-$n.stderr
  grep -h "^summary:" /work/anchor-$n.callgrind
done
for n in 1024 2048; do
  valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/mode-single-$n.callgrind \
    /work/kernel /work/inputs.txt "$n" single > /dev/null 2> /work/mode-single-$n.stderr
  grep -h "^summary:" /work/mode-single-$n.callgrind
done
/work/kernel /work/inputs.txt 8 dump > /work/dump-8.txt 2> /work/native-8.stderr
cat /work/native-8.stderr
