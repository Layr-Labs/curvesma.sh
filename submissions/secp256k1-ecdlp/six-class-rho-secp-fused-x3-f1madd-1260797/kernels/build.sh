#!/bin/sh
# Rebuild the shipped assembly (needs python3 and clang with the AArch64 target).
set -e
# counted loops and the emitter test wrappers: emitted by the generator (byte for byte the shipped files)
python3 gen_loops.py elf > secp_loops.s
python3 gen_loops.py elf tests > secp_tests.s
# both assemble as AArch64 (ADDS with a shifted register is a standard A64 form)
clang -target aarch64-linux-gnu -c secp_loops.s -o secp_loops.o
clang -target aarch64-linux-gnu -c secp_tests.s -o secp_tests.o
# inversion kernels: compiled from the C reference (unchanged from the upstream package)
clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
      -ffreestanding -fno-builtin -fPIC -S -o secp_walk.s secp_walk.c
python3 strip_asm.py secp_walk.s .
