#!/bin/sh
set -e
python3 gen_loops.py elf > secp_loops.s
python3 gen_loops.py elf tests > secp_tests.s
clang -target aarch64-linux-gnu -c secp_loops.s -o secp_loops.o
clang -target aarch64-linux-gnu -c secp_tests.s -o secp_tests.o
clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
      -ffreestanding -fno-builtin -fPIC -S -o secp_walk.s secp_walk.c
python3 strip_asm.py secp_walk.s .
