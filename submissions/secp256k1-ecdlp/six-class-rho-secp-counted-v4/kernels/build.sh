#!/bin/sh
# Rebuild the counted assembly (needs clang) and run the native test (needs a C compiler on 64-bit AArch64 or x86-64).
set -e
clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
      -ffreestanding -fno-builtin -fPIC -S -o secp_walk.s secp_walk.c
python3 strip_asm.py secp_walk.s .
cc -O3 -fno-tree-vectorize -mgeneral-regs-only -ffreestanding -fno-builtin -fPIC -shared -o libsecpw.so secp_walk.c 2>/dev/null \
  || cc -O3 -shared -fPIC -o libsecpw.so secp_walk.c
python3 test_native.py ./libsecpw.so
