#!/bin/sh
# Rebuild the shipped assembly and run the native tests (needs clang, python3, and a 64-bit AArch64 machine
# for the assembly tests: Linux or macOS on ARM).
set -e
# counted loops: emitted by the generator (byte for byte the shipped file)
python3 gen_loops.py elf > secp_loops.s
# counted inversion kernels: compiled
clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
      -ffreestanding -fno-builtin -fPIC -S -o secp_walk.s secp_walk.c
python3 strip_asm.py secp_walk.s .
# native test build
case "$(uname -s)" in Darwin) FLAVOR=macho ;; *) FLAVOR=elf ;; esac
python3 gen_loops.py $FLAVOR > loops_native.s
python3 gen_loops.py $FLAVOR tests > fe_tests_native.s
cc -O3 -shared -fPIC -o libsecpw.so secp_walk.c loops_native.s fe_tests_native.s
python3 test_native.py ./libsecpw.so
