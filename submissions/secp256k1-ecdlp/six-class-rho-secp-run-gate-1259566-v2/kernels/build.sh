#!/bin/sh
# Rebuild the shipped assembly (needs python3 and clang with the AArch64 target).
set -e
# counted loops and the emitter test wrappers: emitted by the generator (byte for byte the shipped files)
python3 gen_loops.py elf > secp_loops.s
python3 gen_loops.py elf tests > secp_tests.s
# both assemble as AArch64 (ADDS with a shifted register is a standard A64 form)
clang -target aarch64-linux-gnu -c secp_loops.s -o secp_loops.o
clang -target aarch64-linux-gnu -c secp_tests.s -o secp_tests.o
# secp_fe.s (fe_mul_f, fe_sqr_f, fe_inv) ships unchanged from b5de673a; it was produced
# there with strip_asm.py from the clang output of the upstream C reference secp_walk.c,
# which is not carried in this archive (256 KiB submission limit). experiment.py parses
# and executes the shipped assembly directly.
