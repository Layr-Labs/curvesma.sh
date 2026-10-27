#!/bin/sh
# Independent measurement of the three cost-model routines; run via ../README.md
# (docker refcount-harness, /blst at f7a88067; /work receives objects and records).
set -eu
echo "toolchain: $(cc --version | head -1) | $(valgrind --version) | debian $(cat /etc/debian_version)"
# COST_MODEL flags for blst's C; upstream AArch64 assembly as shipped.
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
for k in cyc2 fpmul; do
  cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/$k.c /work/blst.o /work/assembly.o -o /work/$k
done
sha256sum /work/blst.o /work/assembly.o /work/cyc2 /work/fpmul
echo "--- cyc/mul via cyc2 (kind 0 = cyclotomic_sqr, kind 1 = fp12_mul) ---"
for kind in 0 1; do
  for off in 0 3 101 1009 4001; do
    for n in 1 65; do
      valgrind --tool=callgrind --collect-atstart=no \
        --callgrind-out-file=/work/v2-$kind-$off-$n.callgrind \
        /work/cyc2 "$n" "$off" "$kind" > /work/v2.stdout 2> /work/v2.stderr
      echo "kind $kind off $off n $n Ir $(grep -h totals: /work/v2-$kind-$off-$n.callgrind | awk '{print $2}') $(cat /work/v2.stdout)"
    done
  done
done
echo "--- blst_fp_mul ---"
for off in 0 3 101 1009 4001; do
  for n in 1 65; do
    valgrind --tool=callgrind --collect-atstart=no \
      --callgrind-out-file=/work/v2-fp-$off-$n.callgrind \
      /work/fpmul "$n" "$off" > /work/v2.stdout 2> /work/v2.stderr
    echo "fp off $off n $n Ir $(grep -h totals: /work/v2-fp-$off-$n.callgrind | awk '{print $2}') $(cat /work/v2.stdout)"
  done
done
objdump -d --no-show-raw-insn /work/cyc2 > /work/cyc2.dis
objdump -d --no-show-raw-insn /work/fpmul > /work/fpmul.dis
echo "--- done ---"
