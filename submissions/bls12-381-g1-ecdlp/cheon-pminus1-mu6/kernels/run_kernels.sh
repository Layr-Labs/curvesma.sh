#!/bin/sh
set -eu
cc --version | head -1
cc -O3 -fno-tree-vectorize -mgeneral-regs-only -o /work/secp_step /k/secp_step.c
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc -O3 -fno-tree-vectorize -mgeneral-regs-only /k/bls_step.c /work/blst.o /work/assembly.o -o /work/bls_step
run() { # name binary args...
  name=$1; shift; bin=$1; shift
  for n in 4096 8192; do
    valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$name-$n.cg /work/$bin $n "$@" > /work/$name-$n.out 2>/dev/null
  done
  a=$(grep '^summary:' /work/$name-4096.cg | cut -d' ' -f2); b=$(grep '^summary:' /work/$name-8192.cg | cut -d' ' -f2)
  echo "$name $a $b slope=$(( (b-a)/4096 )) rem=$(( (b-a)%4096 )) out=$(cat /work/$name-8192.out)"
}
run secp_m0 secp_step 0
run secp_m1 secp_step 1
run secp_m2 secp_step 2
run bls_rho_m0 bls_step 0 0
run bls_rho_m1 bls_step 1 0
run bls_rho_m2 bls_step 2 0
run cheon_m0 bls_step 0 1
run cheon_m1 bls_step 1 1
run cheon_m2 bls_step 2 1
run cheon_m3 bls_step 3 1
