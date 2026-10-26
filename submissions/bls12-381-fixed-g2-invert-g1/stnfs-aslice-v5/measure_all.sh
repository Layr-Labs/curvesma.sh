#!/bin/bash
# Retired-instruction measurements of slice6 (AArch64, /usr/bin/time -l).  Every figure is the difference between two
# runs of the same binary with different amounts of repeated work, divided by the extra work; the repeated work is
# chosen >= ~1e8 instructions so fixed costs (allocation, lattice reduction) and run-to-run noise are negligible.
# Three trials; the MAXIMUM is reported and charged.
cd "$(dirname "$0")"
icount() { /usr/bin/time -l "$@" 2>&1 | awk '/instructions retired/ {print $1}'; }
mx() { python3 -c "import sys; print(max(float(x) for x in sys.argv[1:]))" "$@"; }
echo "== per bucket update, by expected points per coset h (slice 2^26 points; reps scaled to >= 4e6 updates per run)"
for spec in "26 257 77 5" "26 2053 999 7" "26 8209 4321 999" "26 65537 4321 999" "26 524309 33333 123456"; do
  set -- $spec
  h=$(./slice6 icount $1 $2 $3 $4 1 2>/dev/null | sed 's/.*hits_per_rep=\([0-9]*\).*/\1/')
  reps=$(python3 -c "import math; print(max(16, math.ceil(4e6/$h)))")
  vals=""
  for t in 1 2 3; do a=$(icount ./slice6 icount $1 $2 $3 $4 $reps); b=$(icount ./slice6 icount $1 $2 $3 $4 $((3*reps))); vals="$vals $(python3 -c "print(round(($b-$a)/(2*$reps)/$h,2))")"; done
  echo "p=$2 h=$h reps=$reps instr/update trials:$vals -> max $(mx $vals)"
done
echo "== per coset in the sparse regimes (slice 2^27 points; 20000 vs 60000, 100000 cosets)"
for p in 8388617 33554467 134217757 536870923 2147483659; do
  h=$(./slice6 icountc 26 $p 999983 20000 2>/dev/null | sed 's/.*hits\/coset=\([0-9.]*\).*/\1/')
  vals=""
  for t in 1 2 3; do a=$(icount ./slice6 icountc 26 $p 999983 20000); b=$(icount ./slice6 icountc 26 $p 999983 100000); vals="$vals $(python3 -c "print(round(($b-$a)/80000,1))")"; done
  echo "p=$p h=$h instr/coset trials:$vals -> max $(mx $vals)"
done
echo "== per-(slice, pair) setup: Horner a(s) mod p (Barrett), c = a(s) r^-1, coset coordinates + Gram test (1e6 vs 3e6 trials)"
vals=""; for t in 1 2 3; do a=$(icount ./slice6 setup 134217757 999983 1000000); b=$(icount ./slice6 setup 134217757 999983 3000000); vals="$vals $(python3 -c "print(round(($b-$a)/2000000,1))")"; done
echo "setup instr trials:$vals -> max $(mx $vals)"
