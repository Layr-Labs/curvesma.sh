#!/bin/bash
# slice8 retired-instruction measurements (AArch64, /usr/bin/time -l): differences between runs with 3x and 1x the
# repeated work (reps chosen so the repeated work is >= ~1e8 instructions), three trials, maximum reported.
cd "$(dirname "$0")"
icount() { /usr/bin/time -l "$@" 2>&1 | awk '/instructions retired/ {print $1}'; }
mx() { python3 -c "import sys; print(max(float(x) for x in sys.argv[1:]))" "$@"; }
hits() { ./slice8 $1 $2 $3 $4 $5 $6 1 2>/dev/null | sed 's/.*hits_per_rep=\([0-9]*\).*/\1/'; }
echo "== LINE path (p < P_MID): instructions per bucket update by expected points per line coset (slice 2^28 points)"
for spec in "28 0 257 77 5" "28 0 4099 999 7" "28 0 65537 4321 999" "28 0 1048583 777 31337" "28 0 8388617 55555 1234567" "28 0 67108879 999983 4242"; do
  set -- $spec; h=$(hits icount5 $1 $2 $3 $4 $5); reps=$(python3 -c "import math; print(max(4, math.ceil(3e6/max($h,1))))")
  vals=""; for t in 1 2 3; do a=$(icount ./slice8 icount5 $1 $2 $3 $4 $5 $reps); b=$(icount ./slice8 icount5 $1 $2 $3 $4 $5 $((3*reps))); vals="$vals $(python3 -c "print(round(($b-$a)/(2*$reps)/max($h,1),2))")"; done
  lines=$(./slice8 icount5 $1 $2 $3 $4 $5 1 2>&1 >/dev/null | sed 's/.*lines=\([0-9]*\).*/\1/')
  echo "p=$3 hits/rep=$h lines=$lines reps=$reps instr/update:$vals -> max $(mx $vals)"
done
echo "== RECORD path (p >= P_MID): instructions per hit, walk+record (icount6) and walk+record+apply (icount6a)"
for spec in "28 0 257 77 5" "28 0 65537 4321 999" "28 0 8388617 55555 1234567" "28 0 67108879 999983 4242" "28 0 536870923 999983 777" "28 0 2147483659 999983 777"; do
  set -- $spec; h=$(hits icount6 $1 $2 $3 $4 $5); reps=$(python3 -c "import math; print(max(4, math.ceil(3e6/max($h,1))))")
  for m in icount6 icount6a; do
    vals=""; for t in 1 2 3; do a=$(icount ./slice8 $m $1 $2 $3 $4 $5 $reps); b=$(icount ./slice8 $m $1 $2 $3 $4 $5 $((3*reps))); vals="$vals $(python3 -c "print(round(($b-$a)/(2*$reps)/max($h,1),2))")"; done
    echo "p=$3 $m hits/rep=$h reps=$reps instr/hit:$vals -> max $(mx $vals)"
  done
done
echo "== LINE path setup per (pair, line, coset): c_k, exact recentring, walk centre (1e6 vs 3e6 trials)"
vals=""; for t in 1 2 3; do a=$(icount ./slice8 setup 134217757 999983 1000000); b=$(icount ./slice8 setup 134217757 999983 3000000); vals="$vals $(python3 -c "print(round(($b-$a)/2000000,1))")"; done
echo "setup instr trials:$vals -> max $(mx $vals)"
echo "== RECORD path setup per (slice, pair): a(s) mod p (Horner, Barrett), c, exact recentring, walk centre (1e6 vs 3e6 trials)"
vals=""; for t in 1 2 3; do a=$(icount ./slice8 setup6 134217757 999983 1000000); b=$(icount ./slice8 setup6 134217757 999983 3000000); vals="$vals $(python3 -c "print(round(($b-$a)/2000000,1))")"; done
echo "setup6 instr trials:$vals -> max $(mx $vals)"
