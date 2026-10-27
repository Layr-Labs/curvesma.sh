#!/bin/sh
set -e
cd "$(dirname "$0")"
CC=${CC:-gcc}
$CC -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o slice12 slice12.c -lm
$CC -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o nbound12 nbound12.c -lm
$CC -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns -o ecm12 ecm12.c
$CC --version | head -1; valgrind --version
tot() { valgrind --tool=callgrind --callgrind-out-file=/tmp/cg.out "$@" > /tmp/stdout.txt 2>/dev/null
        echo "$(callgrind_annotate /tmp/cg.out 2>/dev/null | awk '/PROGRAM TOTALS/ {gsub(",","",$1); print $1}') $(cat /tmp/stdout.txt)"; }
field() { echo "$1" | tr ' ' '\n' | sed -n "s/^$2=//p"; }
diff2() {   # diff2 <label> <unit-field> <cmd with REPS placeholder> R1 R2
  label=$1; uf=$2; cmd=$3; r1=$4; r2=$5
  o1=$(tot $(echo "$cmd" | sed "s/REPS/$r1/")); o2=$(tot $(echo "$cmd" | sed "s/REPS/$r2/"))
  i1=${o1%% *}; i2=${o2%% *}; u1=$(field "$o1" $uf); u2=$(field "$o2" $uf)
  python3 -c "print('$label: Ir %d -> %d, $uf %s -> %s, Ir per unit %.2f' % ($i1, $i2, $u1, $u2, ($i2-$i1)/($u2-$u1)))"
  echo "$o2" | cut -d' ' -f2- | sed 's/^/    /'
}
LP=26
echo "== LINE path per update, setups included (h = points per line coset = cells / K / p)"
for spec in "257 77 5 2 4" "4099 999 7 4 8" "65537 4321 999 40 80" "1048583 777 31337 400 800"; do set -- $spec
  diff2 "line p=$1" hits "./slice12 bench line $LP 2 $1 $2 $3 REPS" $4 $5; done
echo "== LINE path setup per (pair, line, coset)"
diff2 "setup_line p=65537" units "./slice12 bench setupl 20 2 65537 777 5 REPS" 100000 300000
echo "== RECORD path per hit, per-coset recentring/centre included (h = 2 cells / p)"
for spec in "257 77 5 1 2" "65537 4321 999 20 40" "8388617 55555 1234567 400 800" "67108879 999983 4242 2000 4000" "268435459 999983 777 6000 12000"; do set -- $spec
  diff2 "record p=$1" hits "./slice12 bench rec $LP 2 $1 $2 $3 REPS" $4 $5; done
echo "== APPLY per record (record path + apply minus record path alone, same records)"
for spec in "65537 4321 999 40"; do set -- $spec
  oa=$(tot ./slice12 bench apply $LP 2 $1 $2 $3 $4); orr=$(tot ./slice12 bench rec $LP 2 $1 $2 $3 $4)
  python3 -c "print('apply p=$1: Ir per record %.2f' % ((${oa%% *}-${orr%% *})/$(field "$oa" units)))"; done
echo "== per-(slice, pair) setup6"
for p in 1048583 9223372036854775783; do diff2 "setup6 p=$p" units "./slice12 bench setup6 20 2 $p 999983 1 REPS" 100000 300000; done
echo "== per cell: zero pass (bench zero minus bench build: one extra zero pass)"
oz=$(tot ./slice12 bench zero $LP 2 1048583 777 5 1); ob=$(tot ./slice12 bench build $LP 2 1048583 777 5 1)
python3 -c "print('zero_cells: Ir per cell %.3f' % ((${oz%% *}-${ob%% *})/$(field "$ob" cells)))"
python3 -c "print('build (whole program incl. allocation, layout, first zero pass) Ir per suffix-box entry %.2f, per cell %.3f; rows %s entries %s' % (${ob%% *}/$(field "$ob" entries), ${ob%% *}/$(field "$ob" cells), '$(field "$ob" rows)', '$(field "$ob" entries)'))"
echo "== per cell: scan (bench scan minus bench apply on the same records), block length 16 and 4"
oa=$(tot ./slice12 bench apply $LP 2 65537 4321 1234567 20)
for bl in 16 4; do os=$(tot ./slice12 bench scan $LP 2 65537 4321 $bl 20)
  python3 -c "print('scan block=$bl: Ir per cell %.3f survivors %s' % ((${os%% *}-${oa%% *})/$(field "$os" units), '$(field "$os" survivors)'))"; done
echo "== per record: recover (bench recover minus bench scan)"
os=$(tot ./slice12 bench scan $LP 2 65537 4321 16 20); orc=$(tot ./slice12 bench recover $LP 2 65537 4321 16 20)
python3 -c "print('recover: Ir per record %.3f (records %s)' % ((${orc%% *}-${os%% *})/$(field "$orc" units), '$(field "$orc" units)'))"
echo "== norm-bound blocks (nbound12): per block evaluation (bench minus benchsetup), per-row setup (benchrow), per-a setup; 64/128 production blocks (explore)"
V=${VECTORS:-nbound-vectors.txt}
b1=$(tot ./nbound12 bench $V 1); b3=$(tot ./nbound12 bench $V 3); s1=$(tot ./nbound12 benchsetup $V 1); s3=$(tot ./nbound12 benchsetup $V 3)
r1=$(tot ./nbound12 benchrow $V 1); r5=$(tot ./nbound12 benchrow $V 5); e1=$(tot ./nbound12 explore $V 64 128 1); e3=$(tot ./nbound12 explore $V 64 128 3)
ev16=$(field "$b1" evaluations); ev64=$(./nbound12 explore $V 64 128 | sed 's/evaluations_total=\([0-9]*\).*/\1/')
python3 -c "
b=(${b3%% *}-${b1%% *})/2; s=(${s3%% *}-${s1%% *})/2; r=(${r5%% *}-${r1%% *})/4; e=(${e3%% *}-${e1%% *})/2
print('nbound: Ir per block evaluation %.1f (16-blocks) / %.1f (64/128 blocks); row_setup %.1f Ir per row; setup_a + row_setup %.1f Ir per case' % ((b-s)/$ev16, (e-s)/$ev64, r/2401, s/2401))"
echo "== ecm12 mont_mul (5 limbs) Ir per product"
m1=$(tot ./ecm12 bench 100000); m3=$(tot ./ecm12 bench 300000); python3 -c "print('mont_mul: %.2f Ir per product' % ((${m3%% *}-${m1%% *})/200000))"

