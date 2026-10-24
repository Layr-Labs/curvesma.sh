#!/bin/bash
# v19 evidence driver: audits first, then the production-shaped sweep + the floor census.
set -u
cd "$(dirname "$0")/.."
E=evidence/shards
mkdir -p "$E"
PY="nice -n 10 python3"

# ---- audits (gate the sweep spend only for the SHIP decision, not for CPU scheduling)
( echo '{"mode":"coalesce","pairs":1600,"steps":100,"rules":["v19","tn"]}' | $PY scripts/toy19.py > evidence/coalesce19.json 2> evidence/coalesce19.err ) &
( echo '{"mode":"trap","trials":1600,"steps":200,"tables":8,"rules":["v18","v19"]}' | $PY scripts/sweep19.py > evidence/trap19.json 2> evidence/trap19.err ) &
( echo '{"mode":"tn","samples":1200}' | $PY scripts/toy19.py > evidence/tn19.json 2> evidence/tn19.err ) &

# ---- density sweep: one shard per density (budget max(S_ref, 8*2^d), first report)
for d in 5 7 9 11; do
  case $d in
    5)  T=1600 ;;
    7)  T=1600 ;;
    9)  T=800  ;;
    11) T=200  ;;
  esac
  ( echo "{\"mode\":\"sweep\",\"dpbits\":$d,\"trials\":$T,\"trails\":4}" | \
      $PY scripts/sweep19.py > "$E/sweep_d$d.json" 2> "$E/sweep_d$d.err" ) &
done

# ---- design-budget floor census: 10 independent shards, one configuration
for s in 0 1 2 3 4 5 6 7 8 9; do
  ( echo "{\"mode\":\"census\",\"dpbits\":0,\"R\":128,\"trails\":32,\"trials\":4800,\"seed\":$((100+s))}" | \
      $PY scripts/sweep19.py > "$E/census_$s.json" 2> "$E/census_$s.err" ) &
done
wait
echo "ALL DONE $(ls "$E"/*.json | wc -l)"
