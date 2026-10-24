#!/bin/sh
# run_experiments.sh -- regenerate every number in the v15p package from scratch.
set -e
cd "$(dirname "$0")/.."
clang -O2 -arch arm64 -I scripts -S scripts/walk_v15o.c -o scripts/walk_v15o.s
python3 scripts/count_symbols.py scripts/walk_v15o.s > evidence/counts-v15o.json
python3 scripts/check_layout.py > /dev/null
python3 scripts/ledger_v15o.py > evidence/ledger-run-v15p.json
python3 scripts/label_invariant.py > evidence/label-invariant-run-v15p.txt || true
python3 scripts/gen_cert_v15p.py
echo "regenerated: walk_v15o.s counts layout ledger label-invariant certificate"
