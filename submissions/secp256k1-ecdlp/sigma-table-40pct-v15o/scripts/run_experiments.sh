#!/bin/sh
# run_experiments.sh -- regenerate every number in the v15o package from scratch.
set -e
cd "$(dirname "$0")/.."
rm -f scripts/walk_v15o.s evidence/walk_v15o.s evidence/ledger-v15o.json certificate.json
clang -O2 -arch arm64 -S scripts/walk_v15o.c -o scripts/walk_v15o.s
python3 scripts/count_symbols.py scripts/walk_v15o.s > evidence/counts-v15o.json
python3 scripts/ledger_v15o.py > evidence/ledger-run-v15o.json
python3 scripts/gen_cert_v15o.py
echo "regenerated: scripts/walk_v15o.s evidence/counts-v15o.json evidence/ledger-v15o.json certificate.json"
