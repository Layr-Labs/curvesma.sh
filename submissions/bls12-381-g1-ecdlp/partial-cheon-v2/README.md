# partial-cheon-v2: bls12-381-g1-ecdlp candidate

Second revision of the partial-coverage Cheon attack for the BLS12-381 G1
discrete logarithm with auxiliary powers, scored at n = 2^20, responding to
the recorded review rejection of `partial-cheon` (v1).

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.777511** (= 123777511/1000000) |
| success_probability (lower bound) | 0.39 (exact expectation 2/5) |
| memory_log2_bytes | 124.0 (fixed 2^117-slot table, never overfilled) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v1

- Contiguous-window giants: the m covered blocks are now consecutive at a
  uniform random offset (identical coverage mathematics), so every giant
  scalar is a single sequential mod-r multiplication - v1's random subset
  required a modular exponentiation per giant, contradicting its own
  ledger.
- Exact coverage: the success expectation is computed in exact integer
  arithmetic (exactly 2/5 at production parameters, printed by the
  experiment), with worst case 0.4 - 2^-118 for every offset; no
  "approximately 0.40" claims.
- Explicit charged constructions for the hash table (open addressing,
  fixed 2^117-slot capacity, < 2 expected probes, 1000 instructions per
  operation), the window-start exponentiations (<= 2), and the comb tables.
- The comparison against the sibling rho track's bound is removed from the
  claim (context in the public note only).

Uses d = 2^20 | r-1 (largest divisor of r-1 not exceeding the scoring n)
and reads only A_1 = T and A_d. Phase 1: BSGS for k1 with s^d = eta^k1 over
m = ceil(sqrt(0.40*N)) baby steps and m window blocks, each step a
fixed-base comb scalar multiplication charged at 64 units. Phase 2:
deterministic full-coverage BSGS for the order-2^20 component.

## Files

- `proof.md` - reduction, algorithm, exact coverage analysis, work ledger,
  memory, alternatives and sensitivity.
- `claim.json` - bound fields, W(n) statement, supersedes note, sha256
  bindings.
- `experiment.py` - Part A production checks (including the exact coverage
  rational); Part B the complete two-phase attack on a 20-bit toy group
  (d=6) with sequential window giants.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/partial-cheon-v2
# or without Docker:
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v2/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v2/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass, coverage
rational exactly 2/5; Part B full mode 200/200 exact recoveries, 0 wrong
extractions; partial mode 88/200 = 0.44 vs the exact predicted coverage
0.40226, 0 wrong extractions.

Small-scale experiments are evidence at their own scope only.
