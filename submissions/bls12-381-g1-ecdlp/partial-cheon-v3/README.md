# partial-cheon-v3: bls12-381-g1-ecdlp candidate

Third revision of the partial-coverage Cheon attack for the BLS12-381 G1
discrete logarithm with auxiliary powers, scored at n = 2^20, responding to
the second-round review of v2.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.783097** (= 123783097/1000000) |
| success_probability (lower bound) | 0.39 (exact rational expectation >= 39/100 by integer verification) |
| memory_log2_bytes | 125.1 (2^118-slot table of exactly 2^125 bytes, peak above it) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v2

- The score is computed from the ledger exactly as the review did:
  2m(64+0.192) + 2m(0.5) = 129.384m plus one-offs gives 123.783097; v2's
  123.777499 was not supported by its own table.
- No "exactly 2/5" prose: the coverage expectation is an exact rational
  (extremely close to 0.4, printed by the experiment) verified >= 39/100
  by integer comparison, with worst-single-offset bound 0.4 - 2^-118.
- m is computed by pure integer arithmetic: m = isqrt(floor(2N/5)) + 1
  (v2's int(0.40*N) was float-rounded before the square root).
- Hash table restated honestly: 2^118 slots (load 0.4256, expected probes
  1.74, 1000 instructions per operation, 16-byte indices); memory 125.1
  (v2 claimed exactly 2^124 bytes with 8-byte indices that cannot hold
  117-bit positions).

## Files

- `proof.md` - reduction, algorithm, exact coverage analysis, honest work
  ledger, memory, alternatives and sensitivity.
- `claim.json` - bound fields, W(n) statement, supersedes notes, sha256
  bindings.
- `experiment.py` - Part A production checks (including the exact
  coverage rational and the integer-only m); Part B the complete
  two-phase attack on a 20-bit toy group (d=6) with sequential window
  giants.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/partial-cheon-v3
# or without Docker:
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v3/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v3/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass (coverage
rational printed, >= 39/100 verified in integers); Part B full mode
200/200 exact recoveries, 0 wrong extractions; partial mode 88/200 = 0.44
vs the exact predicted coverage 0.40226, 0 wrong extractions.

Small-scale experiments are evidence at their own scope only.