# partial-cheon-v10: bls12-381-g1-ecdlp candidate

Tenth revision: fixes the seventh-round objections (peak-memory phase
analysis with a blocked comb; deterministic capped factorization with
per-batch enforcement and a corrected per-iteration charge).

v9 is v7 resubmitted (v7 and v8 both failed server-side during a
platform incident before any review; content unchanged from v7, plus a
note table-header fix).

Seventh revision of the partial-coverage Cheon attack, fixing the
sixth-round objections: the factorization bound was an expectation with
a units normalization error; instruction allowances lacked line items;
the toy did not exercise the batched comb construction.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.836032** (= 1934938/15625) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118) |
| memory_log2_bytes | 124.8 (main array + merge buffer) |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v9

- Blocked comb generation (2^24 scalars per block): the live working
  storage during giant generation is O(block) records (<= 2^32 bytes)
  instead of O(m) - the review correctly counted 304m bytes of live
  state (baby table + giant scalars + accumulators + prefix products)
  exceeding the claimed 2^124.8. With blocking, the global peak is
  exactly the sort phase: 256m = 2^124.7676 bytes -> memory_log2_bytes
  = 124.8 unchanged, now derived by phase analysis. The merge buffer is
  freed after the sort.
- Rho caps enforced BEFORE each consumed batch (a completed cycle
  provably uses at most 2^17 iterations - the v9 assert fired only
  after a whole doubling round), hard 12-cycle guard added, and the
  per-iteration charge corrected to 1056 instructions (two 512-instr
  modular multiplications + 32 of loop work - the v9 520-instr
  allowance undercounted the double-mulmod inner loop). Worst case
  12 x 2^17 x 1056 + bounded backtracks = 2^30.63 instructions =
  2^19.63 units, still inside the 2^21 preprocessing budget.
- The birthday-tail probability model is withdrawn entirely: the input
  r-1 is the fixed public parameter, the seed sequence is fixed, and the
  construction is deterministic - measured 111,717 mulmods, 310 gcds,
  9 cycles, reproduced by the hosted run. Success is the clean exact
  statement 0.4 - 2^-118 >= 0.39.
- Toy exercises the blocked path (5 blocks of 64; per-block pending and
  accumulator sizes recorded: max 63 of 64).

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v10/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v10/experiment.py
```

Observed locally: Part A all checks pass, including the executed
capped Brent-rho factorization of the production r-1 (111,717 mulmods,
310 gcds, 9 cycles; matches the supplied list exactly) and the
per-offset minimum >= 39/100 in integers; Part B full mode 200/200,
0 wrong extractions; partial mode 88/200 = 0.44 vs the exact
wrap-inclusive expectation 0.40226, 0 wrong extractions, blocked-comb
equality assertions, sort compare/move upper bounds and the 9-probe toy
budget asserted per run.

Small-scale experiments are evidence at their own scope only.
