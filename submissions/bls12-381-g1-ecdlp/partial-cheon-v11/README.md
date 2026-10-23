# partial-cheon-v11: bls12-381-g1-ecdlp candidate

Eleventh revision: fixes the eighth-round objections (joint rho
worst-case accounting with exact computed figures; streaming comb that
actually demonstrates block-local live storage).

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

## What changed from v10

- Joint rho accounting: backtrack iterations consume the same 2^17
  per-cycle cap (enforced before each step), and each iteration costs
  at most 1536 instructions (512-instr mulmod + up to 1024-instr gcd for
  a backtrack step). Worst case: 12 x 2^17 x 1536 = 2,415,919,104 =
  2^31.169925 instructions = 2^20.204141 units, plus certification and
  comb tables: preprocessing <= 2^20.3 units inside the 2^21 budget.
  The v10 figures (additive sum, mis-normalized 2^19.63) are replaced
  by these exact computed values - the review's 2^19.663572 for the
  main loop alone is confirmed.
- The comb now STREAMS: v10's toy still allocated the accumulator
  array for all m scalars and materialized every giant scalar (the
  review correctly noted the assertion did not establish O(block)
  live storage). v11 generates giant scalars block by block (never an
  m-sized list), allocates each block's accumulator array for exactly
  the block size, searches each block's points immediately and
  discards them; baby scalars stream the same way. Per-block live
  accumulator and scalar array lengths are recorded AND asserted
  (live_accs = live_scalars = 64 at toy scale, 5 blocks).
- All figures re-verified: 88/200 partial (0.44 vs the exact 0.40226
  expectation), 200/200 full, 0 wrong extractions, rho 111,717
  mulmods / 310 gcds / 9 cycles.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v11/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v11/experiment.py
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
