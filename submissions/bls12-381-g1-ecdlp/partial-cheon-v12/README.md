# partial-cheon-v12: bls12-381-g1-ecdlp candidate

Twelfth revision: fixes the ninth-round objections (itemized block
scratch including the inverse-output array; explicit binary-gcd
construction; universal equality assertions matching the claim text).

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

## What changed from v11

- Block scratch itemized to 304 B/record (the review correctly added
  the separately allocated inverse-output array: +48 B) and charged at
  <= 512 B/record x 2^24 = <= 2^33 bytes; no in-place aliasing claimed.
- The rho gcd now has an explicit bounded construction (binary gcd,
  <= 255 rounds x 16 instr = 4096 instr per gcd) and the joint worst
  case is three bounded terms: main 2^30.629357 + batched gcds
  2^25.584963 + backtrack 2^22.754888 = 2^30.678380 instructions =
  2^19.712596 units; preprocessing total <= 2^19.80 units.
- Equality with plain double-and-add is now asserted for EVERY toy
  point (all babies and every giant of every trial) - v11 sampled 8
  giants per block after trial 0, which did not match the universal
  claim text.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v12/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v12/experiment.py
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
