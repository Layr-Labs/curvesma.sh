# partial-cheon-v13: bls12-381-g1-ecdlp candidate

Thirteenth revision, regenerated whole: fixes the tenth-round
objections (every stale 2^32-byte scratch bound replaced consistently
across all five files; the gcd price is now an executed, round-counted
binary-gcd construction; all documents rewritten from a single
coherent source so no generation can leave stale text).

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.836032** (= 1934938/15625) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118; no probabilistic quantity anywhere) |
| memory_log2_bytes | 124.8 (global peak = sort phase, 256m = 2^124.7676 bytes) |
| preprocessing_log2 | 21.0 (included; preprocessing total <= 2^20.01 units) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v12

- Every remaining stale "2^32-byte" block-scratch bound (the review
  found it still in the claim statement, claim limitations and proof
  paragraphs while the README itemized 304 B/record) is replaced by
  the itemized allowance - 304 B/record (scalars 32, accumulators 96,
  prefix 48, inverse-output 48, pending 72, indices 8), charged at
  <= 512 B/record x 2^24 = <= 2^33 bytes - enforced by a consistency
  grep across all five files.
- The gcd price is now an EXECUTED construction: a round-counted
  binary gcd (Stein's algorithm) replaces math.gcd throughout the rho
  factorization; <= 48 instructions per round (4-limb compare 10,
  conditional subtract 10, multiword trailing-zero elimination 12,
  multiword shift 10, control 6), <= 512 rounds per call (each round
  strips >= 1 bit of a < 2^256 total), asserted per call, measured
  maximum 153, verified against math.gcd on 201 sampled pairs
  including (r-1, 254760293^2).
- Joint rho worst case recomputed: main 1,660,944,384 + gcds
  301,989,888 + backtrack 38,535,168 = 2,001,469,440 = 2^30.898412
  instructions = 2^19.932628 units; preprocessing total <= 2^20.01
  units inside the 2^21 budget.
- The proof, claim and manifest are regenerated from scratch (v12's
  accumulated patches had left stale text in several files, which is
  what the review caught); a grep gate now blocks submission on any
  stale token.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v13/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v13/experiment.py
```

Observed locally: Part A all checks pass, including the executed
capped Brent-rho factorization (111,717 mulmods, 310 gcds at max 153
binary-gcd rounds, 9 cycles; matches the supplied list exactly; the
binary gcd matches math.gcd on 201 sampled pairs) and the per-offset
minimum >= 39/100 in integers; Part B full mode 200/200, 0 wrong
extractions; partial mode 88/200 = 0.44 vs the exact wrap-inclusive
expectation 0.40226, 0 wrong extractions, streaming-comb live-array
assertions (live_accs = live_scalars = 64), sort compare/move upper
bounds and the 9-probe toy budget asserted per run.

Small-scale experiments are evidence at their own scope only.
