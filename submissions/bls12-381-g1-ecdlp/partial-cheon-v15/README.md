# partial-cheon-v15: bls12-381-g1-ecdlp candidate

Fifteenth revision: fixes the twelfth-round objections (the offset
sampler's law stated exactly instead of described as uniform; the
evidence prose matched to the declared full-mode run's recorded
7-block shape). No algorithm, ledger or figure changes.

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.836032** (= 1934938/15625) |
| success_probability (lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118; no probabilistic quantity anywhere) |
| memory_log2_bytes | 124.8 (global peak = sort phase, 256m = 2^124.7676 bytes) |
| preprocessing_log2 | 21.0 (included; preprocessing total <= 2^20.01 units) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v14

- The sampler's law is stated exactly: o = floor(U*M/2^237) for U
  uniform in [0, 2^237) assigns each offset floor(2^237/M) or
  ceil(2^237/M) of the 2^237 values - a relative deviation below
  2^-119 from 1/M. The claim no longer describes the offset as
  uniform (exactly uniform only if M divides 2^237, which it does
  not); the per-offset-exact coverage bound - which the review
  agreed is not invalidated by the bias - makes the exact law
  immaterial.
- The evidence description matches the declared full-mode run's
  recorded output: m_toy = M = 418, seven blocks of capacity 64 with
  the final block holding 34 records, live_accs = live_scalars = 64
  as a MAXIMUM over blocks (v14's prose had described the
  partial-mode 5-block shape).
- Everything else unchanged: score 123.836032, memory 124.8, the
  bulk-tz binary gcd, the joint rho worst case, the coverage bound.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v15/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v15/experiment.py
```

Observed locally (declared input): Part A all checks pass, including
the executed capped Brent-rho factorization (111,717 mulmods, 310 gcds
at max 154 binary-gcd rounds, 9 cycles; matches the supplied list
exactly; the bulk-tz binary gcd matches math.gcd on all 207 sampled
and adversarial pairs, gcd(2^100, 1) = 2 rounds) and the per-offset
minimum >= 39/100 in integers; Part B full mode 400/400 exact
recoveries, 0 wrong extractions, streaming-comb live-array assertions
(live_accs = live_scalars = 64), sort compare/move upper bounds and
the 9-probe toy budget asserted per run.

Small-scale experiments are evidence at their own scope only.
