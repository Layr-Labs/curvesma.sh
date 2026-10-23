# automorphism-rho-v4: secp256k1-ecdlp candidate

Fourth revision of the automorphism-aware Pollard rho, fixing the two
third-round objections (epoch vs reporting conflict; self-contradictory
cycle-rate bound).

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **127.622849 (= 127622849/1000000; baseline control 129.000802)** |
| success_probability (lower bound) | 0.39 |
| memory_log2_bytes | 102.0 (hard 2^95-slot table at 2^-32 density) |
| preprocessing_log2 | 25.6 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v3

- Distinguished points at 2^-32 density (reporting interval 2^32), epoch
  E = 2^44: the reporting interval sits 2^12 below the epoch, so a merged
  pair reports before reset with probability >= 1 - 2^-11 (exact
  averaging, subtracted). v3's 2^-64 density with E = 2^36 gave reporting
  probability 2^-28 - the fatal conflict the review identified.
- 16-state look-back: every cycle of length 2..16 is detected
  deterministically at its first return (v3 only looked back 2 and 3).
- Exact enumeration of the algebraic cycle families: zero-sum product-1
  multiplier tuples over G6 for L = 2..8 number exactly 1, 2, 15, 60, 340,
  1680, 9135 (entry rates 2^-7.58 down to 2^-42.5); the order-3 identity
  makes them systematic. Lengths >= 17: trivial tuple bound 2^-80 per
  step, fenced by the epoch (trapped fraction 2^-36). v3's false
  "< 2^-40 for lengths >= 4" claim is withdrawn.
- Waste accounting fixed: escape/trapped steps are charged once (inside
  t); the useful-visit fraction carries only the lost-visit effect.
- Census denominators fixed to per-executed-step (v3 divided by
  trials x budget).

## Reproduce

```sh
echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/automorphism-rho-v4/experiment.py
```

Observed locally: Part A all checks pass (including the exact tuple
enumeration); Part B at the declared input: 85/200 = 0.425 vs 0.4845
ideal, 335 escapes, 0 extraction-verification failures; census: 347
returns in 40,000 steps (rate 0.00868/executed step), none at distances
3..16 above census resolution.

Small-scale experiments are evidence at their own scope only.
