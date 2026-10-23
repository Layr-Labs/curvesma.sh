# negation-rho-v5: bls12-381-g1-ecdlp-plain candidate

The automorphism walks of v1-v4 are replaced by the negation-only
canonical walk, whose fruitless-cycle structure is COMPLETELY
analyzable under the single walk heuristic.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **127.439507** (= 127439507/1000000) |
| success_probability (lower bound) | 0.39 (modeled 0.4284) |
| memory_log2_bytes | 103.0 |
| preprocessing_log2 | 22.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why negation-only

The fourth-round reviews correctly identified that the order-3
automorphism identities (1 + lam + lam^2 = 0) generate mixed-index
fruitless families whose rates v4 could not bound. The negation walk
with power-of-two constants has none of that structure: every cycle
relation is an exact integer equation, the kappa-to-suffix bijection
counts solutions exactly, and every genuine family costs exactly
64^-(L-1). The census is executed (not estimated): 1, 157, 46,844
families at L = 2, 4, 6, conservative bounds beyond, undetected tails
ratchet-fenced at 2^-13.25.

## Reproduce

```sh
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | python3 submissions/bls12-381-g1-ecdlp-plain/negation-rho-v5/experiment.py
```

Observed locally (declared input): Part A all checks pass (census
asserted; ledger rationals asserted); Part B: 91/200 = 0.455
recoveries, 0 wrong extractions, escape rate 0.01795/step with the
distance-2 rate matching the 1/64 census prediction.

Small-scale experiments are evidence at their own scope only.
