# automorphism-rho-v3: secp256k1-ecdlp candidate

Third revision of the automorphism-aware Pollard rho for secp256k1 ECDLP,
responding to the second-round review of v2.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **127.633635** (= 25526727/200000; baseline control 129.000802) |
| success_probability (lower bound) | 0.39 |
| memory_log2_bytes | 70.1 (peak strictly above 2^70, rounded up) |
| preprocessing_log2 | 25.6 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v2

- **Period-3 fruitless cycles.** The review exhibited the order-3 identity
  1+lam+lam^2 = 0 creating period-3 cycles (three same-constant steps with
  multiplier lam or lam^2 return the full state, lazy representation
  included) that the distance-2 detector cannot see. v3 adds a distance-3
  look-back (both families detected at first return, full-state
  fruitless/useful discrimination, one-addition escape) and a per-walker
  epoch limit of 2^36 steps fencing every longer zero-sum family of the
  multiplier ring (pattern-count bound < 2^-40/step for lengths >= 4).
- **Ledger arithmetic.** v2's table summed 644 non-reference instructions
  against a stated total of 641. v3 recomputes: 578 direct + 74.3 amortized
  escapes = 652.3 non-reference instruction-equivalents, score derived from
  that sum.
- **Memory.** Reported as strictly-above-2^70 rounded up: 70.1.
- **Detection and capacity tails.** The post-merge detection wait is now an
  explicitly charged additive allowance (K*2^72 walker-steps = 2^88 units,
  residual failure < 2^-369) instead of a "maximum", and the table-capacity
  Chernoff bound uses the correct ratio (0.065, bound exp(-2^58.2)).

## Files

- `proof.md` - construction, dual-family cycle analysis, epoch fencing,
  success decomposition, full work ledger, memory, heuristics, limitations.
- `claim.json` - bound fields, supersedes notes, sha256 bindings.
- `experiment.py` - Part A production checks (now including the order-3
  identity and the 96-step lazy-representation invariant); Part B the
  complete v3 attack plus a 40,000-step cycle-return census measuring both
  fruitless-family rates.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/automorphism-rho-v3
# or, without Docker:
echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/secp256k1-ecdlp/automorphism-rho-v3/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass (96/96
lazy-representation steps at production size; 1+lam+lam^2 = 0 verified);
Part B at the declared input: 85/200 = 0.425 successes vs 0.4845 ideal,
335 distance-2 escapes, 0 distance-3 escapes, 0 extraction-verification
failures. Census over 40,000 single-walk steps: distance-2 full-state
return rate 1/115 (below the charged 1/64 bound), distance-3 returns 0
(below census resolution, consistent with the ~2^-16.7 algebraic estimate).

Small-scale experiments are evidence at their own scope only.