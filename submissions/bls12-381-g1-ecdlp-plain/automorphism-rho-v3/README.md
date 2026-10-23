# automorphism-rho-v3: bls12-381-g1-ecdlp-plain candidate

Third revision of the automorphism-aware Pollard rho for BLS12-381 G1 ECDLP,
responding to the second-round review of v2 (same objections as the sibling
secp256k1 package).

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **126.982316** (= 31745579/250000) |
| success_probability (lower bound) | 0.39 |
| memory_log2_bytes | 70.1 (peak strictly above 2^70, rounded up) |
| preprocessing_log2 | 25.6 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v2

- Distance-3 look-back added: the order-3 identity 1+lam+lam^2 = 0 mod r
  creates period-3 fruitless cycles (proven algebraically) that the
  distance-2 detector missed; both families are now detected at first
  return with full-state discrimination and one-addition escapes, all
  charged.
- Per-walker epoch limit of 2^36 steps fences every longer zero-sum family
  of the multiplier ring (pattern-count bound < 2^-40/step).
- Ledger recomputed exactly (986 direct + 140.7 amortized escapes = 1126.7
  non-reference instruction-equivalents against C = 2,000); score derived
  from that sum.
- Memory reported as strictly-above-2^70 rounded up (70.1); detection wait
  charged as an explicit additive allowance (2^88 units, residual failure
  < 2^-369); capacity Chernoff bound stated with the correct ratio 0.065.

## Files

- `proof.md` - construction, dual-family cycle analysis, epoch fencing,
  success decomposition, full work ledger, memory, GT/NFS discussion.
- `claim.json` - bound fields, supersedes notes, sha256 bindings.
- `experiment.py` - Part A production checks (including the order-3
  identity and the 96-step lazy-representation invariant); Part B the
  complete v3 attack plus a 40,000-step cycle-return census.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v3
# or, without Docker:
echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v3/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass (96/96
lazy-representation steps at production size, including
[lam]P0 = zeta(P0) and 1+lam+lam^2 = 0); Part B at the declared input:
85/200 = 0.425 successes vs 0.4845 ideal, 335 distance-2 escapes,
0 distance-3 escapes, 0 extraction-verification failures; census:
distance-2 rate 1/115, distance-3 below census resolution.

Small-scale experiments are evidence at their own scope only.