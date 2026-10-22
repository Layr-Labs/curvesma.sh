# partial-cheon: bls12-381-g1-ecdlp candidate

Partial-coverage Cheon attack for the BLS12-381 G1 discrete logarithm with
auxiliary powers, scored at n = 2^20 under `classical-group-work-v2`
(G1 reference price C = 2,000).

## Summary

| Field | Value |
| --- | --- |
| time_log2 = log2(W(2^20)) | **123.777499** |
| success_probability (lower bound) | 0.39 (exact expectation 0.40) |
| memory_log2_bytes | 124.8 |
| preprocessing_log2 | 21.0 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

Uses d = 2^20 | r-1 (largest divisor of r-1 not exceeding the scoring n) and
reads only A_1 = T and A_d = [s^d]P0. Phase 1: partial-coverage BSGS for
k1 with s^d = eta^k1, eta of order N = (r-1)/2^20; m = ceil(sqrt(0.40 N))
baby steps plus m fresh-random giant blocks, each step a fixed-base comb
scalar multiplication charged at 64 units. Phase 2: full-coverage BSGS for
the order-2^20 component. Success probability is exact combinatorics
(expectation 0.40, concentration within 2^-57), not a walk heuristic.
Generic rho on T alone (even with the full order-6 automorphism group)
would cost ~2^127.14; the r+1 Cheon variant is worse (largest divisor of
r+1 <= 2^20 is 484486).

## Files

- `proof.md` - complete reduction, algorithm, exact coverage analysis, work
  ledger, memory, alternatives and sensitivity.
- `claim.json` - claimed bound fields, W(n) statement, sha256 bindings.
- `experiment.py` - stdlib-only experiment: Part A verifies the Cheon
  preconditions exactly at production parameters (2^20 | r-1, verified
  factorization, generator, subgroup orders); Part B runs the complete
  two-phase attack on a 20-bit toy group (d=6) in full and partial modes.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/partial-cheon
# or, without Docker:
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass (zeta = 7
generates Z_r*, eta order N = (r-1)/2^20, zeta_e order 2^20, factorization
verified); Part B full mode 200/200 exact recoveries, 0 wrong extractions;
partial mode 85/200 = 0.425 vs 0.4019 predicted, 0 wrong extractions.

Small-scale experiments are evidence at their own scope only; see proof.md
for the analytical production bound.
