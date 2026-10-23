# automorphism-rho-v2: bls12-381-g1-ecdlp-plain candidate

Second revision of the automorphism-aware Pollard rho for BLS12-381 G1
ECDLP, responding to the recorded review rejections of `automorphism-rho`
(v1): the incorrect O(k^2|G|/N) fruitless-cycle claim, the divisor-based toy
distinguished-point predicate, and the ledger omissions/inconsistency.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **126.976677** (= 126976677/1000000) |
| success_probability (lower bound) | 0.39 |
| memory_log2_bytes | 70.0 (hard 2^63-slot table cap) |
| preprocessing_log2 | 25.6 (included in total; v1 slip corrected) |
| nonuniform_advice_log2_bytes | 0 |

## What changed from v1

- Corrected fruitless-cycle analysis: the negation-map 2-cycle triggers at
  a constant-order rate (~1/192 per step, bounded by 1/64) - measured at
  1/183-1/287 in the toy. It is detected by a 2-step look-back on the full
  walk state, escaped with one charged addition; 3-8-term constant
  constellations are excluded by a charged setup check plus union bounds; a
  5% residual efficiency slack is charged.
- Canonical-state walk with lazy representations: per-step representation
  update is a single mod-r addition plus table lookups; the canonical-frame
  multiplier is one of six values, materialized only at collisions. The
  invariant is verified along 96 production walk steps in Part A.
- Production distinguished-point predicate everywhere (low bits zero); the
  toy measures it via the dbits parameter.
- Ledger completed: representation addition, table/loop work, explicit
  analytic mod-r prices, corrected preprocessing arithmetic
  (2^17 x 384 = 2^25.585), hard memory cap with Chernoff overflow bound,
  detection lag subtracted from the success probability.

## Files

- `proof.md` - corrected construction, cycle analysis, success
  decomposition, full work ledger, memory caps, GT/NFS discussion,
  heuristics, limitations.
- `claim.json` - bound fields, supersedes note, sha256 bindings.
- `experiment.py` - Part A production checks (automorphism, G6 closure,
  constant setup check, 96-step lazy-representation invariant); Part B the
  complete v2 attack on a 20-bit j=0 toy curve.
- `submission.json` - claim statement, exact score fraction, experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v2
# or, without Docker:
echo '{"trials": 200, "dbits": 2, "walkers": 8}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho-v2/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass (96/96
lazy-representation steps at production size, including
[lam]P0 = zeta(P0) with
lam = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001);
Part B at the declared input: 85/200 = 0.425 successes vs 0.4845 ideal
(toy-scale detection lag, analyzed in proof.md), 335 escapes (rate 1/287
per step, inside the predicted 1/64-1/192 band), 0 extraction-verification
failures. dbits sweep (200 trials, 8 walkers): 0 -> 0.545, 2 -> 0.425,
4 -> 0.315.

Small-scale experiments are evidence at their own scope only.
