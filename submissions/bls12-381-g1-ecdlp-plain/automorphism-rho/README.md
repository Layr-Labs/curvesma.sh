# automorphism-rho: bls12-381-g1-ecdlp-plain candidate

Automorphism-aware Pollard rho for BLS12-381 G1 ECDLP, submitted under the
protected `classical-group-work-v2` cost model (G1 reference price C = 2,000).

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **127.143439** |
| success_probability (lower bound) | 0.39 |
| memory_log2_bytes | 70.3 |
| preprocessing_log2 | 25.4 (included in total) |
| nonuniform_advice_log2_bytes | 0 |

The attack walks on the (r-1)/6 classes of BLS12-381 G1 under its order-6
automorphism group (negation plus the j = 0 cube-root map (x, y) -> (zeta*x, y),
whose multiplier lam on <P0> satisfies [lam]P0 = zeta(P0)). It uses 2^16
parallel distinguished-point walkers, a 32-adding walk, a single-shot budget of
1.1*sqrt((r-1)/6) total walker steps, and worst-case per-step charging of
1.8275 reference units. This is the first qualified score for a track whose
selected template is unscored.

## Files

- `proof.md` - complete construction, success analysis, work ledger, memory
  accounting, GT/NFS discussion, heuristics and limitations.
- `claim.json` - claimed bound fields with profile/cost-model sha256 bindings.
- `experiment.py` - stdlib-only experiment: Part A verifies the order-6
  automorphism structure and canonicalization bookkeeping exactly at
  production parameters; Part B runs the full attack on a 20-bit j = 0 toy
  curve and reports the measured single-shot success rate at the claimed
  budget shape.
- `submission.json` - claim statement, exact score fraction and experiment
  declaration.

## Reproduce

```sh
bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/automorphism-rho
# or, without Docker:
echo '{"trials": 200, "dp": 1}' | python3 submissions/bls12-381-g1-ecdlp-plain/automorphism-rho/experiment.py
```

Observed locally (arm64, CPython 3.14): Part A all checks pass, including
[lam]P0 = zeta(P0) with
lam = 0x73eda753299d7d483339d80809a1d804a7780001fffcb7fcfffffffe00000001;
Part B at the declared input: 84/200 successes (0.42) vs ideal birthday
prediction 0.4545, 0 extraction-verification failures, 348 degenerate
reseeds; dp sweep 0/1/2/4 -> 0.42/0.42/0.355/0.33 (toy-scale
distinguished-point lag, analyzed in proof.md).

Small-scale experiments are evidence at their own scope only; see proof.md
for the analytical production bound and its disclosed heuristics.
