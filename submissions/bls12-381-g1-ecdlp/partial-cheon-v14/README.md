# partial-cheon-v14: bls12-381-g1-ecdlp candidate

Fourteenth revision: fixes the eleventh-round objections (bulk-tz
binary gcd with all rounds counted; bounded rejection-free offset
sampler; the declared run moved to deterministic full mode so no
Monte Carlo realization of the success probability is submitted).

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

## What changed from v13

- The binary gcd now performs BULK trailing-zero removal (hardware
  count-trailing-zeros model: ONE constant-cost multiword operation
  per strip; no data-dependent one-bit shift loops anywhere) and
  counts every round - each subtraction+strip pair AND the <= 3
  initial bulk strips - asserted <= 520 rounds per call at <= 50
  instructions per round: <= 26,000 instructions per gcd. The
  review's example gcd(2^100, 1) is an explicit test case costing
  exactly 2 rounds. Measured maxima: 154 rounds over the production
  factorization, 253 over the adversarial power-of-two pairs.
- The window offset is drawn by a BOUNDED REJECTION-FREE sampler:
  U uniform in [0, 2^237), o = (U*M) >> 237 - one 512-instruction
  multiplication and a shift, charged in the ledger, no rejection
  loop, no failure branch - and the per-offset-exact coverage bound
  makes the sampler's law immaterial to the success bound.
- The declared run is now {"trials": 400, "mode": "full"}: 400/400
  exact recoveries, 0 wrong extractions, deterministic. The partial
  regime (the claimed >= 0.39) is established by Part A's exact
  per-offset integer computation, which the same declared run
  executes - not by Monte Carlo. (Partial-mode toy runs realize
  0.375-0.44 across seeds, consistent with the exact expectation
  0.4022613461525126 and its ~0.024 standard error at 400 trials.)
- Joint rho worst case recomputed: main 2^30.629357 + gcds
  2^28.251187 + backtrack 2^25.279320 = 2,021,154,816 = 2^30.912533
  instructions = 2^19.946748 units; preprocessing total <= 2^20.02
  units inside the 2^21 budget.

## Reproduce

```sh
echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v14/experiment.py
echo '{"trials": 200, "mode": "full"}'    | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v14/experiment.py
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
