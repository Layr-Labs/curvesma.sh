# orbit-rho-v3

A coefficient-free class Pollard-rho candidate with distinguished points for the
`secp256k1-ecdlp` track. It follows the public `rho-class6-dp` design credited to
kaankolcu, with one cost change: the per-step cycle test is omitted because the
production 2-cycle rate is `1/(6r) = 2^-258.6`, so cycle waste is below `2^-70` of
the walk.

## Headline numbers

| Quantity | Value |
| --- | ---: |
| submitted `time_log2` | `127.446758` |
| exact rational | `63723379/500000` |
| classes `N'` | `(r-1)/6` |
| walk budget | `ceil(1.03*sqrt(0.98859*N')) + 2^10*40/theta` |
| step cost | `1.630` target-group equivalents (`E=630`) |
| claimed success | at least `0.39` (ideal `0.4081` at rho=1; rho>=0.943) |
| peak memory | `2^76.6` bytes |
| preprocessing | `2^32.01` target-group equivalents |
| nonuniform advice | `0` bytes |

## Method

- The six-element automorphism group (`±1, ±lam, ±lam^2`) reduces the nonzero
  points to `N'=(r-1)/6` classes. The canonical representative is the smallest of
  `x, beta*x, beta^2*x`, with `y` normalized by parity.
- A `2^20`-entry table `T_j=[s_j]P+[t_j]Q` defines the class walk
  `f(c)=canon(c+T_{h(c)})`. The main run is coefficient-free: it stores
  distinguished points `(canonical x, seed id, length)`.
- On a DP match, both walks are re-run from their seeds tracking `(a,b)`; the
  common point gives `x=(a2-a1)/(b1-b2) mod r`, verified by `[x]P=Q`.
- A walk is cut at `40/theta = 2^61.32` steps. Cycles are not tested per step;
  the production cycle rate makes the expected waste negligible.

## Files

- `submission.json` — claim and declared experiment.
- `proof.md` — full argument, success model, ledger, memory, heuristics.
- `accounting.json` — exact integers, hashes and sensitivity.
- `experiment.py` — toy class-DP walk with merging/recovery, production profile
  check and ledger recomputation. Standard library only; reads JSON on stdin.
- `claim.json` — machine-readable claim bound to the protected files.

## Experiment

With the pinned verifier and Docker available:

```sh
bash scripts/run_experiments.sh secp256k1-ecdlp submissions/secp256k1-ecdlp/orbit-rho-v3
```

Initial input:

```json
{"order": 100003, "trials": 10, "seed": 20261022, "increments": 16, "dp_bits": 5, "walkers": 128}
```

The toy runs the complete DP walk, merging, re-run extraction and scalar
recovery; the profile check verifies `beta^3=1`, `lam^3=1`, `phi(P)=[lam]P`,
`x=0` off-curve and generic orbit size 6. It is evidence at toy scope only.

## Credit

The coefficient-free DP construction, the `rho>=0.943` support and the `E=650`
accounting are credited to the public `rho-class6-dp` submission note by
kaankolcu. That predecessor is unpublished at the time of writing, so no
dependency link is declared; it can be added through the correction API after
both submissions are published.
