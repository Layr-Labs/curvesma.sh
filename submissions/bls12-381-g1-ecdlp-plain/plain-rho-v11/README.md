# plain-rho-v11: bls12-381-g1-ecdlp-plain candidate

v6's structured constants are replaced by published pseudorandom
constants; the confinement obstruction found by the v6 review is
measured away, cycle words of length <= 4 are excluded by executed
enumeration, and every ledger allowance is itemized.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **64044349/500000 = 128.088698** |
| success_probability (lower bound) | 0.39 (modeled 0.6053 at Teske's measured factor) |
| memory_log2_bytes | 105.1 (over-report) |
| preprocessing_log2 | 25.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why pseudorandom constants

The v6 review found that with c_j = jQ + 2^j and r = 32Q + rho,
the rescaling z = 32k mod r turns every step into a positive
increment <= 2^36: each trail is confined to width 2^36 * L and
cross-walker collisions are a < 2^-75 event at production scale.
v11 uses c_j = int(sha512("curvesmash/bls12-381-g1-ecdlp-plain/constant/j")) mod r -
published, deterministic, structure-free. Part A measures every
rescaling up to 2^12 (narrowest step window 0.6897 of the
modulus; the same check detects the v6 constants at 1.31e-66) and
the counting heuristic 32*W^31/r^30 excludes confining rescalings
below width 2^240 for all M. Cycle words of length <= 4 are
excluded by enumeration; longer cycle entries occur at the
self-pair birthday scale t^2/(2Kr), with waste beyond a 2^-23
budget fraction of probability <= 2^-10.9 by Markov on the
bounded expectation and reseed work charged. The toy
measures 0.985 of the birthday ideal
with 0 wrong extractions.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only; deterministic; single JSON object on stdout. The
declared run also records the dbits=6 sparse-DP census.
