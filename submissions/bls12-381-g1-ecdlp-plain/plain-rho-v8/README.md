# plain-rho-v8: bls12-381-g1-ecdlp-plain candidate

v6's structured constants are replaced by published pseudorandom
constants; the confinement obstruction found by the v6 review is
measured away, cycle words of length <= 4 are excluded by executed
enumeration, and every ledger allowance is itemized.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **1997276/15625 = 127.825664** |
| success_probability (lower bound) | 0.39 (modeled 0.4066) |
| memory_log2_bytes | 104.1 (over-report) |
| preprocessing_log2 | 23.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why pseudorandom constants

The v6 review found that with c_j = jQ + 2^j and r = 32Q + rho,
the rescaling z = 32k mod r turns every step into a positive
increment <= 2^36: each trail is confined to width 2^36 * L and
cross-walker collisions are a < 2^-75 event at production scale.
v8 uses c_j = int(sha512("curvesmash/bls12-381-g1-ecdlp-plain/constant/j")) mod r -
published, deterministic, structure-free. Part A measures every
rescaling up to 2^12 (narrowest step window 0.6897 of the
modulus; the same check detects the v6 constants at 1.31e-66) and
the counting heuristic 32*W^31/r^30 excludes confining rescalings
below width 2^240 for all M. Cycle words of length <= 4 are
excluded by enumeration; longer cycles are 2^-127 events under
the random-function model, detected free by the table. The toy
measures 1.052 of the birthday ideal
with 0 wrong extractions.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only; deterministic; single JSON object on stdout. The
declared run also records the dbits=6 sparse-DP census.
