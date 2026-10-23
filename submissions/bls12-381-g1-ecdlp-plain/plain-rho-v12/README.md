# plain-rho-v12: bls12-381-g1-ecdlp-plain candidate

Distinguished points go to an append-only arena; duplicates are
found by one final merge sort. Every work term is a hard worst-case
bound; expected values enter only as Markov probabilities inside
the abort total, and the collision-factor tolerance is computed
after all deductions. No external factor is used.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **64044349/500000 = 128.088698** |
| success_probability (lower bound) | 0.39 (modeled 0.6096 at the toy's factor) |
| memory_log2_bytes | 105.1 (over-report) |
| preprocessing_log2 | 25.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why the arena

The v11 review demanded bounded operation counts and a table
construction whose links fit the record scale. v12 deletes the
table: appends cost 20 instructions worst case, the arena holds at
most 2^98 records (next append aborts, Chernoff-subtracted), and
one final sort (2^97.8
units) finds duplicates. No links, no indices, no chains.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only; deterministic; single JSON object on stdout. The
declared run also records the dbits=6 sparse census.
