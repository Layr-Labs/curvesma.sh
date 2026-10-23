# plain-rho-v17: bls12-381-g1-ecdlp-plain candidate

Distinguished points go to an append-only arena; duplicates are
found by one in-place heapsort. Every work term is a hard worst-case
bound; expected values enter only as Markov probabilities inside
the abort total, and the collision-factor tolerance is computed
after all deductions. No external factor is used.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **128090067/1000000 = 128.090067** |
| success_probability (lower bound) | 0.39 (modeled 0.6053 at the toy's factor) |
| memory_log2_bytes | 63.1 (over-report) |
| preprocessing_log2 | 25.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why the arena

The v11 review demanded bounded operation counts and a table
construction whose links fit the record scale. v17 deletes the
table: appends cost 20 instructions worst case (64-bit linear
addressing), the arena holds at most 2^cap records, and the whole
arena fits the 64-bit address space (next append aborts,
Markov-subtracted), and one in-place heapsort plus scan
(2^57.431356 units) finds
duplicates. No links, no indices, no chains.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only; deterministic; single JSON object on stdout. The
declared run also records the dbits=6 sparse census.
