# Peak storage

Memory is reported, not scored, and has no tie-break (COST_MODEL.md). We bound
peak total storage (code, advice, tables, matrix, scratch).

## Dominant stores

| Store | Size (log2 bytes) | Basis |
| --- | ---: | --- |
| Sparse LA matrix (block Wiedemann) | ~ 60.6 | D ~ 2^49.6 rows, weight ~ 2^7 nonzeros, ~16 B/entry |
| Relation buffer (pre-filter) | ~ 63.8 | ~ 2^54.85 relations, few hundred bytes each |
| Factor-base / sieve tables | ~ 55 | ~ B = 2^54.85 primes |
| Code + advice | < 2^30 | no parameter-dependent precomputed table embedded |

Peak ~ 2^63.8 bytes; stated bound **2^66 bytes**. No unbounded table, no embedded
searched advice (nonuniform advice = 0). The LA and relation stores are consistent
with the matrix dimension D = 2^49.6 of proof.md section 5; `experiment.py` reports
log2_peak_bytes and asserts it exceeds 2^50 (consistency with the D ~ 2^50 matrix).
