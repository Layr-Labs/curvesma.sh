# Peak storage

Memory is reported, not scored, and has no tie-break (COST_MODEL.md). `experiment.py`
recomputes the two dominant stores from the same factor-base size and matrix
dimension used in the cost ledger.

## Dominant stores

| Store | Size (log2 bytes) | Basis |
| --- | ---: | --- |
| Relation buffer (pre-filter) | ~66.0 | ~ |FB| = 2^57.97 relations, few hundred bytes each |
| Sparse LA matrix (block Wiedemann) | ~66.3 | D = 2^57.97 rows, weight ~18 nonzeros, ~16 B/entry |
| Factor-base / sieve tables | ~58 | ~ |FB| = 2^57.97 primes |
| Code + advice | < 2^30 | no parameter-dependent precomputed table embedded |

Peak ~ **2^67.1 bytes** (`experiment.py` `storage_log2_bytes`); stated bound
**2^70 bytes**. Both dominant stores are consistent with the matrix dimension
D = 2^57.97 and factor-base size of proof.md section 5. No unbounded table and no
embedded searched advice (nonuniform advice = 0).
