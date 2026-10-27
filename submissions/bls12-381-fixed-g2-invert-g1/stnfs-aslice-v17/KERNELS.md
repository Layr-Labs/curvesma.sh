# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19), by the
differential method: two runs that differ only in the measured unit, divided by
the difference in units, so setup and process costs cancel. The sandbox reads
these sources; it does not compile them.

## Descent sieve: `dsieve16.c`, `measurements16.txt`

For a prime p < 2^70 and a root-pair functional w (12 residues mod p) the kernel
builds the hit lattice {c : w.c = 0 mod p} (determinant p), LLL-reduces it
(long-double decisions, exact __int128 updates) and certifies the result exactly
(every basis vector satisfies the congruence; Bareiss determinant = p). It
enumerates the ball ||c||^2 <= R^2 by Schnorr-Euchner in double precision with
an exact integer ball test at each leaf, sign-canonicalises each hit, ranks it
by the exact lexicographic ball index (unsigned __int128; half counts < 2^84 at
R = 128) and emits (owner = index >> 30, offset, 16 log2 p).

Correctness (`check`, 15 cases, all PASS): set equality with brute force and a
rank bijection for R <= 4 (p from 307 to 2^40); at p = 2^20..2^69 and
R = 10..128 every emitted record is unranked and verified (congruence, ball,
sign, no duplicates) and hit counts match volume/p. Exact ball counts for
R = 4..128 lie inside their volume brackets.

| measured | Ir | charged |
|---|---:|---:|
| enumeration + rank + record, per hit (p = 2^20..2^69, R = 12..192) | 287-718 | 1000 |
| transfer to the owner, per hit | (STORAGE.md) | 32 |
| apply, per record | 9.00 | 9 |
| SIEVE_HIT total | | 1050 |
| scan, per cell | 5.000 | 5 |
| zero, per cell | | 1 |
| LLL + certificate per (p, root pair), p = 2^69 | 2^24.98 | 2^25.1 |

A complete (p, root pair) at p = 2^69 and R = 96, setup plus the enumeration of
its 728 hits, measures 2^25.04 Ir (full-enumeration differential: 47,233.5 Ir
per hit including two setups over 1,456 hits), below PAIR_SETUP alone; the
1,000-per-hit charge is added on top. At the largest production radii (p = 2^69,
R = 160 and 192) `check` verifies all 318,329 and 2,839,153 emitted records,
the per-hit cost is 333.3 and 367.1 Ir, and a complete R = 160 enumeration costs
516.8 Ir per hit including two setups; half counts stay below 2^91 in the
unsigned __int128 rank, and the index tables take 2^31.35 bytes per worker
(STORAGE.md). The `apply`
benchmark runs at R = 4, where all 13.3M half cells have owner 0, so its single
array is exactly that owner's score array. The `realloc` record buffer belongs
to the measurement harness; production nodes stream records to owners
(STORAGE.md).

## Main phase: `slice12.c`, `ecm12.c`, `la_mac.s`, `measurements12.txt`

- `slice12.c`: line path and record path of the relation-collection sieve with
  exact recentring and Barrett setups; per-update costs by hits per coset
  (line: 24.10-310.70; record: 44.90-2231.66), apply 7.00 and recovery 7.003
  per record, layout 12.325, zero 0.625 and scan 8 per cell. These are the
  constants of `model16.py`.
- Certified norm bounds: 5,900 instructions per evaluation (measured 5,862.7),
  4,200 per row, 15,300 per slice; `embedding_bounds.py` is the executable
  reference and [2] checks its bounds contain the exact pool norms.
- `ecm12.c`: five-limb Montgomery multiplication, 367.76 Ir, checked on 20,000
  products; ECM curve counts follow the GMP-ECM tables.
- `la_mac.s`: 33-instruction multiply-accumulate loop with 64-bit column
  indices; [10] recounts it.
