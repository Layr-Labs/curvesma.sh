# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence; `measurements12.txt` is the last section of `kernels_asm.b85` (base85 of lzma
text). The raw `dsieve16` log is summarized in the table below.

## Descent sieve: `dsieve16`

Builds the determinant-p hit lattice of a root-pair functional (p < 2^70), LLL-reduces
and certifies it exactly (congruences, Bareiss determinant = p), enumerates the ball
||c||^2 <= R^2 (Schnorr-Euchner, exact leaf test), sign-canonicalises and ranks each
hit (exact lexicographic ball index) and emits (owner = index >> 40, offset, 16 log2 p)
records; owners stay below 2^62.44 at R = 384 (`tables` refuses half counts >= 2^103).
Its `check` passed: set equality and rank bijection with brute force for R <= 4;
every record unranked and verified at p = 2^20..2^69, R = 10..224 (18,051,918 at
R = 224) and the first 2^24 at R = 256; the split is audited per node ([7]).

| measured | Ir | charged |
|---|---:|---:|
| enumeration + rank + record, per hit (p = 2^20..2^69, R = 12..224) | 302-793 | 1000 |
| transfer to the owner, per hit (four radix passes) | (STORAGE.md) | 84 |
| apply, per record | 9.00 | 9 |
| SIEVE_HIT total | | 1100 |
| scan, per cell | 5.000 | 5 |
| zero, per cell | | 1 |
| certified per-cell norm bound, both sides (DESCENT.md) | ~1,500 | 2,400 |
| LLL + certificate per (p, root pair), p = 2^69 | 2^24.98 | 2^25.1 |

At p = 2^69 the per-hit cost decreases with the radius (417.8, 370.0, 353.9 Ir
at R = 160, 192, 224: fewer enumeration nodes per hit, the rank is 12 lookups at
every radius); the 1000 charge is 2.8x the R = 224 value and covers R up to 384.
A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per hit
including two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold the
index tables (2^34.34 bytes at R = 384, STORAGE.md); `apply` runs at R = 4; the
`realloc` record buffer belongs to the harness, production streams records to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Line and record paths of the relation sieve (per-update costs 24.10-298.78 and
44.67-1711.20 by hits per coset, table below; line path for p < 2^48, record walk above (v3 band strip); apply 7.00, recovery 7.003 per record; layout
12.325, zero 0.625, scan 8 per cell), certified norm bounds (1,508.2 per
evaluation, 4,300 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`; the
1,508.2 is the measured marginal per-block cost, charged exactly in this revision), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the LA kernels of STORAGE.md (group loops 11/33, the v6 five-vector Krylov loops 39/33 of `krylov5.s`, the row tail and reductions of `row_tail.c`, and the 97-instruction multiply into a 9-limb accumulator `bw_full`; listings in `kernels_asm.b85`, counted by [10]): the constants of `model16.py`.

### Slice sieve per-update costs (re-measured, fused exact recentring; 2c2ce2ab)

Coset centres are recentred from the exact adjugate remainder (a_i c - q_i p in
__int128) with per-pair affine forms; all `slice12` checks (line, record, apply,
scan, recover vs brute force, primes up to 2^63) MATCH. Callgrind Ir per unit:
line 24.10 / 31.07 / 87.57 / 298.78 (p = 257, 4099, 65537, 1048583); record
44.67 / 77.42 / 250.54 / 810.82 / 1711.20 (p = 257, 65537, 8388617, 67108879,
268435459); setup6 406.87 / 374.07. These are `LINE_TAB` / `REC_TAB`; the log is
the `measurements12.txt` section of `kernels_asm.b85`.
