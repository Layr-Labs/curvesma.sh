# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the logs are the evidence.

## Descent sieve: `dsieve16`, `measurements16.txt`

For a prime p < 2^70 and a root-pair functional w (12 residues mod p) the kernel
builds the hit lattice {c : w.c = 0 mod p} (determinant p), LLL-reduces it
(long-double decisions, exact __int128 updates) and certifies the result exactly
(every basis vector satisfies the congruence; Bareiss determinant = p). It
enumerates the ball ||c||^2 <= R^2 by Schnorr-Euchner in double precision with
an exact integer ball test at each leaf, sign-canonicalises each hit, ranks it
by the exact lexicographic ball index (unsigned __int128) and emits the record
(u64 owner = index >> 40, u64 offset = index mod 2^40, u16 16 log2 p). The
largest shipped ball (R = 384) has half count 2^102.44, so owners stay below
2^62.44; `tables` refuses any ball with half count >= 2^103.

Correctness (`check`, all PASS on the shipped source): set equality with brute
force and a rank bijection for R <= 4 (p from 307 to 2^40); at p = 2^20..2^69 and
R = 10..224 every emitted record is unranked from (owner, offset) and verified
(congruence, ball, sign, no duplicates) and hit counts match volume/p (all
18,051,918 records at R = 224); at R = 256 the first 2^24 records, with indices
up to 2^95.4 (owners up to 2^55.42). The split is a pure function of the index,
audited for every sieved node ([7], store16.check).

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

## Main phase: `slice12`, `ecm12`, `la_mac.s`, `bw_dot.s`

Line and record paths of the relation sieve (per-update costs 24.10-310.70 and
44.90-2231.66 by hits per coset; apply 7.00, recovery 7.003 per record; layout
12.325, zero 0.625, scan 8 per cell), certified norm bounds (1,550 per
evaluation, 4,300 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the 33-instruction `la_mac.s` loop and the 97-instruction full-field multiply-accumulate loop bw_full in `bw_dot.s` (projections, dense Schirokauer columns, reconstruction) ([10]): the constants of `model16.py`.
