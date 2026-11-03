# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence. The raw `dsieve16` log is omitted for size; the table below summarizes it.

## Descent sieve: `dsieve16`

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
| enumeration + rank + record, per hit (p = 2^20..2^69, R = 12..224) | 302-793 | 800 (R8; plan R >= 192, measured <= 418 there) |
| transfer to the owner, per hit (four radix passes) | (STORAGE.md) | 84 |
| apply, per record | 9.00 | 9 |
| SIEVE_HIT total | | 850 |
| scan, per cell | 5.000 | 5 |
| zero, per cell | | 1 |
| certified per-cell norm bound, both sides (DESCENT.md) | ~1,500 | 2,000 (R8) |
| LLL + certificate per (p, root pair), p = 2^69 | 2^24.98 | 2^25.1 |

At p = 2^69 the per-hit cost decreases with the radius (417.8, 370.0, 353.9 Ir
at R = 160, 192, 224: fewer enumeration nodes per hit, the rank is 12 lookups at
every radius); the 800 charge is 2.3x the R = 224 value and covers R up to 384.
A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per hit
including two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold the
index tables (2^34.34 bytes at R = 384, STORAGE.md); `apply` runs at R = 4; the
`realloc` record buffer belongs to the harness, production streams records to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `la_mac.s`, `bw_dot.s`, `measurements12.txt`

Line and record paths of the relation sieve (per-update costs 24.0984-310.7003 and
44.9007-2231.6591 by hits per coset h (line: h = points per line coset =
cells/K/p; record: h = 2 cells/p), exact raw-count quotients; apply 7.00, recovery 7.003 per record; layout
12.325, zero 0.625, scan 8 per cell), certified norm bounds (1,508.1813 per
evaluation, 4,300 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`; the
1,508.1813 is the exact raw-count quotient of the measured marginal per-block cost, charged exactly in this revision), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the 33-instruction `la_mac.s` loop and the 97-instruction full-field multiply-accumulate loop bw_full in `bw_dot.s` (projections, dense Schirokauer columns, reconstruction) ([10]): the constants of `model16.py`.

R15: line h = S/(K p) >= 2^6.79 at every row and band for p < 2^48
(audit `line_measured`), inside the measured line-tab h-range with the flat
charge above it; secant interpolation overestimates the convex-decreasing
cost in log h (marginal hit + 488.24 setup per coset over h hits).

R16: band-strip uses the `ecm16.py` schedule at (48 bits, 832/576 bits,
99.95%): B1 = 2^13, B2 = 30 B1, 101 curves/round, two rounds/side (audit
`band_strip`); exhaustion divisions ride the charged gcd/redo budget.

R17: grid takes every log2 B1 in 10..36 (was even only); main 69-bit
now B1 = 2^15, 764 curves.

R18: band-strip at 99.9% (92 curves/round); grid adds B2/B1 = 20 (boot B1 = 2^28).
