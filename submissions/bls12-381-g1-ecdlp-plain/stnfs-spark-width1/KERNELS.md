# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence. The large raw `dsieve16` log (`measurements16.txt`) is omitted from
this archive for size; its counts and checks are summarized in the table below.

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

Correctness (`check`, all PASS): set equality with brute force and a rank
bijection for R <= 4 (p from 307 to 2^40); at p = 2^20..2^69, R = 10..224 every
record is unranked and verified (congruence, ball, sign, no duplicates) and hit
counts match volume/p (18,051,918 records at R = 224); at R = 256 the first
2^24 records (indices up to 2^95.4). The owner split is a pure function of the
index, audited for every sieved node ([7], store16.check).

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

At p = 2^69 the per-hit cost falls with the radius (417.8, 370.0, 353.9 Ir at
R = 160, 192, 224); the 1000 charge is 2.8x the R = 224 value and covers R up
to 384. A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per
hit including two setups over 1,396 hits (2^24.99 per setup, below
PAIR_SETUP). Workers hold the index tables (2^34.34 bytes at R = 384,
STORAGE.md); production streams records to owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `la_mac.s`, `bw_dot.s`, `measurements12.txt` (the last three are sections of `kernels_asm.b85`, base85 of lzma)

Line and record paths of the relation sieve (per-update costs 24.10-310.70 and
44.90-2231.66 by hits per coset, charged at the exact logged quotients since R11; apply 7.00, recovery 7.003 per record; layout
12.325, zero 0.625, scan 8 per cell), certified norm bounds (1,508.2 per
evaluation, 4,300 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`; the
1,508.2 is the measured marginal per-block cost, charged exactly in this revision), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the 33-instruction `la_mac.s` loop and the 97-instruction full-field multiply-accumulate loop bw_full in `bw_dot.s` (projections, dense Schirokauer columns, reconstruction) ([10]): the constants of `model16.py`.

**Primes in [2^63, 2^64) (R6).** Prices are measured for p < 2^63 (largest
setup point 2^63 - 25); these primes need double-word residues, so their
record walk, empty cosets, apply and setup6 are charged 2x (`WIDE_REC`,
`WIDE_SETUP`): a stated cushion under H4.

**Line sieve to 2^47 (scores only; band strip).** The line path sieves p in [2^11, 2^47).
Its price is parameterised by h (points per coset), not p: production h spans
2^7.90-2^45, inside the logged 2^1.82-2^13.82 span plus the flat cap above it.
The line implementation is prime-agnostic mod-p arithmetic; AArch64 is
64-bit-native, so the wider prime executes the same instruction count (H4). No
line prime reaches 2^63. No line records; survivor factors in (2^20, 2^47]
come from the band strip, below 2^20 from trial division.
