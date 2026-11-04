# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit); correctness
checks logged alongside. C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) described here; shipped logs are the
evidence; `measurements12.txt` is the last section of `kernels_asm.b85` (base85 lzma
text). The raw `dsieve16` log (`measurements16.txt`) is omitted
for size; counts and checks summarized below.

## Descent sieve: `dsieve16`

For prime p < 2^70 and root-pair functional w (12 residues mod p) the kernel
builds the hit lattice {c : w.c = 0 mod p} (det p), LLL-reduces it
(long-double decisions, exact __int128 updates), certifies exactly
(Bareiss det = p),
enumerates ||c||^2 <= R^2 by Schnorr-Euchner in double precision with
exact ball test per leaf, sign-canonicalises, ranks
by ball index, emits the record
(u64 owner = index >> 40, u64 offset = index mod 2^40, u16 16 log2 p). The
largest shipped ball (R = 384) has half count 2^102.44, owners below
2^62.44; `tables` refuses half count >= 2^103.

Correctness (`check`, all PASS): set equality with brute
force + rank bijection for R <= 4 (p 307 to 2^40); at p = 2^20..2^69,
R = 10..224 every record unranked and verified
(congruence, ball, sign, no duplicates), hit counts match volume/p (all
18,051,918 at R = 224); at R = 256 the first 2^24 records, indices
to 2^95.4. Split is a pure function of the index,
audited per sieved node ([7], store16.check).

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

At p = 2^69 per-hit cost falls with radius (417.8, 370.0, 353.9 Ir
at R = 160, 192, 224); the 1000 charge is 2.8x the R = 224 value, covers R to 384.
A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per hit
with two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold
index tables (2^34.34 bytes at R = 384, STORAGE.md); `apply` at R = 4; records stream to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Line and record paths of the relation sieve (per-update costs 24.10-298.78 and
44.67-1711.20 by hits per coset, table below; line path for p < 2^48, record walk above; apply 7.00, recovery 7.003 per record; layout
12.325, zero 0.625, scan 8 per cell), certified norm bounds (1,508.2 per
evaluation, 4,300 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`; the
1,508.2 is the measured marginal per-block cost, charged exactly), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the LA kernels of STORAGE.md (group loops 11/33, the row tail and reductions of `row_tail.c`, and the 97-instruction multiply into a 9-limb accumulator `bw_full`; listings in `kernels_asm.b85`, counted by [10]): the constants of `model16.py`.

### Slice sieve per-update costs (re-measured, fused exact recentring; 2c2ce2ab)

Coset centres are recentred from the exact adjugate remainder (a_i c - q_i p in
__int128) with per-pair affine forms; all `slice12` checks (line, record, apply,
scan, recover vs brute force, primes up to 2^63) MATCH. Callgrind Ir per unit:
line 24.10 / 31.07 / 87.57 / 298.78 (p = 257, 4099, 65537, 1048583); record
44.67 / 77.42 / 250.54 / 810.82 / 1711.20 (p = 257, 65537, 8388617, 67108879,
268435459); setup6 406.87 / 374.07. These are `LINE_TAB` / `REC_TAB`; the log is
the `measurements12.txt` section of `kernels_asm.b85`.
