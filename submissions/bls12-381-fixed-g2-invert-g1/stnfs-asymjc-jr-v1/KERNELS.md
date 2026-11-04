# Kernels and measurements

Kernel counts are callgrind Ir (record route: qemu instruction count,
STORAGE.md) on linux/arm64 (Debian 12, GCC 12.2 -O3, valgrind 3.19) by the
differential method, with correctness checks logged alongside. The C kernels
(`dsieve16`, `slice12`, `ecm12`, `nbound31`) are described here (H4);
`measurements12.txt` is the last section of `kernels_asm.b85`; the raw log and
C sources are not shipped. Descent charges (SIEVE_HIT 850, BOUND_CELL 2,000)
are the f34f63d9 line's (R8+, same measurements); main-phase constants are
42e80f80's (fused exact recentring).

## Descent sieve: `dsieve16`

For a prime p < 2^70 and a root-pair functional w (12 residues mod p) the kernel
builds the hit lattice {c : w.c = 0 mod p} (det p), LLL-reduces it
(long-double decisions, exact __int128 updates) and certifies it exactly
(congruence per basis vector; Bareiss det = p). It enumerates ||c||^2 <= R^2
by Schnorr-Euchner in double precision with an exact integer ball test per
leaf, sign-canonicalises each hit, ranks it by exact lexicographic ball index
(__int128) and emits (owner = index >> 40, offset = index mod 2^40, 16 log2 p).
The largest shipped ball (R = 384) has half count 2^102.44, so owners stay below
2^62.44; `tables` refuses half count >= 2^103.

Correctness (`check`, all PASS on the shipped source): brute-force set equality
and a rank bijection for R <= 4; every record unranked and verified at p =
2^20..2^69, R = 10..224 (18,051,918 records at R = 224) and the first 2^24 at
R = 256; the split is audited for every sieved node ([7], store16.check).

| measured | Ir | charged |
|---|---:|---:|
| enumeration + rank + record, per hit (p = 2^20..2^69, R = 12..224) | 302-793 | 800 (plan R >= 192, measured <= 418 there) |
| transfer to the owner, per hit (four radix passes) | (STORAGE.md) | 84 |
| apply, per record | 9.00 | 9 |
| SIEVE_HIT total | | 850 |
| scan, per cell | 5.000 | 10 |
| zero, per cell | | 2 |
| certified per-cell norm bound, both sides (DESCENT.md) | ~1,500 | 2,000 |
| LLL + certificate per (p, root pair), p = 2^69 | 2^24.98 | 2^25.1 |

At p = 2^69 per-hit cost falls with radius (417.8/370.0/353.9 Ir at R =
160/192/224); every shipped node sieves at R >= 192 ([7] prints radii), where
the max measured is 418; the 800 charge is 1.9x that and 2.3x the R = 224
value, covering R up to 384. A complete (p, root pair) at p = 2^69, R = 96:
2^24.99/setup, below PAIR_SETUP. The certified per-cell norm bound
(`descent_witness.cell_bound`, 84 linear forms, outward rounding, checked by
[11b]) is ~1,500 instructions/cell, charged 2,000; form roots have certified
enclosures ([7r]).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Line and record paths of the relation sieve (per-update 24.10-298.78 and
44.67-1711.20 by hits/coset, fused exact recentring; line path p < 2^48;
apply 7.00, recovery 7.003/record; layout 12.325, zero 0.625, scan 8/cell),
certified norm bounds (1,508.2/eval, 4,300/row, 17,500/slice;
`embedding_bounds.py`, `nbound31`), five-limb Montgomery multiply (367.76 Ir,
20,000 products checked) and the LA kernels of STORAGE.md (loops 11/33, row
tail and reductions (C source in 42e80f80; `la_mac.s` shipped),
97-instruction `bw_full`; counted by [10]): the `model16.py` constants. Line
h = S/(K p) >= 2^6.79 at every row/band for p < 2^48, inside the measured
range (flat charge above; secant interpolation overestimates the
convex-decreasing cost in log h). The band strip uses the `ecm16.py` schedule
at (48 bits, 832/576 bits, 99.9%), two rounds/side; the main 69.25-bit (f) and
69-bit (g) schedules at 99.9% ([12] prints them); the g side runs only after
the f side has factored completely (early abort, `model16.EA_FACTOR`).
