# Kernels and measurements

Kernel counts are callgrind Ir (the record route is a qemu executed-instruction
count, STORAGE.md) on linux/arm64 (Debian 12, GCC 12.2 -O3, valgrind 3.19) by
the differential method (two runs differing only in the measured unit), with
each kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here, with counts summarised below (H4);
`measurements12.txt` is the last section of `kernels_asm.b85` (base85 of lzma text);
the raw `dsieve16` log and the C sources are not shipped. The descent charges
below (SIEVE_HIT 850, BOUND_CELL 2,000) are those of the f34f63d9 line (R8 and
later, same kernel measurements); the main-phase constants are those of
42e80f80 (fused exact recentring).

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

At p = 2^69 the per-hit cost decreases with the radius (417.8, 370.0, 353.9 Ir
at R = 160, 192, 224: fewer enumeration nodes per hit); every shipped node sieves at R >= 192 (the plan search is
restricted to R >= 192, [7] prints the radii), where the largest measured value
is 418; the 800 charge is 1.9x that value and 2.3x the R = 224 value and covers
R up to 384.
A complete (p, root pair) at p = 2^69, R = 96: 2^24.99 per setup, below PAIR_SETUP.
The certified per-cell norm bound
(`descent_witness.cell_bound`: 84 linear forms in double precision with outward
rounding, [11b] checks it against exact norms) is about 1,500 instructions per
cell (72 + 12 complex evaluations of degree-5 forms with the error term), charged
2,000; the roots entering the forms have certified enclosures ([7r]).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Line and record paths of the relation sieve (per-update costs 24.10-298.78 and
44.67-1711.20 by hits per coset, fused exact recentring; line path for p < 2^48;
apply 7.00, recovery 7.003 per record; layout 12.325, zero 0.625, scan 8 per
cell), certified norm bounds (1,508.2 per evaluation, 4,300 per row, 17,500 per
slice; `embedding_bounds.py`, `nbound31`), the five-limb Montgomery
multiplication (367.76 Ir, checked on 20,000 products) and the LA kernels of
STORAGE.md (group loops 11/33, the row tail and reductions (C source in 42e80f80;
the compiled `la_mac.s` shipped), the 97-instruction `bw_full`; counted by [10]): the constants of `model16.py`. Line h = S/(K p) >=
2^6.79 at every row and band for p < 2^48, inside the measured line-table range
(flat charge above it; the secant interpolation overestimates the
convex-decreasing cost in log h). The survivor band strip uses the `ecm16.py`
schedule at (48 bits, 832/576 bits, 99.9%), two rounds per side; the main
69-bit schedule at 99.9% ([12] prints both).
