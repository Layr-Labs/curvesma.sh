# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence: `measurements12x.txt`, `wide12.c`, `wide12.txt`, `rec7.txt` and `measurements-oh.txt` are sections of
`kernels_asm.b85` (base85 of lzma text). Two sections are omitted for size and shipped byte-identically in
v9 734a9723's archive: `measurements12.txt` (the counts below) and `row_tail.c`. The large raw `dsieve16` log
(`measurements16.txt`) is likewise summarized below.

## Descent sieve: `dsieve16`

Builds the determinant-p hit lattice of a root-pair functional (p < 2^70), LLL-reduces
and certifies it exactly (congruences, Bareiss determinant = p), enumerates the ball
||c||^2 <= R^2 (Schnorr-Euchner, exact leaf test), sign-canonicalises and ranks each
hit (exact lexicographic ball index) and emits (owner = index >> 40, offset, 16 log2 p)
records; owners stay below 2^62.44 at R = 384 (`tables` refuses half counts >= 2^103).
`check`: brute-force set equality and rank bijection at R <= 4; every record unranked/verified at
p = 2^20..2^69, R = 10..224 (18,051,918 at R = 224), first 2^24 at R = 256; split audited per node ([7]).

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

Per-hit cost falls with radius at p = 2^69 (417.8 / 370.0 / 353.9 Ir at R = 160 / 192 / 224; rank = 12 lookups);
the 1000 charge is 2.8x the R = 224 value, covering R <= 384. A complete pair at R = 96 measures 47,853.7 Ir
a hit incl. two setups (2^24.99 each, below PAIR_SETUP). Index tables: 2^34.34 bytes at R = 384 (STORAGE.md);
`apply` at R = 4; records stream to owners in production.

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings; logs `measurements12x.txt`, `wide12.txt`, `measurements-oh.txt`)

Line and record paths of the relation sieve (line path for p < 2^48, record walk above; per-update costs by h, the points
per coset). Overhead kernels (v12oh, slice12oh.c, measurements-oh.txt): scan2 (SWAR, fused zeroing; sums/thresholds < 2^15) at 3.002 + 48.05 leaves/blen; apply/recovery 5.751 measured, 5.76 charged; slice layout per suffix-box entry (below); certified norm bounds
(1,508.2 per evaluation, 4,600 per row, 17,500 per slice; `embedding_bounds.py`, `nbound31`); five-limb Montgomery
multiplication 367.76 Ir; the LA kernels of STORAGE.md (listings counted by [10]): the constants of `model16.py`.

**Record walk below 2^63 (`slice12g`, 2c2ce2ab).** Exact adjugate recentring; all checks MATCH brute force. Record
44.67 / 77.42 / 250.54 / 810.82 / 1711.20 Ir per hit (p = 257 ... 268435459), setup6 406.87 / 374.07: `REC_TAB`.

**Line path (v9): the public `slice12` line walk (0bbca0b2).** `measurements12x.txt` (same container, flags and
differential method) holds (i) its base table on the 2^26 slice, one (p, s) per point: 24.098 / 31.113 / 88.313 / 310.700
at h = 2^13.82 / 2^9.82 / 2^5.82 / 2^1.82; (ii) random-pair groups: each seed draws 4 pairs (prime q in [p, p + p/8], s
uniform) and walks random cosets over all lines of a 2^26-2^32 slice; 3-8 seeds per group. Single pairs scatter (at
2^9.8 random pairs average 35.85 against the base point 31.11), so each entry of `LINE_TAB` is the larger of (i) and the
pooled cost + 2 SE of the group at or below its h, and h >= 2^18.81 is capped at **17.2** >= the bounds of the
three cap groups (p ~ 26 / 84 / 257 on 2^28 / 2^30 / 2^32 slices, pooled 17.04 / 16.62 / 16.29,
falling with p); linear in log2 h between entries (costs are convex in log2 h). [k] recomputes every quotient from raw Ir and hits and checks the table against every
group. The bench diffs (line range, multi-pair mode) are in v9 734a9723's log.

**Primes in [2^63, 2^64) (v9; v8 charged them twice).** Only the long long basis and the k = 64 Barrett reciprocal
break. `wide12.c`, appended to the public slice12 in its own text section (its object code unchanged), adds
`prepare_pair_w` (int128 LLL, abort unless all entries < 2^40, then the unchanged `lll()`; unreachable, since a nonzero
x with sum x_i s^i = 0 mod p has ||x|| > 2^9.5, so Minkowski bounds the last minimum below 2^28) and `setup6w`
(Moller-Granlund reciprocal). `wide12.txt`: reductions exact; 200 wide pairs at three slice sizes certified (congruence,
|det| = p), planted hits, digests, survivors and records MATCH brute force. setup6w 141.00 Ir (within the setup6
charge); the record walk per coset at three wide p is within 2.17 Ir of 611.17 at 2^63 - 25 (same code), so wide primes
pay the public record table `REC_TAB_W` (44.901 / 78.211 / 266.972 / 938.085 / 2231.664) x 1.01; empty cosets and apply x 1.01.
Wide factor-base lattices (<= 5,615,924 Ir per (p, s)) are a 2^70.6 ledger row.

**Slice layout per suffix-box entry (v10; `measurements-oh.txt`, our overhead lane).** v8/v9 charged the `build` bench's
12.325 Ir per cell, which holds only at its 2^26-cell slice: the work is layout per suffix-box entry. Per-function
callgrind of `slice12 bench build` at 2^20, 2^22, 2^24, 2^26 cells (same container): `build_layout` 28.16 / 27.21 / 26.67 / 25.86 Ir per
entry, `zero_cells` 0.625 per cell, remainder < 2^20 Ir per run. Charged: **28.17** per entry (`C_BUILD_E`; K + 2 lines and
2 floor(rb_j s) + 3 per coordinate, `entries_upper`, so every count is raised past rounding), 2^21 per slice and 2^17 per
2^32-byte segment; [10o] checks the log. Production rows have 0.0039 entries per cell (cell-weighted mean over pools
B, C, D), so the layout costs 0.109 Ir per cell instead of 12.325.

**Record walk per canonical a-run (v11; `rec7.c`, `rec7.txt`; our BIGLEVER-1500 lane).** For p in [2^52, 2^63) the
canonical a-values a + t e_0 of one run (`runs16.py`: the maximal t-interval around a on which no listed unit associate is
smaller in (Q_a, lex) and Q_a stays in the region; 11,325 units of log-radius <= 9) share one 7-dimensional lattice per
ideal; `walk7` emits the records of all the run's slices in one pass and `prepare7` keeps two reduced bases per ideal (short
and long runs). [k7]: `rec7.txt` CHECK, record sets equal brute force: "check (record sets, folded, b5 = 0 skipped): ok (0 mismatches over 270 cosets, 2087 records)"; 130 differential callgrind points
(same container; slices 2^60-2^66, runs of 32-1024 a-values, central and off-centre slabs, two lattices) by h' = hits per
run-coset; `REC7_TAB` is the maximum per 1-bit level of log2 h', made non-increasing, read as a step (the level at or below
h'; worst ratio 1.0000). Charged per record: step7(log2 h + log2 g) + 2 Ir, g = run cells / row cells; setup6 once per
(run, ideal); rows whose h' falls below 2^4.99 keep `REC_TAB`. prepare7 <= 3,986,661 Ir per ideal (`rec7.txt` PREP,
p = 2^52..2^62): two bases for 2^59.56 ideals, ledger row 2^68.3.
