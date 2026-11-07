# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3 -fno-tree-vectorize
-fno-tree-loop-distribute-patterns, valgrind 3.19), differential (two runs differing only in the measured unit), checks logged beside counts. Logs: `measurements12.txt`, `wide12.txt`,
`measurements-oh.txt`, `measurements12x.txt`, `rec7.patch`, `rec7.txt`, `slice12oh.c` and `toom3.txt`
are sections of `kernels_asm.b85` (base85 of lzma text); `row_tail.c` (sha256 5f2e0ad5...), `wide12.c` (f7013d87...) and
`slice_mb.h` (57e47e56..., rec7's b-block metric: `config.json` matrix_b in full double precision) are 1b8fa720's, not repeated; the `dsieve16` log is summarized.

## Descent sieve: `dsieve16`

Builds the determinant-p hit lattice of a root-pair functional (p < 2^70), LLL-reduces and certifies it (congruences,
Bareiss determinant = p), enumerates ||c||^2 <= R^2 (Schnorr-Euchner, exact leaf test), ranks each hit (exact
lexicographic ball index) and emits (owner = index >> 40, offset, 16 log2 p) records; owners stay below 2^62.44 at R = 384 (`tables` refuses half counts >= 2^103).
`check` passed: set equality and rank bijection with brute force for R <= 4;
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

At p = 2^69 the per-hit cost falls with R (417.8, 370.0, 353.9 Ir at R = 160, 192, 224): 1000 covers R <= 384; a
(p, root pair) at R = 96: 47,853.7 Ir per hit, two setups over 1,396 hits (2^24.99 per setup < PAIR_SETUP). Index tables take 2^34.34 bytes at R = 384; `apply` measured at R = 4; the bench's `realloc` buffer is harness only (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Per-update costs by hits per coset: line path (p in [2^14, 2^48); 17.2-352.1, `measurements12x.txt`, below), record walk
(44.67-1711.20), the rec7 run walk, apply 5.76, recovery 5.76 per record and scan2 (3.002 per cell + 90.98 per leaf; below); layout 28.17 per suffix-box entry (`measurements-oh.txt`); certified norm bounds 1,508.2 per
evaluation, 4,600 per row (4,300 + 300 box-hull share), 17,500 per slice (`embedding_bounds.py`, `nbound31`); Montgomery
multiplication 367.76 Ir; the LA kernels of STORAGE.md (counted by [10]). (`model16.py` constants).

### Slice sieve per-update costs (fused exact recentring)

Coset centres are recentred exactly (adjugate remainder in __int128, per-pair affine forms); all `slice12` checks
MATCH. Record Ir per hit 44.67 / 77.42 / 250.54 / 810.82 / 1711.20 (p = 257, 65537, 8388617, 67108879, 268435459),
setup6 406.87 / 374.07: `REC_TAB` (exact quotients of the `measurements12.txt` log).

### Primes in [2^63, 2^64) off runs: the WIDE1 path

Only slice12's `long long` basis (p > 2^63 - 1) and the k = 64 Barrett reciprocal break. `wide12.c` (own text section of
slice12 0bbca0b2) adds `prepare_pair_w` (LLL on __int128 entries,
abort unless all < 2^40, then the unchanged `lll()`; for this 6-dim lattice Minkowski puts the last minimum below 2^28)
and `setup6w` (Horner with a 2-by-1 division by p). `wide12.txt`: slice12's .text, main and relocations
BYTE-IDENTICAL; 10^6 reductions exact; 200 wide pairs at three slice sizes: lattices certified, 400/400 planted hits,
records MATCH brute force; setup6w 141.00 Ir; per coset within 2.17 Ir of p = 2^63 - 25. So the per-hit price
(`REC_TAB_W`) x 1.01 (`WIDE_REC`, 2.5x the measured spread) is charged, setup6 at its narrow price, empty cosets 1,600
per pair (2.6x the measured 613). The wide lattices (6 pi(2^64) (p, s) at 2 x 5,615,924 Ir) are the 2^70.6 row; [10w]
recomputes these from the logs.

### Slice layout per suffix-box entry (`measurements-oh.txt`)

`build` bench 27.27 Ir per entry; per-function callgrind at 2^20..2^26 cells:
`build_layout` 28.16 / 27.21 / 26.67 / 25.86 per entry, `zero_cells` 0.625 per cell (fused into scan2), remainder
< 2^20 Ir per run. Entries per cell K L4 / S: 0.452 at the bench, 0.0037-0.0049 in production. Charged 28.17 per entry
(K and each L4 count raised past rounding, `entries_upper`), 2^21 per slice, 2^17 per 2^32-byte segment; [10o].

### Line path: the public slice12 line walk (`measurements12x.txt`)

(i) Base table on the 2^26 slice, one (p, s) per point: 24.098 / 31.113 / 88.313 / 310.700 at h = 2^13.82 / 2^9.82 /
2^5.82 / 2^1.82; (ii) random-pair groups (4 pairs per seed, prime q in [p, p + p/8], random cosets over a
2^26-2^32 slice, 3-8 seeds). Each `LINE_TAB` entry is the larger of (i) and the group's pooled cost + 2 SE at or below its h; h >= 2^18.81 capped at 17.2 (cap groups 17.04 / 16.62 / 16.29) linear in log2 h between entries; [10l].

### Record walk per canonical a-run: `rec7.patch`, `rec7.txt`, `runs16.py`

rec7 extends slice12's record walk to L7(P) = {(b, t) : t = r b(s) mod p} (fourth affine state t; slab pruning by a
Cauchy-Schwarz bound on the t-range of every completion). Source: 1b8fa720's `rec7.c` (sha256 in the patch header) +
`rec7.patch`. For p >= 2^62 bases come from `prepare7w` (__int128 LLL, quad-precision GSO; rows checked against t = r
b(s) mod p; |det| = p; aborts at an entry >= 2^40, as none of the 734 logged ones with p >= 2^62 did; an aborted ideal
takes the slice walk). Walks recentre c = c_hi 2^32 + c_lo through a stored
representative R32 of 2^32 e_6 (float error only moves the centre); c_hi R32 fits 63 bits for reduced entries below
2^28, as in every check and campaign run. Checks (`-DCAPTURE`): brute force on 270 cosets (slices 2^9-2^11, T = 1..33, five slabs, p 97..4099,
stretches 0/16/64/660): 0 mismatches, 2,087 records; WCHECK, 432 pairs with p in [2^62, 2^64) (half above 2^63),
1,728 planted points: all found, every record in its coset and region. Campaign: the
130 CAMP points (p 2^44-2^62.86, 8 on prepare7w; slices 2^60-2^66; T = 32, 64, 128 at stretch 660 and 181, 362, 1,024
at 90; two slabs, two lattices) re-measured, and 42 CAMPW points with 2^63 < p < 2^64 (T = 32..1,024, slices to 2^68,
h' = 2^5.4..2^14.5; 12 below 2^4.99, off the table, not shipped). The model uses rec7 on every slice of a clean run of
32 to 1,024 a-values (`RUN_T`; stretch 660 up to 128, 90 above) in every band of [2^56, 2^64) with h' >= 2^4.99,
reading `REC7_TAB` (bucket maxima over these 184 points and the previous build's 130) as a step function; every point
lies at or below its entry ([10r]). Margin: h' counts full-ellipsoid hits (2 x run cells / p), the model half-slice
cells. Each record pays 2 more (`EXTRA7`), routing stays 10.01, setup6 once per (run, ideal). Bases: at most 3,995,157
(prepare7) and 150,271,737 (prepare7w) instructions; two for each of the 2^59.56 ideals of [2^52, 2^64) (from 2^52:
safe), all at the larger price, x 2^0.05: 2^73.49 <= 2^73.5 units. Runs (`runs16.py`): per listed unit eps (all
u1^i u2^j u3^k of log-radius <= 9: 11,325; the closest included and excluded units lie >= 0.0013 and 0.0011 from the
radius) the exact quadratic in t of Q_a(eps (a + t e_0)) - Q_a(a + t e_0); float roots locate the ends, exact tests
fix them; a run has at most 2101 a-values (Q_a(e_0) = 0.5528759885).

### Batch band strip: Toom-3 (`toom3.txt`)

`toom3` (points 0, 1, -1, 2, inf; Bodrato interpolation in two's complement; schoolbook below 48 limbs) matches
schoolbook on 108 cases (n 47..2000) and Python at n = 20000. Differential Ir at 28 sizes 8..34,992 limbs stay at or
below the recurrence of `model16.batch_cost` (worst 0.9857, flat 0.856-0.896 over the last five levels), the code's own
recursion. `SRT` (C build, guard 64 bits per level = G = L + 2 limbs; exact by proof.md's lemma, as `model16.srt`) returns the exact band parts
of 256 planted survivors against trial division by the 995,846 primes of (2^20, 2^24]. Per wave: 2^39 13-limb leaves, P 2^42.53 limbs.

### Scan, apply, recover: `slice12oh.c`

The public slice12 with three functions replaced: `scan_segment` tests four 16-bit cells per word, ((w | H) - T) & H,
eight cells per step, then four, then single cells, zeroing both sum arrays with its own stores; `apply_batch` is
unrolled x4; `recover_segment` tests eight records per step. The lane test is exact: sums and thresholds lie in
[0, 2^15) (a cell sum is <= 16 x 851). Brute-force checks MATCH at 2^18-2^24 cells, segment shifts 0-3.
Differentials (`measurements-oh.txt`): apply 5.751 and recovery 5.751 per record (p = 257),
charged 5.76. Scan cost per leaf depends on the leaf length L: measured at every L = 1..64 and 256 (65 points), charged 3.002
per cell + 90.98 per leaf, the smallest leaf charge covering every point at 3.002 per cell (L = 15 is the worst); [10s].
