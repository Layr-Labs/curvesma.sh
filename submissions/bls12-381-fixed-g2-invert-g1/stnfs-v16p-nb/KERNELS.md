# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence; `measurements12.txt`, `wide12.c`, `wide12.txt`, `measurements-oh.txt`, `measurements12x.txt`, `rec7.patch`, `slice_mb.h`, `rec7.txt` and `slice12oh.c` are sections of
`kernels_asm.b85` (base85 of lzma text). The large raw `dsieve16` log (`measurements16.txt`) is omitted from
this archive for size; its counts and checks are summarized in the table below.

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

At p = 2^69 the per-hit cost falls with R (417.8, 370.0, 353.9 Ir at R = 160, 192, 224): the 1000 charge covers R <= 384. A complete (p, root pair) at p = 2^69, R = 96: 47,853.7 Ir per hit including
two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold the index tables (2^34.34 bytes at
R = 384); `apply` runs at R = 4; the bench's `realloc` record buffer is harness only, production streams records to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Per-update costs by hits per coset: line path (p in [2^14, 2^48); 17.2-352.1, `measurements12x.txt`, below), record walk
(44.67-1711.20), the rec7 run walk, apply 5.76, recovery 5.76 per record and scan2 (3.002 per cell + 90.98 per leaf; below); layout 28.17 per suffix-box entry (`measurements-oh.txt`); certified norm bounds 1,508.2 per
evaluation, 4,600 per row (4,300 + 300 box-hull share), 17,500 per slice (`embedding_bounds.py`, `nbound31`); Montgomery
multiplication 367.76 Ir; the LA kernels of STORAGE.md (counted by [10]). These are the constants of `model16.py`.

### Slice sieve per-update costs (fused exact recentring)

Coset centres are recentred exactly (adjugate remainder in __int128, per-pair affine forms); all `slice12` checks
MATCH. Record Ir per hit 44.67 / 77.42 / 250.54 / 810.82 / 1711.20 (p = 257, 65537, 8388617, 67108879, 268435459),
setup6 406.87 / 374.07: `REC_TAB` (exact quotients of the `measurements12.txt` log).

### Primes in [2^63, 2^64): the WIDE1 path

Only slice12's `long long` basis (p > 2^63 - 1) and the k = 64 Barrett reciprocal break. `wide12.c`, appended to
slice12 0bbca0b2 (the code base of `measurements12x.txt`) in its own text section, adds `prepare_pair_w` (LLL with
__int128 entries, abort unless all < 2^40, then the unchanged `lll()`; unreachable, Minkowski puts the last minimum
below 2^28) and `setup6w` (Horner with a 2-by-1 division by the invariant p). `wide12.txt` (same container and flags):
slice12's .text, main and relocations BYTE-IDENTICAL; 10^6 reductions exact; 200 wide pairs at three slice sizes:
lattices certified (congruence, |det| = p), 400/400 planted hits, digests, survivors and records MATCH brute force;
setup6w 141.00 Ir; record path per coset at three wide p within 2.17 Ir of p = 2^63 - 25. The walk is 0bbca0b2's,
so its per-hit price (`REC_TAB_W`, the `measurements12x.txt` TAB rec rows) x 1.01 (`WIDE_REC`, 2.5x the measured
spread) is charged; setup6 keeps its narrow charge (`WIDE_SETUP` = 1); empty cosets keep 1,600 per pair (2.6x the
measured 613 Ir). The wide lattices (6 pi(2^64) (p, s) at 2 x 5,615,924 Ir) are the ledger row 2^70.6; [10w]
recomputes these prices from the logs.

### Slice layout per suffix-box entry (`measurements-oh.txt`)

The `build` bench (27.27 Ir per entry; 12.325 per cell only at its 2^26-cell slice) is layout per suffix-box entry.
Per-function callgrind at 2^20..2^26 cells: `build_layout` 28.16 / 27.21 / 26.67 / 25.86 per entry, `zero_cells`
0.625 per cell (a separate pass in the public slice12; scan2 below fuses it), remainder < 2^20 Ir per run. Entries per cell K L4 / S: 0.452 at the bench, 0.0037-0.0049 in
production. Charged: 28.17 per entry (K and each L4 count raised past rounding, `entries_upper`), 2^21 per slice,
2^17 per 2^32-byte segment; [10o] checks the log.

### Line path: the public slice12 line walk (`measurements12x.txt`)

Same container, flags and method: (i) base table on the 2^26 slice, one (p, s) per point: 24.098 / 31.113 / 88.313 /
310.700 at h = 2^13.82 / 2^9.82 / 2^5.82 / 2^1.82; (ii) random-pair groups (4 pairs per seed, prime q in [p, p + p/8],
random cosets over all lines of a 2^26-2^32 slice, 3-8 seeds per group). Each `LINE_TAB` entry is the larger of (i)
and the pooled cost + 2 SE of the group at or below its h; h >= 2^18.81 is capped at 17.2 (cap groups pooled 17.04 /
16.62 / 16.29); linear in log2 h between entries. [10l] recomputes every quotient from raw Ir and hits. The walk is the public code.

### Record walk per canonical a-run: `rec7.patch`, `rec7.txt`, `runs16.py`

rec7 extends slice12's Schnorr-Euchner record walk to L7(P) = {(b, t) : t = r b(s) mod p}: a fourth affine state t and
slab pruning (a Cauchy-Schwarz bound on the t-range of every completion at each level; exact clip of the innermost run);
each record goes to the slice of its t (a table lookup inside the measured walk). The measured source is the previous
build's `rec7.c` (1b8fa720) plus `rec7.patch`;
`slice_mb.h` is the benchmark's b-block metric. For p >= 2^62 the basis comes from `prepare7w` (LLL on __int128 entries
with quad-precision GSO, abort unless every entry < 2^40 (unreachable: reduced entries are ~p^(1/7)), every row checked against t = r b(s) mod p, |det| = p, then
prepare7's tail); every walk recentres c = c_hi 2^32 + c_lo through a stored short representative of 2^32 e_6, so all
recentring arithmetic stays below 2^53. Checks (`-DCAPTURE`): walk7's record sets equal a brute-force scan on 270
cosets (slices 2^9-2^11, T = 1..33, five slab offsets, p 97..4099, stretches 0/16/64/660): 0 mismatches, 2,087
records; WCHECK, 432 pairs with p in [2^62, 2^64) (half above 2^63), 1,728 planted region points: all found, every
record in its coset and region. Campaign (`rec7.txt`: geometry, log2 S, T, slab, p, stretch, lattice, Ir and hits of
two runs, h'): the 130 CAMP points (slices 2^60-2^66; T = 32, 64, 128 at stretch 660 and 181, 362, 1,024 at 90; centre
and off-centre slabs; two lattices) re-measured, plus 54 CAMPW points with 2^63 < p < 2^64 (T = 32..1,024, slices up
to 2^68, h' = 2^-1.0..2^14.5; 42 at h' >= 2^4.99). The model uses rec7 for every slice of a clean run of 32 to 1,024
a-values (`RUN_T`; stretch 660 up to 128 a-values, 90 above) in every band p in [2^56, 2^64) with h' >= 2^4.99 (one
implementation per run and band; contamination charged RUN_ALLOW x D_SUP per slice, proof.md section 3) and reads
`REC7_TAB` as a step function (the entry at or below log2 h'). Its entries are the bucket maxima over these 184 points
and the 130 points of the previous build; every point with h' >= 2^4.99 lies at or below its entry ([10r]). Margin: h'
counts full-ellipsoid hits, 2 x run cells / p, while the model feeds half-slice run cells / p from `runs16.lrs`. Each
rec7 record pays 2 more (`EXTRA7`); routing stays 10.01 per record; setup6 is paid once per (run, ideal). Bases cost at
most 3,995,157 (prepare7, p < 2^62) and 150,271,737 (prepare7w) instructions per ideal (PREP rows); two bases for each
of the 2^59.56 ideals with p in [2^52, 2^64), all at the larger price, x 2^0.05: 2^73.49 <= 2^73.5 units (ledger row).
Runs (`runs16.py`): per listed unit eps (all u1^i u2^j u3^k of log-radius <= 9: 11,325; the closest included and
excluded units lie at least 0.0013 and 0.0011 from the radius) the exact integer quadratic in t of Q_a(eps (a + t e_0))
- Q_a(a + t e_0); float roots locate the ends, exact integer tests fix them; a run has at most 2101 a-values (Q_a(e_0)
= 0.5528759885).

### Scan, apply, recover: `slice12oh.c`

The public slice12 with three functions replaced: `scan_segment` tests four 16-bit cells per word, ((w | H) - T) & H,
eight cells per step, then four, then single cells, zeroing both sum arrays with its own stores; `apply_batch` is
unrolled x4; `recover_segment` tests eight records per step. The lane test is exact: every sum and threshold lies in
[0, 2^15) (a cell sum is at most 16 x 851). Brute-force checks MATCH at 2^18-2^24 cells, segment shifts 0-3.
Differentials (`measurements-oh.txt`, same image and flags): apply 5.751 and recovery 5.751 per record (p = 257),
charged 5.76 (the unmodified functions: 7.000). Scan cost per leaf depends on the leaf length L (a fixed part, 24 per
8-cell step, 22 for a 4-cell step, 14 per single cell): measured at every L = 1..64 and 256 (65 points), charged 3.002
per cell + 90.98 per leaf, the smallest leaf charge covering every point at 3.002 per cell (L = 15 is the worst); [10s].
