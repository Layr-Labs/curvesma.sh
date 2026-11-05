# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence; `measurements12.txt`, `wide12.c`, `wide12.txt`, `measurements-oh.txt` and `measurements12x.txt` are sections of
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

At p = 2^69 the per-hit cost falls with R (417.8, 370.0, 353.9 Ir at R = 160, 192, 224); the 1000
charge is 2.8x the R = 224 value and covers R up to 384.
A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per hit
including two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold the
index tables (2^34.34 bytes at R = 384, STORAGE.md); `apply` runs at R = 4; the
`realloc` record buffer belongs to the harness, production streams records to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings and the `measurements12.txt` log)

Line and record paths of the relation sieve (per-update costs 17.2-352.1 (line, `measurements12x.txt`, below) and
44.67-1711.20 (record) by hits per coset; line path for p < 2^48, record walk above (v3 band strip); apply 7.00, recovery 7.003 per record; layout
28.17 per suffix-box entry (`measurements-oh.txt`, below), zero 0.625 per cell, scan 8 per cell plus 21 per leaf), certified norm bounds (1,508.2 per
evaluation, 4,600 per row (4,300 + 300 box-hull share), 17,500 per slice; `embedding_bounds.py`, `nbound31`; the
1,508.2 is the measured marginal per-block cost, charged exactly in this revision), the
five-limb Montgomery multiplication (367.76 Ir, checked on 20,000 products) and
the LA kernels of STORAGE.md (group loops 11/33, the v6 five-vector Krylov loops 39/33 of `krylov5.s`, the row tail and reductions of `row_tail.c`, and the 97-instruction multiply into a 9-limb accumulator `bw_full`; listings in `kernels_asm.b85`, counted by [10]): the constants of `model16.py`.

### Slice sieve per-update costs (re-measured, fused exact recentring; 2c2ce2ab)

Coset centres are recentred from the exact adjugate remainder (a_i c - q_i p in
__int128) with per-pair affine forms; all `slice12` checks (line, record, apply,
scan, recover vs brute force, primes up to 2^63) MATCH. Callgrind Ir per unit: record
44.67 / 77.42 / 250.54 / 810.82 / 1711.20 (p = 257, 65537, 8388617, 67108879,
268435459); setup6 406.87 / 374.07. These are `REC_TAB`; the log is
the `measurements12.txt` section of `kernels_asm.b85`.

### Primes in [2^63, 2^64): the WIDE1 path

Residues fit one word; only slice12's `long long` basis (p > 2^63 - 1) and the k = 64 Barrett
reciprocal break. `wide12.c` is appended to slice12 0bbca0b2 (the code base of `measurements12x.txt`)
in its own text section: `prepare_pair_w` (per (p, s): LLL with __int128 entries, abort unless all
< 2^40, then the unchanged `lll()`; unreachable, Minkowski puts the last minimum below 2^28) and
`setup6w` (Horner with a 2-by-1 division by the invariant p). `wide12.txt` (same container and
flags): slice12's .text, main and relocations BYTE-IDENTICAL with the code appended; 10^6
reductions exact; 200 wide pairs at three slice sizes: lattices certified (congruence, |det| = p),
400/400 planted hits, digests, survivors and records MATCH brute force; setup6w 141.00 Ir against
setup6 375.07/407.87; record path per coset at three wide p within 2.17 Ir of p = 2^63 - 25.
After the basis preparation the walk is the 0bbca0b2 code, so its per-hit price by hits per coset
(`measurements12x.txt` TAB rec rows, `REC_TAB_W`, above the recentred slice12g prices of `REC_TAB`) is charged
for these primes, x 1.01 (`WIDE_REC`, 2.5x the measured spread); setup6 keeps its narrow charge
(`WIDE_SETUP` = 1); empty cosets keep 1,600 per pair (2.6x the measured 613 Ir). The wide
factor-base lattices (6 pi(2^64) (p, s) at 2 x 5,615,924 Ir) are the ledger row 2^70.6. [10w]
recomputes these prices from the shipped logs.

### Slice layout per suffix-box entry (`measurements-oh.txt`)

The `build` bench (27.27 Ir per entry; 12.325 per cell only at its 2^26-cell slice) is layout per suffix-box entry.
Per-function callgrind at 2^20..2^26 cells: `build_layout` 28.16 / 27.21 / 26.67 / 25.86 per entry, `zero_cells`
0.625 per cell (its own row), remainder < 2^20 Ir per run. Entries per cell K L4 / S: 0.452 at the bench, 0.0039 in
production. Charged: 28.17 per entry (K and each L4 count raised past rounding, `entries_upper`), 2^21 per slice,
2^17 per 2^32-byte segment; [10o] checks the log.

### Line path: the public slice12 line walk (recover-Q v9 734a9723, `measurements12x.txt`)

Same container, flags and differential method: (i) the base table on the 2^26 slice, one (p, s) per point: 24.098 /
31.113 / 88.313 / 310.700 at h = 2^13.82 / 2^9.82 / 2^5.82 / 2^1.82; (ii) random-pair groups: each seed draws 4 pairs
(prime q in [p, p + p/8], s uniform) and walks random cosets over all lines of a 2^26-2^32 slice, 3-8 seeds per group.
Single pairs scatter (random pairs average 35.85 at 2^9.8 against 31.11), so each `LINE_TAB` entry is the larger of
(i) and the pooled cost + 2 SE of the group at or below its h; h >= 2^18.81 is capped at 17.2 (cap groups pooled
17.04 / 16.62 / 16.29); linear in log2 h between entries. [10l] recomputes every quotient from raw Ir and hits and checks the
table against every group and the small-p runs. The bench adds a line-range argument and a multi-pair mode (diff in
recover-Q v9's log); the walk is the public code.
