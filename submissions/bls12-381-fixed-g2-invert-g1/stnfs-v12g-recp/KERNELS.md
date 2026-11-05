# Kernels and measurements

All counts are callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, valgrind 3.19) by the
differential method (two runs differing only in the measured unit), with each
kernel's correctness checks logged alongside its counts. The C kernels (`dsieve16`,
`slice12`, `ecm12`, `nbound31`) are described here; the shipped logs are the
evidence: `measurements12.txt`, `measurements12x.txt`, `wide12.c`, `wide12.txt` and `measurements-oh.txt` are sections of
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

At p = 2^69 the per-hit cost decreases with the radius (417.8, 370.0, 353.9 Ir
at R = 160, 192, 224: fewer enumeration nodes per hit, the rank is 12 lookups at
every radius); the 1000 charge is 2.8x the R = 224 value and covers R up to 384.
A complete (p, root pair) at p = 2^69, R = 96 measures 47,853.7 Ir per hit
including two setups over 1,396 hits (2^24.99 per setup, below PAIR_SETUP). Workers hold the
index tables (2^34.34 bytes at R = 384, STORAGE.md); `apply` runs at R = 4; the
`realloc` record buffer belongs to the harness, production streams records to
owners (STORAGE.md).

## Main phase: `slice12`, `ecm12`, `kernels_asm.b85` (listings; logs `measurements12.txt`, `measurements12x.txt`, `wide12.txt`, `measurements-oh.txt`)

Line and record paths of the relation sieve (line path for p < 2^48, record walk above; per-update costs by h, the points
per coset), apply 5.76, recovery 5.76 per record; slice layout per suffix-box entry; scan 3.002 per cell + 90.98 per leaf with zeroing fused (both below); certified norm bounds
(1,508.2 per evaluation; 4,600 per row = `row_setup` 4,292.1 measured (`nbound31` benchrow, `kernels_asm.b85` section C)
rounded to 4,300, plus 300 for the hull update (at most 4 places x 25 instructions = 100 counted, charged 300); 17,500 per slice; `embedding_bounds.py`, `nbound31`); five-limb Montgomery
multiplication 367.76 Ir; the LA kernels of STORAGE.md (listings counted by [10]): the constants of `model16.py`.

**Record walk below 2^63 (`slice12g`, 2c2ce2ab).** Exact adjugate recentring; all checks MATCH brute force. Record
44.67 / 77.42 / 250.54 / 810.82 / 1711.20 Ir per hit (p = 257 ... 268435459), setup6 406.87 / 374.07: `REC_TAB`.
`C_SETUP6` = 406.87 + 72.0 + 21.13 = 500.00 per (slice, (p, s)) (one coset origin c = -a(s) mod p serves the conjugate ideal
pair (r, -r); `pairs_all` counts NBAR/2 = 1 per prime): the larger measured value plus 93.13 carried from v6's
constant as an unallocated margin (no kernel work is assigned to it; it only raises the charge).

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
entry, `zero_cells` 0.625 per cell, remainder < 2^20 Ir per run. Charged: **28.17** per entry (`C_BUILD_E`, measured; the production entry count is the in-model formula: K + 2 lines and
2 floor(rb_j s) + 3 per coordinate, `entries_upper`, so every count is raised past rounding), 2^21 per slice and 2^17 per
2^32-byte segment; [10o] checks the log. Production rows have 0.0039 entries per cell (cell-weighted mean over pools
B, C, D), so the layout costs 0.109 Ir per cell instead of 12.325.

**Record walk per canonical a-run (`rec7.c`, `rec7.txt`; our BIGLEVER-1500 lane).** For p in [2^52, 2^63) the a-values
a + t e_0 of one run (`runs16.py`: the maximal t-interval around a on which no listed unit associate has smaller Q_a, by
exact integer quadratics in t (associates with a5 = 0 included, so a run can be shorter than the exact canonical interval), and
Q_a stays in the region; 11,325 units of log-radius <= 9) share one 7-dimensional
lattice per ideal; `walk7` emits the records of all the run's slices in one pass and `prepare7` keeps two reduced bases per
ideal: stretch 660 for runs of 32-128 a-values and stretch 90 for 129-1,024 (the campaign measured 32/64/128 at 660 and
181/362/1,024 at 90; 129-180 is in-model). [k7]: `rec7.txt` CHECK, record sets equal brute force: "check (record sets, folded, b5 = 0 skipped): ok (0 mismatches over 270 cosets, 2087 records)". Cells with b5 = 0 lie outside
the collected half-slice (b5 >= 1, proof.md section 3): no walk emits them and no count credits them. 130 differential
callgrind points (same container; slices 2^60-2^66, runs of 32-1,024 a-values, central and off-centre slabs, two
lattices; the measured walk includes record emission into the benchmark's small meta tables, whose masked index costs one
AND that production's full tables do not need) by h' = full-ellipsoid hits per run-coset; `REC7_TAB` is the maximum per 1-bit
level of log2 h', made non-increasing, read as a step (campaign primes 2^44.00-2^62.86, so [2^62.86, 2^63) is read at the same
h', in-model; the level at or below h'; worst ratio 1.0000).
**Schedule: one implementation per (run, band).** A run of 32 to 1,024 a-values (`RUN_T`, the campaign's run lengths)
takes rec7 in a quarter-bit band of [2^52, 2^63) iff log2(run cells) - m >= 4.99 (`REC7_TAB`'s bottom; m = the band's
midpoint log2 p, as for the other tables), run cells = the
half-slice cells of the run's a-values summed exactly (`runs16.lrs`: closed-form integer power sums of the cubed Q_a
quadratic); run length, run cells and band are the same for every row of the run, so all its slices agree, and rec7 then replaces
that band's slice record walk or empty cosets and slice setup6 for every slice of the run. Every other (run, band) (other run
lengths, run cells / p below the table, p < 2^52 or p >= 2^63) runs v10's slice kernel on each slice unchanged: record walk or
empty cosets by the slice's own h, a per-slice schedule of slices that share no work. No per-row choice between rec7 and the
slice walk. The model reads `REC7_TAB` at run cells / p (half-slice cells), one bit below the campaign's full-ellipsoid h' (2 x run
cells / p): a conservative offset, the table being non-increasing. Charged per record: step7 + `EXTRA7` = 2 Ir (a margin on
the measured walk, which already includes emission; it covers no unmeasured work); setup6 at `C_SETUP6` once per (run, (p, s)); per-coset walk costs are inside `REC7_TAB`'s differential.
**Clean runs only.** Production walks a run with rec7 only if every a-value of it is in the canonical list (the exact
canonical test of every eligible a: proof.md section 8, the 2^90 row; sorted by (a_1, ..., a_5, a_0) each run is one
contiguous block, so the check is O(T) per run, < 2^10 instructions per a, and the run boundaries (`runs16`, 11,325 unit quadratics) < 2^24 per a,
both inside that row's 2^41.8 per a). A run containing
an a-value the exact test rejects is contaminated and takes the slice kernel on every slice, so rec7 never walks a
non-canonical slice and every walked slice is a charged row. The model cannot test every run: it charges every run in the
window at the rec7 price and adds to every row `RUN_ALLOW` x `D_SUP` (`contam_sup`, [k7c]). RUN_ALLOW = 0.011 >= the share of
canonical rows whose run is contaminated (0 of 272 exact runs, 95%; proof.md H2). D_SUP = 2^73.067 Ir per slice bounds, for any
slice of the region, the slice-kernel work in [2^52, 2^63) plus those bands' setup6 (record walk at the exact maximum of
h x REC_TAB(h) up to the largest slice, else empty cosets): at least what a row of a contaminated run pays beyond its rec7
charge, however many non-canonical a-values its run has. [k7c] recomputes that work at 2001 slice sizes (maximum 2^73.015). prepare7 <= 3,986,661 Ir per
basis (`rec7.txt` PREP: `rec7.c` draws p in [2^(P2-1), 2^P2), so p in [2^51, 2^62) is measured; x 2^0.05 for
[2^62, 2^63), in-model): two bases for each of the 2^58.59 ideals with p in
[2^52, 2^63), ledger row 2^67.3. rec7 never runs at p >= 2^63 (walk7 is measured below 2^63): [2^63, 2^64) stays on v9's wide
path above (`REC_TAB_W` x 1.01, `setup6w`).

**Scan, apply, recover (v12; `slice12oh.c`, `measurements-oh.txt`; our overhead lane).** The public slice12 with three
functions replaced: scan2 tests four 16-bit cells per 64-bit word, ((w | H) - T) & H, eight cells per step, and zeroes both
sum arrays with its own stores (survivor lanes are rewritten SENTINEL by the scalar path and cleared per survivor before
reuse); apply_batch is unrolled x4; recover_segment tests eight sentinels at once. The lane test is exact because every sum
and threshold is < 2^15 (a cell sum is at most 16 log2|norm| <= 16 x 851). The shipped brute-force checks MATCH at 2^18-2^24
cells with segment shifts 0-3. Scan cost per leaf depends on the leaf length L (24 per 8-cell step, 22 per 4-cell step, 14 per
single cell): measured on the 2^26 slice at every L = 1..64 and 256 (our plain-track lane), charged 3.002 per cell + 90.98 per
leaf, the smallest leaf charge covering every point (L = 15 is the worst; the four block lengths 4/16/64/256 alone would give
48.05). Apply 5.751 and recovery 5.751 per record (p = 257), charged 5.76; the unmodified functions in the same runs give 7.000.
[10s] recomputes every point.
