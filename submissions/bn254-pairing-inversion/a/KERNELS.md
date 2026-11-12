# Kernels

Counts are differential callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns,
Valgrind 3.19). kernels_asm.b85 is base85/lzma text sections: retained listings, sources and measurement logs. Historical run-walk generator gen7.py and patch rec7.ed are omitted; their measurement logs are retained. model16.py holds the charged
constants. No SIMD is used. Machine premise: FEAT_LSE. The line run adds to an aligned 16-bit cell by one staddh (mod 2^16, no flags), one instruction as Valgrind executes and counts it; as ldrh/add/strh (the carried base table) the claim is 0.017 bit higher.

Two kinds of prices are used, and they are kept apart:

- **Measured at BN254 parameters** (this package): the exact-norm kernel for
  this f and g (norms.c, norms.txt), the filtering ratio and the unit share of
  matrix entries (la_bn.py), the descent's laws (DESCENT.md).
- **Carried**: instruction sequences measured for the BLS12-381 instance of the same sieve (same six-dimensional base field degree, same 5-, 6-
  and 7-dimensional lattices, same record format), whose prices are functions of the hit density only and are read here inside the density
  range in which they were measured: the line walk (LINE_TAB), the slice record walk (REC_TAB, empty cosets), setup, apply, routing, scan and
  layout, the FFT products of the batch trees, and the four-limb linear-algebra kernels (r has 254 bits here, four limbs as there). Their
  measurement logs (and the sieve sources) are in the archive as measured. The kernels are basis-blind: an ideal (p, s, rho) is a(s) = rho b(s) mod p with
  x(s) = sum x_j e_j(s), and the basis values e_j(s) (here: the power basis for a, common12.SB for b) are loaded from the factor-base entry.

## Main collection (carried)

Line path: a line c_5 = k is its L4 x W0 box (cell = idx4 W0 + r0 + c_0); a walk level keeps one index, the inner run (x8) is staddh at q, then q += 2 D0.
LINE_TAB (instructions per hit against log2 hits per line): the larger of the base table and pooled groups + 2 SE (measurements12x.txt, cachegrind; base kernel = callgrind on all 281 run pairs), cap
4.91 above 2^18.78; used here for p in [2^14, 2^40) at 2^5.48 to 2^33.29 hits per line (table range 2^1.83 and up). Slice record
walk: exact quotients in REC_TAB (measurements12.txt), setup 406.87 + 72 + 21.13 instructions per (slice or run, ideal pair (p, s, +-rho)); f only: 0.5 pairs per prime, as for empty cosets (proof.md section 3); the slice walk is read at the per-ideal density (half the pair's hits), which overstates its per-hit price; 1,600 per empty
coset; used here only for the cells outside the run rule, inside the table range (table bottom 2^-1.09; below it the
coset is charged as empty). Apply is 5.76 per record (3.63 for a run-walk record: owner path below); cell residue 0.02 per cell; scan (one side) 1.125 per cell plus 42.88 per leaf (leaf
lengths 1..64 and 256; slice12oh.c, measurements-oh.txt); layout 28.57 per suffix entry plus 2^21 per slice and 2^17 per 2^32-byte segment.
Routing: 10.01 per retained record of the slice walk, 1.14 for run-walk records (owner path below).

## Run walk on g ideals: walk7, measured

This package sieves g ideals. Every run-walk record is produced by walk7 (dadon's kernel, accepted on the f-sieved line in bd5e74d1 and 4c961bd2; scalar C, callgrind Ir; source not shipped). It walks the same coset as rec7 below (fixed-depth enumeration in an LLL-reduced basis, c_5 < 0 folded); an ideal reaches it only as (p, s, rho) and the basis values e_j(s). It was measured on ideals of g of this sieve (h(s) = 0, rho^2 = u w(s) mod p, p in [2^40, 2^56.25)) at the production region (theta = 2362/8000, A^2 = 300000). No price of the run walk is carried from f ideals or from another curve any more: model16.K7_TOP = 56.25.

Campaign (rec7bn.txt): 9 density strata x 1,024 walks on runs of 16 to 1024 a-values with l >= 4.99; a draw takes a (run row, band) cell with probability proportional to the ledger's charge and an ideal pair with probability 1/p; seeds fixed before any result was read, every draw kept. kappa = Ir / (H x ref) with ref the reference curve REC7_TAB at the model density l. Charged per stratum, as in the lineage's check: the larger of mean + 2 SE and the hit-weighted sum Ir / sum (H x ref) + 2 SE: [4.99) 0.7723, [8) 0.807, [10) 0.9447, [12) 0.9152, [14) 0.9075, [16) 0.9566, [18) 0.8521, [20) 0.761, [22) 0.6943; [10k] recomputes the table from the K7BN / K7HW lines and requires equality with REC7_K. Largest density drawn 24.53; above it the factor is 1.0. records / H between 0.995 and 0.997. Correctness on g ideals: 120 cosets against the exact enumeration, 144,370 records, 0 mismatches (CHECK7G).

Low densities. Run-rule bands with 3.5 <= l < 4.99 (before: the slice walk; on this sieve only the band [2^56, 2^56.25) reaches them) are walked too, at whole prices per charged hit (low points of REC7_TAB, factor 1): mean + 2 SE of 1,024 walks per stratum, [3.5) 266.07, [4) 232.72, [4.5) 202.71 (K7GL).

Short runs. Runs of 2 to 15 a-values (before: the slice walk) use walk7 with one more stretch class in every band with l >= -0.5, at whole prices per charged hit (REC7S): mean + 2 SE of 1,024 walks per stratum, [-0.5) 997.79, [1) 555.99, [2.5) 350.04, [3.5) 261.85, [4.5) 186.61, [6) 130.34, [8) 93.29, [10) 71.4, [12) 57.3, [14) 47.84, [16) 40.58, [18) 34.21 (K7GS). A stratum's price is its hit-weighted average, not a bound per cell; the lowest stratum is drawn from 24 pool cells and carries 0.02% of the short-run hits.

Owner path. walk7's measured count includes the emission of its records into the chunk pages, so EXTRA7 = 0. The owner of a chunk copies a page (measured 1.1252 instructions per record; + 0.01 hand-off allowance: REC_ROUTE7 = 1.14) and applies it (measured 3.6252: address, shift, one staddh; C_APPLY7 = 3.63). Both loops are branch-free per record; check: 4 pages, 524,288 records, every cell sum equal to the reference (OWNER line).

Bases: five reduced bases per ideal (four stretch classes and one for short runs) at 1.0e9 instructions each (exact-rational preparation, largest measured 5.8e8). Only per-stratum n, sum, sum of squares, max and the hit-weighted sums are shipped, not the per-walk records; the SE treats the draws as independent.

## rec7 (reference curve; its own prices are no longer charged)

REC7_TAB, the reference curve of rec7, remains the reference ref of step7 for l >= 4.99.

rec7 enumerates L7 = {(b, t): t = rho b(s) - a(s) mod p} over one canonical a-run: slab t_lo <= t <= t_hi, the b metric and the run's
parabola (lam A^2 - Q_a(a + t)) as radius. Historical source: rec7.c + rec7.ed and walk7a.S generated by gen7.py (fixed-level DFS and page emission inside the walk). The patch/generator are not included here; retained logs state the brute-force checks and prices. The same binary is run here on BN254 inputs.

rec7's own BN254 campaign (2,304 draws on ideals of f, accepted in the lineage) priced the run walk when f was sieved. It prices nothing here; its lines left the archive for those of the g campaign. Its basis-preparation log (prep7bn.txt, largest 4,093,887) stays as a record.

model16.step7 = ref x REC7_K with ref the reference curve REC7_TAB at the model density l (the hits the model charges: 2 x run cells x
lower/upper / p); EXTRA7 = 0 (owner path above). ref is not a per-walk bound; K bounds every stratum's mean (+ 2 SE).

## Exact norms (measured at BN254)

norms.c (this f and g; limit argument T, 9-limb norms): 27,000 charged; 3,115 pairs (500 uniform to the coordinate caps, 500 uniform on the region box, 1,600 real pool pairs, 512 corners, 3 trivial): largest 26,877 Ir, every output equal to common12.norms (norms.txt).

## g thresholds (table, no kernel)

The thresholds of the sieved side are read from one table of leaves over the b-ball, built once from exact integer norms N(b) (proof.md section 3, Lemma G; gthr.py is the reference). No norm-bound kernel runs per slice; the nb2 logs in the archive are historical. Charged: C_GTAB = 1,024 instructions per cell of the b-ball for the table (per row seven fraction-free 6 x 6 determinants on entries below 2^15, about 2^13 instructions each; per cell six two-limb additions, a leading-zero count with a 2^10-entry logarithm table and the leaf comparison, under 64 instructions) and C_GROW = 472.5 per (slice, row) for cutting the row's leaves to the slice (the row's leaves are stored in order: two searches among at most 2^5 breakpoints and two boundary records, under 200 instructions). The scan kernel is unchanged: 1.125 per cell and 42.88 per leaf, leaves cut into its measured lengths.

## Batch trees (carried)

srtf.c: the SRT lemma of proof.md over GMP 6.2.1 (exact vs mpz_mod and a gcd strip; archive fft.txt). FFTPTS N:charged are differential Ir
+ 1000 per allocator call + realloc copies, x 1.10; all are <= F(N) = 8.5 N log2(N)^2 + 1200 ([5t]). Above 2^26 limbs a product is a
distributed Schoenhage-Strassen multiplication (6 transposes at 84 per record: + 504 N); H-FFT: local work <= F. The trees here are those of
proof.md section 3 with 9-limb leaves (both norms are below 2^576, [5t]); the batch test now runs on N_f: C_GBATCH = 2^23.05 >= the 9-limb tree cost ([5t]).

## Descent and linear algebra

DESCENT.md (dsieve16: hit 1100, scan 5, zero 1, bound 2400 per cell; pair setup 2^25.1) and STORAGE.md (owner_part, _row_red, _chi_mul;
krylov5.s in the archive; the other listings stay with the BLS12-381 package) carry their kernels unchanged: they act on six- and twelve-dimensional integer lattices and on
four-limb residues mod r and do not depend on the curve.

## Code archive

lib.b85 (base85/lzma text, '@@ name' sections; common12.py serves them to import) holds modules used unchanged: counts16.py (sha256 4b52153b2560b46c), ecm16.py (sha256 e75b63b5f03c7d12), filtsim.py (sha256 f2114f534819d39f), ideal_canonical.py (sha256 1cbfd5251b957bf7), lattice_count_bounds.py (sha256 f60b7c292f96055c), unit_orbit.py (sha256 e62732368a5d4755), byte-identical to the BLS12-381 package's; and descent_witness.py (sha256 dfe4c434628540c6), the descent's lattice and norm sampler for this field (DESCENT.md).

## Composite-modulus ECM extraction evidence

stage2-summary.json reports 24 scaled cases (480 curves) at 36/40-bit factor sizes and composite moduli up to 1746 bits. Curves and their unconditional replay use the original composite modulus; extraction takes stage-1 prime-power prefix gcds and stage-2 giant-block product gcds, without per-term gcds. Planted-factor reductions audit detection separately from extraction. 133 of 135 detected factor events were isolated on that curve; 35 of 36 small factors were separated over their 20-curve schedules, with one zero-hit case. These data support H-SEP at the measured scales, not a proof at 132 bits. The model discounts the per-curve order-hit probability by one half and amplifies the discounted probability to the stated per-prime schedule target; the extrapolation is explicit.

The finalmodel4.py shim loads model4core.py from lib.b85; that section contains the complete fixed plan and evaluator source, not an external dependency. The archive also includes the BN toy/child and PGF replay sources.

The package stores model16core.py and ledger16core.py in lib.b85 with their original top-level module names served by import shims. prefix_price.py and model4core.py are embedded there too. No external evaluator is required. The legacy480-curve separation experiment used replay; the stored-prefix lemma and six exact partition comparisons in DESCENT.md transfer the output partition to the new first-pass algorithm.

## ntt3 multiplication prices (dadon's kernel)

ntt3 (dadon's, scalar C; counts in ntt3.txt, checked in [5t]) prices the products of every batch tree through mul_F = min(F, F_ntt): three number-theoretic transforms modulo c 2^55 + 1, c = 131, 174, 197, one limb per point, n = 2^k <= 2^55 points (exact since n (2^64 - 1)^2 < p1 p2 p3). No data-dependent branch, so Ir = 69.75 n k + 90 N + 151.5 k + (119 n + 371, k even; 139.25 n + 271.5, k odd), radix 4 (dadon dd1d38fd) exactly (N output limbs), equal to callgrind on all 49 measured products (k = 2..22), every product exact. Above 2^26 points: + 3 x 504 n (transposes) + 128 n (twiddle pass), allowances, not measured. Here it prices the f batch (9-limb leaves: 2^23.043 <= C_GBATCH 2^23.05 instead of 2^24.22) and, through gbatch_cost, the exact strip of the descent nodes.

The experiment body is served, unchanged, from libextra.b85 (section evidence_core.py); evidence.py is its entry point.
