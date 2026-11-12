# Kernels

Counts are differential callgrind Ir on linux/arm64 (Debian 12, GCC 12.2 -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns,
Valgrind 3.19). kernels_asm.b85 is base85/lzma text sections: retained listings, sources and measurement logs. Historical run-walk generator gen7.py and patch rec7.ed are omitted; their measurement logs are retained. model16.py holds the charged
constants. No SIMD is used. Machine premise: FEAT_LSE. The line run and the run-walk apply add to an aligned 16-bit cell by one staddh (mod 2^16, no flags), one instruction as Valgrind executes and counts it; as ldrh/add/strh (the base tables of measurements12x.txt and linex.txt, apply 3.63) the claim is 0.043 bit higher.

Two kinds of prices are used, and they are kept apart:

- **Measured at BN254 parameters** (this package): the run walk on g-ideal draws at theta 2189/8000, unscaled b form, carried per density stratum (rec7bn.txt, rec7own.txt; below), the exact-norm kernel for
  this f and g (norms.c, norms.txt), stage B on N_g (sbn.txt), the norm-bound kernel on this field's rows (nb2bn.txt), the filtering ratio and the unit share of
  matrix entries (la_bn.py), the descent's laws (DESCENT.md).
- **Carried**: instruction sequences measured for the BLS12-381 instance of the same sieve (same six-dimensional base field degree, same 5-, 6-
  and 7-dimensional lattices, same record format), whose prices are functions of the hit density only and are read here inside the density
  range in which they were measured: the line walk (LINE_TAB), the slice record walk (REC_TAB, empty cosets), setup, apply, routing, scan and
  layout, the FFT products of the batch trees, and the four-limb linear-algebra kernels (r has 254 bits here, four limbs as there). Their
  measurement logs (and the sieve sources) are in the archive as measured. The kernels are basis-blind: an ideal (p, s, rho) is a(s) = rho b(s) mod p with
  x(s) = sum x_j e_j(s), and the basis values e_j(s) (here: the power basis for a, common12.SB for b) are loaded from the factor-base entry.

## Main collection (carried)

Line path: a line c_5 = k is its L4 x W0 box (cell = idx4 W0 + r0 + c_0); a walk level keeps one index, the inner run (x8) is staddh at q, then q += 2 D0.
LINE_TAB (instructions per hit against log2 hits per line): the larger of the base table and pooled groups + 2 SE (measurements12x.txt, cachegrind; base kernel = callgrind on all 281 run pairs) to 2^18.78 (4.91),
then steps from single lines of 2^28/2^30 slices, same kernel and differential (linex.txt; per step the largest run of every family at or above it), cap 3.85 above 2^23.42; used here for p in [2^14, 2^40) at 2^5.48 to 2^33.29 hits per line (table range 2^1.83 and up). Slice record
walk: exact quotients in REC_TAB (measurements12.txt), setup 406.87 + 72 + 21.13 instructions per (slice or run, ideal pair (p, s, +-rho)); g only: 0.5 pairs per prime, as for empty cosets (proof.md section 3); the slice walk is read at the per-ideal density (half the pair's hits), which overstates its per-hit price; 1,600 per empty
coset; used here only for the cells outside the run rule, inside the table range (table bottom 2^-1.09; below it the
coset is charged as empty). Apply is 5.76 per record (run walk: 3.13, below); cell residue 0.02 per cell; scan (sieved side only) 0.53 per cell plus 49.51 per leaf: our scana.S (source kept with its campaign, not shipped, for bytes;
scana.txt: scan2's outputs; a 32-cell group is tested exactly only if its OR meets lane bits floor(log2 T)..14; Ir <= 0.5 L + 49.51 at L = 1..64, 128, 192, 256; + 0.03 for candidates); layout 28.57 per suffix entry plus 2^21 per slice and 2^17 per 2^32-byte segment.
Routing: 10.01 per retained record of the slice walk, 1.14 for run-walk records (below).

## Run walk (rec7): kernel carried, price measured at BN254

rec7 (gbotrel's page-emission walk7a: rec7.c + rec7.ed + gen7.py, generator not shipped) and walk7 enumerate the same coset L7 = {(b, t): t = rho b(s) - a(s) mod p} of one canonical a-run (slab t_lo..t_hi, b metric, the run's parabola as radius) and emit the same record stream.

Run-walk prices: walk7 on g ideals. walk7 (dadon's kernel, 258122fe; scalar C, callgrind Ir; source not shipped) is measured on the g ideals of this sieve (h(s) = 0, rho^2 = u w(s) mod p) at theta = 2189/8000, on the unscaled b form; the production region here (theta = 2280/8000, b form / 0.93^2) is read with the same per-density-stratum prices, as in kaankolcu's a676d4c9 / 8785cc1c (prices are per density stratum, and a walk's cost is set by its coset density and ideal, which the strata condition on, not by the band of Q_a(a); own pool of 50,703 canonical rows): no f-ideal table and no G7 remain. 1,024 walks per stratum; a draw takes a (row, band) with probability proportional to the ledger's charge and an ideal pair with probability 1/p; seeds fixed first, every draw kept. p < 2^56, l >= 4.99 (rec7bn.txt): kappa = Ir / (H x ref), charged the larger of mean + 2 SE and the hit-weighted sum Ir / sum (H x ref) + 2 SE: [4.99) 0.721, [8) 0.7227, [10) 0.8305, [12) 0.8025, [14) 0.8019, [16) 0.8614, [18) 0.7674, [20) 0.6848, [22) 0.6255; 1.0 above l = 24.5 ([10k] recomputes, requires equality). [2^56, 2^57) (rec7own.txt, K7TOP), mean + 2 SE: at least mean + 2 SE (the larger of it and the hit-weighted value): [4.99) 0.6591, [6.5) 0.6793, [8.5) 0.6447. 2.5 <= l < 4.99 (K7LOW), per charged hit: [2.5) 306.75, [3.5) 258.35, [4) 223.79, [4.5) 194.79 ([10w]). Runs of 2 to 15 a-values (RUN_S; before: slice walk) use walk7 with one more stretch class where l >= -0.5, per charged hit (REC7S, K7HS; l from -0.5, 1, 2.5, 3.5, 4.5, 6, 8, ..., 18): 1011.38, 545.41, 350.26, 250.57, 179.26, 119.78, 82.68, 62.66, 50.44, 42.05, 35.92, 30.49; a stratum's price is its hit-weighted average, not a per-cell bound. Check: 120 g-ideal cosets against the exact enumeration, 144,370 records, 0 mismatches. Two builds since this version: walk7 was rewritten for speed (inner three levels in one function, hand-written emission loop) in two builds of one source that differ only in how a run whose c_5 changes sign is emitted (record by record in build j, split by one division in build i). Same enumeration, same records: the record hash of both builds equals the first build's on all 53,248 campaign walks, and the C path of the same sources has 0 mismatches against brute force (400 toy cosets) and the exact enumeration (308 production cosets). Both were measured on the same walk lists; each stratum is charged at the cheaper build's own bound (j for l < 18, the top octave, low densities and short runs; i for l >= 18 and the last short-run stratum; the two lowest short-run strata keep the first build's price). The other build's charged values are in rec7bn.txt. With these tables walk7 is the cheaper kernel in every cell class of step7. Only per-stratum sums are shipped; the SE treats the draws as independent. The f-ideal campaign lines of rec7 and walk7 left the archive: they price nothing.

model16.step7 = ref x REC7_K, ref the reference curve REC7_TAB at the model density l (REC7_K56 on [2^56, 2^57); below l = 4.99 the factor is 1 and REC7_TAB holds whole prices). Owner path: walk7's count includes the emission into the chunk pages (EXTRA7 = 0); the owner copies a page (measured 1.1252 per record + 0.01 hand-off: REC_ROUTE7 = 1.14) and applies it with one staddh per record (measured 3.12513, OWNERL: C_APPLY7 = 3.13); 524,288 records, every cell sum equal to the reference (OWNER). Bases: five per ideal at 1.0e9 (largest measured 5.8e8).

rec7 on g ideals (kshitij-hash's; rec7g.txt, [10r]): same region as walk7's measurement, draw rule and Ir method (Ir(2 cosets) - Ir(1 coset)); 0 mismatches against brute force (2,160 toy g cosets) and an exact enumeration (98 production cosets). REC7_KR (7 strata, p < 2^56, l < 19.96) and REC7_KR56 ([2^56, 2^57)) are charged as REC7_K is; rec7 emits up to 1.5% more records than charged hits, so each rec7 hit pays EXTRA_R = 2.11 >= 2 x 1.0151 + 4.27 x 0.0151 = 2.0947 (rounded up). Each cell class (band, stratum) runs the cheaper measured kernel, min(walk7, rec7 + EXTRA_R) (model16.step7); four rec7 bases per ideal at 5.0e6 are added to the bases row.

## Sieved side g

The kernels take an ideal (p, s, rho) of g (rho^2 = u w(s) mod p) as one of f. The line and slice lattices {b : b(s) = c mod p} do not depend on rho (only the coset does): their tables apply unchanged. The run-walk lattice L7 = {(b, t): t = rho b(s) mod p} (determinant p, rho a root of an even polynomial for either side) does; the earlier campaigns drew ideals of f. Paired check (not replayed): 400 random prime ideals (p, s), p in [2^40, 2^57), one root of f and one of g each, L7 LLL-reduced in the metric Q_b + Q_a(1) t^2: mean Gram-Schmidt log2-length differences f - g 0.010, 0.023, 0.012, -0.002, 0.003, -0.026, -0.019 (paired SE 0.012-0.022), all within 2 SE. Since this version the run walk is measured on g ideals (above); G7 = 1.

## Exact norms (measured at BN254)

norms.c (this f and g; limit argument T, 9-limb norms): 27,000 charged; 3,115 pairs (500 uniform to the coordinate caps, 500 uniform on the region box, 1,600 real pool pairs, 512 corners, 3 trivial): largest 26,877 Ir, every output equal to common12.norms (norms.txt).

## Norm bounds (two-level thresholds)

The norm-bound kernel (nb2) is the BLS12-381 package's binary, run on this field's rows: 256 production vectors, every output MATCHes the Python reference; 1187.5 instructions per evaluation and 155.4 per row (mean + 2 SE; nb2bn.txt). It evaluates both norms; the thresholds use its g output. Where a place's interval meets 0 the measured binary returns no bound for that side; the evaluator keeps its certified upper product (formed before the sign test) with lower bound 0; this branch occurs at 398 of 1,068 g evaluations on 300 S rows and is charged +16 per evaluation. Charged: these BN254 values rounded up, 1187.6 per evaluation and 155.4 per row ([10k] checks both against nb2bn.txt), plus 16 per evaluation for the zero branch (above): 1203.6; the slice term was not measured here and stays the BLS12-381 package's 18,500 x NB_ALLOW 1.25 = 23125.

## Stage B (measured at BN254, on N_g)

sbn.c (our code, sha256 in sbn.txt; source not shipped, for bytes, as for our accepted c27/c28 stage-B kernels): the K2 leaf of our
BLS12-381 c27 with 9 limbs plus one mpn_tdiv_qr for y = P14 mod N_g (P14: the primes below 2^14). On the exact N_g of the 9,856 primitive canonical rows of pool S
with Q_a <= (2189/8000) A^2 (the unscaled b form; production N_g here are smaller) and on 264 planted inputs every output equals the exact value (sbn.txt). Ir per call = (Ir(3 calls) - Ir(1 call))/2 on
1,024 of them: mean + 2 SE 71,153, largest 73,442. C_SB = 80,000 >= largest + 200 (hand-over and window test: allowance) per stage-A survivor ([5s]).

## Batch trees (carried)

srtf.c: the SRT lemma of proof.md over GMP 6.2.1 (exact vs mpz_mod and a gcd strip; archive fft.txt). FFTPTS N:charged are differential Ir
+ 1000 per allocator call + realloc copies, x 1.10; all are <= F(N) = 8.5 N log2(N)^2 + 1200 ([5t]). Above 2^26 limbs a product is a
distributed Schoenhage-Strassen multiplication (6 transposes at 84 per record: + 504 N); H-FFT: local work <= F. The trees here are those of
proof.md section 3 with 9-limb leaves for N_f (below 2^576; N_f < 2^521.1 on the region, [5t]); C_GBATCH = 2^23 >= the 9-limb tree cost of the f batch (waves of 2^49) (ntt3 prices through mul_F; 2^24.22 with F alone).

## Descent and linear algebra

DESCENT.md (dsieve16: hit 1100, scan 5, zero 1, bound 2400 per cell; pair setup 2^25.1) and STORAGE.md (owner_part, _row_red, _chi_mul;
the listings, krylov5.s included, stay with the BLS12-381 package) carry their kernels unchanged: they act on six- and twelve-dimensional integer lattices and on
four-limb residues mod r and do not depend on the curve.

Dense-band tables ([6g], tab.txt). Per Krylov step and 5-vector slot, sums of all multisets of g orbit columns of norm <= 2^28, 2^18.25, 2^14 (g = 2, 3, 4), terms y, chi y, 2r - y, chi(2r - y); a row stores one pointer per group, added by unit_run5 (bit-exact). Pairs (< 4r) unreduced, g >= 3 reduced to [0, r): [6b]'s mass bound holds. Rows group greedily, densest band first; adds are exact Poisson sums over the [6] rows (la_bn). Band edges are absolute norms on the rows' per-bit (log2 q) composition; the ratio falls with scale. Row re-encoding (< 2^47 instructions) is in filtering.

## Code archive

lib.b85 and libextra.b85 (base85/lzma '@@ name' sections served to import by common12.py) hold every module: ecm16, filtsim (with the column-partition clamp of our c29q2, sha256 7691b5e3d0d5339c), ideal_canonical, lattice_count_bounds, unit_orbit (byte-identical to the BLS12-381 package's), descent_witness, model16core, ledger16core, model4core, prefix_price, evidence_core (body of evidence.py).

## Composite-modulus ECM extraction evidence

stage2-summary.json: 24 scaled cases (480 curves); DESCENT.md states what they support (H-SEP at the measured scales) and the explicit extrapolation.

Shims (finalmodel4.py, model16.py, ledger16.py) load the embedded modules (complete plan and evaluator; no external evaluator). The legacy 480-curve separation experiment used replay; the stored-prefix lemma and six exact partition comparisons in DESCENT.md transfer the output partition to the new first-pass algorithm.

ntt3 (dadon's; ntt3.txt, checked line by line in [5t]): three NTTs mod c 2^55 + 1 (c = 131, 174, 197), one limb per point, exact for n <= 2^55; radix 4 (dadon 8cc015a0): Ir = ntt_F (model16core), equal to callgrind on all 49 measured products, 2 <= k <= 22 (37 listed); above 2^26 points + 1640 n, allowances. mul_F = min(F, F_ntt) prices every batch tree here (f batch, descent strip).
