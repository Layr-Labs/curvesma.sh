# Kernels

Counts are differential callgrind Ir on linux/arm64, Debian12, GCC12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, Valgrind3.19. kernels_asm.b85 is base85/lzma text sections: listings,
sources, patches and logs. model16.py holds charged constants. No SIMD is used. Machine premise: FEAT_LSE. Line run and run-walk apply add to an aligned 16-bit cell by one staddh (mod 2^16, no flags), one instruction as Valgrind counts.

## Main collection

Basis. Coordinates use common12's basis e_j of O_h (e_0 = 1; [0c]); an ideal (p, s, r) is a(s) = r b(s) mod p,
x(s) = sum x_j e_j(s): setups load e_1(s)..e_5(s) from the factor-base entry (measured: s^j by four power steps).
Walks, records and bounds are basis-blind.

slice12 recentres by exact adjugate remainders (measurements12.txt: checks MATCH). Record walk prices are exact quotients in REC_TAB;
setup6 is406.87+72+21.13 instructions, empty cosets1600. Line path: (cat linewalk.ed; echo w) | ed -s slice12.c (084b480c; no comments, .text as measured). A line
b5 = k is its L4 x W0 box (cell = idx4 W0 + r0 + b0): a walk level keeps one index, the inner run (x8) is staddh at q,
then q += 2 D0. Reference binary64 operations in order: same nodes and hits in order; in-order traces vs c9's
walk (each hit's (k, idx4, b0), each node's bits) (and c13's) match on 48 runs; brute force MATCH; log hits as c9's. LINE_TAB: the larger of base TAB and pooled
groups+2SE (measurements12x.txt; cachegrind, = callgrind on the old walk's 281 run pairs); cap4.91 above2^18.78. Apply is5.76 per slice-walk record. CELLRES (no hits): total minus layout, scan (with leaf lists) and zero is0.002 per cell(2^24->2^26), charged0.02 (CELL_ROUTE). Scan (f only) is1.125 per
cell plus42.88 per leaf, covering every leaf length1..64 and256; slice12oh.c
and measurements-oh.txt. Cells hold H + s (H = 2^15, s, T < 2^15): w - T has bit 15 iff s >= T. Layout is28.57 per suffix entry plus2^21 per slice and2^17
per2^32-byte segment; logs include allocation/zeroing. Norms (norms.c, norms.txt): common12.norms by 11 products mod h and the norm through x -> -1/x (fixed limbs, no division), exact for |coordinates| <= 2^13 ([5t]), with L and the compare (limit T = 16 SA, + 96 on class Y): 2,600 pairs, 0 mismatches, <= 48,231 Ir (C_NORMS 48,300). Bounds (f only): one evaluation is nbound31f's, bit-identical to the reference (measurements12.txt: 1154.6 on single rows). nb2 (C, not shipped: bytes) = model16.two_level around it (both levels, merge, scan-leaf cut) and embedding_bounds.hull; row contexts are inputs. Campaign (setting, law, size and rule fixed before the first draw): draw i by random.Random("c26-nb2:" + str(i)), a priced S then G region row with probability proportional to weight(qa) x evaluations / box cells / fk; unit = its box-block; Ir = (3 calls - 1 call)/2. On all 512 the hulls, leaves (interval, bounds, scope), evaluations and scan pieces equal the Python reference. Per evaluation mean 1173.6, sd 41.1, + 2 SE 1177.3, largest 1301.7 (255,728 evaluations): C_NB_EVAL = max(1154.6, 1178). Hulls per row 154.9, + 2 SE 155.10: C_NB_ROW = ceil(222.15 + 155.10) = 378 (nbrow.txt: 222.15 per row, 1254 per slice, charged 18500; nbrow.c not shipped). nb2.txt: n, sum, sum of squares, maximum of ceil(10 x); [10s] recomputes both prices.

For2^63<=p<2^64, wide12.c prepares a __int128 lattice (entries guarded at2^40),
then the unchanged walk; wide12.txt: identical .text, exact reductions, 400/400 planted
hits, MATCH. setup6w=141.00; REC_TAB_W x1.01, setup6 and1600 per empty coset cover these;
four preparations (A and coset) at5615924 per ideal fit the2^71.6 row.

## rec7: fixed-level DFS

rec7 enumerates L7={(b,t):t=r*b(s) mod p} over one canonical a-run, with a slab
t_lo<=t<=t_hi and the b metric in slice_mb.h. Source is original rec7.c
(sha25604a2ab937a56eebf)
plus rec7.ed ((cat rec7.ed; echo w) | ed -s rec7.c; wide basis preparation; walk7 calls walk7a.S from gen7.py). Prepare7w uses __int128 rows,
quad-precision GSO, a2^40 guard, congruence/determinant checks. A failed basis takes the slice walk.
Two-limb recentring is the inherited rec7 range assumption (entries below2^28).

Refinement lemma. walk7a runs the reference DFS levels 6..1 as fixed loops with the reference's binary64
operations on the same operands in order, so the same nodes and leaves. Integer bookkeeping is exact: slab clip
l0=max(l0,-floor(v/ad)), h0=min(h0,floor((W-v)/ad)), v=sgd t-TLs; runs split at k=0 into folded one-sign
sub-runs; so the same records in the same order. Page emission: the walk runs STORAGE.md's scatter itself. With {row, tails}
at entry (t << 14) + k of the run's line table (line k < 2^14 of a-value t; 8 Ir per entry, layout: carried over, not re-measured at this pitch): cell = row[i4] + b0, one post-indexed stp of (cell mod 2^32, log)
at tails[cell >> 30], tail stored back, hand-off iff its bit 20 is set (11 instructions; 9 where the entry address is a stepped register). Lemma: cell < 2^60, so the
index is < 2^30; a page is the lower half of a 2^21-aligned window, a tail a multiple of 8 in [page, page + 2^20], so stores stay in the page, full exactly when
bit 20 is set; a chunk gets its records in emission order. The stream equals the previous emitter's: page fingerprints and counts on all 5,120 draws, record by record (capture hook) on strata 0-7. Checks, 0 mismatches: 270 brute-force cosets, wide 432 pairs
(1728/1728 planted), in-order trace vs c14's walk (TCHECK; 32-record pages, a table per 4 lines); eight wrong
variants fail it.

Price: production campaign of x-walk-asm (ours; walk7a's cosets and records; rec7.txt XA, K7 lines; strata, sizes, estimator and seeds fixed before the first draw). A draw is one walk production makes: (row, quarter-octave band, ideal) with probability proportional to the reference charge weight(qa) x ln(band) x ref(l, T) x 1/p: a pool S or G region row the model prices with rec7; a degree-one prime ideal of the band (1/p); walked on the row's own run: sub-run t = 2u of the row's class (X: all b; Y: b_0 = 2k + 1 in the P2 basis), the config's b metric in those coordinates, the run's parabola lam_c A^2 - kap Q_a(a + t e_0) as centre and radius, residues r e_k(s)/2 mod p, coset -a(s)/2, stretch 660 up to 128 a-values, 90 above. Draw i of stratum k: random.Random("xwa-bls31:k:i"); 512 draws in each of 10 strata of the model density l (edges 4.99, 8, 10, ..., 24); all 5,120 kept. kappa = Ir / (H x ref): Ir = instructions of 2 cosets - 1 (callgrind); H = 2 x 0.999595 vb fk sum_u (lam_c - kap qa(2u))^3 / p, the hits the model charges for the walk (measured hits/H 1.0003..1.0008 per stratum); ref = REC7_TAB's interpolant, the envelope of c21's 544 bench points (x REC7_R2 = 1.29 above 1,024 a-values). REC7_K = mean + 2 SE per stratum (0.4612..0.6009; sd 0.03..0.08); model16.step7 = ref x REC7_K. [10r] recomputes every K from each stratum's count, sum and sum of squares of the values rounded up to 0.001. ref bounds every stratum's mean (K <= 1, checked); no per-walk bound is claimed (0 of the 5,120 walks exceed it; largest 0.813). Not modelled: the half-step centre offset of b_0 = 2k + 1. K7 counts whole walks with their page emission: no further per-record term. Run-walk records are 8 bytes (cell mod 2^32, 16 log2 p), stored by stp w7, w18, [x1], #8; campaign, trace and checks are this walk's. Owner (lean8.txt; MATCH vs C): page copy 1.0626 per record (16 ldp + 16 stp + subs + b.ne per 32), apply 3.1252 (per record ldp w, w; add; staddh; subs + b.hs per 16; 16-bit sums, the chunk's two low id bits in the array base). REC_ROUTE7 = 1.07 (copy + STORAGE.md's 0.0009 hand-off), C_APPLY7 = 3.13. Pair walk: f is
even in x ([0b]), so (b,t) is a hit of (p,s,r) iff (-b,t) is one of (p,s,-r): one walk per pair over both signs of b5
(the bench) emits both ideals' hits, 2 x coset cells/p. model16 reads ref at l = lg+1+log2(lower/upper) - log2 p ([3]), one
setup6 per pair. Basis preparation is unchanged, <=3995157 narrow and
150271737 wide instructions; six bases (2 for A, 2 per run coset) per2^59.56 ideals from2^48 to2^64, all at
the larger price times2^0.05, fit2^75.09 units.

## Coset walks (P2)

On b_0 = 2k + 1 or t = 2u the kernels run unchanged on the lattice in (k, b_1..b_5) or u (constants halved mod p, one
region axis halved); the target (c - 1)/2 adds 9.0 Ir to setup6 (415.87), charged 10. COSL lines: the line walk on b_0 = 2k at 2^25, below LINE_TAB ([10r]).

## Descent

dsieve16 certifies congruences/determinant, enumerates, leaf-tests and ranks hits
(checks throughR256, DESCENT.md). Hit1100 (walk/rank1000, transfer84, apply9); scan5, zero1,
bound2400 per cell; pair setup2^25.1.

## Stage B (K1, K2)

K1 = the FFT SRT below with P24 (1,077,871 primes, 24,197,143 bits, checked), root by two exact divisions (k1.c, sha256 602bb96d7ee63748). On 6 waves of real pool norms every remainder equals mpz_mod; charged as FFTPTS, per leaf: 758,037 (2^14 leaves), 748,999 (2^15; 3 waves, + 2 SE 749,038), 870,463 (2^16). K2 (k2.c, sha256 ce36ea41483bd5f1, plain C): 2-part, Montgomery reduction of y, 10 Montgomery squarings, divstep gcd, exact division, both logs; on 1,024 real norms every output equals Python's; Ir = (3 calls - 1 call)/2 per leaf: mean 70,433, sd 2,330, + 2 SE 70,579, largest 78,471. C_SB = 749,038 + 70,579 + 200 (25-word hand-over, window test; allowance) <= 820,000. Sources not shipped (bytes); sb.txt: n, sum, sum of squares, maximum; [10k].

## FFT SRT: band strip and g batch

srtf.c (sha256 ac0437e59fb4865a): the proof.md SRT lemma over GMP 6.2.1 (Debian arm64, static; same toolchain).
Exact vs mpz_mod and a gcd strip (archive fft.txt SRTWAVES). FFTPTS N:charged (rounded up): differential
Ir + 1000 per allocator call + realloc copies, x 1.10 (NEON copies <= 8.1%); all <= F(N) = 8.5 N log2(N)^2 + 1200
([10b], audit); charged/F 0.85 at 2^20, 0.46 at 2^26. Above 2^26 limbs (P_g 2^64.03, root ~2^65.2) a product is a
distributed Schoenhage-Strassen multiplication (pieces of ~sqrt(4N) limbs multiplied locally; 6 transposes at the
descent's 84 per record: + 504 N, model16.fft_F). H-FFT: local work <= F; the largest measured c = charged/(N log2 N
log2 log2 N) is 33.31 (2^20): 13.1k per limb at 2^65.2 vs F's 36.1k.

_combine is hand-written (two adds/adcs chains, 32 instructions; idea: yukonautominer, Provenance); 2,000,000 inputs equal to the compiled one.
