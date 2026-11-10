# Kernels

Counts are differential callgrind Ir on linux/arm64, Debian12, GCC12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, Valgrind3.19. kernels_asm.b85 is base85/lzma text sections: listings,
sources, patches and logs. model16.py holds charged constants. No SIMD is used.

## Main collection

Basis. Coordinates use common12's basis e_j of O_h (e_0 = 1; [0c]); an ideal (p, s, r) is a(s) = r b(s) mod p,
x(s) = sum x_j e_j(s): setups load e_1(s)..e_5(s) from the factor-base entry (measured: s^j by four power steps).
Walks, records and bounds are basis-blind.

slice12 recentres by exact adjugate remainders (measurements12.txt: checks MATCH). Record walk prices are exact quotients in REC_TAB;
setup6 is406.87+72+21.13 instructions, empty cosets1600. Line path: (cat linewalk.ed; echo w) | ed -s slice12.c (084b480c; no comments, .text as measured). A line
b5 = k is its L4 x W0 box (cell = idx4 W0 + r0 + b0): a walk level keeps one index, the inner run (x8) is ldrh/add/strh at
q + jD0. Reference binary64 operations in order: same nodes and hits in order; in-order traces vs c9's
walk (each hit's (k, idx4, b0), each node's bits) (and c13's) match on 48 runs; brute force MATCH; log hits as c9's. LINE_TAB: the larger of base TAB and pooled
groups+2SE (measurements12x.txt); cap6.06 above2^18.78. Apply is5.76 per record. CELLRES (no hits): total minus layout, scan (with leaf lists) and zero is0.002 per cell(2^24->2^26), charged0.5 (CELL_ROUTE). Scan (f only) is1.125 per
cell plus42.88 per leaf, covering every leaf length1..64 and256; slice12oh.c
and measurements-oh.txt. Cells hold H + s (H = 2^15, s, T < 2^15): w - T has bit 15 iff s >= T. Layout is28.57 per suffix entry plus2^21 per slice and2^17
per2^32-byte segment; logs include allocation/zeroing. Bounds (f bound only, bit-identical; measurements12.txt) cost1154.6 per block,
535 per row and18500 per slice (nbrow.c/.txt: 222.2 per row, 1254 per slice); warm bisection changes only charged counts.

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
of the line from a per-slice table (8 Ir per entry, layout): cell = row[i4] + b0, one post-indexed stp of (cell, ideal)
at tails[cell >> 30], tail stored back, hand-off iff it is 0 mod 2^20 (12 instructions). Lemma: cell < 2^60, so the
index is < 2^30; a tail is a multiple of 16 inside its 2^20-aligned page, so stores stay in the page, full exactly at
0 mod 2^20; a chunk gets its records in emission order. Checks, 0 mismatches: 270 brute-force cosets, wide 432 pairs
(1728/1728 planted), in-order trace vs c14's walk (TCHECK; 16-record pages, a table per 4 lines); six wrong
variants fail it.

rec7.txt prices walk7a at the same130 narrow+42 in-range wide points
plus all 96 CAMPH points (p 2^47.03..2^57, h' 2^23..2^30); REC7_TAB is the monotone envelope of the points with T <= 1,024; runs above 1,024 a-values use REC7_TAB2, that of all points with groups G7 and W2048 (T = 2,048; 28 + 4 + 6 wide); both are
rounded up to0.01, a step function in log2 hits (bottom4.99; [10r] checks every point). EXTRA7=2 remains; run-walk routing is3.01 (STORAGE.md). Pair walk: f is
even in x ([0b]), so (b,t) is a hit of (p,s,r) iff (-b,t) is one of (p,s,-r): one walk per pair over both signs of b5
(the bench) emits both ideals' hits, 2 x coset cells/p. model16 reads REC7_TAB at lg+1+log2(lower/upper) ([3]), one
setup6 per pair. Basis preparation is unchanged, <=3995157 narrow and
150271737 wide instructions; six bases (2 for A, 2 per run coset) per2^59.56 ideals from2^48 to2^64, all at
the larger price times2^0.05, fit2^75.09 units.

## Coset walks (P2)

On b_0 = 2k + 1 or t = 2u the kernels run unchanged on the lattice in (k, b_1..b_5) or u (constants halved mod p, one
region axis halved); the target (c - 1)/2 adds 9.0 Ir to setup6 (416.28), charged 10. COS lines: the 122 narrow CAMP
points walked as the two run cosets (X: t = 2u; Y: also b_0 = 2k); REC7_TAB is the envelope over CAMP, CAMPW, CAMPH and
COS. COSL lines: the line walk on b_0 = 2k at 2^25, below LINE_TAB. [10r] checks both.

## Descent

dsieve16 certifies congruences/determinant, enumerates, leaf-tests and ranks hits
(checks throughR256, DESCENT.md). Hit1100 (walk/rank1000, transfer84, apply9); scan5, zero1,
bound2400 per cell; pair setup2^25.1.

## FFT SRT: band strip and g batch

srtf.c (sha256 ac0437e59fb4865a): the proof.md SRT lemma over GMP 6.2.1 (Debian arm64, static; same toolchain).
Exact vs mpz_mod and a gcd strip (archive fft.txt SRTWAVES). FFTPTS N:charged (rounded up): differential
Ir + 1000 per allocator call + realloc copies, x 1.10 (NEON copies <= 8.1%); all <= F(N) = 8.5 N log2(N)^2 + 1200
([10b], audit); charged/F 0.85 at 2^20, 0.46 at 2^26. Above 2^26 limbs (P_g 2^63.53, root ~2^65.2) a product is a
distributed Schoenhage-Strassen multiplication (pieces of ~sqrt(4N) limbs multiplied locally; 6 transposes at the
descent's 84 per record: + 504 N, model16.fft_F). H-FFT: local work <= F; the largest measured c = charged/(N log2 N
log2 log2 N) is 33.31 (2^20): 13.1k per limb at 2^65.2 vs F's 36.1k.
