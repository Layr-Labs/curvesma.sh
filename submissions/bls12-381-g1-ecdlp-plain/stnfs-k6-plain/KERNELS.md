# Kernels
Counts are differential callgrind Ir on linux/arm64, Debian12, GCC12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, Valgrind3.19. These count
instructions, not cycles. `kernels_asm.b85` is base85/lzma text sections: listings,
sources, patches and logs. `model16.py` holds charged constants. No SIMD is used.
## Main collection
The slice12 kernel uses exact adjugate-remainder recentring. `measurements12.txt`
records all checks MATCH. Record walk prices are exact quotients in REC_TAB;
setup6 is406.87+72+21.13 instructions, empty cosets1600. Line path: slice12 +
`linewalk.patch` (fixed five-level walk, pointer inner run unrolled x4, n=5 setup;
comments dropped, .text byte-identical to the measured source, sha256 bae54e70).
Refinement: the reference's binary64 operations (explicit fma, fcvtas=llround,
fsqrt, fmaxnm) on the same operands in order, hence the same cells in order;
0 mismatches over 6260 cosets (bitwise setups, 7.86M trace records, brute force). LINE_TAB: the larger of base TAB and pooled
groups+2SE (`measurements12x.txt`); cap8.69 above2^18.78. Apply/recover are5.76 per record. Whole-program CELLRES runs (no hits): total minus layout, scan and zero is0.096 per cell(2^24->2^26), charged0.5 (CELL_ROUTE). Scan (f array only) is1.376 per
cell plus37.62 per leaf, covering every leaf length1..64 and256; `slice12oh.c`
and `measurements-oh.txt`. Packed16-bit sums and thresholds are below2^15, since
each sum is<=16*851. Layout is28.17 per suffix entry plus2^21 per slice and2^17
per2^32-byte segment; logs include allocation/zeroing. Bounds (f bound only, bit-identical; `measurements12.txt`) cost1154.6 per block,
4600 per row and17500 per slice; warm bisection changes only charged counts.
For2^63<=p<2^64, `wide12.c` prepares a __int128 lattice (entries guarded at2^40),
then the unchanged walk; `wide12.txt`: identical .text, exact reductions, 400/400 planted
hits, MATCH. setup6w=141.00; REC_TAB_W x1.01, setup6 and1600 per empty coset cover these;
two preparations at5615924 per ideal fit the2^70.6 row.
## rec7: fixed-level DFS
rec7 enumerates L7={(b,t):t=r*b(s) mod p} over one canonical a-run, with a slab
t_lo<=t<=t_hi and the b metric in `slice_mb.h`. Source is original rec7.c
(sha25604a2ab937a56eebf69425f2738d7a5b889f9900f4e3e113459946a4b6852a153)
plus `rec7.patch`. The patch includes the wide basis preparation and this
fixed-level walk. Prepare7w uses __int128 rows, quad-precision GSO, a2^40 guard,
congruence/determinant checks, then unchanged setup. A failed basis takes the
slice walk. Two-limb recentring uses the stored representative of2^32*e6;
c_hi*R32 fits63 bits for reduced entries below2^28, as in all logged checks
and campaign runs (the inherited rec7 range assumption).
**Refinement lemma.** The walk (fixed five-level loops of the reference DFS, then
rec7fast: run sign folded once per run of >=3 points with k never 0, record stores
as fixed asm blocks, fsqrt/fmaxnm/fcvtas=llround single instructions, branch-free
slab clip) takes each level's integer z in the reference's ascending order with the
same binary64 operations on the same operands (explicit fma), so the same leaf
intervals, slab clips, ring_wrap calls and record stores (segment, slot, cell, id)
occur in the same order, with the same hit count. It inherits the reference's
recentring, GSO/bound precision and integer-range conditions (|z| < 2^53); BND_EPS and integer certificates are unchanged.
Checks, 0 mismatches: 270 brute-force cosets (2087 records), wide 432 pairs
(1728/1728 planted), a trace against the shipped walk on all 38004 cosets of the
184 campaign points (bounds, binary64 states, leaves, every record in order), equal
final rings, and on 192 cosets of the 96 CAMPH points (TCHECK); two wrong variants fail the trace on their first coset.
`rec7.txt` prices rec7fast at the same130 narrow+42 in-range wide points
(12 wide points below the table excluded) plus all 96 CAMPH points (p 2^47.03..2^57, h' 2^23..2^30; pool runs
have h'<=2^28.01 at p>=2^48); REC7_TAB is their monotone envelope
rounded up to0.01, a step function in log2 hits (bottom4.99; [10r] checks every point). EXTRA7=2, routing10.01, setup6 per(run,ideal) remain; bench hits are
full-ellipsoid (twice the model's), a margin. Basis preparation is unchanged, <=3995157 narrow and
150271737 wide instructions; two bases per2^59.56 ideals from2^48 to2^64, all at
the larger price times2^0.05, fit2^73.5 units. Runs are defined by `runs16.py`.
## Descent, arithmetic and band strip
`dsieve16` certifies congruences/determinant, enumerates, leaf-tests and ranks hits
(checks throughR256, DESCENT.md). Hit1100 (walk/rank1000, transfer84, apply9); scan5, zero1,
bound2400 per cell; pair setup2^25.1; owner/index bounds atR384 2^62.44/2^34.34 bytes.
Montgomery five-limb multiplication costs 367.76 Ir.
## FFT SRT: band strip and g batch
`srtf.c` (sha256 ac0437e59fb4865a): the proof.md SRT lemma over GMP 6.2.1 (Debian arm64, static; same toolchain).
Exact vs mpz_mod and a gcd strip (P = primorial(B), B 1e5..3e8, K 1..2^20, planted/limb-edge z, Newton offsets
-50..+50, G = 1, 2): archive `fft.txt` `SRTWAVES` K:|P| bits:mismatches. `FFTPTS` N:charged (rounded up): differential
Ir + 1000 per allocator call + realloc copies, x 1.10 (NEON copies <= 8.1%); all <= F(N) = 8.5 N log2(N)^2 + 1200
([10b], audit); charged/F 0.85 at 2^20, 0.46 at 2^26. Above 2^26 limbs (P_g 2^63.53, root ~2^65.2) a product is a
distributed Schoenhage-Strassen multiplication (pieces of ~sqrt(4N) limbs multiplied locally; 6 transposes at the
descent's 84 per record: + 504 N, `model16.fft_F`). H-FFT: local work <= F; the largest measured c = charged/(N log2 N
log2 log2 N) is 33.31 (2^20): 13.1k per limb at 2^65.2 vs F's 36.1k. `gbatch_cost`: g (K = 2^61) 2^24.973 <= C_GBATCH
2^25; band (K = 2^39) 2^23.921 <= C_BATCH 2^24.
