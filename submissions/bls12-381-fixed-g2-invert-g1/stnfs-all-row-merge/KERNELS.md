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
The linewalk uses the same binary64 operations/cells in order (0mismatches over6260cosets/7.86Mrecords); LINE_TAB is the larger of base and group+2SE, cap8.69. Apply/recover are5.76 per record. Whole-program CELLRES runs (no hits): total minus layout, scan and zero is0.096 per cell(2^24->2^26), charged0.5 (CELL_ROUTE, was8). Scan (f array only) is1.376 per
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

The fixed-level rec7 walk preserves the reference's binary64 operation order, slab clips, record order and ring state under the inherited range/precision conditions. Raw archive traces cover270brute-force cosets,432wide pairs and38004campaign cosets with zero mismatches; two deliberately wrong variants fail immediately. It inherits the original recentering/GSO certificates and integer-range assumptions.

`rec7.txt` prices rec7fast at the same130 narrow+42 in-range wide points
(12 wide points below the table excluded); REC7_TAB is their monotone envelope
rounded up to0.01, a step function in log2 hits (bottom4.99; H4 extrapolates by
hits/coset). EXTRA7=2, routing10.01, setup6 per(run,ideal) remain; bench hits are
full-ellipsoid (twice the model's), a margin. Basis preparation is unchanged, <=3995157 narrow and
150271737 wide instructions; two bases per2^59.56 ideals from2^52 to2^64, all at
the larger price times2^0.05, fit2^73.5 units. Runs are defined by `runs16.py`.

## Descent, arithmetic and band strip

`dsieve16` certifies congruences/determinant, enumerates, leaf-tests and ranks hits
(checks throughR256, DESCENT.md). Hit1100 (walk/rank1000, transfer84, apply9); scan5, zero1,
bound2400 per cell; pair setup2^25.1; owner/index bounds atR384 2^62.44/2^34.34 bytes.

Montgomery five-limb multiplication costs 367.76 Ir.

## FFT SRT: band strip and g batch

`srtf.c` (sha256 ac0437e59fb4865a): the supplied model16.py srt construction over GMP 6.2.1 (Debian arm64, static; same toolchain).
Exact vs mpz_mod and a gcd strip (P = primorial(B), B 1e5..3e8, K 1..2^20, planted/limb-edge z, Newton offsets
-50..+50, G = 1, 2): archive `fft.txt` `SRTWAVES` K:|P| bits:mismatches. `FFTPTS` N:charged (rounded up): differential
Ir + 1000 per allocator call + realloc copies, x 1.10 (NEON copies <= 8.1%); all <= F(N) = 8.5 N log2(N)^2 + 1200
([10b], audit); charged/F 0.85 at 2^20, 0.46 at 2^26. Above2^26limbs H-FFT uses distributed Schoenhage–Strassen pieces and six charged transposes (+504N). The measured local envelope and extrapolation remain a heuristic; logs and every FFT/SRT point are in the archive. `gbatch_cost`: g (K = 2^61) 2^24.973 <= C_GBATCH
2^25; band (K = 2^39) 2^23.921 <= C_BATCH 2^24.

## Additional credited rec7 accounting transfer

Credit gbotrel ee8ff14b-9960-4045-80fd-eb0025102a5a(stnfs-k6) and d9002d43-c521-4e63-a641-72414554bba8 for the P7=48 extension and corrected setup refund. Latest native rec7 records have268points(130CAMP,42CAMPW,96CAMPH); the expanded exact REC7_TAB includes the newly measured48..55-bit range. The configuration, candidate predicates, retained yield sums and descent law are unchanged; this changes the charged record walk/setup work. The pool work sums/squares, allowance1.00852 and bothS/G anchor computations are recomputed from the new costs.

For a run withslice sizes S_i and R=sum S_i, the once-per-run setup refund is c(1/R-1/S_i) percell. Since dp already counts the f-side ideals, this refund is OUTSIDE the SIEVE_SIDES=1/2 scaling. Exact summed work is sum_i S_i[c/S_i+c(1/R-1/S_i)]=c. The prior double-halved expression summed to c(T+1)/2. No equal-slice-size premise is needed. The new u7 term is separate from the other side-scaled rec7 savings, and all wider bases/tables are included. No success or matrix theorem is imported from the prior accepted score.

The funded .0355 premerge row envelope is intentionally wider than the earlier.03330 envelope. H3 mean is.0354 withvariance1/U, giving10^8/Ulo rather than10^10/Ulo. The exact ledger includes seven owners:6sends at33,6combinations at64,finish92,fold180,reduction183,2control,7headers at10,routing40 and14sign groups at14 give1345, charged1350. Local56-bit offsets and all accumulator/storage bounds remain valid.

## Exact convolution in the model certificate

The full62node/boot replay uses exact binary-rational nonnegative convolution. Scale each float's exact numerator/denominator to a common power of2, then pack coefficients into base2^(8w), where8w>=bit(maxA)+bit(maxB)+ceil(log2(minlength)). Each product coefficient is strictly below that base, so no coefficient carries into the next slot. Integer multiplication/unpacking is exact; int/int division rounds each output once. No coefficient, mesh cell or subtree is removed. The previous direct sum used at most2*384 roundings per coefficient. Across<=17 convolution/division stages, gamma_(17*(2*384+1))<2^-39; the unchanged2^-30 numerical guard dominates. Retention-tail accumulation preserves its original term order exactly.

The certificate counts464nonzero convolutions, maximum5233input slots,384summands,241denominator bits,256slot bits and operands/products<2^21bits. Thus fewer than2^9multiplications on fewer than2^15 64-bit limbs, with fewer than2^16 Python30-bit digits per operand. Even charging256instructions per schoolbook digit product gives<2^49instructions. Allow2^50including packing/control: below the existing C*2^50 construction row, not free arbitrary-width arithmetic. Every quantity and bound is recomputed, with no stored descent-output substitution.
