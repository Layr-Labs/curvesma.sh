# Kernels

Counts are differential callgrind Ir on linux/arm64, Debian12, GCC12.2 -O3
-fno-tree-vectorize -fno-tree-loop-distribute-patterns, Valgrind3.19. These count
instructions, not cycles. `kernels_asm.b85` is base85/lzma text sections: listings,
sources, patches and logs. `model16.py` holds charged constants. No SIMD is used.

## Main collection

The slice12 kernel uses exact adjugate-remainder recentring. `measurements12.txt`
records all checks MATCH. Record walk prices are exact quotients in REC_TAB;
setup6 is406.87+72+21.13 instructions, empty cosets1600. Line prices LINE_TAB are
the larger of the base measurement and pooled random-pair groups plus two SE,
interpolated in log2 hits: `measurements12x.txt` includes every group and checks.
Its high-density cap is17.2. Apply/recover are5.76 per record. Scan is3.002 per
cell plus90.98 per leaf, covering every leaf length1..64 and256; `slice12oh.c`
and `measurements-oh.txt`. Packed16-bit sums and thresholds are below2^15, since
each sum is<=16*851. Layout is28.17 per suffix entry plus2^21 per slice and2^17
per2^32-byte segment; logs include allocation/zeroing. Bounds cost1508.2 per block,
4600 per row and17500 per slice; warm bisection changes only charged counts.

For2^63<=p<2^64, `wide12.c` prepares a __int128 lattice, guards entries at2^40,
then uses unchanged LLL/walk. `wide12.txt`: common .text, main and relocations
byte-identical;10^6 reductions exact;200 pairs certified,400/400 planted hits,
brute-force MATCH. setup6w=141.00, per-coset spread<=2.17. REC_TAB_W times1.01,
ordinary setup6 and1600 per empty coset cover these. Six prime-root lattices,
two preparations at5615924 per ideal, fit the2^70.6 ledger row.

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

**Refinement lemma.** A reference DFS frame at level k consists of its centre
C[k], remaining squared radius, slab centre TP[k], and four integer affine
states. Its bounds are exactly LEVEL_BOUNDS(k). The fixed-level code initializes
the same frame, takes its integer z in the same ascending order, and computes
the same remainder, slab centre, child centres and affine states. At entry to a
child, these equal the reference values: affine states start at base+(lo-1)*d
and advance by d before processing lo. Induction down k=6..1 proves that the
same leaf intervals are visited in the same order. At k=1 the slab clipping,
sign fold, cell/routing index and record emission are unchanged. Thus the new
walk is a control-flow specialization, not new pruning or an approximation.
It inherits the reference's recentring, GSO/bound precision and integer-range
conditions; all loop endpoints and affine intermediates stay in signed64-bit
ranges. Production endpoints are below LLONG_MAX, so the final loop increment
is defined. No altered tolerance or reassociation is part of the algorithm:
expression order and BND_EPS remain the reference's. Integer congruences and
certificates are unchanged.

Checks:270 brute-force cosets at slice sizes2^9..2^11,T=1..33,five slabs,p97..4099,
four stretches:0 mismatches,2087 records. Wide check:432 pairs in[2^62,2^64),
1728 planted points found, every record valid. A dual-walk trace build compared
each level's bounds, binary64 centre, remainder and slab centre in traversal
order: all184 campaign points and70599411 hits,0 discrepancies. Tests are scoped, not a production run.

`rec7.txt` prices this build at the same130 narrow+42 in-range wide points
(12 additional wide points below the table excluded). REC7_TAB is their
monotone cost envelope, rounded upward to0.01, a step function in log2 hits,
bottom4.99. H4 extrapolates that the cost depends on hits/coset at production
scale. EXTRA7=2, routing10.01 and setup6 once per(run,ideal) remain. Bench
full-ellipsoid hits are twice the model's half-slice hits, an existing margin.
Runs outside32..1024 and contaminated runs take the slice walk with unchanged
contamination allowance. Basis preparation is unchanged, <=3995157 narrow and
150271737 wide instructions; two bases per2^59.56 ideals from2^52 to2^64, all at
the larger price times2^0.05, fit2^73.5 units. Runs are defined by `runs16.py`.

## Descent, arithmetic and band strip

`dsieve16` certifies congruences/determinant, enumerates, leaf-tests and ranks
hits. DESCENT.md and logs give checks throughR256. Hit price1100 includes1000
walk/rank,84 transfer and9 apply (rounded); scan5,zero1,bound2400 per cell;
pair setup2^25.1. Owner/index bounds atR384 are2^62.44/2^34.34 bytes; H4
extrapolates costs toR384.

Montgomery five-limb multiplication costs367.76 Ir. `toom3.txt` has108
schoolbook matches at47..2000 limbs and a Python match at20000. Differential
measurements at28 sizes8..34992 fit `model16.batch_cost`; worst ratio0.9857.
H-BT extrapolates that recurrence beyond measured sizes. SRT correctness and
precision are the exact lemma in proof.md, not H-BT. Its256 planted832-bit
survivors match trial division by995846 band primes; guardG=L+2, minimumG=1
also tested. Band splitting remains separately charged in proof.md. SRT,
descent and survival probabilities did not change with DFS specialization.
