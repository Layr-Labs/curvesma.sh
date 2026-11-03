# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels, the descent plan and
H1-H5 are unchanged; kernel prices, filtering basis, block sizes, charges and
reserve are revised below. Full R1-R4 text in submission 5f39af0 and git.

## R1 (from 27202622): stronger main-phase ECM buys a smaller region

ECM_FIND 0.99 -> 0.999 (403 curves/side, B1 = 2^16, B2 = 30 B1; +5 instr/cell)
gains 0.066 yield bits, spent on theta = 87/200 -> 207/500 (1378 canonical
rows, radius 931, count products 0.973023183/1.026085506), slacks (300, 228).
Rels 2^64.1311, margin 0.3012; mean norms 747.5/498.9; E 2^117.9122, meter
2^118.635, success >= 0.391949.

## R2 (be446c45): direct chunk scatter replaces line buffer + radix sort

Records scatter at emission to the open page of chunk = cell >> 30 (7-instr
STORAGE.md hot path, `store16.scatter` == `store16.route`): 36 -> 21 per record
(7 + flush 4 + directory 2 + dispatch 8; R3 measures the tail). Per cell
450.8 -> 414.6; relcol 2^121.4690 -> 2^121.3481.

## R3: the whole record route, measured

qemu-aarch64 counts the complete route (bench source in 5f39af0): 10.0009 per
record beyond source 7 + apply 9 (R3 lines in measurements12.txt). Charged 10.1
(R5 reprices to 10.01, R7 to 10.005). Per cell 414.6 -> 388.2; allowance 1.0597;
relcol <= 2^121.2533.

## R4: explicit mod-r reductions in linear algebra

Per-row finishing made explicit (group sums < 2^319, one C_RED = 775.52
reduction of v < 2^518): 1,211.52 <= 16 x 101. Projections reduced once each,
generator at 97 + C_RED: generator 2^106.4172 -> 2^109.5864; LA <= 2^119.36845.

Per-revision ledger columns (R1-R17) live in git history; R18 drops the table
for the package cap. Current ledger: source-core.md section 8.

Yield, margin, descent and success are unchanged by R3 and R4. R5-R9 columns in 339cb2b7 and git.

## R5: exact prices, 10.01 route/finishing, 23x4 blocks, 0.05% reserve

Exact raw-count prices (`prices_exact`). Route 10.01, finishing 13/entry
(C_ENTRY 85; 33/97 static). Blocks (23, 4). Empty branch dead. Reserve 0.05%
(residual 2.99e-5). route_bench.c, _bw_dot, _row_dot dropped for the cap.

## R6: refiltered linear algebra

Old N/U basis (0.0471) irreproducible (0.0497-0.0551). 33 fresh 2^12-2^18 sims; per-scale N/U max (mean) in git (means fall monotonically). Five-seed 2^17 basis; PAD 1.2 -> N/U 0.0662 (covers 28 runs >= 2^13). W_ROW 100.7 (worst 100.61).

## R7: budgeted transfers, 10.005 route, 22x4 blocks

Transfers: <= 28 groups/row of 5 limbs over depth-<= 5 tree; 12/hop, 8/group,
20 dispatch: <= 1,924/row, charged 20/entry (4.7% margin). C_ENTRY 73; blocks
(22, 4). Route 10.005 (measured 10.0015).

## R8: tighter descent kernels, lower meter

SIEVE_HIT 1100->850 (covers 793 max; radii>=192 cost<=418).
BOUND_CELL 2400 -> 2000 (33% over the ~1500 estimate). E[two logs] 117.9122
-> 117.66092; meter 118.635 -> 118.38379 at the same Markov gap, so success
is unchanged (all other success terms identical). New total (2001/2000) x
2^121.65114 + 2^80/20000 < 2^121.651861167243655.

## R9: rigorous loss, one-SE rels allowance, meter 118.54

loss hoists Pss_hi at s+ (DESCENT.md; values same). Rels 64.131->64.063 (one pool SE below; 16% tail disclosed). E->117.81630, meter 118.54, success >=0.3917; total < 2^121.668894167243655.

## R10: corrected loss proof (numbers unchanged)

R9's hoisting had monotonicity backwards (review). R10 bounds loss by num/den (DESCENT.md; matches R9 to 6e-05). Numbers stand.

## R11: 2^19 refilter + retention allowances

Eight 2^19 runs (worst 0.0491); PAD 1.25->0.0614, N 2^59.4237 (LA 118.570161, gen 108.828524). Prune/_DOWN/banding per DESCENT.md. Core <= 2^121.640393, total < 2^121.641114167243655.

## R12 (muse-spark-1.3-contributor): descent plan re-optimised at current inputs

Coord search over all entries+boot (mu>=16,phi<=0.05,psi<=0.8): 18 move, boot same. E[two logs]
117.81629554 -> 117.76547297 (-0.05082); meter 118.54 -> 118.48918 at the same
0.72371 gap (Markov 0.60554); boot phi 3.7997e-04 -> 3.2307e-04, success
0.39172 -> 0.39187. Core <= 2^121.634559, total < 2^121.635280167243655.

## R13 (muse-spark-1.3-contributor): retention allowances union-bounded

The review refuted R11's own-mass `_DOWN` weights and above-s banding. R13
exact t,k<=19 conv to S<=19 + union rest 4.9e-41, _DOWN 1e-50, _REM 1.0 (DESCENT.md). Loss 8.6e-7->1.0e-6; E 7th-dec; rest same.

## R14: sieve to 2^64, smaller region/slacks, retention numerical fix

lx2 63->64, theta 207/500->19/50, slacks (300,228)->(293,222): relcol
121.07769, E 117.80618, meter 118.52989, margin 0.3007, success 0.39170,
LOG2_RELS 64.050. K-tail geometric majorant (rest 6.9e-41). Core <= 2^121.507523, total < 2^121.508244167243655.

## R15: line/record split 2^20 -> 2^48 (charged-cost crossover)

Line (streaming, no route/apply/recover) carried p < 2^20; record (emit +
route + apply + recover) carried [2^20, 2^64). Both tabs price by hits/coset h: line h>=2^6.79 for p<2^48 (audit line_measured), in measured range, flat above; secant overestimates convex-decreasing cost. P_MID 20 -> 48: per-cell work 408.2 -> 328.4;
yield/margin/descent/LA/success untouched; SE re-derived (1.724%/3.019%/
3.477%), allowance 1.0676 -> 1.0696. Core <= 2^121.281215, total <
2^121.281936167243655.

## R16: survivor band-strip recovers the R15 line-credited (2^20, 2^48] primes

R15's line path credits band primes unrecorded; R16 recovers them on survivors: two fixed ECM rounds/side (<= 2^48-bit at 99.95%: B1 = 2^14, 62 curves; 2^38.1359), factors divided out. Unstripped cofactors <= 832/576 bits; strip failure <= 41/28 x (5e-4)^2/side (union), counted lost (yield x0.99998275). C_surv 2^39.0895->2^39.6899; per-cell 328.4->336.0; SE re-derived, allowance 1.0696->1.0700; success same. Core <= 2^121.30450, total < 2^121.305221167243655.

## R17: descent plan re-optimised at current inputs; ECM grid to odd log2 B1

Coord search as R12 at R14-R16 inputs + pairwise probes: 8 move, boot same, E 117.80618->117.78173. The `ecm16.py`
B1 grid takes all log2 B1 in 10..36 (was even only), moving the main 69-bit
schedule to B1 = 2^15 (764 curves; 2^38.0758/2^38.0742) and the band 48-bit
schedule to B1 = 2^13 (101 curves; 2^36.4645/2^35.4050). C_surv 2^39.6899 ->
2^39.6451, per-cell 336.0 -> 335.3, relcol 120.79404; E 117.78162, meter
118.50533 at the same Markov gap (0.60554); boot phi 3.8889e-04 -> 3.6910e-04,
success 0.391700 -> 0.391750.
Also corrects two stale DESCENT.md figures (boot trials, lattices).
Core <= 2^121.29889, total < 2^121.299611167243655.

## R18: excess 1.2, new filter basis, re-optimized descent/band/grid

Excess 1.2; 2^19 basis (N/U <= 0.0556) + trend (worst 0.0644 >= 2^15);
N/U 0.0695, N 2^59.6045, LA 118.931853. 16 owners, depth-4 tree.
Theta 719/2000, slacks (293,224); rels 64.01555, margin 0.30119;
relcol 120.61756. Descent 63.925, 744.8/498.9; Markov co-opt (21+boot):
E 118.03397, meter 118.76810, success 0.391913. Band 99.9%; B2/B1=20 grid.
Core <= 2^121.28540, total < 2^121.286121167243655.

## R19: filter pad 1.2 on a 16-seed exact-means basis

Sixteen 2^19 sims at exact means (worst N/U 0.0552, max coef
49) + 2^12-2^20 trend (worst 0.0644 >= 2^15). Charge 1.2x worst
= 0.0663 (R6 precedent); N 2^59.6045 -> 2^59.5360, LA -> 118.794879,
gen -> 109.009417. R18 layout kept. Rest untouched.
Core <= 2^121.25958, total < 2^121.260301167243655.

## R21: triples plan re-search

(R20 UDUP 0.4 rejected on H5 support; stratified relation d/r 0.44/0.28
min -- 0.25 stands. Base reverted to R19.) Coordinate descent with
(R,n,kA) triples at UDUP 0.25 (28 moves): E 118.0192, meter 118.75286.
Rest untouched.
Core <= 2^121.25689, total < 2^121.257611167243655.

## R22: storage occupancy recount, blocks (21, 4)

N 2^59.5360 < 12x2^56: blocks 12..15 empty, 24 sums/row (not 32).
Finishing 1,163.52 <= 12x100.7 (3.9%); transfers 1,364 <= 15x100.7 (10.7%).
C_ENTRY 73->67; blocks (22,4)->(21,4); LA 118.714572, gen 108.85662.
Core <= 2^121.24261, total < 2^121.243331167243655.

## R23: filter basis at relation-weighted means, pad 1.35

Sim inputs nf_mean->nf_rel (744.8405/498.9650): filter eats relations,
not cells (+19 f-bits overstated). 16 2^19 sims (worst 0.0463,
maxcoef 50) + trend (2^15 n=30: 0.0602; 2^20: 0.0445).
Charge 1.35x = 0.0626 (covers all, +3.9%).
N 2^59.5360 -> 2^59.4527, LA 118.548101, gen 108.722616, BW (21,4).
12-block count re-verified (11 occupied). [6] to nf_rel.
Core <= 2^121.21506, total < 2^121.215781167243655.

## R24: occupancy recount to 11 blocks, C_ENTRY 66

N 2^59.4527 < 11x2^56 = 2^59.4595: 22 sums/row (not 24). Finishing
1,139.52 <= 12x100.7 (6.0%); transfers 1,252 <= 14x100.7 (12.6%; 13/entry
= 4.5%, under R7's 4.7%, rejected). C_ENTRY 67->66; BW (21,4).
LA -> 118.534187, gen -> 108.722616.
Core <= 2^121.21288, total < 2^121.213601167243655.

## Evidence

`audit_p.py` recomputes every changed quantity (proof.md, Evidence).
