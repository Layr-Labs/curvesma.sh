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

| core ledger row | R1-R9 | R1-R11 |
| --- | ---: | ---: |
| relation collection | <= 121.25242 | <= 121.25242 |
| linear algebra; generator | 118.788193; 109.004035 | 118.570161; 108.828524 |
| two-log meter | 118.54 | 118.54 |
| fixed rows | unchanged | unchanged |
| core ledger | <= 2^121.668173 (unrounded 2^121.66817191) | <= 2^121.640393 (unrounded 2^121.64039187) |

Yield, margin, descent and success are unchanged by R3 and R4. R5-R9 columns in 339cb2b7 and git.

## R5: exact prices, 10.01 route/finishing, 23x4 blocks, 0.05% reserve

Exact raw-count prices (`prices_exact`). Route 10.01, finishing 13/entry
(C_ENTRY 85; 33/97 static). Blocks (23, 4). Empty branch dead. Reserve 0.05%
(residual 2.99e-5). route_bench.c, _bw_dot, _row_dot dropped for the cap.

## R6: refiltered linear algebra

Old N/U basis (0.0471) irreproducible (0.0497-0.0551). 33 fresh sims at
2^12-2^18; N/U max (mean): 0.0764 (0.0733), 0.0631 (0.0588), 0.0631 (0.0580),
0.0576 (0.0551), 0.0595 (0.0546), 0.0551 (0.0523), 0.0509 (0.0501); means fall
monotonically. FILT_TABLE is the five-seed 2^17 basis; PAD 1.2 -> N/U 0.0662
(covers 28 runs >= 2^13). W_ROW 100.7 (worst 100.61). Rerun: shipped filtsim.

## R7: budgeted transfers, 10.005 route, 22x4 blocks

Transfers: <= 28 groups/row of 5 limbs over depth-<= 5 tree; 12/hop, 8/group,
20 dispatch: <= 1,924/row, charged 20/entry (4.7% margin). C_ENTRY 73; blocks
(22, 4). Route 10.005 (measured 10.0015).

## R8: tighter descent kernels, lower meter

SIEVE_HIT 1100 -> 850 (enum/record 1000 -> 800: covers the 793 measured max;
plan radii are all >= 192 where costs are <= 418 and falling in R).
BOUND_CELL 2400 -> 2000 (33% over the ~1500 estimate). E[two logs] 117.9122
-> 117.66092; meter 118.635 -> 118.38379 at the same Markov gap, so success
is unchanged (all other success terms identical). New total (2001/2000) x
2^121.65114 + 2^80/20000 < 2^121.651861167243655.

## R9: rigorous loss, one-SE rels allowance, meter 118.54

`retain16.loss` hoists Pss_hi at s+ out of the tail sum (DESCENT.md proof;
values unchanged, flat hi). Rels input 64.131 -> 64.063: one pool SE below
the central estimate (0.51 SE in the safe factors + 0.5 SE explicit; 16%
tail disclosed). E -> 117.81630, phi -> 3.7998e-04, counter 2^-9.50;
meter 118.54, Markov 0.60555, success >= 0.3917; total < 2^121.668894167243655.

## R10: corrected loss proof (numbers unchanged)

R9's hoisting proof had the smoothness monotonicity backwards (review caught it:
Pss(m-x) increases in x). R10 bounds loss by num/den with the original tail sum
and den = Pss_lo(m)(1-P_up(X>s)) (DESCENT.md proof; values match R9 to 6e-05).
All R9 numbers stand.

## R11 (muse-spark-1.3-contributor): 2^19 refilter + retention allowances

LA: FILT_TABLE becomes eight fresh 2^19 runs (worst N/U 0.0491, weights <=
100.03, maxcoef <= 39); PAD 1.25 -> N/U 0.0614, N = 2^59.4237, covering the R6
sims at 2^15-2^18 and all 16 new 2^18/2^19 runs (worst 0.0595; means fall
monotonically 0.0733 -> 0.0480 over 2^12-2^19). W_ROW 100.7 kept; blocks (22,4)
re-verified optimal at the new N. LA 118.788193 -> 118.570161, generator
109.004035 -> 108.828524. Retention: R10 proof kept; the prune now deletes
sub-threshold atoms and counts them in _DOWN at their own mass (no moves, no
teleport), local valuations tabulate exactly to 400 bits with the ~1e-61
remainder counted at 1.0, and the omitted large-block mass is resolved into
bands at hi(band top) under the same block model (DESCENT.md). Losses match R10
to 4 figures; E, success and the meter are unchanged by the hardening. Core <=
2^121.640393, total < 2^121.641114167243655. Also fixes stale doc numbers
(STORAGE 10.005, KERNELS 800, model docstring, margin/relcol/E[N] figures).

## Evidence

`audit_p.py` recomputes every changed quantity (proof.md, Evidence).
