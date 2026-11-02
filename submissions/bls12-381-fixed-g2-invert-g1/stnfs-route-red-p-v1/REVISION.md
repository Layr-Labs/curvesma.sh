# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels, the descent plan and
H1-H5 are unchanged; kernel prices, filtering basis, block sizes, charges and
reserve are revised below.

## R1 (from 27202622): stronger main-phase ECM buys a smaller region

The yield model multiplies each relation by ECM_FIND^(lambda_f + lambda_g) for
primes in (2^63, 2^69]. Raising ECM_FIND from 0.99 to 0.999 (403 instead of 269
curves per side, B1 = 2^16, B2 = 30 B1, wheel 30030; per-survivor charge
2^38.51 -> 2^39.09, about 5 instructions per cell) gains about 0.066 bits of
yield, spent on a smaller sub-region theta = 87/200 -> 207/500 (1378 canonical
rows, 1285 primitive; smallest a-radius 931, certified count products
0.973023183 / 1.026085506), slacks (298, 229) -> (300, 228). Retained relations
2^64.1311, margin 0.3012 over 1.3 U (pool SE 0.137 bits); relation-weighted
mean norms 747.5/498.9 (descent inputs); two-log expectation 2^117.9122 with
the unchanged plan, meter 2^118.635, Markov 0.60589, boot failure 1.7864e-04,
success >= 0.391949. The in-package 2^12 filtering rerun gives N/U = 0.0715,
below the charged 1.61 x 0.0471 = 0.0758.

## R2 (be446c45, same account): direct chunk scatter replaces line buffer + radix sort

Source routing per retained record: line-buffer append 6, two 16-bit radix
passes 16, page flush 4, directory 2, sentinel and dispatch 8 = 36, whose only
purpose is to group records by chunk = cell >> 30 in stable order. Writing each
record at emission to the open page of its chunk through a 2^29-entry tail table
gives the same order (`store16.scatter` == `store16.route`, checked) with no line
buffer or radix passes: the 7-instruction STORAGE.md listing (audit check
`route_measured`), priced 7 + 4 + 2 + 8 = 21 per record in be446c45 and replaced
by the measured route in R3. REC_TAB stays charged in full, so any buffer write
inside it is counted twice; memory is reported separately (open pages <= 2^49
bytes per line).

Effect: 450.8 -> 414.6 charged instructions per cell; relation collection
2^121.4690 -> 2^121.3481.

## R3 (this package): the whole record route, measured

R2 kept the inherited flush 4 + directory 2 + sentinel/dispatch 8 per record.
The R3 bench program (source in base package 5f39af0) implements the complete route (scatter at emission, page hand-off
to the chunk directory, a full 16-byte copy of every page into the owner's
receive area, owner page walk into apply) and qemu-aarch64 counts its executed
instructions at n and 2n records (measurements12.txt, R3 lines): route
**10.0009** per record beyond the 7-instruction source loop and the 9.000 apply,
of which 7 is the R2 hot path, 3 the copy and 0.0009 page terms (37.7 per page).
Charged **REC_ROUTE = 10.1** (STORAGE.md). Per cell 414.6 -> 388.2; work SE
2.983% (conditional 1.463%, incidence 2.599%), allowance 1.0597 >= 1.05966;
relation collection 2^121.1696 x 1.0597 = 2^121.25329 <= 2^121.2533.

## R4 (this package): explicit mod-r reductions in linear algebra

STORAGE.md (Finishing) gives the per-row computation inside the 16-per-entry
finishing charge with no per-row weight bound: group sums < 2^319, six-limb sign
accumulators, one reduction of v < 2^518 including the nine-limb dense
accumulator by two measured MONT5 products, C_RED = 775.52; per row
1,211.52 <= 16 x 101. Projection inner products are reduced once each and every
generator multiplication is charged 97 + C_RED: generator 2^106.4172 ->
2^109.5864; linear algebra stays <= 2^119.36845.

| core ledger row | R1-R5 | R1-R6 | R1-R7 | R1-R8 | R1-R9 |
| --- | ---: | ---: | ---: | ---: | ---: |
| relation collection | <= 121.252465 | <= 121.252465 | <= 121.25242 | <= 121.25242 | <= 121.25242 |
| linear algebra; generator | 119.334642; 109.470557 | 118.935684; 109.15179 | 118.788193; 109.004035 | 118.788193; 109.004035 | 118.788193; 109.004035 |
| two-log meter | unchanged | unchanged | unchanged | 118.38379 | 118.54 |
| fixed rows | unchanged | unchanged | unchanged | unchanged | unchanged |
| core ledger | 2^121.766352 | 2^121.700199 | 2^121.679361 | 2^121.65114 | 2^121.66817245 (<= 2^121.668173; unrounded 2^121.66817191) |

Yield, margin, descent and success are unchanged by R3 and R4.

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
tail disclosed). E 117.66092 -> 117.81630, boot phi 1.7864e-04 -> 3.7998e-04,
counter 2^-9.70 -> 2^-9.50; meter 118.38379 -> 118.54, Markov 0.60555,
success >= 0.3917; total < 2^121.668894167243655.

## Evidence

`audit_p.py` recomputes every changed quantity (proof.md, Evidence).
