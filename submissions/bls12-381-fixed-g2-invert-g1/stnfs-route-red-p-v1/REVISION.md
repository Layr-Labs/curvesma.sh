# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels and their measured prices,
filtering, block Wiedemann block sizes, the descent plan and H1-H5 are unchanged.

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

| core ledger row | R1 + R2 | R1-R4 | R1-R5 |
|---|---:|---:| ---: |
| relation collection | <= 121.3481 | <= 121.2533 | <= 121.252465 |
| linear algebra; generator | 119.36845; 106.4172 | 119.36845; 109.5864 | 119.334642; 109.470557 |
| two-log meter | 118.635 | unchanged | unchanged |
| fixed rows | 90, 88, 80, 72, 50, 11 | unchanged | unchanged |
| core ledger | <= 2^121.8398 | 2^121.773284 (<= 2^121.7733; unrounded 2^121.7732726) | 2^121.76635122 (<= 2^121.766352; unrounded 2^121.76635083) |

Yield, margin, descent and success are unchanged by R3 and R4.

## R5 (this package): exact prices, repriced route/finishing, 23x4 blocks, 0.05% reserve

Line/record/norm-bound prices are the exact raw-count quotients (`prices_exact`;
e.g. 3828670/123056 = 31.113233). Route 10.1 -> 10.01 (measured 10.0015, 0.1%
margin). LA finishing at the explicit R4 budget (1,211.52/row) as 13/entry
(1,313, 8.4% margin): C_ENTRY 85; the 33/97 loops are static counts of the
shipped assembly. Blocks re-minimised to (23, 4): the minimum of the per-N^2
expression over m 1..128, n 1..40 (recompute from model16 in three lines).
Empty-coset branch never triggers (min log2h 2.32 vs floor -1.09;
`empty_never`). Reserve 1% -> 0.05%: residual rounding is bounded by half an
ulp per unshipped constant weighted by core share (2.99e-5, `reserve_residual`),
covered tenfold. New total: (2001/2000) x 2^121.766352 + 2^80/20000 <
2^121.767073167243655. Yield, margin, descent, success unchanged.
route_bench.c (R3 bench source) is dropped to fit the package cap; its counts
(R3 lines) and recomputation (`route_measured`) stay, and the source is
byte-identical in declared base 5f39af0. bw_dot.s keeps only the charged
_bw_full loop (the uncharged _bw_dot function is dropped for the same cap).

## Evidence

`audit_p.py` recomputes every changed quantity (proof.md, Evidence).
