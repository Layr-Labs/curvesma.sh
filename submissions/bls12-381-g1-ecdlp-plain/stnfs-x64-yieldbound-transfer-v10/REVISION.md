# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, kernel prices, filtering, descent plan and H1-H5
are unchanged.

## R1 (from 27202622): stronger main-phase ECM buys a smaller region

ECM_FIND 0.99 -> 0.999 above X2 (403 curves, was 269) gains about 0.066 bits, spent on
theta 87/200 -> 207/500 and slacks (298, 229) -> (300, 228); descent re-derived
at the new inputs (Markov 0.60589, success >= 0.391949); 2^12 filtering rerun
N/U 0.0715 < charged 0.0758.

## R2 (recover-P 43d603a4/be446c45): direct chunk scatter

Each record is written at emission to its chunk's open page through a 2^29-entry
tail table (same stable order: `store16.scatter` == `store16.route`; 7-instruction
hot path, `route_hot_path`), replacing the line buffer and two radix passes:
36 -> 21 per record, per cell 450.8 -> 414.6, relation collection 2^121.3481.

## R3 (this package): unit entries and a combined row tail in the matrix product

88 -> **44** per sparse entry (STORAGE.md, Linear algebra): +-1 entries in an
11-instruction loop (`_group_run`, new here; unit share 0.90-0.95 measured with
multiplier tracking, charged 0.88, H4u) and recover-Q 06aac4a1's combined row
tail (2,090 per row). Block Wiedemann (18, 4): linear algebra 2^119.36845 ->
**2^118.74271**, generator 2^106.4172 -> 2^105.5088.

## R4 (recover-P thread): measured end-to-end record route

Scatter 7 + page copy 13/4 + page hand-off < 64/2^16 = 10.26 per record (STORAGE.md),
charged **10.3** (was 21; recover-P measured 10.0009 under qemu-aarch64).
Relation collection <= 2^121.2551.

## R5 (withdrawn): survivor-gated recovery

Every record is re-read at 7.003, as in the source.

## R6 (from the recover-Q thread): record walk to 2^64

Primes in [2^63, 2^64) move from ECM to the record walk; `slice12` is measured
only below 2^63, so their record-walk hits, empty cosets, apply and setup6 are
charged **2x** (WIDE_REC = WIDE_SETUP = 2, KERNELS.md). At theta 19/50 and
slacks (291, 225): relation collection <= 2^121.1492.

## R7 (from the recover-Q thread): line sieve from 2^11

`P_MIN` = 2048, `LX0` = 11, config `prime_min` 2048, passed to both `Retention`
objects; `threshold_model` counts p < prime_min as the unsieved small group
against the slacks. At slacks (297, 233): relation collection <= 2^121.1117.

## R8 (this package): relation input to the descent at a yield lower bound

The descent's missing-log input LOG2_RELS = 64.131 must hold at a linear 2-SE
lower confidence bound on the yield: rels x (1 - 2 x relative pool SE) >= 2^64.131
(`model16.YIELD_SE_K`, audit `yield_lower_bound_covers_descent`). 300-point grid:
theta 19/50 -> 56/125
(1722 canonical rows, a-radius >= 969, count products 0.974069635 /
1.025050907), slacks (299, 237) (cofactors 307/245 bits). Relations 2^64.4074,
pool SE 0.1250 bits, lower bound 2^64.1328, margin over 1.3 U 0.5775. Per cell
423.3 over 2^127.122 cells; work SE 2.663%, allowance 1.0533; relation
collection 2^121.4850 x 1.0533 = 2^121.55989 <= 2^121.5599. A guard continues
collection until 2^64.131 relations exist; the charged region covers this
unless the yield is more than 2 SE below its estimate. Descent at mean norms 749.0/499.0: two-log
2^117.90781, boot 1.7558e-04, E[N] 2^17.593; meter 2^118.635 -> **2^118.627**
(Markov 0.60744), success >= 0.390390.

| core ledger row | R1-R3 | R1-R4 | + R6 | + R6, R7 | + R6-R8 |
|---|---:|---:|---:|---:|---:|
| relation collection | <= 121.3481 | <= 121.2551 | <= 121.1492 | <= 121.1117 | <= 121.5599 |
| linear algebra | <= 118.7428 | <= 118.7428 | <= 118.7428 | <= 118.7428 | <= 118.7428 |
| core ledger | <= 121.7452 | <= 121.6752 | <= 121.5968 | <= 121.5694 | **<= 121.908** (rows sum 121.907975; unrounded 121.9079584) |

Generator 105.5088 and fixed rows unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
