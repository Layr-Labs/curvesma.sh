# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels and their measured prices,
filtering, the descent plan and H1-H5 are unchanged.

## R1 (from 27202622): stronger main-phase ECM buys a smaller region

ECM_FIND 0.99 -> 0.999 for primes above X2 (403 instead of 269 curves per side;
per survivor 2^38.51 -> 2^39.09) gains about 0.066 bits of yield, spent on
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

Scatter 7 + owner's page copy 13 per 4 records (STORAGE.md listings) + page
hand-off and directory walk < 64 per 2^16-record page = 10.26, charged
**10.3** (was 21); recover-P's executed-instruction count (gcc 13, qemu-aarch64)
gives 10.0009. Per cell 414.6 -> 388.7; allowance 1.0597; relation collection
<= 2^121.2551.

## R5 (withdrawn): survivor-gated recovery

Every record is re-read at 7.003, as in the source.

## R6 (from the recover-Q thread): record walk to 2^64, region 19/50

Primes in [2^63, 2^64) move from ECM to the record walk; `slice12` is measured
only below 2^63, so their record-walk hits, empty cosets, apply and setup6 are
charged **2x** (WIDE_REC = WIDE_SETUP = 2). The yield buys theta 207/500 ->
19/50 (1040 canonical rows, a-radius >= 892, count products 0.971857579 /
1.027240120): relations 2^64.1311, margin 0.3013, mean norms 746.0/499.0,
relation collection <= 2^121.1492 at slacks (291, 225). Descent at the new
norms: two-log 2^117.91615, Markov 0.60758, boot 1.8144e-04, E[N] 2^17.60,
counter 2^-9.69, success >= 0.390235.

## R7 (from the recover-Q thread): line sieve from 2^11

The line sieve starts at 2^11 instead of 2^8 (`P_MIN` = 2048, `LX0` = 11,
config `prime_min` 2048, passed to both `Retention` objects). `threshold_model`
treats p < prime_min as the unsieved small group counted against the
slacks, so the cheaper line path costs yield, recovered by slacks
(291, 225) -> (297, 233) (784-point grid; cofactors 305/241 bits).
Relations 2^64.1311, margin 0.3013, mean norms unchanged (746.0/498.9, so the
descent is unchanged). Per cell 417.9; work SE 3.377%, allowance 1.0676;
relation collection 2^121.0173 x 1.0676 = 2^121.11165 <= 2^121.1117.

| core ledger row | R1 + R2 | R1-R3 | R1-R4 | + R6 | + R6, R7 |
|---|---:|---:|---:|---:|---:|
| relation collection | <= 121.3481 | <= 121.3481 | <= 121.2551 | <= 121.1492 | <= 121.1117 |
| linear algebra | <= 119.36845 | <= 118.7428 | <= 118.7428 | <= 118.7428 | <= 118.7428 |
| core ledger | <= 121.8398 | <= 121.7452 | <= 121.6752 | <= 121.5968 | **<= 121.5694** (rows sum 121.569314; unrounded 121.5692664) |

## Relation count used by the descent (guard and sensitivity)

The descent's missing-log model takes LOG2_RELS = 64.131, the floor of the
pool point estimate 2^64.13109 (pool SE 0.149 bits). The algorithm enforces it:
relation collection continues over further sub-region slices until 2^64.131
retained relations exist, so the descent input holds by construction (H2'). The
work charge uses the point-estimate yield; a yield 1 SE (2 SE) lower would need
2^0.149 (2^0.298) times the collection work, a core ledger of about 2^121.68
(2^121.79). Disclosed under H2'; not covered by the allowance, which bounds the
work estimator only. The grid optimum lies on a flat ridge (neighbours
within 0.01 bits), so grid selection moves the estimate far less than its SE.

Generator 105.5088, meter 118.635 and fixed rows unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
