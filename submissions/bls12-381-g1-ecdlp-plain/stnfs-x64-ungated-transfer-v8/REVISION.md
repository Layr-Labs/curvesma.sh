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

Gating recovery on a chunk's survivors (our 07742818, then a per-survivor form)
needs the record count of a chunk conditional on a survivor, which is not
measured. This package re-reads every record at 7.003 (no gating), as the
source does.

## R6 (from the recover-Q thread): record walk to 2^64, region 19/50

X2 = 2^63 -> 2^64: primes in [2^63, 2^64) move from ECM to the record walk.
Their residues need double-word arithmetic and `slice12` is measured only
below 2^63, so every record-walk hit, empty coset, apply and setup6 for p in
[2^63, 2^64) is charged at **2x** (WIDE_REC = WIDE_SETUP = 2). The extra yield
buys theta 207/500 -> 19/50 (1040 canonical rows, 965 primitive, a-radius
>= 892, count products 0.971857579 / 1.027240120) and slacks (300, 228) ->
(291, 225) (best of a 1,080-point grid over theta and slacks). Relations 2^64.1311, margin 0.3013,
mean norms 746.0/499.0. Per cell 428.9 over 2^126.692 cells; work SE 3.354%,
allowance 1.0671; relation collection 2^121.0554 x 1.0671 = 2^121.14914 <=
2^121.1492. Descent at the new norms: two-log
2^117.91615, Markov 0.60758, boot 1.8144e-04, E[N] 2^17.60, node counter
2^-9.69, success >= 0.390235.

| core ledger row | R1 + R2 | R1-R3 | R1-R4 | R1-R4 + R6 |
|---|---:|---:|---:|---:|
| relation collection | <= 121.3481 | <= 121.3481 | <= 121.2551 | <= 121.1492 |
| linear algebra | <= 119.36845 | <= 118.7428 | <= 118.7428 | <= 118.7428 |
| core ledger | <= 121.8398 | <= 121.7452 | <= 121.6752 | **<= 121.5968** (upper-rounded rows sum 121.596717; unrounded 121.5966599) |

Generator 105.5088, meter 118.635 and fixed rows unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
