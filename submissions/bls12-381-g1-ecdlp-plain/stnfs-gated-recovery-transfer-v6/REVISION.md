# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels and their measured prices,
filtering, the descent plan and H1-H5 are unchanged.

## R1 (from 27202622): stronger main-phase ECM buys a smaller region

ECM_FIND 0.99 -> 0.999 for primes in (2^63, 2^69] (403 instead of 269 curves
per side; per-survivor 2^38.51 -> 2^39.09, about 5 instructions per cell) gains
about 0.066 bits of yield (factor ECM_FIND^(lambda_f + lambda_g)), spent on
theta = 87/200 -> 207/500 (1378 rows, 1285 primitive, a-radius >= 931, count
products 0.973023183 / 1.026085506) and slacks (298, 229) -> (300, 228).
Relations 2^64.1311, margin 0.3012; mean norms 747.5/498.9; two-log
2^117.9122 (unchanged plan), meter 2^118.635, success >= 0.391949.

## R2 (from recover-P 43d603a4/be446c45): direct chunk scatter

The source routed each record through a line buffer and two 16-bit radix passes
(6 + 16) to group records by chunk = cell >> 30 in stable order, then flush 4,
directory 2, dispatch 8 = 36. Writing each record at emission to its chunk's open
page through a 2^29-entry tail table gives the same stable order
(`store16.scatter` == `store16.route`) with a 7-instruction hot path (STORAGE.md;
`route_hot_path`): 7 + 4 + 2 + 8 = 21, per cell 450.8 -> 414.6, relation
collection 2^121.3481. R4 replaces the remaining 4 + 2 + 8.

## R3 (this package): unit entries and a combined row tail in the matrix product

The source charges 88 per sparse entry: 33 (`mac_run` byte multiply-accumulate)
+ block selection 4 + sign 3 + transfers 32 + finishing and headers 16.

- **Unit entries (this package).** Rows are stored per (block, sign) group as
  a unit list and a general list; `_group_run` (clang 18, `la_mac.s`) adds a +-1
  entry in **11** instructions and keeps the 33-instruction loop for the rest.
  The unit share of filtered row weight was measured here with a copy of
  `filtsim.merge` that tracks each row's relation multipliers: 0.897-0.923 at
  2^12, 0.933-0.939 at 2^14-2^15 and 0.946-0.947 at 2^17 (rising with scale;
  `audit_g1.py {"unit_entries": [k, seed]}`); charged at 0.88 (H4u).
- **Row tail (recover-Q 06aac4a1).** Partial sums are combined up the 20-owner
  tree: send 22, combine 51, finish 119 (Barrett), **2,090 per row** (was 48
  per entry); recompiled here, `_finish` tested on 20,000 exact pairs.

Entry price 4 + 3 + 2090/101 + 0.88 x 11 + 0.12 x 33 + 168/101 = 43.00, charged
**44**; block Wiedemann re-minimised at (m, n) = (18, 4). Linear algebra
2^119.36845 -> **2^118.74271**, generator 2^106.4172 -> 2^105.5088. N, filtering,
the descent and success are unchanged.

## R4 (from the recover-P thread): measured end-to-end record route

R2 kept flush 4 + directory 2 + dispatch 8 per record. Measured route: scatter 7
+ owner page copy 13 per 4 records (STORAGE.md listing) + page terms < 64 per
2^16-record page = 10.26, charged **10.3** (recover-P's executed count, gcc 13
under qemu-aarch64: 10.0009). Per cell 414.6 -> 388.7, collection 2^121.2551.

## R5 (from the recover-Q thread): survivor-gated recovery

Recovery (7.003 per record) re-reads a chunk's records to rebuild survivors'
factorizations, but a 2^30-cell chunk needs it only if its scan found a
survivor, with probability <= 2^30 s_f s_g (the model's per-side survivor upper
bounds; worst row 2^-0.03). Charge 7.003 x min(1, 2^30 s_f s_g) per record + 16
per chunk + 8 per released 2^16-record page. Per cell 388.7 -> 372.2; work SE
3.000%, allowance 1.0600; relation collection 2^121.1085 x 1.0600 = 2^121.19252
<= 2^121.1926.

| core ledger row | R1 + R2 | R1-R3 | R1-R4 | R1-R5 |
|---|---:|---:|---:|---:|
| relation collection | <= 121.3481 | <= 121.3481 | <= 121.2551 | <= 121.1926 |
| linear algebra | <= 119.36845 | <= 118.7428 | <= 118.7428 | <= 118.7428 |
| core ledger | <= 121.8398 | <= 121.7452 | <= 121.6752 | **<= 121.6287** (unrounded 121.6285982) |

Generator 105.5088 (R3), two-log meter 118.635 and fixed rows are unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
