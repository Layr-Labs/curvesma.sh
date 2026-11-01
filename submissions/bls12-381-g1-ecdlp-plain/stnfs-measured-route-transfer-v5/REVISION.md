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
2^117.9122 (unchanged plan), meter 2^118.635, Markov 0.60589, boot 1.7864e-04,
success >= 0.391949; 2^12 filtering rerun N/U 0.0715 < charged 0.0758.

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
  tree, one 10-word message per owner per row: send 22, combine 51, finish 119
  (Barrett), **2,090 per row** in all (was 48 per entry). Recompiled here from
  the same C, with `_finish` tested against exact integers (20,000 pairs).

Entry price 4 + 3 + 2090/101 + 0.88 x 11 + 0.12 x 33 + 168/101 = 43.00, charged
**44**; block Wiedemann re-minimised at (m, n) = (18, 4). Linear algebra
2^119.36845 -> **2^118.74271**, generator 2^106.4172 -> 2^105.5088. N, filtering,
the descent and success are unchanged.

## R4 (from the recover-P thread): measured end-to-end record route

The R2 price kept the source's flush 4, directory 2 and dispatch 8 per record.
The whole route was then measured: scatter at emission 7, the owner's plain
copy of each page into its receive area 13 instructions per 4 records
(STORAGE.md listing), and page hand-off plus directory walk below 64
instructions per 2^16-record page. That is 10.26 per record, charged **10.3**
(was 21); recover-P's executed-instruction count (gcc 13, qemu-aarch64)
gives 10.0009. Per cell 414.6 -> 388.7; work SE 2.982%, allowance 1.0597;
relation collection 2^121.1714 x 1.0597 = 2^121.25509 <= 2^121.2551.

| core ledger row | R1 | R1 + R2 | R1-R3 | R1-R4 |
|---|---:|---:|---:|---:|
| relation collection | <= 121.4690 | <= 121.3481 | <= 121.3481 | <= 121.2551 |
| linear algebra | <= 119.36845 | <= 119.36845 | <= 118.7428 | <= 118.7428 |
| generator | 106.4172 | 106.4172 | 105.5088 | 105.5088 |
| core ledger | <= 121.9268 | <= 121.8398 | <= 121.7452 | **<= 121.6752** (unrounded 121.6751060) |

Two-log meter 118.635 and fixed rows 90, 88, 80, 72, 50, 11 are unchanged.

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
