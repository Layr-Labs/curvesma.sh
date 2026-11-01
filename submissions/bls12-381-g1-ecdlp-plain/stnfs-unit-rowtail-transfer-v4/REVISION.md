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

Source routing per retained record: line-buffer append 6, two 16-bit radix
passes 16, page flush 4, directory 2, sentinel and dispatch 8 = 36, only to
group records by chunk = cell >> 30 in stable order. Writing each record at
emission to the open page of its chunk through a 2^29-entry tail table gives the
same stable order (`store16.scatter` == `store16.route`, checked) with no line
buffer or radix pass. Hot path: 7 AArch64 instructions (STORAGE.md listing;
audit check `route_hot_path`); a full 1 MiB page still pays 4 per record. New
price 7 + 4 + 2 + 8 = **21** (REC_ROUTE). The measured record-path cost is
still charged in full. Effect: 450.8 -> 414.6 charged instructions per cell;
work SE 2.961%, allowance 1.0593; relation collection 2^121.2649 x 1.0593 =
2^121.34803 <= 2^121.3481.

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

| core ledger row | R1 | R1 + R2 | R1 + R2 + R3 |
|---|---:|---:|---:|
| relation collection | <= 121.4690 | <= 121.3481 | <= 121.3481 |
| linear algebra | <= 119.36845 | <= 119.36845 | <= 118.7428 |
| generator | 106.4172 | 106.4172 | 105.5088 |
| two-log meter | 118.635 | 118.635 | 118.635 |
| fixed rows | 90, 88, 80, 72, 50, 11 | same | same |
| core ledger | <= 121.9268 | <= 121.8398 | **<= 121.7452** (unrounded 121.7451299) |

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
