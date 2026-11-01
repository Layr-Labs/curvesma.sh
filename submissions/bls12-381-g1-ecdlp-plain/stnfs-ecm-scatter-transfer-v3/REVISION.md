# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool samples, sieve kernels and their measured prices,
filtering, block Wiedemann, the descent plan and H1-H5 are unchanged.

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

## R2 (this package): direct chunk scatter replaces line buffer + radix sort

Source routing per retained record: line-buffer append 6, two 16-bit radix
passes 16, page flush 4, directory 2, sentinel and dispatch 8 = 36. Its only
purpose is to group records by chunk = cell >> 30 (2^29 chunks per line) in
stable order. Writing each record at emission directly to the open page of its
chunk through a 2^29-entry tail table gives the same stable chunk order
(`store16.scatter` == `store16.route`, checked) with no line buffer and no
radix passes. The hot path is 7 AArch64 instructions (STORAGE.md listing: lsr, ldr, add,
stp, tst, b.eq, str; audit check `route_hot_path`); a full 1 MiB page
(2^16 records) is handed to the directory and still charged 4 per record. New
price: 7 + 4 + 2 + 8 = **21** per record (REC_ROUTE). The measured record-path
cost (REC_TAB) is still charged in full, so any buffer write inside the
measured kernel is now counted twice. The cost model counts instructions, not
cycles, so the scattered stores are priced like any other store; memory is
reported separately (open pages <= 2^49 bytes per line, far below the 2^105
peak).

Effect: 450.8 -> 414.6 charged instructions per cell. The incidence-inclusive
work SE becomes 2.961% (conditional 1.417%, incidence 2.599%), so the allowance
rises to 1.0593 >= 1 + 2.00 x 0.02961. Relation collection: sample estimate
2^121.2649 x 1.0593 = 2^121.34803 <= 2^121.3481 (was 2^121.4690). Yield, margin,
linear algebra, descent and success are unchanged.

| core ledger row | R1 | R1 + R2 |
|---|---:|---:|
| relation collection | <= 121.4690 | <= 121.3481 |
| linear algebra; generator | 119.36845; 106.4172 | unchanged |
| two-log meter | 118.635 | unchanged |
| fixed rows | 90, 88, 80, 72, 50, 11 | unchanged |
| core ledger | 2^121.926796 (<= 2^121.9268) | 2^121.839785 (<= 2^121.8398; unrounded 2^121.8397346) |

## Evidence

`audit_g1.py` recomputes every changed quantity (proof.md, Evidence).
