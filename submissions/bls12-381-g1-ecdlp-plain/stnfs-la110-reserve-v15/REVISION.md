# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, sieve kernel prices, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

- **R1** (our 27202622): main-phase ECM_FIND 0.99 -> 0.999 (403 curves) buys a smaller region.
- **R2** (recover-P be446c45): records scattered at emission to their chunk's page (`store16.scatter` == `store16.route`).
- **R3**: replaced by R10. **R5**: withdrawn (every record is re-read at 7.003).
- **R4** (recover-P thread): measured record route 7 + 13/4 + < 0.01 = 10.26, charged **10.3** (STORAGE.md).
- **R6** (recover-Q thread): record walk to 2^64; primes in [2^63, 2^64) charged **2x** on record walk, empty cosets, apply and setup6 (KERNELS.md).
- **R7** (recover-Q thread): line sieve from 2^11; p < 2^11 is the unsieved small group of `threshold_model`.
- **R8**: the descent input holds at a 1-SE lower bound on the yield, rels x (1 - relative pool SE) >= 2^LOG2_RELS, with the fixed-pool estimator over all 20,000 draws (audit `yield_lower_bound_covers_descent`); a guard continues collection until 2^LOG2_RELS relations exist.
- **R9** (recover-Q thread): six descent-plan entries re-optimised for fewer relations (90g, 95g,
  f-side 100, 110, 120, 150); jointly optimised **LOG2_RELS 63.97**, theta 190/500, slacks
  (297, 233): collection <= 2^121.1117, meter **2^118.996**, success >= 0.390320
  (source-core.md sections 3, 6, 7).

## R10 (from our recover-Q 711418d3, AI-valid): fully priced linear algebra, pad 1.10

The linear-algebra accounting of 711418d3, unchanged (STORAGE.md, source-core.md section 4):
every reduction mod r is a compiled loop-free kernel (C source `row_tail.c` and the
listings are sections of `kernels_asm.b85`; audit input {"source": "row_tail.c"} prints any section): the row tail, which reduces the sparse pair and the dense 9-limb
accumulator (fold9 + Barrett reduce6), at C_ROW 2,750 per row and product (counted 2,708);
each lone 9-limb reduction at C_RED9 420 (counted 385: Horner elements, projection entries,
generator multiplications); six-limb sums valid for any row; 21 per sparse entry (4 + 3 +
0.88 x 11 + 0.12 x 33 = 20.64, H4u). The matrix pad drops 1.61 -> **1.10** on the measured N/U
trend (`SCALE_NU`: max over seeds 0.0850, 0.0635, 0.0563, 0.0516, 0.0469, 0.0461, 0.0425 at
2^12..2^18, 31 runs, checked non-increasing; charged 1.10 x 0.0471 = 0.05182, above every value
from 2^15 up; that production N/U stays at or below these is part of H3). N = 2^59.1811,
blocks (19, 4): linear algebra 2^118.7428 -> **2^117.7327**, generator 2^105.5088 ->
2^107.2082 (every multiplication now carries a reduction).

| core ledger row | 3eccc6c9 (R1-R9) | + R10 |
|---|---:|---:|
| relation collection | <= 121.1117 | <= 121.1117 |
| linear algebra | <= 118.7428 | **<= 117.7327** |
| generator | <= 105.5088 | <= 107.2082 |
| two-log meter | 118.996 | 118.996 |
| core ledger | <= 121.6221 | **<= 121.5198** (rows 121.5197816; unrounded 121.5197467) |

## G1 wrapper reserve (proof.md)

The whole-core reserve drops from 101/100 to **1001/1000**: `audit_g1.py` `rounding_ratios` shows every measured
kernel price is within 0.08% of its charged value (exact Ir ratios at most 1.000104; worst printed-rounding case
1.0008). G1: 124.85609 -> **124.84318**.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`; proof.md,
Evidence); `evidence16.py` [10] recounts the kernels.
