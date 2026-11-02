# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, sieve kernel prices, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

- **R1** (our 27202622): ECM_FIND 0.99 -> 0.999 above X2 (403 curves) buys a smaller region.
- **R2** (recover-P 43d603a4/be446c45): records written at emission to their chunk's page
  through a 2^29-entry tail table (`store16.scatter` == `store16.route`): 36 -> 21 per record.
- **R3**: replaced by R10.
- **R4** (recover-P thread): measured route 7 + 13/4 + < 64/2^16 = 10.26 per record, charged
  **10.3** (recover-P measured 10.0009 under qemu-aarch64).
- **R5** withdrawn: every record is re-read at 7.003, as in the source.
- **R6** (recover-Q thread): primes in [2^63, 2^64) move from ECM to the record walk; measured
  prices stop below 2^63, so their record-walk hits, empty cosets, apply and setup6 are charged
  **2x** (WIDE_REC, WIDE_SETUP; KERNELS.md).
- **R7** (recover-Q thread): line sieve from 2^11 (`P_MIN` 2048, `LX0` 11, config `prime_min`;
  `threshold_model` counts p < prime_min against the slacks).
- **R8**: the descent's missing-log input must hold at a linear 1-SE lower bound on the yield,
  rels x (1 - relative pool SE) >= 2^LOG2_RELS (`YIELD_SE_K` = 1, `DESCENT_RELS_MIN`; audit
  `yield_lower_bound_covers_descent`), with the fixed-pool estimator over all N = 20,000 draws,
  relative variance (Q/S^2 - 1/N) N/(N - 1) (`N_POOL`). A guard continues collection until
  2^LOG2_RELS relations exist. 2-SE versions: c53a007e and 30d365cb (both AI-valid).
- **R9** (recover-Q thread): six descent-plan entries re-optimised for fewer relations (old iA
  units): 90g (75,42) -> (75,43), 95g (80,24) -> (80,25), f-side n 100: 83 -> 84, 110: 85 -> 86,
  120: 94 -> 95, 150: 110 -> 109; boot unchanged. Jointly optimised: **LOG2_RELS 63.97**, theta
  190/500, slacks (297, 233); relations 2^64.1311, pool SE 0.1495 bits, lower bound 2^63.9733,
  margin 0.3012; 417.9 per cell over 2^126.692 cells; work SE 3.377%, allowance 1.0676;
  collection <= 2^121.1117; descent at 745.9/498.9: two-log 2^118.27320, boot 7.746e-04, E[N]
  2^18.0001, counter 2^-9.29; meter **2^118.996** (Markov 0.60593), success >= 0.390320.

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

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`; proof.md,
Evidence); `evidence16.py` [10] recounts the kernels.
