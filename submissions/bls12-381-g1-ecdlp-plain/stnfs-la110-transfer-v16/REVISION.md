# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

- **R1**: ECM_FIND 0.99 -> 0.999 above X2 (403 curves) buys a smaller region.
- **R2**: chunk-scatter record route 36 -> 21 per record.
- **R3**: replaced by R10. **R4**: measured route 10.26, charged **10.3**.
- **R5** withdrawn: every record re-read at 7.003.
- **R6**: primes in [2^63, 2^64) to the record walk, charged **2x**.
- **R7**: line sieve from 2^11 (`P_MIN` 2048).
- **R8**: descent relation input at a linear 1-SE yield lower bound + collection guard.
- **R9**: six descent-plan entries re-optimised; jointly **LOG2_RELS 63.97**, theta
  190/500, slacks (297, 233); relations 2^64.1311, margin 0.3012; work SE 3.377%,
  allowance 1.0676; collection <= 2^121.1117; two-log 2^118.27320, meter **2^118.996**
  (Markov 0.60593), success >= 0.390320.

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

## R11 (as in our aa2a2ba5, 124.843970069062): exact price ceilings, 0.01 table shift, reserve 1001/1000

The nine line/record prices become the exact logged differential quotients
(`LINE_TAB`/`REC_TAB` Fractions; audit `exact_prices` re-parses the Ir pairs
from measurements12.txt and checks equality, plus the quoted decimals match to
0.01; setup6 406.87 and norm-bound 1508.2 verified as ceilings). Both tables
are checked strictly decreasing and convex in log2h (exact rational slopes;
audit `price_shape`), and `tab()` charges 0.01 below its argument, covering
two-decimal boundary rounding twofold under the H4 shape; production log2h is
checked in the covered region (line 34.91+, always flat; record 1.57+, above
the last segment; audit `tab_shift_ranges`). Collection with every
quoted-decimal constant at +half-ulp moves the ledger 1.000036x (audit
`residual_bound`), so the reserve drops 101/100 -> **1001/1000**. Yield,
descent, success and all other parameters are untouched. Collection <=
121.1117 -> **<= 121.1127** (shift +0.0010); core ledger <= 121.5198 ->
**<= 121.5206** (unrounded 121.5204978); G1 124.85609 -> **124.843970069062**.

## R12 (this package): line sieve to 2^48.75

The line/record path split moves P_MID 20 -> **48.75** (0.25-grid, so the
record-walk quadrature still hits the 2^63 wide threshold exactly; audit pins
the value and the grid). Primes below 2^48.75 are line-sieved at the exact
h-parameterised price (flat 24.098 up to h = 2^13.82, then the measured
interpolation down to h = 2^1.82); only primes above move by record walk, so
apply/recovery/routing shrink with the record count. Production line log2h
stays in the measured span (min 6.16, audit `tab_shift_ranges`); the line
implementation is prime-agnostic 64-bit arithmetic (KERNELS.md), and no line
prime reaches the 2^63 wide cushion. Yield, descent, success and all other
parameters are untouched (P_MID enters only the cost split). Collection <=
121.1127 -> **<= 120.8043** (337.7 per cell; work SE 3.471%, allowance
1.06943); core ledger <= 121.5206 -> **<= 121.2945** (unrounded 121.2943751);
G1 124.843970069062 -> **124.617870069062**.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`; proof.md,
Evidence); `evidence16.py` [10] recounts the kernels.
