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

The 711418d3 linear-algebra accounting, unchanged: every reduction mod r is a
compiled loop-free kernel (`row_tail.c` + listings in `kernels_asm.b85`): row
tail at C_ROW 2,750/row (counted 2,708); lone 9-limb reductions at C_RED9 420
(counted 385); 21 per sparse entry (20.64, H4u). Pad 1.61 -> **1.10** on the
measured N/U trend (`SCALE_NU` max 0.0850-0.0425 at 2^12..2^18, 31 runs,
non-increasing; charged 0.05182, above every value from 2^15 up; H3). N =
2^59.1811,
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

The nine line/record prices become the exact logged quotients (audit
`exact_prices` re-parses the Ir pairs; setup6/norm-bound verified ceilings).
Both tables are checked decreasing and convex (exact slopes; `price_shape`),
and `tab()` charges 0.01 below its argument, covering boundary rounding
twofold; production log2h is checked covered (`tab_shift_ranges`). Collection
at +half-ulp on every quoted decimal moves the ledger 1.000036x
(`residual_bound`), so the reserve drops 101/100 -> **1001/1000**. Nothing else
changes. Collection <= 121.1117 -> **<= 121.1127**; ledger <= 121.5198 ->
**<= 121.5206**; G1 124.85609 -> **124.843970069062**.

## R12 (submitted as 494e4872, rejected: unpriced line-factor recovery)

The line/record split moved P_MID 20 -> 48.75, but line-sieved primes in
[2^20, 2^48.75) were neither trial-divided nor recorded, leaving their
cofactor removal unpriced. R12b below completes the revision.

## R12b (this package): line sieve to 2^47.5 with priced line-hit recovery

The split moves P_MID 20 -> **47.5** (0.25-grid, so the record-walk quadrature
still hits the 2^63 wide threshold exactly; audit pins the value and the grid).
Primes below 2^47.5 are line-sieved at the exact h-parameterised price (flat
24.098 up to h = 2^13.82, then the measured interpolation down to h = 2^1.82);
only primes above move by record walk, so apply/recovery/routing shrink with
the record count. Production line log2h stays in the measured span (min 7.41,
audit `tab_shift_ranges`); the line implementation is prime-agnostic 64-bit
arithmetic (KERNELS.md), and no line prime reaches the 2^63 wide cushion. A
second line pass recovers line hits on survivors (STORAGE.md): 6 per multiple
(5 counted), measured setups, <= 2^43 primes (Rosser-Schoenfeld), 1 flag store
per cell; a brute-force toy verifies completeness (audit `line_recov_loop`,
`line_recov_nlp`, `line_recov_toy`). Yield, descent, success and all other
parameters are untouched. Collection <= 121.1127 -> **<= 120.8569** (350.3 per
cell; work SE 3.452%, allowance 1.06905); core ledger <= 121.5206 ->
**<= 121.3321** (unrounded 121.3320234); G1 124.843970069062 ->
**124.655470069062**.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`; proof.md,
Evidence); `evidence16.py` [10] recounts the kernels.
