# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged. Each entry
names its submission and verdict; superseded text is in those submissions.

- **R1-R9** (AI-valid 3eccc6c9, 124.9584): ECM_FIND 0.999; chunk-scatter route 10.3;
  R5 withdrawn; [2^63, 2^64) walk 2x; line sieve from 2^11; R8 yield guard; R9 plan,
  theta 190/500, slacks (297, 233).
- **R10** (from our recover-Q 711418d3, AI-valid): fully priced LA with loop-free
  kernels (row tail 2,750, red9 420, 21/entry under H4u); pad 1.61 -> 1.10 on the
  measured N/U trend (H3). Ledger <= 121.5198.
- **R11** (aa2a2ba5, failed on packaging size only): exact price ceilings re-parsed
  from measurements12.txt, convex tables, lookups 0.01 below their argument,
  reserve 101/100 -> 1001/1000.
- **R12, R12b** (494e4872, 43f5ea7e: rejected, line-prime recovery unpriced /
  unjustified) and **R13** (16253187, cancelled): superseded by R12c; R13's SLACK
  96 -> 60 kept.
- **R12c** (AI-valid 0f2da92, 124.754970069062): recovery reran the first-pass coset
  enumerator with flag checks; P_MID 45.5 on the 0.25 grid.
- **R14** (v20tdiv, 124.697570069062): trial division of every prime below 2^30 per
  survivor (C_TDIV 2^35.95); the rerun covered [2^30, 2^45.5).
- **R15** (v21, 124.687770069062): the first pass records line hits in [2^30, 2^45.5)
  through the chunk-scatter path; the rerun is deleted (STORAGE.md). Work SE 3.444%;
  collection <= 120.9128 (364.1 per cell).
- **R16** (AI-valid 356b5f33, 124.587070069062): boot + twelve plan entries
  re-searched at slack 60; two-log 2^117.5873, cap 2^118.3128; stated residual
  multipliers corrected and audit-checked.
- **R17** (AI-valid d605c287, 124.548670069062): merge weight 100 -> 300, blocks
  (26, 5): N/U 0.03225 at 2^17 (charged 0.03547), row weight 304, H4u 0.889+;
  LA 2^117.2095, ledger <= 121.2253.

Evidence for R10-R17: `audit_g1.py` (`la_prices`, `matrix_pad`, `exact_prices`,
`price_shape`, `tab_shift_ranges`, `residual_bound`, `residual_ledger`,
`no_rerun`, `tdiv_prime_bound`), `evidence16.py` [10].

## R18 (this package): descent pr integrated over measured norm sub-cells, plan re-searched, 2-SE relation input, modelled success >= 0.50, row tail over the 7 occupied owners

Changes on R17, each recomputed by `audit_g1.py` and `evidence16.py`:

1. **Descent P integrated over measured norm sub-cells** (DESCENT.md): the per-cell
   semismooth probability is the equal-weight sum over 24 cf-quantile sub-cells
   (own cg maxima) and 10 boot deciles instead of one evaluation per sextile at
   its largest norms, same +23/+6 (+4) shifts; a lower bound on the cell average
   under the unchanged norm heuristic; boost capped at 4x per node and 1.5x at
   the boot. On R17's plan alone: two-log 2^117.5873 -> 2^116.5576.
2. **Plan re-searched** at the integrated P with singles, pairs, (R, n, kA)
   triples and boot moves, for the smallest metered cap that holds modelled
   success 0.45 with the boot's per-relation failure <= 0.6732 (the base's
   validated level); the meter is then set for 0.50: two-log -> **2^115.1678**, cap **2^116.1806** (Markov 0.49557), boot
   (215, 368, 14, 4), 44 entries moved; min mu 31.0, max psi 0.536, max phi 0.0362.
3. **Relation input at 2 SEs**: the descent's occurrence model uses LOG2_RELS =
   63.79 <= the retained count times (1 - 2 relative pool SEs) = 2^63.7959
   (was 1 SE, 63.97); `YIELD_SE_K` = 2. Missing-log rates, descent work and
   success are all evaluated at this input.
4. **Row tail over the occupied owners** (STORAGE.md): N < 7 x 2^56, so the
   row tree spans 7 block owners (6 sends and combines, headers 10 x 7): row
   tail 2,708 -> **1,317**, C_ROW 2,750 -> **1,350** (audit `la_occupancy`,
   `la_prices`); block Wiedemann re-scanned over m in 2..64, n in 2..32:
   **(23, 4)**. LA 2^117.2095 -> **2^117.0498**, generator 2^107.0654 ->
   **2^106.9498**.

Success: 1 - 2 x 1.256e-03 - 0.49557 - 2^-9.58 - 2^-10.70 - 2/r >= **0.50**
(model 0.5000103; was 0.390320; the claim stays 0.39, the gap is margin). Collection, yield margin, filtering, prices,
UDUP, kernels and the wrapper are untouched. Ledger <= 121.2253 -> **<= 121.058691**
(unrounded 121.0586517, rows 121.0586896); G1 124.548670069062 -> **124.382061069061269**.
Sub-cell tables, their fresh-lattice reproduction and the relation-weighted
unit-duplicate campaign: DESCENT.md (live checks `subbins`, `udup_relw`).
Lineage: our recover-Q v14 line (sub-cell integration), recover-P R22
(occupied-block recount, f34f63d9), recover-Q 57598408 (unit-duplicate method).
