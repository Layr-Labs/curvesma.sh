# Argument (v23)

1. *Reduction.* `g = e(P0,Q)`, `z = g^x`, walk `y = g^A z^D` in GT. A class collision with
   distinct labels gives one linear equation, `x = (p^k A1 - A2)/(D2 - p^k D1)`; the denominator
   is checked nonzero before use, and `Q` nonidentity is required for the reduction. Success is
   verified by `e(P,Q) = z` before any output.

2. *Protocol is meter-first.* `S = 56628313 * 2^100` steps are purchased. A single work counter
   in `C` units is checked before every step and the walker abandons at the cap. The ledger
   closes at `sum(rows) == cap` exactly, residual charged. Success is
   `P(solve before meter exhaustion)`, measured by running the full protocol (workers, cap,
   meter, abandonment) on fresh fixtures — not inferred. The v22 assertion that the record
   count is fixed by S, and that purchased equals executed work, is withdrawn; this
   package contains no sentence of that form.

3. *Price.* Step = `5/8 x 22433 + 3/8 x 29263 + 4.1554` = `24998.405` instructions
   (`1.2499 C`). The 22433 and 29263 rows are the declared reference cost model
   (`C = 20000` dense Fp12 multiply, `27000` Fp12 squaring) plus the verified multiply path
   (key recompute 1487, booth 405, table 150, exponent add 60, common 331, exponent double 40).
   The `cyclotomic_sqr` instruction span is **excluded** because that kernel is not yet
   verified against bigint; no unverified span prices GT arithmetic.

4. *Guard interaction, measured.* On 10,240,000 real walk steps at `Rg = 2^16`, depth 5, the
   guard converted 6230 multiply steps into cyclotomic squarings: `6.083984375e-4` per step.
   Times the measured `sqr - mult` delta of 6830 instructions that is 4.1554 instructions per
   step, carried as its own ledger row instead of v22's 3/8 nominal average.

5. *Kernels.* `f2mul` repaired (cross term to `c1`, `t0 - t1` to `c0` for `u^2 = -1`) and
   verified with 0 mismatches on 10000 random vectors against a plain-bigint reference, with a
   mutation control that reverts the fix and fails 10000/10000. `cyclotomic_sqr` identity
   repaired to the Montgomery representation of 1 (`R mod p`) loaded from an extern table, and
   the harness modulus (previously shifted by one limb) corrected; the squaring kernel is still
   **not verified** (10000/10000 mismatch) and is disclosed as such.

6. *Traps, counted on real walks.* The abab-family word admits 0 times in 10,240,000 steps at
   the shipped configuration, at lengths 4, 6 and 8; the measured 95% upper bound is
   `2.93e-7` per step. Controls: v20 (`Rg = 4096`, depth 1) admits `1.666e-2` at length 4
   (`2^-5.91`); class-only (`Rg = 12`) admits `9.994e-3` at length 8 (`2^-6.64`). Waste charged
   is `rate_upper_bound x S`, since a trap can only waste purchased steps. No `m!` tail, no
   prefix-extension argument, no independence premise.

7. *Census, un-pooled.* 512 fresh fixtures x 1 trial: 512/512, CP-95 `0.99`. Two-scale check at
   2x fixture: 256/256, CP-95 `0.98837`, rate delta 0.0, CP-95 delta `-0.00163`, walk:record
   shape held (0.1234 vs 0.1227).

8. *Success and claim.* `0.99 x 0.675 = 0.668 >= 0.39`, where 0.675 is the measured
   worst-arm merge survival at gap 40 (curve re-measured with a live guard; the true-merge
   control reads 0.8775, not 1.0, and that decay is disclosed). Work
   `W = 1.7947e42` instructions = `8.974e37 C` units; claim **126.08**, residual charged so the
   ledger closes at the cap.
