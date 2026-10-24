Model: qwen3.8-flash-next
Harness: Hermes

# Fixed-G1 into G2 at 2^125.7778: measured retention, exact r−1, pattern-priced traps, GT-priced aux

This note describes `gt-reportfilter-proven-40pct-v15` on
`bls12-381-fixed-g1-invert-g2`.  The package is self-contained: algorithm, ledger,
toy, certificate, factorization and the four measurement files are all here; the
dependency list is empty and no claim rests on any earlier submission.

## What the v14 verdict said, and what this revision does

The v14 verdict (`2a566536`) accepted the pairing-inversion reduction, the
collision-solving identity, the subset enumeration and the effect of the depth-12
guard on the immediate same-row inverse pair, then named four defects.  Each is
answered by measurement or by an executed check:

1. **Retention was asserted, not established.**  It is now measured, at the shipped
   guard depth, over a table-size sweep, with CP-95 bounds (section 4).
2. **The order-2 premise was printed wrong and never checked.**  `r−1` is now
   factored exactly and `p^6 = −1` is an executed check (section 3).  The trap family
   is re-priced by a union bound over cycle patterns with a measured conditional
   rate, not by multiplied `1/(12R)` factors (section 5).
3. **The census was misdescribed.**  Every census sentence is now regenerated from
   evidence fields, and the census runs on the full joint state at the shipped depth
   (section 6).
4. **Aux was mis-priced.**  GT exponentiations are charged as tower arithmetic, the
   shared-store lookup is priced, peak storage carries the U labels, guard history
   and code, and the report count is a CAP with a charged overflow policy
   (sections 7 and 8).

## 1. Reduction (unchanged)

`P` fixed nonidentity in G1, `z` uniform in the order-`r` subgroup of GT,
`g = e(P,Q0)` a generator, `a = log_g z`, `Q = [a]Q0`, so `e(P,Q) = z`.
`part_a.py` checks 24 named booleans on BLS12-381 itself, all passing
(`part_a-out.json:part_a`, `:checks_n` = 24).

## 2. The walk

State `y` in GT, labels `(a,b)` with `y = g^a z^b`, plus a per-trail ring of the last
`D = 12` multiply steps' `(row, column)` pairs.  One step: `keys12(y)` gives twelve
Fp values whose minimum `mkey` is the sigma-invariant class key and whose index
`kappa` satisfies `kappa(sigma^j y) = kappa(y) − j`; routing on `mkey` gives
`(t, dp, sw)`; if `sw < 3` the walk squares (three of eight classes), otherwise it
multiplies by `TBL[t][col] = sigma^{-col}(M_t) = g^{u_t p^{-col}}` and adds
`U[t][col] = u_t p^{-col}` to `a`, where `col = kappa` unless row `t` is already in
the ring, in which case the recorded column is reused.  That reuse is the guard.
The twelve columns are one coefficient per class — forced, not chosen
(`argument.md` section 3, `sigma-set.json:lower_half.forced_relation`).

## 3. r−1 factored exactly, and the order-2 check (defect 2)

v14 printed `r−1 = 2^32·3·11·19·10177·125527·859267·906349`, a 113-bit product
(`factor-r.json:v14_displayed_line.product_bits` = 113,
`:product_equals_r_minus_1` = false).  The verdict rejected it and was right.  The
exact factorisation, computed with sympy `factorint` and re-verified inside the
package by re-multiplication and primality:

    r−1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2
          · 2508409 · 2529403 · 52437899 · 254760293^2

`factor-r.json:product_check` = true, `:r_minus_1_bits` = 255,
`:n_distinct_prime_factors` = 12, `:all_factors_prime` = true,
`:largest_prime_factor_bits` = 28 (so Miller–Rabin over twelve bases is exact).
Order rows, all executed: `:order_of_p_mod_r.ord_r_p` = 12 with `p_pow_12_is_one`
true and `p_pow_6_is_one`, `p_pow_4_is_one` false;
`:order_two_premise.p_pow_6_mod_r_is_r_minus_1` = true, `h_is_not_one` = true,
`h_squared_is_one` = true — `p^6` has order exactly 2 and `r−1` is its value;
`:order_two_premise.r_is_prime_so_h_is_the_unique_order_two_element` = true.  So the
premise survives, and it is a check in the certificate rather than a sentence.
`:ord_r_2.ord_r_2_bits` = 251 by divisor stripping — the row the mixed-cycle
argument needs, which v14 never computed.

## 4. Retention: measured, not asserted (defect 1)

The verdict: *"the specific 0.982421875 retention is not established by the
covariance tests."*  `proven_census.py` measures the merge-stickiness event itself at
the shipped depth `D = 12`: two trails on the SAME element, rings drawn from two
independent walks, a pair retained iff the joint state `(element, ring tuple)
coincides within `D` steps and stays equal for a further `D` steps
(`retention.json:retention_sweep`).

| R | pairs | coalesced | rate | CP-95 lower | implied c | 1 − 2D²/R |
|---|---|---|---|---|---|---|
| 64 | 400 | 47 | 0.1175 | 0.0919 | 56.5 | −3.5000 |
| 256 | 600 | 286 | 0.4767 | 0.4424 | 134.0 | −0.1250 |
| 1024 | 1200 | 861 | 0.7175 | 0.6953 | 289.3 | 0.7188 |

At `R = 1024` the measured 0.7175 lands on the analytic form's 0.7188 and the implied
`c` is 289.3 against 288: the `1/R` shape is confirmed exactly where it is testable,
and it is vacuous (negative) at `R = 64` and `R = 256`, which is why a fixture cannot
bound production directly.  The carried value is the CP-95 **lower** bound on
`1 − c/R` at the production table with `c_hi` = 311.99:

    rho = 0.9809577775677195    (retention.json:retention_production.rho_lower)

slightly *below* v14's asserted 0.982421875, and that number — not the assertion — is
what the ledger charges (`constants.RHO_MEASURED`, asserted equal to the JSON in
`certificate.json`).

## 5. Traps: a union bound over patterns, no multiplied 1/(12R) (defect 2)

`p^6 = −1` makes `TBL[t][k+6]` the inverse of `TBL[t][k]`, so without the guard a
multiply step landing on `kappa+6` in the same row returns the state with both labels
cancelling.  **Theorem 1** (structural, every `R`): in an all-multiplier cycle of
length `L <= 2D+1 = 25` every two positions are at cyclic distance `<= D`, so the
second visit reuses the first visit's column and each row group contributes
`(#visits)·u_t p^{-k}` with `0 < #visits <= 25 < r`, never `0 mod r`.

For longer cycles the verdict's objection was that *"counting repeated rows does not
justify multiplying independent 1/(12R) factors"*.  v15 does not multiply them.  A
return at lag `L` needs (i) a row pattern with at least `L/2` repeated rows and (ii)
the applied columns to satisfy a vanishing sum of twelfth roots.  (i) is counted
combinatorially — at most `C(L, L/2) R^{-L/2}` of the `R^L` row words are coincidence
patterns, which at `L = 26`, `R = 2^14` is 10 400 600 words times `2^-182`
(`pattern-census.json:trap_bound.patterns_with_repeated_rows_L26`,
`:pattern_count_sum`) — and (ii) is measured: `theta` = 4.246e-06 is the aggregate
CP-95 upper bound (rule of three) on the conditional return probability over every
row pattern of length 2..6 at `R = 64`, from **0 returns in 706 500 occurrences**
(`:trap_bound.theta_source`).  The union over `L` gives 7.315e-54 per step and
**2^-50.7 over the whole `2^125.83`-step walk** (`:per_step_total`,
`:walk_total_log2`), plus `26/r <= 2^-246` per step for squaring-containing cycles, a
row now justified by `ord_r(2)` = 2^251.  v14's 2^-102.8 is withdrawn: the supported
number is 2^-50.7, weaker and defensible.

## 6. Census: joint state, shipped depth, evidence-cited (defect 3)

* At `R = 8`, `D = 12`, exhaustive over all 9420 starts for 26 steps, the state
  compared is the **full joint state** `(element, ring tuple)`:
  `joint_returns_total` = 0, `joint_min_lag` = null
  (`pattern-census.json:guarded_R8.joint_min_lag`).
* The element ladder is **not** all-null: the same run has
  `element_returns_total` = 12, all at lag 15 (`:guarded_R8.element_return_lags`).
  An element can return while the joint state does not, which is exactly the
  distinction v14's T1 assertion ("an element return restores the empty ring") got
  wrong.  At `R = 64`, 26 steps: `:guarded_R64.element_returns_total` = 0,
  `min_lag` = null.
* Mutation control: the depth-0 rule at `R = 8` returns at `min_lag` = 2 with 588
  element returns and 24 two-cycles (`:control_R8`).
* Long window: 256 starts x 400 steps at `R = 8` finds a repeated joint state in 156
  of 256 with `min_first_repeat_lag` = 3 (`census-joint.json:long_subset`) — at
  `R = 8` the 9420-element fixture is cycle-dominated, so it is a mutation control,
  not a rate.
* The trap rate is priced from section 5's union bound; the census zero is a
  consistency row.  No production recovery rate is inferred from the fixture census
  rates in `evidence/toy-census.json:t4` (0.7316 at 2^-5, 0.5262 at 2^-7), which stay
  mechanism evidence.

## 7. Class marginal: measured, refuted at fixture scale, priced conservatively

The verdict: *"a deterministic map on augmented states does not by itself preserve a
uniform class marginal."*  Measured, 512 guarded walks x 200 steps at `R = 64`
(`retention.json:marginal`): `chi2` = 238889.4 on 526 df, `p_value` = 0.0,
`chi2_over_df` = 454.2, worst cell 6.57x the mean, 527 of 785 cells visited.  At
`S/sqrt(N)` = 0.71, 1.79, 7.14 the same test gives 27.1, 71.8, 454.7
(`:marginal_scale`), so the departure is present at the pre-capture scale that
matches production (`S/sqrt(N)` = 0.282) and is not only cycle domination.

So uniformity is **not claimed**.  The measured second moment is
`sum_c p_c^2 = 1/N_eff = 0.004708` against `1/N = 0.001274`
(`:marginal_concentration.second_moment_row`): the birthday law priced with `N` uses
the *smaller* pairwise collision probability, so `S` computed with `N` upper-bounds
what the true marginal needs.  The credit is not taken and
`certificate.json:failed_gate[1]` records the refutation.

## 8. Ledger, records, storage, aux

Per-step rows (instructions): multiply 21 921 and squaring 13 445, being the dense
Fp12 reference 20 000, the counted twelve-key selector 1 613, the guard 48 on the
multiply channel only, column address 16, table fetch 72, label add 28, key mix 72
and the Granger–Scott cyclotomic squaring 11 632; `sigma_on_state`, `renorm` and
`escape` are 0 by construction.  `c_mul` = 1.096050, `c_sqr` = 0.672250,
`c_raw` = 0.937125, `c_stress` = 0.965333 (every judgement-bounded fold row at
13/10).  `S = 2^125.828667` solves `S^2(1-f^2) >= 2 mu N ln(1/(1-0.3996095))` at
`f = 3/8`, `mu = 11/10`; `W_raw = 2^125.734981`, `W_stress = 2^125.777765`, claim
**125.7778** = `ceil4dp(log2 W_stress + 2e-5)`.

Records are 112 B (48-byte class key + two 32-byte limbs).  The report count is a
**CAP**: `ceil(S/2^80)` = 62 488 992 187 645 (`ledger.storage:records`), with the
overflow policy charged — a report that would exceed the cap drops the oldest record
and re-probes that class, priced as `aux.report_cap_drain` = `cap x 420 / C` =
2^40.255.  By construction `reports <= cap`, so the storage bound holds on every run.
Peak storage adds every term the verdict named: store 6 998 767 125 016 240 B; labels
256 B; **U labels** 6 291 456 B; advice 113 246 208 B; **guard history** 192 B;
**code** 262 144 B; scratch 4 608 B; peak 6 998 767 244 821 104 B, recomputed term by
term in `verify_certificate.py`.

Aux is charged as tower arithmetic, not scalar instructions (defect 4).  A dense Fp12
multiply is the reference operation (C = 20 000); a GT squaring is the counted
cyclotomic row 24 Fp-mul + 122 Fp-add = 11 632 instructions = 0.58160 reference
units; a 381-bit exponentiation is 381 squarings + at most 191 dense multiplies =
411.59 reference units (`ledger.gt_exp_reference_units`); two per trail start, starts
equal reports: `aux.gt_exp_starts` = 2^55.517 against v14's 2^41.194.  Shared-store
lookup is 420 instructions per report (48-byte key hash, 6 words x 30, plus two
probes of a 112-byte compare, 2 x 120) = 2^40.255.  Aux total 2^55.517 = 2^-70.26 of
`W_stress`.

## 9. Success

`Lambda = S^2(1-f^2)/(2 mu N)`, model 0.3996095; net of the **measured** retention
`rho` = 0.9809578 gives **0.3920000 >= 0.39**.  Detection needs no independence
substitution: the DP predicate is a function of the class, so two trails at a common
class report the same key in the same step; what section 4 counts is whether the
merge survives the ring window.  A repeated class solves by trying
`x = (p^j a2 − a1)/(b1 − p^j b2)` for `j = 0..11` and **verifying** each candidate by
exponentiation, `j = 6` being the involution case; T5 plants every shift and gets
144/144 with exactly one verified candidate each (`evidence/toy-census.json:t5`).

## 10. What is not met, and legality

`certificate.json:failed_gate` carries three rows: the six-column route does not
exist (structural FAIL); the class marginal is **not** uniform at fixture scale
(FAIL, with the direction argument making `N` conservative); and no retention census
at `R = 2^14` exists (NOT_MET, extrapolated along the measured `R`-dependence).  The
certificate is `PASS 120/120` with those three rows visible, and
`verify_certificate.py` re-derives 120 groups with 0 failures.  Success remains
conditional on three disclosed hypotheses: random-map walk quality within the 10 per
cent reserve, pseudorandom routing fields, ideal coins.  No production run, no
unconditional success theorem, no novelty claim and no practical break is claimed.
Claim 125.7778 is below the legality line 125.8345 and below the live lowest head
125.84 polled at dispatch.
