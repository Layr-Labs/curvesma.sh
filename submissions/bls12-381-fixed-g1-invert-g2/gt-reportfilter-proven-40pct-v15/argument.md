# argument.md — derivation, itemisation and scope for `gt-reportfilter-proven-40pct-v15`

Every number in this file is emitted by `ledger.py`, `gen_cert.py`, `factor_r.py`,
`proven_census.py` or `sigma_set.py`, and `certificate.json` recomputes it from the
raw evidence JSONs.  Each claim carries a `file:field` citation.

The v14 verdict (submission `2a566536`) accepted the reduction, the collision
identity, the subset enumeration and the depth-12 guard's effect on the immediate
same-row inverse pair, and named four defects.  Sections 2, 4a, 4b, 6, 7 and 8 answer
them; section 9 records what is still not met.

## 1. Reduction (unchanged; re-checked on BLS12-381)

`P` fixed nonidentity in G1, `z` uniform in the order-`r` subgroup of GT,
`g = e(P,Q0)`, `a = log_g z`, `Q = [a]Q0`, so `e(P,Q) = z`.  `part_a.py` checks 24
named booleans on the real curve (`part_a-out.json:part_a`, `:checks_n` = 24),
including `pairing_nondegenerate`, `planted_ePQ_equals_z`, `planted_wrong_Q_fails`,
`ord_r_p_is_twelve` and `p_pow6_is_minus_one`.

## 2. r−1, its factorisation, and the order-2 premise (v14 DEFECT 2)

The verdict: *"the displayed factorization in argument.md is far too small to equal
the 255-bit r−1, and no corresponding order check appears in the supplied tests."*
Both halves are repaired in `factor_r.py`, output `evidence/factor-r.json`.

    r−1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
          * 2508409 * 2529403 * 52437899 * 254760293^2

* `factor-r.json:product_check` = true, `:r_minus_1_bits` = 255 — the product is
  re-multiplied and compared to `r-1` by `assert` inside the script.
* `:all_factors_prime` = true, `:largest_prime_factor_bits` = 28, so Miller–Rabin
  over the first twelve prime bases is exact on every factor.
* `:v14_displayed_line.product_bits` = 113 and
  `:v14_displayed_line.product_equals_r_minus_1` = false: the number v14 printed is
  refuted by re-multiplication, exactly as the verdict said.
* `:order_of_p_mod_r.ord_r_p` = 12, with `p_pow_12_is_one` true while
  `p_pow_6_is_one` and `p_pow_4_is_one` are false: the order is established by the
  divisor test, not by a name.
* **The order-2 premise survives and is now checked**:
  `:order_two_premise.p_pow_6_mod_r_is_r_minus_1` = true, `h_is_not_one` = true,
  `h_squared_is_one` = true, so `p^6` has order exactly 2, and `r` is prime so `r−1`
  is the unique order-2 element of `F_r^*`
  (`:order_two_premise.r_is_prime_so_h_is_the_unique_order_two_element`).
* `:ord_r_2.ord_r_2_bits` = 251 by divisor stripping (a linear scan cannot
  terminate).  This is the row the mixed-cycle argument needs; v14 never computed it.

## 3. The multiplier family the quotient forces, and the route that does not exist

`mu(v p^j) = p^j mu(v)` has the single-term solution `mu(v) = c(C) p^{-kappa(v)}`, one
coefficient per class, twelve columns per row.  `sigma_set.py` scans all 2^12 subsets
of `Z/12`: the `k -> k+1`-closed ones are exactly `{}` and `Z/12`
(`sigma-set.json:enumeration.closed`); `Z/12` contains the inverse pairs `{k,k+6}`;
the lower half is inversion-free but not closed, forcing `m^2 = 1`
(`:lower_half.forced_relation`).  So no nonempty sigma-stable inversion-free column
table exists, the guard is shipped by necessity, and
`certificate.json:failed_gate[0]` records that as a FAIL row rather than a choice.

## 4a. The guard, the short cycles, and the pricing that replaces 1/(12R)

**Theorem 1 (structural, every R).**  In an absorbing all-multiplier cycle of length
`L <= 2D+1 = 25` every two positions are at cyclic distance `<= D`, so the second
visit reuses the first visit's column; each row group then contributes
`(#visits) * u_t p^{-k}` with `0 < #visits <= 25 < r`, never `0 mod r`.  No absorbing
all-multiplier cycle of length `<= 25` exists.

**Longer cycles: a union bound over PATTERNS, not a product of independence
factors.**  A return at lag `L` needs (i) a row pattern with at least `L/2` repeated
rows and (ii) the applied columns to satisfy a vanishing sum of twelfth roots.
(i) is *counted*: at most `C(L, L/2) R^{-L/2}` of the `R^L` row words are coincidence
patterns — `L = 26`, `R = 2^14`: 10 400 600 words times `2^-182`
(`pattern-census.json:trap_bound.patterns_with_repeated_rows_L26`,
`:pattern_count_sum`).  (ii) is *measured*: `theta` = 4.246e-06, the aggregate CP-95
upper bound (rule of three, `3/N`) on the conditional return probability over every
row pattern of length 2..6 at `R = 64`, from **0 returns in 706 500 occurrences**
(`:trap_bound.theta_source`).  The union over `L` gives 7.315e-54 per step and
**2^-50.7 over the whole `2^125.83`-step walk** (`:trap_bound.per_step_total`,
`:walk_total_log2`), plus a mixed-cycle row `26/r <= 2^-246` per step justified by
`ord_r(2)` = 2^251 (`factor-r.json:ord_r_2`).

No `1/(12R)` factor is multiplied anywhere, no independence of sampled table entries
is assumed, and no cancellation across rows is claimed away.  v14's 2^-102.8 claim is
withdrawn and replaced by the weaker, supported 2^-50.7.

## 4b. The two things the theory cannot supply — measured (v14 DEFECT 1)

**(a) Class marginal.**  The verdict: *"a deterministic map on augmented states does
not by itself preserve a uniform class marginal."*  `proven_census.py` measures it:
512 guarded walks x 200 steps at `R = 64`, class key `min{keys12(y)}`
(`retention.json:marginal`).  The result **refutes uniformity at fixture scale**:
`chi2` = 238889.4 on 526 df, `p_value` = 0.0, `chi2_over_df` = 454.2, worst cell
6.57x the mean, 527 of 785 cells visited.  At `S/sqrt(N)` = 0.71, 1.79 and 7.14 the
same test gives `chi2_over_df` = 27.1, 71.8 and 454.7 (`:marginal_scale`), so the
departure is present at the pre-capture scale matching production
(`S/sqrt(N)` = 0.282) and is not only cycle domination.

Uniformity is therefore **not** claimed.  What is used is the direction of the error:
the measured second moment is `sum_c p_c^2 = 1/N_eff = 0.004708` against
`1/N = 0.001274` (`:marginal_concentration.second_moment_row`), so the birthday law
priced with `N` uses the *smaller* pairwise collision probability and `S` computed
with `N` upper-bounds what the true marginal needs.  The credit is not taken.
`certificate.json:failed_gate[1]` records the refutation rather than hiding it.

**(b) Retention, measured directly.**  The verdict: *"the specific 0.982421875
retention is not established by the covariance tests."*  v15 measures the
merge-stickiness event at the shipped depth `D = 12`: two trails on the SAME element
with rings drawn from two independent walks, retained iff the joint state
`(element, ring tuple)` coincides within `D` steps and stays equal for a further `D`
steps (`retention.json:retention_sweep`).

| R | pairs | coalesced | rate | CP-95 lower | implied c | 1 − 2D²/R |
|---|---|---|---|---|---|---|
| 64 | 400 | 47 | 0.1175 | 0.0919 | 56.5 | −3.5000 |
| 256 | 600 | 286 | 0.4767 | 0.4424 | 134.0 | −0.1250 |
| 1024 | 1200 | 861 | 0.7175 | 0.6953 | 289.3 | 0.7188 |

At `R = 1024` the measured 0.7175 sits on the analytic 0.7188 and the implied `c` is
289.3 against the form's 288, so the `1/R` shape is confirmed where it is testable
and is vacuous (negative) at `R = 64` and `R = 256`.  The carried value is the CP-95
lower bound on `1 − c/R` at `R = 2^14` with `c_hi` = 311.99 (max over the sweep of the
CP-95 upper failure rate times `R`):

    rho = 0.9809577775677195      (retention.json:retention_production.rho_lower)

against v14's asserted 0.982421875 — slightly worse, measured, and it is the number
the ledger uses (`constants.RHO_MEASURED`, asserted equal to the JSON in
`certificate.json`).  No census at `R = 2^14` exists here and none is claimed
(`failed_gate[2]:production_scale_retention_census`).

## 5. Success, with the measured retention inside it

`Lambda = S^2 (1-f^2) / (2 mu N)` with `f = 3/8` (routing share; one number in the
price, the budget and the model) and `mu = 11/10` (random-map reserve; it only
inflates `S`, and 4b(a) shows the measured marginal departs in the opposite
direction).  Target `0.392 / rho` = 0.39960945207240445 (`constants.TARGET_PAD`), so
`S = 2^125.828667`, model success 0.3996095, and

    net success = rho * 0.3996095 = 0.3920000  >=  0.39.

## 6. Records, the report CAP, and peak storage (v14 DEFECT 4)

A record is `(mkey, a, b)` = 112 B (`constants.RECORD_BYTES` = 48 + 2x32,
`REPORT_WORDS` = 14).

* **The report count is a CAP, not an expectation.**  `records = ceil(S/2^80)` =
  62 488 992 187 645 (`ledger.storage:records`); capacity is the cap and the
  overflow policy is stated and charged — a report that would exceed the cap drops
  the oldest record and re-probes that class, priced as `aux.report_cap_drain` =
  `cap x 420 / C` = 2^40.255.  Since `reports <= cap`, the bound holds on every run.
  v14's `floor(S/2^80)` is gone.
* **Peak storage adds every term the verdict named**: store 6 998 767 125 016 240 B;
  labels 256 B; **U labels** 6 291 456 B (`storage:u_label_bytes` = 2^14 x 12 x 32);
  advice 113 246 208 B; **guard history** 192 B (`:guard_history_bytes` = 4 trails x
  12 slots x 4 B); **code** 262 144 B (`:code_bytes`); scratch 4 608 B; peak
  6 998 767 244 821 104 B (`:peak_bytes`), recomputed term by term in
  `verify_certificate.py`.

## 7. Auxiliary work: GT arithmetic, not scalar instructions (v14 DEFECT 4)

v14 charged `reports x 2 x 381 / C` scalar instructions and called them two GT
exponentiations.  v15 charges the tower:

* a dense Fp12 multiplication **is** the reference operation (C = 20 000), so it
  costs exactly 1;
* a GT squaring is the counted Granger–Scott cyclotomic row,
  `K.cyclo_sqr_instructions()` = 24 Fp-mul + 122 Fp-add = 11 632 instructions
  (`= ledger.rows()['cyclo_sqr']`, one definition, asserted in the certificate),
  i.e. 0.58160 reference units — 24 Fp multiplications, not 54, because the
  cyclotomic form removes the Fp2 multiplications;
* a 381-bit exponentiation is 381 squarings plus at most 191 dense multiplications =
  411.59 reference units (`ledger.gt_exp_reference_units`);
* two per trail start, starts equal reports: `aux.gt_exp_starts` = 2^55.517
  (`ledger.aux`), up from v14's 2^41.194.

Shared-store lookup/index is priced at 420 instructions per report (48-byte key hash,
6 words x 30, plus two probe steps of a 112-byte compare, 2 x 120):
`aux.store_lookup` = 2^40.255.  Aux total 2^55.517 = 2^-70.26 of `W_stress`, so the
repair costs 0.0014 bits of claim, not zero.

## 8. Census, regenerated from the evidence (v14 DEFECT 3)

Every census sentence is a field of `evidence/pattern-census.json` or
`evidence/census-joint.json`, and `gen_cert.py` re-reads those fields.

* **Full joint state, shipped depth.**  `proven_census.pattern_census` compares
  `(element, ring tuple)` — not the element alone and not "an element reached from an
  empty ring".  At `R = 8`, `D = 12`, exhaustive over all 9420 starts for 26 steps:
  `joint_returns_total` = 0, `joint_min_lag` = null
  (`pattern-census.json:guarded_R8.joint_min_lag`).
* **The element ladder is NOT all-null and v15 says so.**  The same run has
  `element_returns_total` = 12, all at lag 15 (`:guarded_R8.element_return_lags`): an
  element can return while the joint state does not — exactly the distinction v14's
  T1 assertion ("an element return restores the empty ring") got wrong.  At
  `R = 64`, 26 steps, `:guarded_R64.element_returns_total` = 0, `min_lag` = null.
* **Mutation control.**  The depth-0 rule at `R = 8` returns at `min_lag` = 2 with
  588 element returns and 24 two-cycles (`:control_R8`), reproducing the family the
  v13 verdict exhibited.
* **Long window, disclosed.**  A 256-start x 400-step window at `R = 8` finds a
  repeated joint state in 156 of 256 starts with `min_first_repeat_lag` = 3
  (`census-joint.json:long_subset`): at `R = 8` the 9420-element fixture is
  cycle-dominated, so it is used as a *mutation control*, not as a rate.
* **Trap rate.**  Priced from the union bound of 4a (2^-50.7 over the walk), with the
  census-observed zero as a consistency row, not as the source of the number.  Scope:
  fixture `r = 9421`, `R = 8/64`, cycle-dominated; no production recovery rate is
  inferred from the census rates in `evidence/toy-census.json:t4`, which remain
  mechanism evidence only.

## 9. What is still not met

`certificate.json:failed_gate` carries three rows: the six-column route does not
exist (`no_inversion_free_sigma_stable_column_table`, FAIL, structural); the class
marginal is **not** uniform at fixture scale
(`class_marginal_uniformity_refuted_at_fixture`, FAIL, with 4b(a)'s direction
argument making `N` conservative); and no retention census at `R = 2^14` exists
(`production_scale_retention_census`, NOT_MET, extrapolated along the measured
`R`-dependence of 4b(b)).  Success remains conditional on the three disclosed
hypotheses in the claim text: random-map walk quality within the 10 per cent reserve,
pseudorandom routing fields, ideal coins.  No practical break is claimed.

## 10. Legality

Claim 125.7778 against the live lowest head 125.84 and the legality line 125.8345
(`ledger.claim_log2`, `certificate.json:claim`).  Nothing in this package depends on
the identity or validity of any other submission.
