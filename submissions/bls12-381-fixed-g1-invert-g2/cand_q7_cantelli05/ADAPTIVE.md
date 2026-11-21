# `ADAPTIVE.md` — Selected-Slot Conditional `H5` Premise & Live Strictly-Descending `k = 1` `Q(zeta_15)` Campaign

All checks below run live inside `audit_cm.py` (pure Python stdlib, `checks_ok = true`, `< 6.5s` locally).

## 1. Explicit Scope of the Selected-Slot `H5` Premise (No Unconditional Marginal Deduction Claimed)

Because strict `k = 1` descent takes the **first** surviving relation at each parent, the law of an encountered child ideal `I` at slot `(b, s)` is **selection-conditioned** (not an independent generic marginal draw). Accordingly, we **explicitly declare** the selected-slot conditional bound as our exploratory STNFS `H5` heuristic:
- **`H5` (Selected-Slot Conditional Emptiness Bound):** For every reached first-choice slot `(b, s)`, `P(M(I) = 0 | reached first-choice slot (b, s)) <= p_0(b, s)`, where one-sided Cantelli `p_0(b, s) = max(0.0001, (1 + 0.5 * mu)^(-1))` (`Var(M(I) | reached first-choice slot) <= max(2 * mu, 0.0001 * mu^2)` (pure one-sided Cantelli inequality, ZERO `NB(r > 1)` parametric tail!)).
- **Why pure `k = 1` Boole (and NOT naive PGF substitution):** `verify_threshold_counterexample` inside `q15_k1_strict_h5_live` exhibits an exact algebraic witness (`0.10*t + 0.90*t^100 <= (1 + 30(1-t))^(-1)` for all `t in [0, 0.8]`, yet `0.10 > G(0.10) = 0.03571`) proving naive PGF composition is invalid for general non-i.i.d. laws. Our descent therefore uses **ONLY** Boole's union bound (`Phi_tree = sum_v p_0(v) <= 0.60`) across slots, wrapped in genuinely i.i.d. fresh-`e` boot retries (`Wald` identity).
- **Honest empirical scope:** None of the toy-scale or single-step live tests below prove `H5` at 381-bit production scale; they stress-test strictly-descending selected-slot behavior over `O_H = Z[zeta_15]` (`deg-(24,24)` multi-level campaign with 100% faithful selection & direct modeled `p_0` denominator) alongside exact `deg-(72,24)` `(9,3)` single-node dispersion/`hi100` checks and 381-bit `U_BLS` encountered-child lattice checks.

## 2. Live 100%-Faithful Strictly-Descending `k = 1` `Q(zeta_15)` Campaign (`q15_k1_strict_h5_live`)

Executed **100% LIVE** inside `audit_cm.py` (`w7_verifiers.py`, zero precomputed control tables) on the degree-`(24,24)` `Q(zeta_15)` pair `F(X) = X^3 - (2+w)w^2`, `G(X) = X^3 - 2 w^2` (exact same `O_H = Z[zeta_15]` and `16D` Bareiss-LLL `99/100`, `B_FB = 8192`):
- **100% Faithful Strictly-Descending Selection (`zero_omitted_primes = true`):** `_strict_fac` strips primes `<= 8192` and **REJECTS** any relation having ANY residual prime `> 8192` that is `>= q_lim`, not `1 mod 15`, or lacks an explicit `_ch_root` witness. Empty child lists occur **ONLY** when both norms are `100% B_FB`-smooth, and every queued child explicitly asserts `8192 < cp < q`.
- **Direct Modeled Claimed `p_0(mu_model)` Denominator (`direct_modeled_p0_denominator = true`):** At every reached first-survivor child node, `mu_model = 0.50 * mu_raw` (matching production `udup = 0.50`) is computed directly from the smooth + 1-LP Dickman integral and plugged directly into the submission's exact claimed bound one-sided Cantelli `p_0(b, s) = max(0.0001, (1 + 0.5 * mu)^(-1))`:
  - **Boundary regime (`36` live `17`-bit roots, `n_pts = 32`):** `O/E = 11 / 21.400 = 0.5140 <= 0.85` (`ne = 62 = 51 + 11`, `depth = 4`, `obs/mod mu = 1.118`), `3/3` bit-exact multi-level witnesses (`depth = 4`).
  - **Low-`p_0` regime (`16` live `17`-bit roots, `n_pts = 160`):** `O/E = 0 / 5.268 = 0.0000 <= 0.50` (`ne = 31 = 31 + 0`, `depth = 4`, `obs/mod mu = 1.319`).
  - **Disjoint matched control (`36` independent `15–16`-bit ideals, zero overlap with campaign):** `O/E = 3 / 11.914 = 0.2518 <= 0.85`.

## 3. Exact Degree-`(72,24)` `(9,3)` Dispersion & High-`mu` Zero-Free Check (`q15_93_b16_live`)

- Across `120` exact degree-`(72,24)` `(9,3)` ideals (`40 g + 40 f + 40 enc`, `U_toy = -5`), evaluated at unit-thinned `0.5 * mu_raw`: measured excess dispersion `k_hat <= 0.1640` (`fano < 1.5`), empirical empty rates strictly below the claimed Cantelli `1/(1 + 0.5*mu)` bound on all three strata.
- Across `100` higher-mean `(9,3)` ideals (`hi100`: `40 enc + 30 g + 30 f`): **`0/100` zeros**, `min_M >= 2`, `k_hat < 0.25`.
- Live direct 16-vector `Z`-basis Bareiss-LLL regeneration on `camp93["f"][0]` confirms `>= 5` exact coprime relations.

## 4. Production `U_BLS` Encountered-Child Lattice & Boot Checks (`q15_enc_bias_live`, `prod_sq24_live`, `voronoi_and_boot_live`)

- **Encountered vs Matched Random (`120` vs `120` ideals, ACTUAL `U_BLS = -0xd201000000010000`):** `f`-side offset diff `enc - rand = -0.131 +/- 0.910` bits (`KS p = 0.8760`), `g`-side diff `-0.562 +/- 0.462` bits (`KS p = 0.1223`), `0/120` above `qbins`, `8/8` live-reverified; smoothness-weighted `ESS = 19.4` (`diff = -1.31 +/- 3.04b`, `0/120` discards).
- **`prod_sq24_live` & `prod_cond_q15_live`:** `24 + 12` live 381-bit single-special-`q` Bareiss-LLL lattices (`det == q`, `q` in `73..258` bits, `0` pre-supplied vectors): all 6 `f`/`g` norm sextiles strictly below `qbins`, `|mean_rho| < 0.04`, `>= 10x` modeled yield headroom.
- **`voronoi_and_boot_live`:** Live 16D Bareiss-LLL (`99/100`) on `1` independent random `G_T` target (`136` pairs, `s_boot = 2^21`, `0` pre-supplied LLL vectors): first-vector norm `2542.56b` and deterministic Lovász bound `2591.9b < m_sv = 2592.0b`, `gcd_gt1e4 == 0/136`.
