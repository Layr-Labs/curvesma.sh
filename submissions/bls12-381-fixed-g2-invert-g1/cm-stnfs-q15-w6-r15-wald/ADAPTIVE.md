# `ADAPTIVE.md` — Selected-Slot Conditional `H5` Premise & Live Strict-`k = 1` `(9,3)` Campaign

All checks below run live inside `audit_cm.py` (pure Python stdlib, `checks_ok = true`).

## 1. Explicit Scope of the Selected-Slot `H5` Premise (No Unconditional Marginal Deduction Claimed)

Because strict `k = 1` descent takes the **first** surviving relation at each parent, the law of an encountered child ideal `I` at slot `(b, s)` is **selection-conditioned** (not an independent generic marginal draw). Accordingly, we **explicitly declare** the selected-slot conditional bound as our exploratory STNFS `H5` heuristic:
- **`H5` (Selected-Slot Conditional Emptiness Bound):** For every reached first-choice slot `(b, s)`, `P(M(I) = 0 | reached first-choice slot (b, s)) <= p_0(b, s)`, where `NB(mu, r = 1.5)` `p_0(b, s) = max(5e-05, (1 + mu / 1.5)^(-1.5))` (excess dispersion `<= 1/1.5 = 0.6667` (`3x..20x` more conservative than measured `k_hat <= 0.164` across `120` `(9,3)` ideals in `q15_93_b16_live`)).
- **Why pure `k = 1` Boole (and NOT naive PGF substitution):** `verify_threshold_counterexample` inside `q15_k1_strict_h5_live` exhibits an exact algebraic witness (`0.10*t + 0.90*t^100 <= (1 + 30(1-t))^(-1)` for all `t in [0, 0.8]`, yet `0.10 > G(0.10) = 0.03571`) proving naive PGF composition is invalid for general non-i.i.d. laws. Our descent therefore uses **ONLY** Boole's union bound (`Phi_tree = sum_v p_0(v) <= 0.60`) across slots, wrapped in genuinely i.i.d. fresh-`e` boot retries (`Wald` identity).
- **Honest empirical scope:** The live campaigns below do **not** prove `H5` at 381-bit production scale; they stress-test the exact selected-slot conditional behavior directly on our **exact degree-`(72,24)` `(9,3)` `Q(zeta_15)` polynomials** and on 381-bit production `U_BLS` encountered-child lattices.

## 2. Live Strict-`k = 1` Zero-Retry First-Survivor Multi-Level `(9,3)` Descent (`q15_k1_strict_h5_live`)

Executed **LIVE** inside `audit_cm.py` on our **exact degree-`(72,24)` `(9,3)` `Q(zeta_15)` polynomials** (`f(X) = X^9 + (w+1)X^6 + X^3 - w`, `g(X) = X^3 + 5 w^2`, `U_toy = -5`, `16D` Bareiss-LLL `99/100`, `LP = 2^24`, `B_stop = 2048`, `N_pts = 448`, taking ONLY the first surviving relation at every node — ZERO sibling retries):
- **Interleaved Disjoint Control Calibration (`287` fresh ideals, bits `12..17`):** Control primes (`even` index in every bit-band) and test roots (`odd` index in every bit-band) are 100% disjoint and distribution-matched within every bit-band, calibrating exact per-`(side, bit)` control empty probability `p_cal(s, b)`.
- **Live `36`-Root Strict-`k = 1` Multi-Level Execution:** Runs `_desc_k1_slot` live on `36` fresh `(9,3)` roots:
  - Exact partition `nonempty == completed + slot_zeros` (`14 == 13 + 1`),
  - **`36` reached first-survivor children** across **depth `<= 3`**,
  - **Direct per-selected-slot `O/E_p0 = 1 / 6.650 = 0.1504 <= 0.85`** (no adverse first-survivor conditioning penalty!),
  - Replays **`3/3` bit-exact multi-level `k = 1` witnesses** (2 completed at depth `2` and `3`, 1 empty-slot termination).
- **Disjoint `221`-Root Offline Companion (`q15_k1_h5.b85`):** Exact partition `54 == 48 + 6`, `75` reached first-survivor children across depth `<= 3`, **direct per-selected-slot `O/E_p0 = 6 / 20.060 = 0.2991 <= 0.85`**.

## 3. Toy Degree-`(72,24)` `(9,3)` Dispersion & High-`mu` Zero-Free Check (`q15_93_b16_live`)

- Across `120` degree-`(72,24)` `(9,3)` ideals (`40 g + 40 f + 40 enc`, `U_toy = -5`), evaluated at the model's unit-thinned mean `0.5 * mu_raw`: measured excess dispersion `k_hat <= 0.1640` (`fano < 1.5`), empirical empty rates strictly below the claimed `p_0` bound on all three strata.
- Across `100` higher-mean `(9,3)` ideals (`hi100`: `40 enc + 30 g + 30 f`): **`0/100` zeros**, `min_M >= 2`, `k_hat < 0.25`.
- Live direct 16-vector `Z`-basis Bareiss-LLL regeneration on `camp93["f"][0]` confirms `>= 5` exact coprime relations.

## 4. Production `U_BLS` Encountered-Child Lattice & Boot Checks (`q15_enc_bias_live`, `prod_sq24_live`, `voronoi_and_boot_live`)

- **Encountered vs Matched Random (`120` vs `120` ideals, ACTUAL `U_BLS = -0xd201000000010000`):** `f`-side offset diff `enc - rand = -0.131 +/- 0.910` bits (`KS p = 0.8760`), `g`-side diff `-0.562 +/- 0.462` bits (`KS p = 0.1223`), `0/120` above `qbins`, `8/8` live-reverified; smoothness-weighted `ESS = 19.4` (`diff = -1.31 +/- 3.04b`, `0/120` discards).
- **`prod_sq24_live` & `prod_cond_q15_live`:** `24 + 12` live 381-bit single-special-`q` Bareiss-LLL lattices (`det == q`, `q` in `73..258` bits, `0` pre-supplied vectors): all 6 `f`/`g` norm sextiles strictly below `qbins`, `|mean_rho| < 0.04`, `>= 10x` modeled yield headroom.
- **`voronoi_and_boot_live`:** Live 16D Bareiss-LLL (`99/100`) on `1` independent random `G_T` target (`136` pairs, `s_boot = 2^21`, `0` pre-supplied LLL vectors): first-vector norm `2542.56b` and deterministic Lovász bound `2591.9b < m_sv = 2592.0b`, `gcd_gt1e4 == 0/136`.
