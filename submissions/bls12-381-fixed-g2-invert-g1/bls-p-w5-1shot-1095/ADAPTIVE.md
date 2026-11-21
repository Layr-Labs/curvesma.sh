# `ADAPTIVE.md` — Supporting Checks for 1-Shot Marginal `NB(mu, r = 2.0)` Descent (`H5`)

All checks below run live inside `audit_cm.py` (pure Python stdlib, `checks_ok = true`). Because Wave-5 descent is **1-shot (`k = 1`, `K_boot = 1`)** with a **pure Boole union bound** over `6,314` deterministic slots, **no sibling-retry law, no conditional-after-failure independence, no `K_boot > 1` power law, and no `0.001` discard tail assumption** are used anywhere. The *only* `H5` premise needed is the **marginal** per-node empty probability `Pr[M_v = 0] <= (1 + mu_v / 2.0)^(-2.0)` for each encountered node `v`.

## 1. Toy Degree-`(72,24)` `(9,3)` Marginal Overdispersion & PGF Check (`q15_93_b16_live`)

- Across `120` degree-`(72,24)` `(9,3)` ideals (`40 g + 40 f + 40 enc`, `U_toy = -5`), evaluated at the **model's** unit-thinned mean `0.5 * mu_raw` against the exact claimed `NB(mu, r = 2.0)` (`1/r = 0.5000`):
  - `g` stratum (`n = 40`): measured excess dispersion `k_hat = 0.1640` (**`3.0x` below `1/r = 0.5000`**), `emp0 = 0.1750 < mod0 = 0.4891`, `worst_pgf_r < 0.85`.
  - `f` stratum (`n = 40`): `k_hat = 0.0328` (**`15.2x` below `1/r`**), `emp0 = 0.0000 < mod0 = 0.3380`, `worst_pgf_r < 0.69`.
  - `enc` stratum (`n = 40`): `k_hat = 0.0375` (**`13.3x` below `1/r`**), `emp0 = 0.2000 < mod0 = 0.5155`, `worst_pgf_r < 0.81`.
- Across `100` higher-mean `(9,3)` ideals (`hi100`: `40 enc + 30 g + 30 f`): `0/100` zeros, `min_M >= 2`, `k_hat < 0.25`.
- Live direct 16-vector `Z`-basis Bareiss-LLL regeneration on `camp93["f"][0]` confirms `>= 5` exact coprime relations.

## 2. Production `U_BLS` Encountered-Child Norm Bias Check (`q15_enc_bias_live`)

- **Encountered vs Matched Random (`120` vs `120` ideals, ACTUAL `U_BLS = -0xd201000000010000`):**
  - `120` genuine encountered child prime ideals dividing actual cells of production `Q(zeta_15)` special-`q` parent balls vs `120` matched independent random ideals of equal side & bit-size (all reduced with exact Bareiss-LLL `99/100`, `det == p`):
  - `f`-side offset diff `enc - rand = -0.131 +/- 0.910` bits (`KS p = 0.8760`), `g`-side diff `-0.562 +/- 0.462` bits (`KS p = 0.1223`), `max_cf = 79.9 < 131.91`, `max_cg = 23.15 < 46.50`, **`0/120` above `qbins`**.
  - `audit_cm.py` live-re-executes `8/8` witnesses (`4 enc + 4 rand`) with exact Bareiss-LLL and ACTUAL `U_BLS`.
- **Smoothness-Weighted Encountered Law (`n = 120`, `ESS = 19.4`):**
  - Weighting each encountered child by its parent cell's modeled smoothness probability gives unweighted vs weighted `f`-headroom `81.08b` vs `79.77b` (`diff = -1.31 +/- 3.04b`, within `2 SE`), `g`-headroom `36.32b` vs `35.68b`, and **`0/120` above `qbins`**.

## 3. Production-Scale 381-Bit Lattice & Boot Checks (`prod_sq24_live`, `prod_cond_q15_live`, `voronoi_and_boot_live`)

- **`prod_sq24_live` & `prod_cond_q15_live`:** `24 + 12` live 381-bit single-special-`q` Bareiss-LLL lattices (`det == q`, `q` in `73..258` bits, `0` pre-supplied vectors): all 6 `f`/`g` norm sextiles strictly below `qbins`, `|mean_rho| < 0.04`, `>= 10x` modeled yield headroom, and `Decimal(50)` Voronoi full unit group (`9` generators).
- **`voronoi_and_boot_live`:** Live 16D Bareiss-LLL (`99/100`) on `1` independent random `G_T` target (`136` pairs, `s_boot = 2^21`, `0` pre-supplied LLL vectors): first-vector norm `2542.56b` and deterministic bound `2591.9b < m_sv = 2592.0b`, `gcd_gt1e4 == 0/136`. Boot semismoothness (`pr_b = 2^-39.5120`) and marginal `NB(mu, r = 2.0)` per-node emptiness are explicitly declared as exploratory `H1`/`H5` heuristics (supported, not proved at 381-bit scale).
