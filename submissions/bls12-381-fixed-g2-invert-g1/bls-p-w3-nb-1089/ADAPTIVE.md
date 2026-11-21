# `ADAPTIVE.md` — Supporting Empirical & Production-Scale Checks for `Q(zeta_15)` Descent (`H5`)

All checks below run live inside `audit_cm.py` (pure Python stdlib, `checks_ok = true`). We explicitly distinguish **small-scale structural / toy empirical support** from **production-scale 381-bit lattice checks** and **declared exploratory heuristics (`H1/H5`)**.

1. **Toy `(9,3)` Negative-Binomial PGF & Overdispersion Check (`q15_93_b16_live`):**
   - Evaluates `120` degree-`(72, 24)` `(9,3)` ideals (`40` `g`-side, `40` `f`-side, `40` encountered children `enc`, `U_toy = -5`) against `NB(0.5 * mu_raw, r = 1.0)` (at the **model's** unit-thinned mean `0.5 * mu_raw`, not only the observed mean; note `NB(mu, r = 0.5)` is strictly heavier-tailed than `NB(mu, r = 1.0)`).
   - Across all three groups: `fano < 1.5`, excess overdispersion `k_hat < 0.25 << 1.0 = 1/r`, empirical zero rate `emp0 < nb_r1_empty < nb_r05_empty`, and worst-case probability generating function ratio `wpgf_r1 < 0.90` over `t in [0, 0.80]`.
   - Across `100` higher-mean `(9,3)` ideals (`hi100`: `40 enc + 30 g + 30 f`): `0/100` zeros, `min_M >= 2`, `k_hat < 0.25`.
   - Live direct 16-vector `Z`-basis Bareiss-LLL regeneration on `camp93["f"][0]` confirms `>= 5` exact coprime relations.

2. **Live 12-Root Unblinded Sibling-Retry Tree Sanity Check (`q15_nb_tree_live`):**
   - Runs full-depth recursive single-special-`q` sibling-retry descent live on `12` fresh `q ~ 2^16` prime ideals (`6 g + 6 f`, `q == 1 mod 15`, `B_fb = 8192`) with an explicit node counter and ancestor-cycle check.
   - Confirms `>= 11/12` roots complete down to `B_fb` reaching depth `>= 2` (`>= 8` depth-1 and `>= 8` depth-2+ nodes). This is a small-scale algorithmic/structural sanity check, NOT a proof of production `H5`.

3. **Production-Scale 381-Bit Single-Special-`q` Lattice Checks (`prod_sq24_live` & `prod_cond_q15_live`):**
   - **`prod_sq24_live`:** Live Bareiss-LLL (`99/100`) on `24` production 381-bit single-special-`q` lattices (`12 g + 12 f`, `q` spanning `73..258` bits, `det(R) == q`, `0` pre-supplied vectors). Verifies all 6 `f`/`g` norm sextiles strictly below `qbins`, `sv_emp <= sv_plan`, modeled smoothness yield `>= 10x` charged `pr`, and pairwise small-prime co-divisibility `|mean_rho| < 0.04`.
   - **`prod_cond_q15_live`:** Live Bareiss-LLL (`99/100`) on `12` production 381-bit encountered-child special-`q` lattices (`det(R) == q`, exact parent divisibility `12/12`, all 6 sextiles below `qbins`, `>= 10x` modeled yield headroom), plus pure `Decimal(50)` Voronoi fundamental domain (`9` generators of full `O_H^*`).

4. **Production-Scale 381-Bit `G_T` Boot Lattice Check (`voronoi_and_boot_live`):**
   - Runs live 16D Bareiss-LLL (`99/100`) on **`1` independent random `G_T` targets (`136` pairwise comparisons)** in default mode (`boot_witness.json`, `s_boot = 2^21`, `0` pre-supplied LLL vectors), checking exact determinant `p^8 * s_boot^16`, Lovász conditions, norm deciles below `m_sv = 1525.8`, and `gcd_gt1e4 == 0 / 136`.
   - The boot semismoothness probability (`pr_b = 2^-39.5120`) and sibling-retry independence are explicitly declared as exploratory `H1`/`H5` heuristics (supported, not proved at 381-bit scale).
