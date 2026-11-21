# `ADAPTIVE.md` — Three-Part `[7e]` Conditioning & Exact-Degree `(9,3)` Sibling-Retry Verification (`DISC = 0.001`)

All checks below run live inside `audit_cm.py` (pure Python stdlib, `checks_ok = true`). We explicitly separate **empirical / structural conditioning evidence** from **production-scale 381-bit lattice checks** and **declared exploratory heuristics (`H1/H5`)**.

## 1. Per-Node Certified Yield Discard (`DISC = 0.001`)

During each special-`q` node's bound pass, the algorithm sums the modeled cell yield `ss_f * ss_g` over the ball `S_16(R_n)` (with `q` removed). Any node whose certified mean yield is below its bin's charged `pr` is immediately discarded as a failure before cofactorization. Consequently, **every used node's own `H1` yield is provably `>=` its charged `pr`**, and adaptive child selection can act only through the discard rate (`DISC = 0.001`, explicitly added inside `phi(node) = max((1 + mu*(1-psi)/r)^(-r) + DISC, 0.02)`) and the conditional failure law measured in `[7e]` below.

## 2. Three-Part `[7e]` Conditioning & Sibling-Retry Campaign (`q15_cond7e_live`)

- **Part (A) — Production `U_BLS` Encountered vs Matched Random Children (`120` vs `120` ideals, ACTUAL degree-`(72,24)` polynomials):**
  - `120` genuine encountered child prime ideals dividing actual cells of production `Q(zeta_15)` special-`q` parent balls vs `120` matched independent random ideals of equal side & bit-size (all reduced with exact Bareiss-LLL `99/100`, `det == p`, evaluated with ACTUAL `U_BLS = -0xd201000000010000`):
  - `f`-side offset diff `enc - rand = -0.131 +/- 0.910` bits (`KS p = 0.8760`), `g`-side diff `-0.562 +/- 0.462` bits (`KS p = 0.1223`), `max_cf = 79.9 < 131.91`, `max_cg = 23.15 < 46.50`, **`discards = 0/120`**.
  - `audit_cm.py` live-re-executes `8/8` witnesses (`4 enc + 4 rand`) with exact Bareiss-LLL and ACTUAL `U_BLS`.

- **Part (B) — Exact-Degree `(72,24)` `(9,3)` One-Level Sibling-Retry Campaign WITH REAL SIBLING FAILURES:**
  - Across `33` `(9,3)` parents yielding `76` distinct `O_H`-coprime sibling relations and `111` encountered children (explicitly excluding the parent's own relation from child scans to prevent cycle double-counting!) vs `111` size-matched random ideals:
    * Encountered vs random child retained counts: **`M_enc = 1.21 +/- 0.11` (`emp0 = 0.3153`)** vs **`M_rand = 1.38 +/- 0.14` (`emp0 = 0.3423`)**, **`KS p = 0.9262`** (no measurable selection bias).
    * **Direct conditioning after earlier sibling failures:** Unconditional relation failure **`P(fail) = 0.3947` (`30/76`)** vs conditional failure after earlier siblings failed **`P(fail | prev failed) = 0.2667` (`4/15`)** — no degradation after earlier sibling failures!
    * **Cross-relation sibling failure clustering:** Observed joint pair failure `P(both fail) = 0.1273` vs independent product `P(fail)^2 = 0.1558`, sibling pair Pearson correlation **`corr_sib = -0.0784`** (no positive clustering across distinct sibling relations).
  - **Live `(9,3)` sibling retries inside `audit_cm.py`:** `q15_cond7e_live` live-re-executes **`6/6` actual `(9,3)` parents** that fail their first relation and recover via subsequent distinct sibling relations (**`live_93_sibling_retries >= 6`** live on the exact degree-`(72,24)` `(9,3)` polynomials!).

- **Part (C) — Production `U_BLS` Smoothness-Weighted Encountered Law (`n = 120`, `ESS = 19.4`):**
  - Weighting each encountered child by its parent cell's modeled smoothness probability gives unweighted vs weighted `f`-headroom `81.08b` vs `79.77b` (`diff = -1.31 +/- 3.04b`, well within `2 SE`), `g`-headroom `36.32b` vs `35.68b`, and **`discards = 0/120`**.

## 3. Additional Supporting Checks (`q15_93_b16_live`, `prod_sq24_live`, `prod_cond_q15_live`, `voronoi_and_boot_live`)

- **`q15_93_b16_live`:** `120` `(9,3)` ideals (`40 g + 40 f + 40 enc`) verify `NB(0.5 * mu_raw, r = 1.0)` PGF domination (`wpgf_r1 < 0.90`, `k_hat < 0.25 << 1.0`, `emp0 < nb_r1_empty < nb_r05_empty`) plus `100` `hi100` ideals (`0/100` zeros) and live direct-16-vector Bareiss-LLL regeneration.
- **`prod_sq24_live` & `prod_cond_q15_live`:** `24 + 12` live 381-bit single-special-`q` Bareiss-LLL lattices (`det == q`, `q` in `73..258` bits, `0` pre-supplied vectors): all 6 `f`/`g` norm sextiles strictly below `qbins`, `|mean_rho| < 0.04`, `>= 10x` modeled yield headroom, and `Decimal(50)` Voronoi full unit group (`9` generators).
- **`voronoi_and_boot_live`:** Live 16D Bareiss-LLL (`99/100`) on `1` independent random `G_T` target (`136` pairs, `s_boot = 2^21`, `0` pre-supplied LLL vectors): first-vector norm `2542.56b` and deterministic bound `2591.9b < m_sv = 2592.0b`, `gcd_gt1e4 == 0/136`. Boot semismoothness (`pr_b = 2^-39.5120`) and sibling-retry independence are explicitly declared as exploratory `H1`/`H5` heuristics (supported, not proved at 381-bit scale).
