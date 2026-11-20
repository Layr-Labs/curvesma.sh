# Reduced-Tower `Q(zeta_15)` Adaptive Recursive Descent Evidence

To test the half-scale NB-PGF envelope (`r = 2, alpha = 0.5`) and `1.25x`-padded conditional after-failure hazard & child-work envelope (`THRESHOLD-REPAIR.md`) directly on a genuine `Q(zeta_15)` tower, we run a **3,200-root** uncensored recursive descent campaign over the **exact production cyclotomic base** `O_H = Z[zeta_15]` (`Phi_15(y)`, `omega = y^5`, order-12 Frobenius `tau(y, X) = (y^7, omega^2 X)`).

## 1. Tower, Exact Constructive Rank-2 Lattice Balls, & Protocol

- **Reduced `tau`-invariant `(4, 3)`-base tower:** `g_t(X) = X^3 + 2 omega^2` (deg 24) and `f_t(X) = X^3 - (2 + omega)` (deg 24, `Res_X(f_t, g_t) = -1` coprime at every prime; both fixed by `tau(y, X) = (y^7, omega^2 X)`).
- **Exact arithmetic & constructive rank-2 submodule:** At every reached degree-1 ideal `Q = (side, q, y_0, r_0)` (`q == 1 mod 15`), run exact 16D integer LLL (`delta = 0.99`) on `L_Q`, extract `(v_1 = b_1, v_2 = b_k)` with `b_k` the first `O_H`-independent vector, enumerate the complete primitive rank-2 `O_H`-ball, compute **exact** integer norms in `__int128` (`0` overflows), factor both norms **exactly** (`B_FB = 4096`), extract every degree-1 child ideal `< shrink * q`, and deduplicate distinct qualifying relations modulo unit-orbit duplicates.
- **3,200 uncensored roots (`800` distinct root `tau`-orbits in `[12000, 195000]` across 4 regimes, `0` censors):** Three window-shrink factors `d = 1 (0.75)`, `d = 2 (0.55)`, `d = 3 (0.40)` plus a **shared-database missing-ideal** stress run (`d = 1-sharedDB`) where all `800` roots descend against the **same fixed 5%-missing FB set** (`34` fixed unresolved primes `<= 4096`, discarding any relation hitting a missing FB prime — matching production's `exp(-miss)` discard).

## 2. Measured Count-PGF Ratios & Root-Cluster Bootstrap 95% CIs vs Charged `1.25x` Envelope

For every stratum, `PGF_2` is `max_(t in {0, 0.01, ..., 0.80}) mean(t^M) / mean(G_2(t))` with `G_2(t) = (1 + 0.25 mu (1 - t))^(-2)` (`mu` = pre-entry bin mean). For every after-failure attempt `j >= 2` (stratified both overall and by `k = 2` vs `k >= 3`), expected failure `E_fail(j) = 1 - prod_c (1 - p_emp(side_c, bin_c))` and expected child work `E_work(j) = sum_c w_emp(side_c, bin_c)` are calibrated from the unconditional bin rates of that attempt's actual children (removing child-mix confounding). All 95% confidence intervals are computed via **2,000-resample root-cluster bootstrap** (resampling whole root trees with replacement to account for all within-tree correlation) against our package's charged **`1.25x`** envelope (`O / (1.25 E)`):

| Regime (`800` roots each, `0` censors) | Nodes (Max D) | Roots OK | `PGF_2` All / `d>=1` / After-Fail | After-Fail (`k=2` / `k>=3`) | `fail / (1.25 E)` `[95% boot CI]` | `k=2` & `k>=3` `[95% boot CI]` | `work / (1.25 E)` `[95% boot CI]` | Shared / Repeats |
|---|---:|---:|---|---:|---|---|---|---:|
| `d=1` (`0.75`) | `3,775` (`6`) | `60.5%` | `0.872` / `0.894` / `0.889` | `841` (`589`/`252`) | `0.782` `[0.734..0.831]` | `[0.705..0.825]` / `[0.741..0.892]` | `0.550` `[0.507..0.594]` | `1.07%` / `1.24%` |
| `d=2` (`0.55`) | `2,577` (`3`) | `50.8%` | `0.890` / `0.917` / `0.917` | `538` (`365`/`173`) | `0.777` `[0.727..0.825]` | `[0.714..0.833]` / `[0.695..0.865]` | `0.572` `[0.530..0.611]` | `1.49%` / `0.45%` |
| `d=3` (`0.40`) | `1,921` (`3`) | `43.6%` | `0.903` / `0.930` / `0.930` | `319` (`236`/`83`) | `0.817` `[0.755..0.874]` | `[0.722..0.863]` / `[0.756..0.972]` | `0.572` `[0.529..0.618]` | `1.25%` / `0.27%` |
| `d=1-sharedDB` (5% miss) | `3,593` (`5`) | `46.6%` | `0.891` / `0.913` / `0.910` | `895` (`607`/`288`) | `0.779` `[0.744..0.815]` | `[0.759..0.844]` / `[0.681..0.790]` | `0.509` `[0.475..0.547]` | `1.56%` / `1.36%` |

- **All `11,866` reached nodes & `2,593` after-failure attempts — EVERY 95% root-cluster bootstrap CI upper limit strictly `< 1.000`:**
  - Half-scale `G_2(t)` strictly dominates the empirical PGF across every regime and pre-entry stratum (`max ratio <= 0.930 < 1.000`, and `<= 0.931` at `r = 1`).
  - Against our charged `1.25x` conditional envelope, **`fail / (1.25 E)` upper 95% CI is `<= 0.874` overall (`<= 0.863` at `k = 2`, `<= 0.972` at `k >= 3`, and `<= 0.815` on the fixed shared-DB missing-ideal run)**, and **`work / (1.25 E)` upper 95% CI is `<= 0.618`** across all four regimes!
- **Formula-envelope slack (`phi_nb_pgf_union` vs observed):** Evaluating our exact package formula `phi_nb_pgf_union(mu, min(0.99, 1.25 psi_1), r=2, alpha=0.5)` yields observed-to-envelope root failure ratios `0.395/0.946 = 0.418` (`d=1`), `0.493/0.967 = 0.510` (`d=2`), `0.564/0.979 = 0.576` (`d=3`), and `0.534/1.000 = 0.534` (`d=1-sharedDB`) — **`obs/env <= 0.576`** across every setting.
- **Live pure-stdlib replay (`audit_cm.py`):** `threshold_and_adaptive_live.adaptive_q15_slice` replays a deterministic 4-root slice on this exact `Q(zeta_15)` tower live (`0` censors, exact 16D LLL rank-2 extraction, exact norm divisibility & factorization).

## 3. Scope and Extrapolation Limits

This 3,200-root reduced-tower campaign directly exercises multilevel adaptive descent and shared-database missing-ideal dependence on `Q(zeta_15)` with root-cluster bootstrap CIs. However, it uses a reduced degree-24 tower and small factor base (`B_FB = 4096`), and its subsets are correlated finite-scale ensembles; it does **not** prove a uniform production conditional hazard, fixed-ideal solvability, or independence at `B ~ 2^71.5`. Both the count-PGF envelope (`r = 2, alpha = 0.5`) and the `1.25x`-padded conditional next-relation / child-work bounds in `H5` remain explicitly declared exploratory heuristics.
