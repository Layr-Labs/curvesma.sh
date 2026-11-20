# Reduced-Tower `Q(zeta_15)` Adaptive Recursive Descent Evidence

To test the half-scale NB-PGF envelope (`r = 2, alpha = 0.5`) and conditional after-failure hazard (`THRESHOLD-REPAIR.md`) directly on a genuine `Q(zeta_15)` tower, we run an uncensored recursive descent campaign over the **exact production cyclotomic base** `O_H = Z[zeta_15]` (`Phi_15(y)`, `omega = y^5`, order-12 Frobenius `tau(y, X) = (y^7, omega^2 X)`).

## 1. Tower, Exact Constructive Rank-2 Lattice Balls, & Protocol

- **Reduced `tau`-invariant `(4, 3)`-base tower:** `g_t(X) = X^3 + 2 omega^2` (deg 24) and `f_t(X) = X^3 - (2 + omega)` (deg 24, `Res_X(f_t, g_t) = -1` coprime at every prime; both fixed by `tau(y, X) = (y^7, omega^2 X)`).
- **Exact arithmetic & constructive rank-2 submodule:** At every reached degree-1 ideal `Q = (side, q, y_0, r_0)` (`q == 1 mod 15`), run exact 16D integer LLL (`delta = 0.99`) on `L_Q`, extract `(v_1 = b_1, v_2 = b_k)` with `b_k` the first `O_H`-independent vector, enumerate the complete primitive rank-2 `O_H`-ball, compute **exact** integer norms in `__int128` (`0` overflows), factor both norms **exactly** (`B_FB = 4096`), extract every degree-1 child ideal `< shrink * q`, and deduplicate distinct qualifying relations modulo unit-orbit duplicates.
- **1,600 uncensored roots (`400` distinct root `tau`-orbits in `[12000, 131000]` across 4 regimes, `0` censors):** Three window-shrink factors `d = 1 (0.75)`, `d = 2 (0.55)`, `d = 3 (0.40)` plus a **shared-database missing-ideal** stress run (`d = 1-sharedDB`) where all `400` roots descend against the **same fixed 5%-missing FB set** (`34` fixed unresolved primes `<= 4096`).

## 2. Measured Count-PGF Ratios & Child-Bin Row-Specific Conditional `O/E`

For every stratum, `PGF_2` is `max_(t in {0, 0.01, ..., 0.80}) mean(t^M) / mean(G_2(t))` with `G_2(t) = (1 + 0.25 mu (1 - t))^(-2)` (`mu` = pre-entry bin mean). For every after-failure attempt `j >= 2` with actual children `{c}`, expected failure `E_fail(j) = 1 - prod_c (1 - p_emp(side_c, bin_c))` and expected child work `E_work(j) = sum_c w_emp(side_c, bin_c)` are calibrated from the unconditional bin rates of those actual children (removing child-mix confounding).

| Regime (`400` roots each, `0` censors) | Nodes (Max D) | Roots OK | `PGF_2` All / `d>=1` / After-Fail | After-Fail Att | `fail O/E` `[95% CI]` | `work O/E` `[95% CI]` | Shared Consec | Tree Repeats |
|---|---:|---:|---|---:|---|---|---:|---:|
| `d=1` (`shrink=0.75`) | `1,341` (`5`) | `56.5%` | `0.880` / `0.901` / `0.895` | `250` | `0.984` `[0.881..1.086]` | `0.805` `[0.705..0.905]` | `2/250` (`0.80%`) | `7/941` (`0.74%`) |
| `d=2` (`shrink=0.55`) | `945` (`3`) | `48.5%` | `0.897` / `0.924` / `0.934` | `148` | `1.105` `[0.985..1.225]` | `0.776` `[0.651..0.901]` | `4/148` (`2.70%`) | `2/545` (`0.37%`) |
| `d=3` (`shrink=0.40`) | `734` (`2`) | `42.8%` | `0.912` / `0.934` / `0.950` | `80` | `1.137` `[0.950..1.324]` | `0.842` `[0.657..1.026]` | `1/80` (`1.25%`) | `1/334` (`0.30%`) |
| `d=1-sharedDB` (5% miss) | `1,272` (`5`) | `45.0%` | `0.897` / `0.917` / `0.921` | `253` | `1.016` `[0.936..1.096]` | `0.734` `[0.644..0.825]` | `4/253` (`1.58%`) | `8/872` (`0.92%`) |

- **All `4,292` reached nodes & `731` after-failure attempts:** Half-scale `G_2(t)` strictly dominates the empirical PGF across every regime and pre-entry stratum (`max ratio <= 0.950 < 1.000`, and `<= 0.949` at `r = 1`). Every row-specific `fail O/E` 95% CI brackets `1.00` (including the fixed shared-DB missing-ideal run at `1.016 [0.936..1.096]`), and every `child_work O/E` is strictly below `0.85`.
- **Live pure-stdlib replay (`audit_cm.py`):** `threshold_and_adaptive_live.adaptive_q15_slice` replays a deterministic 4-root slice on this exact `Q(zeta_15)` tower live (`0` censors, exact 16D LLL rank-2 extraction, exact norm divisibility & factorization).

## 3. Scope and Extrapolation Limits

This reduced-tower campaign directly exercises multilevel adaptive descent and shared-database missing-ideal dependence on `Q(zeta_15)`, providing stronger evidence than pooled single-level counts. However, it uses a reduced degree-24 tower and small factor base (`B_FB = 4096`), and its subsets are correlated finite-scale ensembles; it does **not** prove a uniform production conditional hazard, fixed-ideal solvability, or independence at `B ~ 2^71.5`. Both the count-PGF envelope (`r = 2, alpha = 0.5`) and the conditional next-relation / child-work bounds in `H5` remain explicitly declared exploratory heuristics.
