# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, conservative descent main base `MAIN_LB = 2^71.29` (`= min(B_f, B_g)`), `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (below `lr_rethi - 0.009 - 2 SE = 63.4846`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes** (`rho_dB(u, x) = rho(u) - gamma rho(u - 1) / ln x`), `BOOT_RELS = 48`, `RET_BOOT = 0.50` (deterministic 50% norm-certificate discard filter), `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

**100% Live Semismoothness Execution (Zero Precomputed Probabilities or Costs):** `cm_data.json` stores ONLY the 78 structural node 9-tuples `[lo, hi, side, n, iA32, kA, kB, R, thin]`, the empirical `q`-lattice offsets, and the boot deciles (`m_max = 2600.0`b). Inside Docker, `audit_cm.py` invokes `cm_models.eval_descent_live` to compute EVERY node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `psi`, `W`, `phi`, `EN`, and boot `pr_log2`, `miss`, `psi`, `phi_boot`, `two_log`, `E[N]` 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing the ledger at the coarser default `MESH = 96/192` lower sum (`104.980946 <= 104.981`) is strictly more conservative than `MESH = 384/384` (`104.848708`, runnable on demand via `{"descent_mesh": 384}`).

## 1. Fresh-Coin Target Randomization, 16D Guillevic Boot (`det = p^4`), & 50% Norm-Certificate Filter

For each of the two discrete logarithms (`log g` with `g = e(P0, Q)` and `log z`), draw an independent uniform exponent `e in {1, ..., r - 1}` (fresh algorithm coins; cost `< 2^40` instructions inside the wrapper cap) and boot on `z' = z^e` (resp. `g' = g^e`). Because `G_T` has prime order `r`, BOTH `g'` and `z'` are **exact uniform on `G_T \ {1}`** from fresh coins alone, even when `g = e(P0, Q)` is fixed.

### Algebraic Non-Degeneracy Theorem (`100%` of `G_T \ {1}`)

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)` (since `tau(y, X) = (y^7, omega^2 X)` satisfies `tau^6(y, X) = (y^4, X)`):
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r | Phi_12(p)`, so every nonzero element of `F_(p^4)^* * F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`) and `c_6 = c_(6,0) + c_(6,1) X + c_(6,2) X^2` (`c_(6,i) in F_(p^2)`). Vanishing of the `X^2`-coefficient of `c_6 z'`, namely `c_(6,0) z_2 + c_(6,1) z_1 + c_(6,2) z_0 = 0` in `F_(p^4) ~= F_(p^2)^2`, is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore **always** has a nonzero solution `c_6 in F_(p^6)^*`.
- Writing `c_6 z' = w_0 + w_1 X`: if `w_1 = 0` or `w_0 = 0`, then `c_6 z' in F_(p^4)^*` or `F_(p^4)^* X`, forcing `z' in F_(p^4)^* F_(p^6)^*`, impossible for `z' in G_T \ {1}` of prime order `r`. Setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16-dimensional linear lattice `L_Guillevic = { (a, b) in Z[zeta_15]^2 : a - gamma b == 0 (mod P) } subset Z^16`. Because `O_H / P ~= F_(p^4)` imposes 4 `F_p`-linear congruences, **`det(L_Guillevic) = p^4`** (`~ 2^1524`), and evaluating **`a - b X`** gives **`Psi_g(a - b X) = (a mod P) - (b mod P) X = (b mod P)(gamma - X) = -(b mod P) w_1^(-1) c_6 z'`**. Since `-(b mod P) w_1^(-1) in F_(p^4)^*` and `c_6 in F_(p^6)^*` both have order coprime to `r`, `Psi_g(a - b X)` has exact discrete log `e log_r(z) mod r` (and `|N_(K_H/Q)(b)| ~= 2^680 << p^4` guarantees `b != 0 mod P`). Its exact relative norm in `K_g = K_H[X]/(X^3 - omega^2 u)` is **`N_(K_g/K_H)(a - b X) = Res_X(a - b X, X^3 - omega^2 u) = a^3 - omega^2 u b^3`** (shortest-vector `g`-norm `~ p^(20/3) ~= 2538` bits).

### Free Rank-2 `Z[zeta_15]`-Module Enumeration, Live 4-Target Witness, & Deterministic 50% Norm Filter

Because `K_H = Q(zeta_15)` is a PID (`h = 1`), `L_Guillevic = Z[zeta_15] v_1 + Z[zeta_15] v_2` is a free rank-2 `Z[zeta_15]`-module. Under LLL with the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)`, enumerating `v = 1 * v_1 + beta * v_2` for `beta in B_8(1)` (`17` vectors: `0` and `+-e_i`) yields **zero unit duplicates (`UDUP_BOOT = 1.000`)** and **100% base-primitivity (`gcd(1, beta) = 1`)**.

- **Live 4-target BLS12-381 `w'=omega^2` witness (`boot_witness.json`, verified LIVE in `audit_cm.py`):** Four real reduced BLS12-381 `w' = omega^2` rank-2 boot bases `(gam, v_1, v_2)` are bundled inside the package. `audit_cm.py` (`voronoi_and_boot_live`) verifies LIVE on all 4 bases: exact `F_(p^4)` congruence `a_i == gamma b_i (mod P)` (yielding `a - b X == b (gamma - X) mod P`), exact determinant `|N_(K_H/Q)(a_1 b_2 - a_2 b_1)| == p^4`, all `4 x 17 = 68` exact 2,550-bit relative norms `|N_(K_H/Q)(a^3 - u omega^2 b^3)|` of `a - b X` (mean `2551.75`b, max `2563.01`b `< 2600`b, pairwise `gcd > 10^4` rate `3/544`), AND **`4/4` full passes of the exact 10-decile discard predicate** (`ball_sorted[idx_j] <= boot_sub_deciles[j]` for every decile `j = 0..9`).
- **Deterministic 50% boot-ball norm-certificate filter (`RET_BOOT = 0.50`, `m_max = 2600.0`b):** To eliminate any tail sensitivity beyond a finite empirical sample, after LLL on each fresh boot target the algorithm computes all 17 exact norms of `a - b X` in `B_8(1)` and **deterministically discards** any boot ball whose sorted deciles exceed `boot_sub_deciles` (top decile `2600.0`b, `> 24` bits above the 32-target sample maximum `2575.73`b and `> 36` bits above the live 4-target witness max `2563.01`b). Charging `RET_BOOT = 0.50` (`2x` balls!) requires only that `>= 50%` of uniform rerandomized targets satisfy this conservative `2600`b envelope (supported by `32/32` targets in `[11]` and `4/4` live decile passes in `boot_witness.json`, `kappa_hat = 0.7076 << BALL_PAIR = 2^24`).
- **Explicit Multi-Factor Boot ECM Extraction Loop & Multiplicative Detection (`cm_models.eval_descent_live`):** After trial division (`TDIV`) strips all primes `<= 2^20`, each `<= 2600`-bit boot norm can contain at most `k_lx = floor(2600/20) + 1 = 131` prime factors in `(2^20, LX = 2^70]`, at most `kA = 13` prime factors in window A `(2^70, n_A = 2^234.53125]`, and at most `1` prime factor in window B `(n_A, n_1 = 2^265]`. The extraction loop runs: (1) up to `131` iterations of a `0.999`-detection ECM schedule targeting `<= 2^70` (`ecm_instr(70, 2600, 0.999)`) + 40-round MR (`PRIME`), dividing out each found factor; (2) up to `kA = 13` iterations of a `0.999`-detection ECM schedule targeting `<= 2^234.53125` (`ecm_instr(n_A, 2600, 0.999)`) + 40-round MR (`PRIME`), dividing out each found factor; (3) final 40-round MR verifying the remaining cofactor is `<= 2^265` and prime (or `1`), plus `IDEAL`. Every single boot candidate is conservatively charged the FULL worst-case `131 * (ecm(70) + MR) + 13 * (ecm(n_A) + MR) + MR` work (crediting zero early exit) with exact multiplicative detection probability **`det_b = 0.999^(131 + 13 + 1) = 0.999^145 = 0.8650`**!
- **Certified Boot Parameters (100% LIVE in `audit_cm.py`):** `n_1 = 265`, `n_A = 234.53125`, `k_A = 13`, `R_b = 1` (`17` cells/ball), `m_max = 2600.0`b, `BOOT_RELS = 48`, `RET_BOOT = 0.50`:
  - At default `MESH = 96/192` (strict lower sum charged in ledger): `pr = 2^-37.901262`, `miss = 0.225854`, `psi = 0.686131`, **`phi_boot = psi^48 = 1.4045e-08 <= 2e-08`**, boot own `2^117.4695` + subtree `2^117.0353` = `2^118.2687` instr/log -> **`two_log = 2^104.980946 <= 2^104.981`** core units.
  - At on-demand `MESH = 384/384` (`{"descent_mesh": 384}`): `pr = 2^-37.798554`, `miss = 0.225854`, `psi = 0.679840`, **`phi_boot = 9.0265e-09`**, **`two_log = 2^104.848708`** core units.

## 2. 16D Rank-2 Canonical Special-`q` Nodes & `[7e]` Child Conditioning

Every degree-1 prime ideal `q` above `LX = 2^70` generates a free rank-2 `Z[zeta_15]`-module `L_q = Z[zeta_15] v_1 + Z[zeta_15] v_2` of determinant `q`. Enumerating `v = alpha v_1 + beta v_2` with `alpha` unit-canonical in `B_8(R)` modulo `O_H^* = mu_30 x Z^3` (canonical pass fraction `p_can/30 = 0.008700 = 2^(-6.845)`) and `beta in B_8(R)` eliminates unit duplicates (`UDUP_NODE = 0.95` charged) with base-primitive density `1/zeta_K(2) = 0.9785`.
- Measured across 56 real BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`, 3,584 exact canonical pairs): worst-case normalized offsets `(c_f, c_g) = (+131.66, +46.15)`b (stored in `cm_data.json`: `qbins`, `qbins_sub`, `qbins_lo`).
- **`[7e]` Child-Conditioning Check (BLS12-381 `w' = omega^2`):** Across `22` real children dividing short parent norms vs `24` random split ideals (`q ~ 2^76`), random SV normalized `(+5.58, -13.19)`b vs encountered `(+8.98, -10.84)`b: shift `(+3.40, +2.35)`b, well inside the `(+10, +4)`b tolerance and far below the charged upper tails.
- **Node multi-factor ECM extraction loop:** Identically on every special-`q` node (both `ecm` and `sieve` survivors), after `TDIV` strips primes `<= 2^20`, the loop runs up to `k_lx_f = floor(nf_max/20)+1` and `k_lx_g = floor(ng_max/20)+1` full `ecm(LX, norm_max, 0.999) + MR` iterations for factors `<= LX`, followed by up to `kA` full `ecm(nA, T, 0.999) + MR` iterations per side for window-A factors, plus final cofactor MR (`<= 2^n`). Every survivor is charged the FULL worst-case schedule count and exact multiplicative detection **`det_n = 0.999^(k_lx_f + k_lx_g + 2*(kA + kB))`**.
- **Over-dispersed Negative-Binomial convexity (`r_sh in {0.5, 1.0}` vs Poisson):** Modeling each node's smooth-relation count `K` as Negative-Binomial with shape `r_sh = 0.5` (`CV = sqrt(2)`, single-window nodes) or `r_sh = 1.0` (`CV = 1`, geometric, two-window nodes) is **strictly more heavy-tailed / conservative** than Poisson (`r_sh -> inf`) for every mean `mu > 0` and child-failure probability `psi in (0, 1)`. Indeed, for any `x > 0`, `f(r) = (1 + x/r)^(-r)` satisfies `(d/dr) ln f(r) = -ln(1 + x/r) + x/(r + x) < 0` (since `ln(1 + t) > t/(1 + t)` for `t = x/r > 0`), so `f(r)` is strictly decreasing in `r`; hence BOTH zero-relation failure `P(K = 0) = (1 + mu/r)^(-r) > e^(-mu)` AND the subtree-failure MGF `E[psi^K] = (1 + mu(1 - psi)/r)^(-r) > e^(-mu(1 - psi))` strictly upper-bound the Poisson law, providing an explicit over-dispersed envelope for basis-to-basis relation-count variance.

## 3. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

Every number below is computed 100% LIVE from `plan_struct` by `cm_models.eval_descent_live` inside `audit_cm.py` (columns `own`, `W`, `EN` are `log2` instructions / `log2` expected nodes):

| Bin | S | How | R | n | nA | kA | kB | mu | psi | phi | Nf/Ng | own | W | EN |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| [70,76) | f | sieve | 48 | 70 | 70.0 | 0 | 0 | 336.9 | 0.000 | 0.0395 | 793/800 | 99.95 | 99.95 | 0.00 |
| [70,76) | g | sieve | 48 | 70 | 70.0 | 0 | 0 | 347.0 | 0.000 | 0.0389 | 863/730 | 99.95 | 99.95 | 0.00 |
| [76,81) | f | sieve | 40 | 76 | 76.0 | 1 | 0 | 454.7 | 0.078 | 0.0355 | 791/801 | 95.88 | 101.11 | 1.66 |
| [76,81) | g | sieve | 40 | 76 | 76.0 | 1 | 0 | 468.1 | 0.078 | 0.0350 | 867/725 | 95.88 | 101.11 | 1.66 |
| [81,86) | f | sieve | 48 | 76 | 76.0 | 1 | 0 | 609.8 | 0.078 | 0.0308 | 827/815 | 99.96 | 101.62 | 1.66 |
| [81,86) | g | sieve | 48 | 76 | 76.0 | 1 | 0 | 597.8 | 0.078 | 0.0311 | 908/734 | 99.96 | 101.62 | 1.66 |
| [86,91) | f | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.0 | 0.157 | 0.0431 | 803/809 | 96.09 | 102.22 | 2.52 |
| [86,91) | g | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.0 | 0.157 | 0.0431 | 889/723 | 96.17 | 102.22 | 2.52 |
| [91,96) | f | sieve | 40 | 76 | 76.0 | 2 | 0 | 56.0 | 0.157 | 0.0218 | 843/824 | 97.67 | 102.26 | 2.52 |
| [91,96) | g | sieve | 40 | 76 | 76.0 | 2 | 0 | 52.3 | 0.157 | 0.0232 | 934/733 | 97.64 | 102.26 | 2.52 |
| [96,101) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 73.4 | 0.157 | 0.0200 | 880/838 | 100.31 | 102.54 | 2.52 |
| [96,101) | g | sieve | 48 | 76 | 76.0 | 2 | 0 | 63.9 | 0.157 | 0.0200 | 976/742 | 100.30 | 102.54 | 2.52 |
| [101,106) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 24.0 | 0.165 | 0.0486 | 897/845 | 100.26 | 102.72 | 2.71 |
| [101,106) | g | sieve | 48 | 81 | 75.5 | 1 | 1 | 62.8 | 0.157 | 0.0200 | 998/744 | 100.13 | 103.09 | 3.45 |
| [106,111) | f | sieve | 48 | 81 | 75.5 | 1 | 1 | 24.8 | 0.165 | 0.0471 | 915/853 | 100.12 | 103.22 | 3.56 |
| [106,111) | g | sieve | 48 | 81 | 81.0 | 2 | 0 | 43.8 | 0.157 | 0.0274 | 1021/747 | 100.08 | 103.50 | 4.00 |
| [111,116) | f | sieve | 48 | 81 | 75.5 | 2 | 1 | 28.2 | 0.235 | 0.0454 | 932/860 | 105.64 | 105.92 | 3.86 |
| [111,116) | g | sieve | 48 | 81 | 75.5 | 2 | 1 | 26.6 | 0.244 | 0.0484 | 1043/749 | 105.09 | 105.52 | 3.97 |
| [116,121) | f | sieve | 48 | 86 | 86.0 | 2 | 0 | 28.0 | 0.165 | 0.0421 | 950/868 | 100.03 | 104.05 | 4.09 |
| [116,121) | g | sieve | 48 | 96 | 83.0 | 1 | 1 | 61.4 | 0.165 | 0.0201 | 1066/752 | 100.01 | 104.31 | 4.48 |
| [121,126) | f | sieve | 48 | 96 | 83.0 | 1 | 1 | 28.4 | 0.173 | 0.0418 | 967/875 | 100.01 | 104.38 | 4.55 |
| [121,126) | g | sieve | 48 | 101 | 85.5 | 1 | 1 | 44.4 | 0.165 | 0.0272 | 1088/754 | 100.00 | 104.48 | 4.48 |
| [126,131) | f | sieve | 48 | 96 | 96.0 | 2 | 0 | 29.9 | 0.181 | 0.0403 | 985/883 | 100.03 | 104.66 | 4.90 |
| [126,131) | g | sieve | 48 | 116 | 100.2 | 2 | 0 | 39.1 | 0.172 | 0.0310 | 1111/757 | 100.03 | 104.87 | 4.85 |
| [131,136) | f | sieve | 48 | 96 | 83.0 | 2 | 1 | 44.7 | 0.243 | 0.0297 | 1002/890 | 103.32 | 105.28 | 5.04 |
| [131,136) | g | sieve | 48 | 96 | 83.0 | 2 | 1 | 32.9 | 0.253 | 0.0401 | 1133/759 | 102.53 | 105.17 | 5.10 |
| [136,141) | f | sieve | 48 | 101 | 85.5 | 2 | 1 | 52.3 | 0.243 | 0.0256 | 1020/898 | 103.05 | 105.32 | 5.04 |
| [136,141) | g | sieve | 48 | 101 | 85.5 | 2 | 1 | 35.7 | 0.253 | 0.0371 | 1156/762 | 102.26 | 105.24 | 5.11 |
| [141,146) | f | sieve | 48 | 96 | 96.0 | 3 | 0 | 35.3 | 0.259 | 0.0378 | 1037/905 | 104.43 | 105.91 | 5.57 |
| [141,146) | g | sieve | 56 | 101 | 85.5 | 2 | 1 | 34.1 | 0.253 | 0.0388 | 1194/769 | 104.07 | 105.64 | 5.11 |
| [146,151) | f | sieve | 48 | 106 | 93.6 | 2 | 1 | 39.5 | 0.264 | 0.0343 | 1055/913 | 103.40 | 105.83 | 5.81 |
| [146,151) | g | sieve | 48 | 101 | 101.0 | 3 | 0 | 30.8 | 0.269 | 0.0436 | 1201/767 | 103.69 | 105.95 | 5.62 |
| [151,156) | f | sieve | 56 | 106 | 93.6 | 2 | 1 | 40.6 | 0.264 | 0.0334 | 1088/925 | 104.72 | 106.19 | 5.81 |
| [151,156) | g | sieve | 56 | 106 | 98.1 | 2 | 1 | 32.0 | 0.275 | 0.0423 | 1239/774 | 104.51 | 106.26 | 5.86 |
| [156,161) | f | sieve | 48 | 116 | 100.2 | 2 | 1 | 42.5 | 0.269 | 0.0322 | 1090/928 | 103.69 | 107.57 | 6.21 |
| [156,161) | g | sieve | 56 | 106 | 106.0 | 3 | 0 | 35.3 | 0.286 | 0.0391 | 1262/777 | 105.21 | 106.67 | 6.24 |
| [161,166) | f | sieve | 56 | 116 | 100.2 | 2 | 1 | 43.4 | 0.269 | 0.0316 | 1123/940 | 104.85 | 107.69 | 6.21 |
| [161,166) | g | sieve | 56 | 116 | 105.9 | 2 | 1 | 34.7 | 0.291 | 0.0400 | 1284/779 | 104.69 | 107.78 | 6.56 |
| [166,171) | f | sieve | 56 | 126 | 98.0 | 2 | 1 | 31.9 | 0.269 | 0.0421 | 1141/948 | 104.50 | 107.64 | 6.57 |
| [166,171) | g | sieve | 56 | 131 | 100.5 | 2 | 1 | 38.6 | 0.281 | 0.0357 | 1307/782 | 103.97 | 107.62 | 6.86 |
| [171,176) | f | sieve | 56 | 126 | 106.8 | 2 | 1 | 37.7 | 0.280 | 0.0365 | 1158/955 | 105.29 | 107.93 | 7.15 |
| [171,176) | g | sieve | 56 | 131 | 110.0 | 2 | 1 | 34.9 | 0.280 | 0.0393 | 1329/784 | 104.59 | 107.84 | 7.32 |
| [176,181) | f | sieve | 56 | 141 | 105.5 | 2 | 1 | 52.9 | 0.280 | 0.0266 | 1176/963 | 105.05 | 107.80 | 7.16 |
| [176,181) | g | sieve | 48 | 131 | 100.5 | 3 | 1 | 61.2 | 0.367 | 0.0261 | 1336/782 | 106.80 | 108.39 | 7.23 |
| [181,186) | f | sieve | 56 | 146 | 108.0 | 2 | 1 | 38.3 | 0.280 | 0.0360 | 1193/970 | 105.37 | 107.98 | 7.59 |
| [181,186) | g | sieve | 56 | 131 | 100.5 | 3 | 1 | 58.7 | 0.367 | 0.0272 | 1374/789 | 107.73 | 108.77 | 7.23 |
| [186,191) | f | sieve | 56 | 151 | 110.5 | 2 | 1 | 36.1 | 0.292 | 0.0386 | 1211/978 | 105.50 | 108.12 | 7.84 |
| [186,191) | g | sieve | 48 | 141 | 105.5 | 3 | 1 | 47.0 | 0.384 | 0.0344 | 1381/787 | 106.40 | 108.37 | 7.62 |
| [191,196) | f | sieve | 56 | 161 | 115.5 | 2 | 1 | 68.1 | 0.302 | 0.0216 | 1228/985 | 106.36 | 109.66 | 8.25 |
| [191,196) | g | sieve | 56 | 141 | 105.5 | 3 | 1 | 44.5 | 0.385 | 0.0362 | 1419/794 | 107.30 | 108.66 | 7.62 |
| [196,201) | f | sieve | 56 | 166 | 118.0 | 2 | 1 | 47.1 | 0.302 | 0.0305 | 1246/993 | 106.86 | 110.01 | 8.48 |
| [196,201) | g | sieve | 56 | 151 | 110.5 | 3 | 1 | 88.4 | 0.372 | 0.0200 | 1442/797 | 108.01 | 109.11 | 8.24 |
| [201,206) | f | sieve | 64 | 171 | 120.5 | 2 | 1 | 97.3 | 0.303 | 0.0200 | 1277/1005 | 108.57 | 110.34 | 8.71 |
| [201,206) | g | sieve | 48 | 146 | 119.9 | 3 | 1 | 48.1 | 0.401 | 0.0345 | 1448/794 | 107.66 | 109.85 | 8.39 |
| [206,211) | f | sieve | 56 | 181 | 125.5 | 2 | 1 | 57.1 | 0.303 | 0.0255 | 1281/1008 | 108.82 | 110.64 | 9.15 |
| [206,211) | g | sieve | 56 | 161 | 115.5 | 3 | 1 | 75.0 | 0.401 | 0.0228 | 1487/802 | 107.38 | 110.21 | 8.66 |
| [211,216) | f | sieve | 64 | 181 | 125.5 | 2 | 1 | 44.4 | 0.303 | 0.0323 | 1312/1020 | 109.44 | 110.84 | 9.15 |
| [211,216) | g | sieve | 56 | 166 | 118.0 | 3 | 1 | 57.7 | 0.401 | 0.0291 | 1509/804 | 107.36 | 110.46 | 8.91 |
| [216,221) | f | ecm | 16 | 171 | 136.3 | 5 | 1 | 17520.9 | 0.593 | 0.0200 | 1186/980 | 119.09 | 119.10 | 10.40 |
| [216,221) | g | sieve | 64 | 171 | 120.5 | 3 | 1 | 92.8 | 0.388 | 0.0200 | 1546/812 | 109.11 | 110.79 | 9.07 |
| [221,226) | f | ecm | 16 | 176 | 139.6 | 5 | 1 | 27540.6 | 0.593 | 0.0200 | 1203/987 | 119.93 | 119.93 | 10.62 |
| [221,226) | g | sieve | 56 | 181 | 125.5 | 3 | 1 | 97.7 | 0.402 | 0.0200 | 1554/809 | 108.16 | 110.81 | 9.53 |
| [226,231) | f | ecm | 16 | 171 | 148.9 | 5 | 1 | 78046.1 | 0.593 | 0.0200 | 1221/995 | 121.45 | 121.45 | 10.82 |
| [226,231) | g | sieve | 64 | 181 | 125.5 | 3 | 1 | 70.8 | 0.402 | 0.0241 | 1591/817 | 108.83 | 110.94 | 9.53 |
| [231,236) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 78110.0 | 0.679 | 0.0200 | 1238/1002 | 120.53 | 120.54 | 11.10 |
| [231,236) | g | sieve | 56 | 176 | 139.6 | 3 | 1 | 108.5 | 0.402 | 0.0200 | 1599/814 | 110.07 | 111.23 | 9.72 |
| [236,241) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 43482.4 | 0.690 | 0.0200 | 1256/1010 | 120.53 | 120.54 | 11.16 |
| [236,241) | g | sieve | 40 | 191 | 149.4 | 3 | 1 | 57.0 | 0.402 | 0.0295 | 1587/805 | 108.89 | 111.16 | 10.28 |
| [241,246) | f | ecm | 16 | 171 | 148.9 | 6 | 1 | 119970.4 | 0.679 | 0.0200 | 1273/1017 | 122.21 | 122.21 | 11.36 |
| [241,246) | g | sieve | 48 | 191 | 149.4 | 3 | 1 | 72.1 | 0.402 | 0.0237 | 1628/814 | 109.74 | 111.38 | 10.28 |
| [246,251) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 348110.8 | 0.679 | 0.0200 | 1291/1025 | 122.29 | 122.29 | 11.69 |
| [246,251) | g | sieve | 56 | 191 | 149.4 | 3 | 1 | 67.1 | 0.403 | 0.0253 | 1667/822 | 110.52 | 111.68 | 10.28 |
| [251,256) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 151583.7 | 0.679 | 0.0200 | 1308/1032 | 122.29 | 122.29 | 11.69 |
| [251,256) | g | sieve | 56 | 196 | 152.7 | 3 | 1 | 56.2 | 0.403 | 0.0299 | 1689/824 | 110.71 | 112.14 | 10.47 |
| [256,261) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 87107.0 | 0.691 | 0.0200 | 1326/1040 | 122.29 | 122.29 | 11.74 |
| [256,261) | g | sieve | 64 | 201 | 156.0 | 3 | 1 | 142.1 | 0.403 | 0.0200 | 1726/832 | 112.01 | 112.88 | 10.71 |
| [261,266) | f | ecm | 16 | 186 | 160.6 | 6 | 1 | 367373.4 | 0.679 | 0.0200 | 1343/1047 | 124.43 | 124.43 | 11.93 |
| [261,266) | g | sieve | 64 | 206 | 159.2 | 3 | 1 | 115.8 | 0.403 | 0.0200 | 1748/834 | 112.48 | 113.41 | 10.98 |

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work (charged at stricter `96/192` lower sum):** `2^104.980946` <= **`2^104.981`** core units (`C = 20,000`; `2^104.848708` at `MESH = 384/384`), `-3.883` bits vs `c36` (`108.864`) even with full de Bruijn enabled everywhere, `m_max = 2600`b, `RET_BOOT = 0.50`, and full worst-case multi-factor ECM schedules!
- **Metered cap:** `t = 2^115.6` core units (reached = failure). By Markov's inequality, `P(cap) <= 2^(104.981 - 115.6) = 6.358614e-04` (`6.994476e-04` at `+10%` headroom).
- **Expected node count (`E[N]` computed live in `audit_cm.py`):** `E[N] = 2^16.3392` (`2^16.3022` at `MESH = 384`). Hard node abort `N_MAX = 2^25.8` gives `P(counter abort) <= 2^(16.3392 - 25.8) <= 2^-9.4608`.
- **Repeat-ideal bound (conservative over `tau`-orbits):** The smallest bin `[70, 76)` holds `2 (li(2^76) - li(2^70)) / 12 >= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits across both sides. Bounding all `N_MAX` nodes into `2^66.66` orbits with `2^6` `1/p`-weighting slack gives `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-10.06) <= 2^-10.05`.
