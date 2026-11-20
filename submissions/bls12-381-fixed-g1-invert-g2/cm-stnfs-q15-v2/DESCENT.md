# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, conservative descent main base `MAIN_LB = 2^71.29` (`= min(B_f, B_g)`), `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (below `lr_rethi - 0.009 - 2 SE = 63.4846`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes** (`rho_dB(u, x) = rho(u) - gamma rho(u - 1) / ln x`), `BOOT_RELS = 48`, `RET_BOOT = 0.50` (deterministic 50% norm-certificate discard filter), `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

**100% Live Semismoothness Execution (Zero Precomputed Probabilities or Costs):** `cm_data.json` stores ONLY the 78 structural node 9-tuples `[lo, hi, side, n, iA32, kA, kB, R, thin]`, the empirical `q`-lattice offsets, and the boot deciles (`m_max = 2600.0`b). Inside Docker, `audit_cm.py` invokes `cm_models.eval_descent_live` to compute EVERY node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `psi`, `W`, `phi`, `EN`, and boot `pr_log2`, `miss`, `psi`, `phi_boot`, `two_log`, `E[N]` 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing the ledger at the coarser default `MESH = 96/192` lower sum (`102.907298 <= 102.908`) is strictly more conservative than `MESH = 384/384` (`102.799364`, runnable on demand via `{"descent_mesh": 384}`).

## 1. Fresh-Coin Target Randomization, 16D Guillevic Boot (`det = p^4`), & 50% Norm-Certificate Filter

For each of the two discrete logarithms (`log g` with `g = e(P, Q0)` and `log z`), draw an independent uniform exponent `e in {1, ..., r - 1}` (fresh algorithm coins; cost `< 2^40` instructions inside the wrapper cap) and boot on `z' = z^e` (resp. `g' = g^e`). Because `G_T` has prime order `r`, BOTH `g'` and `z'` are **exact uniform on `G_T \ {1}`** from fresh coins alone, even when `g = e(P, Q0)` is fixed.

### Algebraic Non-Degeneracy Theorem (`100%` of `G_T \ {1}`)

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)` (since `tau(y, X) = (y^7, omega^2 X)` satisfies `tau^6(y, X) = (y^4, X)`):
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r | Phi_12(p)`, so every nonzero element of `F_(p^4)^* * F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`) and `c_6 = c_(6,0) + c_(6,1) X + c_(6,2) X^2` (`c_(6,i) in F_(p^2)`). Vanishing of the `X^2`-coefficient of `c_6 z'`, namely `c_(6,0) z_2 + c_(6,1) z_1 + c_(6,2) z_0 = 0` in `F_(p^4) ~= F_(p^2)^2`, is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore **always** has a nonzero solution `c_6 in F_(p^6)^*`.
- Writing `c_6 z' = w_0 + w_1 X`: if `w_1 = 0` or `w_0 = 0`, then `c_6 z' in F_(p^4)^*` or `F_(p^4)^* X`, forcing `z' in F_(p^4)^* F_(p^6)^*`, impossible for `z' in G_T \ {1}` of prime order `r`. Hence `gamma = w_0 / w_1 in F_(p^4) \ F_(p^2)` is defined for **100% of nonidentity targets**.
- Define the 16-dimensional linear lattice `L_Guillevic = { (a, b) in Z[zeta_15]^2 : a - gamma b == 0 (mod P) } subset Z^16`. Because `O_H / P ~= F_(p^4)` imposes 4 `F_p`-linear congruences, **`det(L_Guillevic) = p^4`** (`~ 2^1524`), giving shortest-vector `g`-norm `N_g(a + b X) ~ p^(20/3) ~= 2538` bits. For every short vector in the boot ball, `|N_(K_H/Q)(b)| ~= 2^680 << p^4 ~= 2^1524`, so `b != 0 (mod P)` is guaranteed and `Psi_g(a + b X) = (b mod P) w_1^(-1) c_6 z'` has exact discrete log `e log_r(z) mod r`.

### Free Rank-2 `Z[zeta_15]`-Module Enumeration, Live 4-Target Witness, & Deterministic 50% Norm Filter

Because `K_H = Q(zeta_15)` is a PID (`h = 1`), `L_Guillevic = Z[zeta_15] v_1 + Z[zeta_15] v_2` is a free rank-2 `Z[zeta_15]`-module. Under LLL with the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)`, enumerating `v = 1 * v_1 + beta * v_2` for `beta in B_8(1)` (`17` vectors: `0` and `+-e_i`) yields **zero unit duplicates (`UDUP_BOOT = 1.000`)** and **100% base-primitivity (`gcd(1, beta) = 1`)**.

- **Live 4-target BLS12-381 `w'=omega^2` witness (`boot_witness.json`, verified LIVE in `audit_cm.py`):** Four real reduced BLS12-381 `w' = omega^2` rank-2 boot bases `(gam, v_1, v_2)` are bundled inside the package. `audit_cm.py` (`voronoi_and_boot_live`) verifies LIVE on all 4 bases: exact `F_(p^4)` congruence `a_i - gamma b_i == 0 (mod P)`, exact determinant `|N_(K_H/Q)(a_1 b_2 - a_2 b_1)| == p^4`, and all `4 x 17 = 68` exact 2,550-bit `g`-norms (mean `2551.75`b, max `2563.01`b `< 2600`b, pairwise `gcd > 10^4` rate `3/544`).
- **Deterministic 50% boot-ball norm-certificate filter (`RET_BOOT = 0.50`, `m_max = 2600.0`b):** To eliminate any tail sensitivity beyond a finite empirical sample, after LLL on each fresh boot target the algorithm computes all 17 exact norms in `B_8(1)` and **deterministically discards** any boot ball whose max norm exceeds `m_max = 2600.0`b or whose sorted deciles exceed `boot_sub_deciles` (top decile `2600.0`b, `> 24` bits above the 32-target sample maximum `2575.73`b and `> 36` bits above the live 4-target witness max `2563.01`b). Charging `RET_BOOT = 0.50` (`2x` balls!) requires only that `>= 50%` of uniform rerandomized targets satisfy this conservative `2600`b envelope (supported by `32/32` targets in `[11]`, SV mean `2540.82`b, ball mean `2552.54`b, ball max `2575.73`b, `kappa_hat = 0.7076 << BALL_PAIR = 2^24`).
- **Certified Boot Parameters (100% LIVE in `audit_cm.py`):** `n_1 = 265`, `n_A = 234.53125`, `k_A = 13`, `R_b = 1` (`17` cells/ball), `m_max = 2600.0`b, `BOOT_RELS = 48`, `RET_BOOT = 0.50`:
  - At default `MESH = 96/192` (strict lower sum charged in ledger): `pr = 2^-37.898475`, `miss = 0.225854`, `psi = 0.681829`, **`phi_boot = psi^48 = 1.0385e-08 <= 2e-08`**, boot own `2^113.7617` + subtree `2^115.8996` = `2^116.1950` instr/log -> **`two_log = 2^102.907298 <= 2^102.908`** core units.
  - At on-demand `MESH = 384/384` (`{"descent_mesh": 384}`): `pr = 2^-37.795803`, `miss = 0.225854`, `psi = 0.675557`, **`phi_boot = 6.6644e-09`**, **`two_log = 2^102.799364`** core units.

## 2. 16D Rank-2 Canonical Special-`q` Nodes & `[7e]` Child Conditioning

Every degree-1 prime ideal `q` above `LX = 2^70` generates a free rank-2 `Z[zeta_15]`-module `L_q = Z[zeta_15] v_1 + Z[zeta_15] v_2` of determinant `q`. Enumerating `v = alpha v_1 + beta v_2` with `alpha` unit-canonical in `B_8(R)` modulo `O_H^* = mu_30 x Z^3` (canonical pass fraction `p_can/30 = 0.008700 = 2^(-6.845)`) and `beta in B_8(R)` eliminates unit duplicates (`UDUP_NODE = 0.95` charged) with base-primitive density `1/zeta_K(2) = 0.9785`.
- Measured across 56 real BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`, 3,584 exact canonical pairs): worst-case normalized offsets `(c_f, c_g) = (+131.66, +46.15)`b (stored in `cm_data.json`: `qbins`, `qbins_sub`, `qbins_lo`).
- **`[7e]` Child-Conditioning Check (BLS12-381 `w' = omega^2`):** Across `22` real children dividing short parent norms vs `24` random split ideals (`q ~ 2^76`), random SV normalized `(+5.58, -13.19)`b vs encountered `(+8.98, -10.84)`b: shift `(+3.40, +2.35)`b, well inside the `(+10, +4)`b tolerance and far below the charged upper tails.

## 3. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

Every number below is computed 100% LIVE from `plan_struct` by `cm_models.eval_descent_live` inside `audit_cm.py` (columns `own`, `W`, `EN` are `log2` instructions / `log2` expected nodes):

| Bin | S | How | R | n | nA | kA | kB | mu | psi | phi | Nf/Ng | own | W | EN |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| [70,76) | f | sieve | 48 | 70 | 70.0 | 0 | 0 | 336.3 | 0.000 | 0.0395 | 793/800 | 99.95 | 99.95 | 0.00 |
| [70,76) | g | sieve | 48 | 70 | 70.0 | 0 | 0 | 346.3 | 0.000 | 0.0390 | 863/730 | 99.95 | 99.95 | 0.00 |
| [76,81) | f | sieve | 40 | 76 | 76.0 | 1 | 0 | 458.2 | 0.079 | 0.0354 | 791/801 | 95.82 | 101.11 | 1.66 |
| [76,81) | g | sieve | 40 | 76 | 76.0 | 1 | 0 | 471.7 | 0.079 | 0.0349 | 867/725 | 95.82 | 101.11 | 1.66 |
| [81,86) | f | sieve | 48 | 76 | 76.0 | 1 | 0 | 615.6 | 0.079 | 0.0307 | 827/815 | 99.95 | 101.62 | 1.66 |
| [81,86) | g | sieve | 48 | 76 | 76.0 | 1 | 0 | 602.0 | 0.079 | 0.0310 | 908/734 | 99.95 | 101.62 | 1.66 |
| [86,91) | f | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.2 | 0.157 | 0.0427 | 803/809 | 92.69 | 102.20 | 2.52 |
| [86,91) | g | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.2 | 0.157 | 0.0428 | 889/723 | 92.71 | 102.20 | 2.52 |
| [91,96) | f | sieve | 40 | 76 | 76.0 | 2 | 0 | 56.4 | 0.157 | 0.0216 | 843/824 | 95.93 | 102.22 | 2.52 |
| [91,96) | g | sieve | 40 | 76 | 76.0 | 2 | 0 | 52.6 | 0.157 | 0.0230 | 934/733 | 95.92 | 102.22 | 2.52 |
| [96,101) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 73.9 | 0.157 | 0.0200 | 880/838 | 99.96 | 102.48 | 2.52 |
| [96,101) | g | sieve | 48 | 76 | 76.0 | 2 | 0 | 64.4 | 0.157 | 0.0200 | 976/742 | 99.96 | 102.47 | 2.52 |
| [101,106) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 24.1 | 0.165 | 0.0483 | 897/845 | 99.96 | 102.67 | 2.71 |
| [101,106) | g | sieve | 48 | 81 | 75.5 | 1 | 1 | 63.2 | 0.157 | 0.0200 | 998/744 | 99.96 | 103.07 | 3.45 |
| [106,111) | f | sieve | 48 | 81 | 75.5 | 1 | 1 | 24.9 | 0.165 | 0.0468 | 915/853 | 99.96 | 103.20 | 3.56 |
| [106,111) | g | sieve | 48 | 81 | 81.0 | 2 | 0 | 44.1 | 0.157 | 0.0272 | 1021/747 | 99.96 | 103.48 | 4.00 |
| [111,116) | f | sieve | 48 | 81 | 75.5 | 2 | 1 | 28.4 | 0.236 | 0.0451 | 932/860 | 101.28 | 103.71 | 3.86 |
| [111,116) | g | sieve | 48 | 81 | 75.5 | 2 | 1 | 26.8 | 0.244 | 0.0481 | 1043/749 | 100.92 | 103.77 | 3.97 |
| [116,121) | f | sieve | 48 | 86 | 86.0 | 2 | 0 | 28.2 | 0.165 | 0.0418 | 950/868 | 99.95 | 104.04 | 4.09 |
| [116,121) | g | sieve | 48 | 96 | 83.0 | 1 | 1 | 61.8 | 0.164 | 0.0200 | 1066/752 | 99.95 | 104.28 | 4.48 |
| [121,126) | f | sieve | 48 | 96 | 83.0 | 1 | 1 | 28.7 | 0.173 | 0.0415 | 967/875 | 99.95 | 104.35 | 4.55 |
| [121,126) | g | sieve | 48 | 101 | 85.5 | 1 | 1 | 44.7 | 0.164 | 0.0270 | 1088/754 | 99.95 | 104.43 | 4.48 |
| [126,131) | f | sieve | 48 | 96 | 96.0 | 2 | 0 | 30.1 | 0.180 | 0.0399 | 985/883 | 99.96 | 104.61 | 4.90 |
| [126,131) | g | sieve | 48 | 116 | 100.2 | 2 | 0 | 39.4 | 0.171 | 0.0307 | 1111/757 | 99.96 | 104.80 | 4.84 |
| [131,136) | f | sieve | 48 | 96 | 83.0 | 2 | 1 | 45.0 | 0.243 | 0.0295 | 1002/890 | 100.43 | 104.90 | 5.04 |
| [131,136) | g | sieve | 48 | 96 | 83.0 | 2 | 1 | 33.1 | 0.252 | 0.0398 | 1133/759 | 100.21 | 104.95 | 5.10 |
| [136,141) | f | sieve | 48 | 101 | 85.5 | 2 | 1 | 52.6 | 0.243 | 0.0255 | 1020/898 | 100.41 | 105.01 | 5.04 |
| [136,141) | g | sieve | 48 | 101 | 85.5 | 2 | 1 | 35.9 | 0.253 | 0.0369 | 1156/762 | 100.18 | 105.06 | 5.10 |
| [141,146) | f | sieve | 48 | 96 | 96.0 | 3 | 0 | 35.6 | 0.257 | 0.0375 | 1037/905 | 101.64 | 105.34 | 5.57 |
| [141,146) | g | sieve | 56 | 101 | 85.5 | 2 | 1 | 34.3 | 0.253 | 0.0385 | 1194/769 | 103.53 | 105.46 | 5.11 |
| [146,151) | f | sieve | 48 | 106 | 93.6 | 2 | 1 | 39.7 | 0.262 | 0.0340 | 1055/913 | 101.06 | 105.56 | 5.80 |
| [146,151) | g | sieve | 48 | 101 | 101.0 | 3 | 0 | 31.0 | 0.267 | 0.0432 | 1201/767 | 101.60 | 105.63 | 5.62 |
| [151,156) | f | sieve | 56 | 106 | 93.6 | 2 | 1 | 40.9 | 0.262 | 0.0331 | 1088/925 | 103.70 | 105.86 | 5.80 |
| [151,156) | g | sieve | 56 | 106 | 98.1 | 2 | 1 | 32.2 | 0.273 | 0.0420 | 1239/774 | 103.71 | 106.02 | 5.85 |
| [156,161) | f | sieve | 48 | 116 | 100.2 | 2 | 1 | 42.8 | 0.267 | 0.0319 | 1090/928 | 101.72 | 106.13 | 6.21 |
| [156,161) | g | sieve | 56 | 106 | 106.0 | 3 | 0 | 35.5 | 0.284 | 0.0388 | 1262/777 | 104.07 | 106.31 | 6.24 |
| [161,166) | f | sieve | 56 | 116 | 100.2 | 2 | 1 | 43.7 | 0.267 | 0.0313 | 1123/940 | 103.88 | 106.35 | 6.21 |
| [161,166) | g | sieve | 56 | 116 | 105.9 | 2 | 1 | 34.9 | 0.290 | 0.0397 | 1284/779 | 103.91 | 106.57 | 6.56 |
| [166,171) | f | sieve | 56 | 126 | 98.0 | 2 | 1 | 32.1 | 0.267 | 0.0418 | 1141/948 | 103.71 | 106.66 | 6.57 |
| [166,171) | g | sieve | 56 | 131 | 100.5 | 2 | 1 | 38.9 | 0.279 | 0.0354 | 1307/782 | 103.60 | 106.87 | 6.85 |
| [171,176) | f | sieve | 56 | 126 | 106.8 | 2 | 1 | 37.9 | 0.279 | 0.0363 | 1158/955 | 104.26 | 107.07 | 7.14 |
| [171,176) | g | sieve | 56 | 131 | 110.0 | 2 | 1 | 35.1 | 0.279 | 0.0390 | 1329/784 | 103.95 | 107.17 | 7.32 |
| [176,181) | f | sieve | 56 | 141 | 105.5 | 2 | 1 | 53.2 | 0.279 | 0.0264 | 1176/963 | 104.07 | 107.22 | 7.16 |
| [176,181) | g | sieve | 48 | 131 | 100.5 | 3 | 1 | 61.5 | 0.365 | 0.0259 | 1336/782 | 104.41 | 107.30 | 7.22 |
| [181,186) | f | sieve | 56 | 146 | 108.0 | 2 | 1 | 38.5 | 0.279 | 0.0358 | 1193/970 | 104.39 | 107.60 | 7.58 |
| [181,186) | g | sieve | 56 | 131 | 100.5 | 3 | 1 | 58.9 | 0.365 | 0.0270 | 1374/789 | 105.54 | 107.52 | 7.22 |
| [186,191) | f | sieve | 56 | 151 | 110.5 | 2 | 1 | 36.3 | 0.290 | 0.0383 | 1211/978 | 104.54 | 107.77 | 7.83 |
| [186,191) | g | sieve | 48 | 141 | 105.5 | 3 | 1 | 47.2 | 0.382 | 0.0341 | 1381/787 | 104.30 | 107.63 | 7.62 |
| [191,196) | f | sieve | 56 | 161 | 115.5 | 2 | 1 | 68.5 | 0.300 | 0.0215 | 1228/985 | 105.32 | 108.37 | 8.24 |
| [191,196) | g | sieve | 56 | 141 | 105.5 | 3 | 1 | 44.7 | 0.382 | 0.0360 | 1419/794 | 105.44 | 107.80 | 7.62 |
| [196,201) | f | sieve | 56 | 166 | 118.0 | 2 | 1 | 47.3 | 0.301 | 0.0303 | 1246/993 | 105.81 | 108.69 | 8.47 |
| [196,201) | g | sieve | 56 | 151 | 110.5 | 3 | 1 | 88.7 | 0.370 | 0.0200 | 1442/797 | 106.28 | 108.36 | 8.23 |
| [201,206) | f | sieve | 64 | 171 | 120.5 | 2 | 1 | 97.7 | 0.301 | 0.0200 | 1277/1005 | 107.80 | 109.31 | 8.70 |
| [201,206) | g | sieve | 48 | 146 | 119.9 | 3 | 1 | 48.2 | 0.398 | 0.0343 | 1448/794 | 106.00 | 108.59 | 8.39 |
| [206,211) | f | sieve | 56 | 181 | 125.5 | 2 | 1 | 57.3 | 0.301 | 0.0253 | 1281/1008 | 107.78 | 109.61 | 9.14 |
| [206,211) | g | sieve | 56 | 161 | 115.5 | 3 | 1 | 75.3 | 0.398 | 0.0226 | 1487/802 | 105.83 | 108.78 | 8.66 |
| [211,216) | f | sieve | 64 | 181 | 125.5 | 2 | 1 | 44.6 | 0.301 | 0.0321 | 1312/1020 | 108.56 | 109.87 | 9.14 |
| [211,216) | g | sieve | 56 | 166 | 118.0 | 3 | 1 | 57.9 | 0.399 | 0.0289 | 1509/804 | 105.85 | 109.06 | 8.90 |
| [216,221) | f | ecm | 16 | 171 | 136.3 | 5 | 1 | 17531.9 | 0.589 | 0.0200 | 1186/980 | 116.77 | 116.79 | 10.38 |
| [216,221) | g | sieve | 64 | 171 | 120.5 | 3 | 1 | 93.0 | 0.385 | 0.0200 | 1546/812 | 107.89 | 109.58 | 9.06 |
| [221,226) | f | ecm | 16 | 176 | 139.6 | 5 | 1 | 27571.6 | 0.589 | 0.0200 | 1203/987 | 117.61 | 117.62 | 10.60 |
| [221,226) | g | sieve | 56 | 181 | 125.5 | 3 | 1 | 97.9 | 0.399 | 0.0200 | 1554/809 | 106.63 | 109.68 | 9.52 |
| [226,231) | f | ecm | 16 | 171 | 148.9 | 5 | 1 | 78075.8 | 0.589 | 0.0200 | 1221/995 | 119.13 | 119.13 | 10.81 |
| [226,231) | g | sieve | 64 | 181 | 125.5 | 3 | 1 | 70.9 | 0.399 | 0.0239 | 1591/817 | 107.72 | 109.86 | 9.52 |
| [231,236) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 78032.9 | 0.675 | 0.0200 | 1238/1002 | 117.95 | 117.96 | 11.08 |
| [231,236) | g | sieve | 56 | 176 | 139.6 | 3 | 1 | 108.7 | 0.400 | 0.0200 | 1599/814 | 108.50 | 110.18 | 9.71 |
| [236,241) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 43426.8 | 0.686 | 0.0200 | 1256/1010 | 117.95 | 117.96 | 11.14 |
| [236,241) | g | sieve | 40 | 191 | 149.4 | 3 | 1 | 57.2 | 0.399 | 0.0293 | 1587/805 | 107.30 | 110.38 | 10.27 |
| [241,246) | f | ecm | 16 | 171 | 148.9 | 6 | 1 | 119715.4 | 0.675 | 0.0200 | 1273/1017 | 119.62 | 119.63 | 11.33 |
| [241,246) | g | sieve | 48 | 191 | 149.4 | 3 | 1 | 72.2 | 0.400 | 0.0236 | 1628/814 | 108.15 | 110.52 | 10.27 |
| [246,251) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 347545.4 | 0.675 | 0.0200 | 1291/1025 | 119.71 | 119.71 | 11.66 |
| [246,251) | g | sieve | 56 | 191 | 149.4 | 3 | 1 | 67.2 | 0.400 | 0.0252 | 1667/822 | 108.95 | 110.71 | 10.27 |
| [251,256) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 151291.4 | 0.675 | 0.0200 | 1308/1032 | 119.71 | 119.71 | 11.66 |
| [251,256) | g | sieve | 56 | 196 | 152.7 | 3 | 1 | 56.3 | 0.400 | 0.0298 | 1689/824 | 109.14 | 111.04 | 10.46 |
| [256,261) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 86912.5 | 0.687 | 0.0200 | 1326/1040 | 119.71 | 119.71 | 11.72 |
| [256,261) | g | sieve | 64 | 201 | 156.0 | 3 | 1 | 142.2 | 0.401 | 0.0200 | 1726/832 | 110.48 | 111.68 | 10.70 |
| [261,266) | f | ecm | 16 | 186 | 160.6 | 6 | 1 | 366571.7 | 0.675 | 0.0200 | 1343/1047 | 121.84 | 121.84 | 11.90 |
| [261,266) | g | sieve | 64 | 206 | 159.2 | 3 | 1 | 115.8 | 0.401 | 0.0200 | 1748/834 | 110.94 | 112.11 | 10.97 |

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work (charged at stricter `96/192` lower sum):** `2^102.907298` <= **`2^102.908`** core units (`C = 20,000`; `2^102.799364` at `MESH = 384/384`), `-5.956` bits below `c36` (`108.864`) even with full de Bruijn enabled everywhere, `m_max = 2600`b, and `RET_BOOT = 0.50`.
- **Metered cap:** `t = 2^115.6` core units (reached = failure). By Markov's inequality, `P(cap) <= 2^(102.908 - 115.6) = 1.511219e-04` (`1.662340e-04` at `+10%` headroom).
- **Expected node count (`E[N]` computed live in `audit_cm.py`):** `E[N] = 2^16.3104` (`2^16.2740` at `MESH = 384`). Hard node abort `N_MAX = 2^25.8` gives `P(counter abort) <= 2^(16.3104 - 25.8) <= 2^-9.4896`.
- **Repeat-ideal bound (conservative over `tau`-orbits):** The smallest bin `[70, 76)` holds `2 (li(2^76) - li(2^70)) / 12 >= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits across both sides. Bounding all `N_MAX` nodes into `2^66.66` orbits with `2^6` `1/p`-weighting slack gives `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-10.06) <= 2^-10.05`.
