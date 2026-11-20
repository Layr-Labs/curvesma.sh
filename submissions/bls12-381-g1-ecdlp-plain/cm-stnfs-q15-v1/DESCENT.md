# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, conservative descent main base `MAIN_LB = 2^71.29` (`= min(B_f, B_g)`), `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (below `log2_rels - 0.009 - 2 SE = 63.4881`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes** (`rho_dB(u, x) = rho(u) - gamma rho(u - 1) / ln x`), `MESH = 384`, `BOOT_RELS = 36`, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

## 1. Fresh-Coin Target Randomization & 16D Guillevic Rank-2 Torus Boot (`det = p^4`)

For each of the two discrete logarithms (`log g` with `g = e(P0, Q0)` and `log z`), draw an independent uniform exponent `e in {1, ..., r - 1}` (fresh algorithm coins; cost `< 2^40` instructions inside the wrapper cap) and boot on `z' = z^e` (resp. `g' = g^e`). Because `G_T` has prime order `r`, BOTH `g'` and `z'` are **exact uniform on `G_T \ {1}`** from fresh coins alone, even when `g = e(P0, Q0)` is fixed.

### Algebraic Non-Degeneracy Theorem (`100%` of `G_T \ {1}`)

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)` (since `tau(y, X) = (y^7, omega^2 X)` satisfies `tau^6(y, X) = (y^4, X)`):
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r | Phi_12(p)`, so every nonzero element of `F_(p^4)^* * F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`) and `c_6 = c_(6,0) + c_(6,1) X + c_(6,2) X^2` (`c_(6,i) in F_(p^2)`). Vanishing of the `X^2`-coefficient of `c_6 z'`, namely `c_(6,0) z_2 + c_(6,1) z_1 + c_(6,2) z_0 = 0` in `F_(p^4) ~= F_(p^2)^2`, is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore **always** has a nonzero solution `c_6 in F_(p^6)^*`.
- Writing `c_6 z' = w_0 + w_1 X`: if `w_1 = 0` or `w_0 = 0`, then `c_6 z' in F_(p^4)^*` or `F_(p^4)^* X`, forcing `z' in F_(p^4)^* F_(p^6)^*`, impossible for `z' in G_T \ {1}` of prime order `r`. Hence `gamma = w_0 / w_1 in F_(p^4) \ F_(p^2)` is defined for **100% of nonidentity targets**.
- Define the 16-dimensional linear lattice
  `L_Guillevic = { (a, b) in Z[zeta_15]^2 : a - gamma b == 0 (mod P) } subset Z^16`.
  Because `O_H / P ~= F_(p^4)` imposes exactly 4 `F_p`-linear congruences, **`det(L_Guillevic) = p^4`** (`~ 2^1524`), giving shortest-vector `g`-norm `N_g(a + b X) ~ p^(20/3) ~= 2538` bits (vs `2850` bits for a naive 24D cubic boot). For every short vector in the boot ball, `|N_(K_H/Q)(b)| ~= 2^680 << p^4 ~= 2^1524`, so `b != 0 (mod P)` is provably guaranteed, and `Psi_g(a + b X) = (b mod P) w_1^(-1) c_6 z'` has exact discrete log `e log_r(z) mod r`.

### Free Rank-2 `Z[zeta_15]`-Module Enumeration (`UDUP_BOOT = 1.0`, `Prim = 1.0`)

Because `K_H = Q(zeta_15)` has class number `h(K_H) = 1` (PID) and `L_Guillevic` is a `Z[zeta_15]`-submodule of `Z[zeta_15]^2` of rank 2, `L_Guillevic = Z[zeta_15] v_1 + Z[zeta_15] v_2` is a **free rank-2 `Z[zeta_15]`-module**. Under LLL with the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)`, vectors `v_1, v_2` with `|det_Z(span(v_1, v_2))| = p^4` are extracted cleanly on 100% of tested BLS12-381 targets (`16/16` in `verify_rank2_module.gp` and `32/32` fresh `w'=omega^2` targets in `[11]`).

Let `B_8(R)` denote the 8D lattice ball `{ x in Z^8 : sum_i x_i^2 <= R^2 }` in the LLL-reduced basis of `O_H` (`B_8(1)` has exactly `17` vectors: `0` and `+-e_i`, `i = 1..8`). Enumerating `v = 1 * v_1 + beta * v_2` with `beta in B_8(R_b)` (`R_b = 1`, `17` cells/ball):
- **Zero unit duplicates (`UDUP_BOOT = 1.000`):** any unit relation `eps (v_1 + beta v_2) = v_1 + beta' v_2` forces `eps = 1` and `beta' = beta` by basis uniqueness — completely eliminating the `15x` `mu_30`-duplicate trap of raw 16D balls.
- **100% base-primitivity (`gcd_(Z[zeta_15])(1, beta) = (1)`).**
- **Boot renewal & `[11]` live `BALL_PAIR` check:** Across `32` fresh BLS12-381 `G_T` targets (`w' = omega^2`, `17` cells/ball = `544` exact norms): SV mean `2540.82`b (max `2553.16`b), `R_b = 1` ball mean `2552.54`b (max `2575.73`b, strictly below the charged top boot decile **`m_max = 2580.00`b**!), pairwise `gcd > 10^4` collisions `3/4352` (`0.069%`), live pairwise smoothness ratio `kappa_hat = 0.7076` against charged `BALL_PAIR = 2^24` (`> 2^24.4 x` headroom).
- **Certified Boot Parameters (`MESH = 384`, de Bruijn ON, verified live in `audit_cm.py`):** `n_1 = 260`, `n_A = 230.3125`, `k_A = 13` (window A `(LX, n_A]` permits up to `k_A` primes; window B `(n_A, n_1]` permits `1` prime), `R_b = 1` (`17` cells/ball), `m_max = 2580.00`b, `miss = 0.224003`, `psi = 0.675693` (`retry = 3.0835`), `BOOT_RELS = 36` => **`phi_boot = psi^36 = 7.4299e-07 <= 7.5e-07`**. Boot own work `2^113.2128` instructions + retried subtree `2^115.5195` instructions = `2^115.7851` instructions per log.

## 2. 16D Rank-2 Canonical Special-`q` Nodes & `[7e]` Child Conditioning

Every degree-1 prime ideal `q` above `LX = 2^70` generates a free rank-2 `Z[zeta_15]`-module `L_q = Z[zeta_15] v_1 + Z[zeta_15] v_2` of determinant `q`. Enumerating `v = alpha * v_1 + beta * v_2` with `alpha` unit-canonical in `B_8(R)` (modulo `O_H^* = mu_30 x Z^3`, measured canonical pass fraction `p_can/30 = 0.008700 = 2^(-6.845)`) and `beta in B_8(R)`:
- Eliminates all `mu_30` and free-unit duplicates (`UDUP = 1.000`, conservatively charged `UDUP_NODE = 0.95`), with base-primitive density `1 / zeta_K(2) = 0.9785`.
- At each bin `[lo, hi)` and side `S in {f, g}`, a smooth child relation has all prime factors `<= 2^n`, with at most `k_A` primes in window A `(LX, n_A]` and at most `k_B in {0, 1}` primes in window B `(n_A, n]`.
- Measured across 56 real BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`, 3,584 exact canonical pairs): normalized worst-case shifted offsets `(c_f, c_g) = (+131.66, +46.15)` bits (tighter than raw 16D balls `(+164.19, +51.57)`).
- **`[7e]` Child-Conditioning Check (BLS12-381 `w' = omega^2`):** Comparing `22` real children dividing short parent norms against `24` fresh random split ideals (`q ~ 2^76`): random SV normalized `(c_f, c_g) = (+5.58, -13.19)`b vs encountered children `(+8.98, -10.84)`b, shift `delta = (+3.40, +2.35)`b — well inside `c36`'s `(+10, +4)` tolerance and far below the charged upper tails.

## 3. Complete Bottom-Up Certified Plan (All `39` Bins / `78` Bin-Side Nodes, `MESH = 384`, Full de Bruijn ON)

Every row below is stored in `cm_data.json` and executed LIVE bottom-up by `audit_cm.py` (recomputing `psi`, `lazy_fraction`, `own`, `W`, `phi`, `EN`, and the final boot derivation from scratch; columns `own`, `W`, `EN` are `log2` instructions / `log2` expected nodes):

| Bin | S | How | R | n | nA | kA | kB | mu | psi | phi | Nf/Ng | own | W | EN |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| [70,76) | f | sieve | 48 | 70 | 70.0 | 0 | 0 | 335.6 | 0.000 | 0.0396 | 793/800 | 99.95 | 99.95 | 0.00 |
| [70,76) | g | sieve | 48 | 70 | 70.0 | 0 | 0 | 345.6 | 0.000 | 0.0390 | 863/730 | 99.95 | 99.95 | 0.00 |
| [76,81) | f | sieve | 40 | 76 | 76.0 | 1 | 0 | 459.1 | 0.079 | 0.0354 | 791/801 | 95.82 | 101.11 | 1.66 |
| [76,81) | g | sieve | 40 | 76 | 76.0 | 1 | 0 | 472.6 | 0.079 | 0.0349 | 867/725 | 95.82 | 101.11 | 1.66 |
| [81,86) | f | sieve | 48 | 76 | 76.0 | 1 | 0 | 616.8 | 0.079 | 0.0306 | 827/815 | 99.95 | 101.62 | 1.66 |
| [81,86) | g | sieve | 48 | 76 | 76.0 | 1 | 0 | 603.2 | 0.079 | 0.0310 | 908/734 | 99.95 | 101.62 | 1.66 |
| [86,91) | f | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.3 | 0.157 | 0.0426 | 803/809 | 92.71 | 102.20 | 2.52 |
| [86,91) | g | sieve | 32 | 76 | 76.0 | 2 | 0 | 27.3 | 0.157 | 0.0427 | 889/723 | 92.73 | 102.20 | 2.52 |
| [91,96) | f | sieve | 40 | 76 | 76.0 | 2 | 0 | 56.6 | 0.157 | 0.0215 | 843/824 | 95.93 | 102.22 | 2.52 |
| [91,96) | g | sieve | 40 | 76 | 76.0 | 2 | 0 | 52.8 | 0.157 | 0.0230 | 934/733 | 95.92 | 102.22 | 2.52 |
| [96,101) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 74.2 | 0.157 | 0.0200 | 880/838 | 99.96 | 102.48 | 2.52 |
| [96,101) | g | sieve | 48 | 76 | 76.0 | 2 | 0 | 64.6 | 0.157 | 0.0200 | 976/742 | 99.96 | 102.48 | 2.52 |
| [101,106) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 24.2 | 0.165 | 0.0481 | 897/845 | 99.96 | 102.67 | 2.71 |
| [101,106) | g | sieve | 48 | 81 | 75.5 | 1 | 1 | 63.7 | 0.157 | 0.0200 | 998/744 | 99.96 | 103.07 | 3.45 |
| [106,111) | f | sieve | 48 | 81 | 75.5 | 1 | 1 | 25.1 | 0.165 | 0.0465 | 915/853 | 99.96 | 103.20 | 3.56 |
| [106,111) | g | sieve | 48 | 81 | 81.0 | 2 | 0 | 44.5 | 0.157 | 0.0270 | 1021/747 | 99.96 | 103.48 | 4.00 |
| [111,116) | f | sieve | 48 | 81 | 75.5 | 2 | 1 | 28.6 | 0.236 | 0.0447 | 932/860 | 101.30 | 103.72 | 3.86 |
| [111,116) | g | sieve | 48 | 81 | 75.5 | 2 | 1 | 27.1 | 0.245 | 0.0476 | 1043/749 | 100.94 | 103.77 | 3.97 |
| [116,121) | f | sieve | 48 | 86 | 86.0 | 2 | 0 | 28.5 | 0.166 | 0.0413 | 950/868 | 99.95 | 104.04 | 4.09 |
| [116,121) | g | sieve | 48 | 96 | 83.0 | 1 | 1 | 63.0 | 0.164 | 0.0200 | 1066/752 | 99.95 | 104.28 | 4.48 |
| [121,126) | f | sieve | 48 | 96 | 83.0 | 1 | 1 | 29.2 | 0.173 | 0.0407 | 967/875 | 99.95 | 104.35 | 4.55 |
| [121,126) | g | sieve | 48 | 101 | 85.5 | 1 | 1 | 45.8 | 0.164 | 0.0264 | 1088/754 | 99.95 | 104.43 | 4.48 |
| [126,131) | f | sieve | 48 | 96 | 96.0 | 2 | 0 | 30.7 | 0.179 | 0.0391 | 985/883 | 99.96 | 104.61 | 4.90 |
| [126,131) | g | sieve | 48 | 116 | 100.2 | 2 | 0 | 40.8 | 0.171 | 0.0297 | 1111/757 | 99.96 | 104.80 | 4.84 |
| [131,136) | f | sieve | 48 | 96 | 83.0 | 2 | 1 | 46.3 | 0.242 | 0.0287 | 1002/890 | 100.43 | 104.90 | 5.04 |
| [131,136) | g | sieve | 48 | 96 | 83.0 | 2 | 1 | 34.1 | 0.252 | 0.0388 | 1133/759 | 100.21 | 104.95 | 5.10 |
| [136,141) | f | sieve | 48 | 101 | 85.5 | 2 | 1 | 54.5 | 0.242 | 0.0246 | 1020/898 | 100.40 | 105.01 | 5.04 |
| [136,141) | g | sieve | 48 | 101 | 85.5 | 2 | 1 | 37.2 | 0.252 | 0.0357 | 1156/762 | 100.18 | 105.06 | 5.10 |
| [141,146) | f | sieve | 48 | 96 | 96.0 | 3 | 0 | 36.7 | 0.256 | 0.0364 | 1037/905 | 101.63 | 105.34 | 5.56 |
| [141,146) | g | sieve | 56 | 101 | 85.5 | 2 | 1 | 35.5 | 0.253 | 0.0373 | 1194/769 | 103.53 | 105.46 | 5.11 |
| [146,151) | f | sieve | 48 | 106 | 93.6 | 2 | 1 | 41.4 | 0.261 | 0.0327 | 1055/913 | 101.04 | 105.56 | 5.80 |
| [146,151) | g | sieve | 48 | 101 | 101.0 | 3 | 0 | 32.1 | 0.266 | 0.0417 | 1201/767 | 101.58 | 105.63 | 5.62 |
| [151,156) | f | sieve | 56 | 106 | 93.6 | 2 | 1 | 42.6 | 0.261 | 0.0318 | 1088/925 | 103.69 | 105.86 | 5.80 |
| [151,156) | g | sieve | 56 | 106 | 98.1 | 2 | 1 | 33.5 | 0.272 | 0.0403 | 1239/774 | 103.71 | 106.02 | 5.85 |
| [156,161) | f | sieve | 48 | 116 | 100.2 | 2 | 1 | 45.2 | 0.266 | 0.0303 | 1090/928 | 101.69 | 106.13 | 6.21 |
| [156,161) | g | sieve | 56 | 106 | 106.0 | 3 | 0 | 37.1 | 0.283 | 0.0373 | 1262/777 | 104.05 | 106.31 | 6.23 |
| [161,166) | f | sieve | 56 | 116 | 100.2 | 2 | 1 | 46.1 | 0.266 | 0.0297 | 1123/940 | 103.87 | 106.35 | 6.21 |
| [161,166) | g | sieve | 56 | 116 | 105.9 | 2 | 1 | 36.9 | 0.288 | 0.0377 | 1284/779 | 103.90 | 106.57 | 6.56 |
| [166,171) | f | sieve | 56 | 126 | 98.0 | 2 | 1 | 34.3 | 0.266 | 0.0392 | 1141/948 | 103.70 | 106.65 | 6.56 |
| [166,171) | g | sieve | 56 | 131 | 100.5 | 2 | 1 | 41.8 | 0.278 | 0.0331 | 1307/782 | 103.59 | 106.86 | 6.85 |
| [171,176) | f | sieve | 56 | 126 | 106.8 | 2 | 1 | 40.5 | 0.277 | 0.0340 | 1158/955 | 104.24 | 107.06 | 7.14 |
| [171,176) | g | sieve | 56 | 131 | 110.0 | 2 | 1 | 37.7 | 0.277 | 0.0364 | 1329/784 | 103.93 | 107.17 | 7.31 |
| [176,181) | f | sieve | 56 | 141 | 105.5 | 2 | 1 | 58.0 | 0.277 | 0.0243 | 1176/963 | 104.04 | 107.21 | 7.15 |
| [176,181) | g | sieve | 48 | 131 | 100.5 | 3 | 1 | 67.4 | 0.363 | 0.0238 | 1336/782 | 104.32 | 107.28 | 7.22 |
| [181,186) | f | sieve | 56 | 146 | 108.0 | 2 | 1 | 42.2 | 0.277 | 0.0327 | 1193/970 | 104.35 | 107.60 | 7.58 |
| [181,186) | g | sieve | 56 | 131 | 100.5 | 3 | 1 | 64.7 | 0.364 | 0.0247 | 1374/789 | 105.47 | 107.49 | 7.22 |
| [186,191) | f | sieve | 56 | 151 | 110.5 | 2 | 1 | 40.1 | 0.289 | 0.0349 | 1211/978 | 104.49 | 107.76 | 7.83 |
| [186,191) | g | sieve | 48 | 141 | 105.5 | 3 | 1 | 52.6 | 0.380 | 0.0308 | 1381/787 | 104.20 | 107.62 | 7.61 |
| [191,196) | f | sieve | 56 | 161 | 115.5 | 2 | 1 | 76.6 | 0.299 | 0.0200 | 1228/985 | 105.23 | 108.36 | 8.24 |
| [191,196) | g | sieve | 56 | 141 | 105.5 | 3 | 1 | 49.8 | 0.381 | 0.0324 | 1419/794 | 105.36 | 107.78 | 7.61 |
| [196,201) | f | sieve | 56 | 166 | 118.0 | 2 | 1 | 53.3 | 0.299 | 0.0271 | 1246/993 | 105.72 | 108.67 | 8.47 |
| [196,201) | g | sieve | 56 | 151 | 110.5 | 3 | 1 | 100.4 | 0.368 | 0.0200 | 1442/797 | 106.16 | 108.33 | 8.22 |
| [201,206) | f | sieve | 64 | 171 | 120.5 | 2 | 1 | 110.8 | 0.299 | 0.0200 | 1277/1005 | 107.72 | 109.28 | 8.69 |
| [201,206) | g | sieve | 48 | 146 | 119.9 | 3 | 1 | 54.1 | 0.396 | 0.0307 | 1448/794 | 105.88 | 108.56 | 8.38 |
| [206,211) | f | sieve | 56 | 181 | 125.5 | 2 | 1 | 65.8 | 0.299 | 0.0222 | 1281/1008 | 107.64 | 109.56 | 9.14 |
| [206,211) | g | sieve | 56 | 161 | 115.5 | 3 | 1 | 86.6 | 0.396 | 0.0200 | 1487/802 | 105.97 | 108.80 | 8.65 |
| [211,216) | f | sieve | 64 | 181 | 125.5 | 2 | 1 | 51.2 | 0.299 | 0.0281 | 1312/1020 | 108.46 | 109.82 | 9.14 |
| [211,216) | g | sieve | 56 | 166 | 118.0 | 3 | 1 | 67.1 | 0.396 | 0.0251 | 1509/804 | 105.72 | 109.04 | 8.90 |
| [216,221) | f | ecm | 16 | 171 | 136.3 | 5 | 1 | 20839.3 | 0.585 | 0.0200 | 1186/980 | 116.77 | 116.79 | 10.37 |
| [216,221) | g | sieve | 64 | 171 | 120.5 | 3 | 1 | 109.0 | 0.383 | 0.0200 | 1546/812 | 107.78 | 109.54 | 9.06 |
| [221,226) | f | ecm | 16 | 176 | 139.6 | 5 | 1 | 33002.4 | 0.586 | 0.0200 | 1203/987 | 117.61 | 117.62 | 10.59 |
| [221,226) | g | sieve | 56 | 181 | 125.5 | 3 | 1 | 116.0 | 0.397 | 0.0200 | 1554/809 | 106.46 | 109.65 | 9.51 |
| [226,231) | f | ecm | 16 | 171 | 148.9 | 5 | 1 | 92017.0 | 0.586 | 0.0200 | 1221/995 | 119.13 | 119.13 | 10.79 |
| [226,231) | g | sieve | 64 | 181 | 125.5 | 3 | 1 | 84.6 | 0.397 | 0.0202 | 1591/816 | 107.80 | 109.87 | 9.51 |
| [231,236) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 93330.3 | 0.670 | 0.0200 | 1238/1002 | 117.95 | 117.96 | 11.06 |
| [231,236) | g | sieve | 56 | 176 | 139.6 | 3 | 1 | 127.6 | 0.397 | 0.0200 | 1599/814 | 108.32 | 110.11 | 9.71 |
| [236,241) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 52165.9 | 0.682 | 0.0200 | 1256/1010 | 117.95 | 117.96 | 11.12 |
| [236,241) | g | sieve | 40 | 191 | 149.4 | 3 | 1 | 67.8 | 0.397 | 0.0249 | 1587/805 | 107.11 | 110.35 | 10.26 |
| [241,246) | f | ecm | 16 | 171 | 148.9 | 6 | 1 | 140995.3 | 0.670 | 0.0200 | 1273/1017 | 119.62 | 119.63 | 11.31 |
| [241,246) | g | sieve | 48 | 191 | 149.4 | 3 | 1 | 85.8 | 0.398 | 0.0200 | 1628/814 | 108.29 | 110.53 | 10.26 |
| [246,251) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 422147.8 | 0.670 | 0.0200 | 1291/1025 | 119.71 | 119.71 | 11.64 |
| [246,251) | g | sieve | 56 | 191 | 149.4 | 3 | 1 | 80.4 | 0.398 | 0.0212 | 1667/822 | 109.08 | 110.74 | 10.26 |
| [251,256) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 184618.8 | 0.670 | 0.0200 | 1308/1032 | 119.71 | 119.71 | 11.64 |
| [251,256) | g | sieve | 56 | 196 | 152.7 | 3 | 1 | 67.8 | 0.398 | 0.0249 | 1689/824 | 108.94 | 110.98 | 10.45 |
| [256,261) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 106565.6 | 0.683 | 0.0200 | 1326/1040 | 119.71 | 119.71 | 11.70 |
| [256,261) | g | sieve | 64 | 201 | 156.0 | 3 | 1 | 174.4 | 0.398 | 0.0200 | 1726/832 | 110.25 | 111.57 | 10.69 |
| [261,266) | f | ecm | 16 | 186 | 160.6 | 6 | 1 | 438902.5 | 0.670 | 0.0200 | 1343/1047 | 121.84 | 121.84 | 11.88 |
| [261,266) | g | sieve | 64 | 206 | 159.2 | 3 | 1 | 144.0 | 0.399 | 0.0200 | 1748/834 | 110.69 | 111.99 | 10.96 |

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work:** `2^102.497353` <= **`2^102.498`** core units (`C = 20,000`), `-6.366` bits below `c36` (`108.864`) even with full de Bruijn enabled everywhere and `m_max = 2580.00`b.
- **Metered cap:** `t = 2^115.6` core units (reached = failure). By Markov's inequality, `P(cap) <= 2^(102.498 - 115.6) = 1.137378e-04` (`1.251116e-04` at `+10%` headroom).
- **Expected node count (`counts` recurrence with unit node cost, verified live in `audit_cm.py`):** `E[N] = 2^16.0696` (`68773.7` nodes for both logs). Setting hard node abort `N_MAX = 2^25.8` gives `P(counter abort) <= 2^(16.0696 - 25.8) <= 2^-9.7304`.
- **Repeat-ideal bound (conservative over `tau`-orbits):** The smallest bin `[70, 76)` holds `2 (li(2^76) - li(2^70)) / 12 >= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits across both sides (even the sub-bin `[70, 75)` holds `>= 2^66.66`). Conservatively bounding all `N_MAX` nodes into `2^66.66` orbits with `2^6` `1/p`-weighting slack gives `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-10.06) <= 2^-10.05`.
