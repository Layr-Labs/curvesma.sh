# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `MAIN_LB = 2^71.29`, `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (`<= lr_rethi - 0.009 - 2 SE = 63.4846`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes**, `BOOT_RELS = 48`, `RET_BOOT = 0.50`, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, the empirical `q`-lattice offsets, and the boot deciles (`m_max = 2600.0`b). `audit_cm.py` invokes `cm_models.eval_descent_live` to compute every node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `psi`, `phi`, `W`, `EN`, and boot work/failure 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`105.614078 <= 105.615`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`105.597846`).

## 1. Fresh-Coin Randomization, Constructive 16D LLL Rank-2 Submodule, & Boot

For each discrete logarithm (`log g` with `g = e(P0, Q0)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy & Constructive 16D Integer LLL Rank-2 Submodule Extraction

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, neither `w_0` nor `w_1` can vanish for `z' in G_T \ {1}`; setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). For any `(a, b) in L_Guillevic` with `b != 0 (mod P)`, evaluating **`a - b X`** gives **`Psi_g(a - b X) = (b mod P)(gamma - X) = -(b mod P) w_1^(-1) c_6 z'`**, which has exact discrete log `e log_r(z) mod r` and exact relative norm **`N_(K_g/K_H)(a - b X) = Res_X(a - b X, X^3 - omega^2 u) = a^3 - omega^2 u b^3`**.
- **Constructive extraction algorithm (NOT assuming LLL generates all of `L_Guillevic`):**
  1. Form the 16 x 16 HNF basis of `L_Guillevic` (4 `F_p`-linear congruences `a == gamma b mod (p, h_pi)`).
  2. Run ordinary 16D integer LLL under the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)` (charged `LLL_NODE = 2^38` Ir per ball; `< 2^28` Ir actual).
  3. Set `v_1 = b_1` (shortest LLL vector) and `v_2 = b_k` (the first LLL basis vector `K_H`-linearly independent of `v_1`, i.e., `det_H(v_1, b_k) = a_1 b_k - b_1 a_k != 0` in `O_H`; since `rank_Z(L_Guillevic) = 16 > 8 = rank_Z(O_H v_1)`, such a `k in {2..9}` always exists).
  4. Then **`M = Z[zeta_15] v_1 + Z[zeta_15] v_2 <= L_Guillevic`** is **provably a free rank-2 `Z[zeta_15]`-submodule**: every vector `v_1 + beta v_2` (`beta in B_8(1)`, `17` vectors `0, +-e_i`) lies in `L_Guillevic`, has zero unit duplicates (`UDUP_BOOT = 1.0`), and is base-primitive (`gcd(1, beta) = 1`).
- **Live verification (`audit_cm.py`: `voronoi_and_boot_live`):** (a) Runs the full constructive 16D integer LLL + `(b_1, b_k)` extractor live from an unreduced 16 x 16 HNF basis on split-base `p_0 = 7207 == 7 (mod 15)` (`k = 2`, `index = 1`, `a_i == gamma b_i mod P`); (b) verifies 4 bundled reduced 381-bit BLS12-381 `w'=omega^2` bases (`boot_witness.json`): `a_i == gamma b_i (mod P)`, `|N(a_1 b_2 - a_2 b_1)| == p^4` (`index = 1`), all `68/68` exact relative norms `|N(a^3 - u omega^2 b^3)| < 2600`b (mean `2551.75`b, max `2563.01`b), AND **`4/4` passes of the exact 10-decile discard predicate** (`ball_sorted[idx_j] <= boot_sub_deciles[j]`).
- **Deterministic 50% norm-certificate filter (`RET_BOOT = 0.50`, `m_max = 2600.0`b):** After LLL and `(v_1, v_2)` extraction, the algorithm computes all 17 exact norms of `a - b X` and discards any ball whose sorted deciles exceed `boot_sub_deciles` (charging `RET_BOOT = 0.50` = `2x` balls; supported by `32/32` targets in `[11]` with max `2575.73`b, `kappa_hat = 0.7076 << BALL_PAIR = 2^24`, and `4/4` live witnesses).
- **Explicit multi-factor boot ECM extraction loop & limb-scaled MR:** After `TDIV` strips primes `<= 2^20`, each `<= 2600`b norm has at most `k_lx = 131` primes in `(2^20, 2^70]`, `kA = 13` primes in `(2^70, n_A = 2^234.53125]`, and `1` prime in `(n_A, n_1 = 2^265]`. Every boot candidate is conservatively charged ALL `131 * (ecm(70, 2600, 0.999) + mr(2600)) + 13 * (ecm(n_A, 2600, 0.999) + mr(2600)) + mr(2600) + IDEAL` instructions (zero early-exit credit; `mr_instr(2600) = 5.14e9` Ir limb-scaled 40-round MR) with multiplicative detection **`det_b = 0.999^145 = 0.8650`**.
- **Certified Boot (`audit_cm.py` live):** `n_1 = 265`, `n_A = 234.53125`, `k_A = 13`, `BOOT_RELS = 48`, `RET_BOOT = 0.50`:
  - Default `MESH = 96/192`: `pr = 2^-37.905005`, `miss = 0.225854`, `psi = 0.284517`, **`phi_boot = psi^48 = 6.2690e-27 <= 1.0e-12`**, boot own `2^116.2845` + subtree `2^118.6451` -> **`two_log = 2^105.614078 <= 2^105.615`** core units.
  - On-demand `MESH = 384/384` (`{"descent_mesh": 384}`): **`two_log = 2^105.597846`**, `phi_boot = 6.2690e-27`.

## 2. Special-`q` Nodes, Exact `P(M < k) + psi^k` Union Bound, & Full Survivor Work

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1, v_2)` and enumerates `alpha v_1 + beta v_2` (`alpha` unit-canonical in `B_8(R)`, `UDUP_NODE = 0.95` charged). Measured across 56 BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`): worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b (`cm_data.json`). `[7e]` child-conditioning shift across `22` encountered vs `24` random ideals: `(+3.40, +2.35)`b inside `(+10, +4)`b.

- **Exact `min_k (P(M < k) + psi^k)` union bound (NO `PGF(psi)` substitution!) & full survivor work:** Under `H5` (node relation count `M ~ Poisson(mu)` and explicit conditional child-failure hazard `P(cand j fails | prior j-1 failed) <= psi`), each node's failure probability is bounded by the elementary union bound **`phi_node = max(phi_floor, disc + min_{k >= 1} (P(M < k) + psi^k))`** (evaluated live in `cm_models.phi_poisson_union`). Every sieve node is charged **100% full survivor work (`1.0 x survivors * surv_cost`)** (zero lazy-ordering premise!), and every survivor is charged the FULL worst-case multi-factor ECM extraction loop (`k_lx_f`, `k_lx_g` schedules at `<= LX`, `2 kA` schedules at `<= nA`, limb-scaled `mr_instr(nb)` throughout, and `det_n = 0.999^(k_lx_f + k_lx_g + 2*(kA + kB))`).

## 3. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

| Bin | S | How | R | n | nA | kA | kB | mu | psi | phi | Nf/Ng | own | W | EN |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| [70,76) | f | sieve | 48 | 70 | 70.0 | 0 | 0 | 336.9 | 0.000 | 0.0200 | 793/800 | 99.95 | 99.95 | 0.00 |
| [70,76) | g | sieve | 48 | 70 | 70.0 | 0 | 0 | 347.0 | 0.000 | 0.0200 | 863/730 | 99.95 | 99.95 | 0.00 |
| [76,81) | f | sieve | 40 | 76 | 76.0 | 1 | 0 | 448.4 | 0.040 | 0.0200 | 791/801 | 96.34 | 101.07 | 1.62 |
| [76,81) | g | sieve | 40 | 76 | 76.0 | 1 | 0 | 461.7 | 0.040 | 0.0200 | 867/725 | 96.38 | 101.07 | 1.62 |
| [81,86) | f | sieve | 48 | 76 | 76.0 | 1 | 0 | 601.5 | 0.040 | 0.0200 | 827/815 | 100.02 | 101.60 | 1.62 |
| [81,86) | g | sieve | 48 | 76 | 76.0 | 1 | 0 | 587.7 | 0.040 | 0.0200 | 908/734 | 100.03 | 101.60 | 1.62 |
| [86,91) | f | sieve | 32 | 76 | 76.0 | 2 | 0 | 26.6 | 0.080 | 0.0200 | 803/809 | 98.79 | 102.21 | 2.42 |
| [86,91) | g | sieve | 32 | 76 | 76.0 | 2 | 0 | 26.6 | 0.080 | 0.0200 | 889/723 | 98.87 | 102.22 | 2.42 |
| [91,96) | f | sieve | 40 | 76 | 76.0 | 2 | 0 | 55.1 | 0.080 | 0.0200 | 843/824 | 100.79 | 102.57 | 2.42 |
| [91,96) | g | sieve | 40 | 76 | 76.0 | 2 | 0 | 51.6 | 0.080 | 0.0200 | 934/733 | 100.69 | 102.54 | 2.42 |
| [96,101) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 72.4 | 0.080 | 0.0200 | 880/838 | 102.30 | 103.19 | 2.42 |
| [96,101) | g | sieve | 48 | 76 | 76.0 | 2 | 0 | 62.9 | 0.080 | 0.0200 | 976/742 | 102.14 | 103.10 | 2.42 |
| [101,106) | f | sieve | 48 | 76 | 76.0 | 2 | 0 | 23.6 | 0.084 | 0.0200 | 897/845 | 101.28 | 102.87 | 2.60 |
| [101,106) | g | sieve | 48 | 81 | 75.5 | 1 | 1 | 61.9 | 0.080 | 0.0200 | 998/744 | 101.39 | 103.21 | 3.30 |
| [106,111) | f | sieve | 48 | 81 | 75.5 | 1 | 1 | 24.4 | 0.084 | 0.0200 | 915/853 | 100.78 | 103.18 | 3.41 |
| [106,111) | g | sieve | 48 | 81 | 81.0 | 2 | 0 | 43.2 | 0.080 | 0.0200 | 1021/747 | 100.91 | 103.46 | 3.85 |
| [111,116) | f | sieve | 48 | 81 | 75.5 | 2 | 1 | 27.8 | 0.120 | 0.0200 | 932/860 | 108.30 | 108.34 | 3.65 |
| [111,116) | g | sieve | 48 | 81 | 75.5 | 2 | 1 | 26.2 | 0.124 | 0.0200 | 1043/749 | 107.66 | 107.73 | 3.74 |
| [116,121) | f | sieve | 48 | 86 | 86.0 | 2 | 0 | 27.6 | 0.084 | 0.0200 | 950/868 | 100.44 | 103.94 | 3.93 |
| [116,121) | g | sieve | 48 | 96 | 83.0 | 1 | 1 | 60.5 | 0.080 | 0.0200 | 1066/752 | 100.53 | 104.38 | 4.27 |
| [121,126) | f | sieve | 48 | 96 | 83.0 | 1 | 1 | 28.0 | 0.084 | 0.0200 | 967/875 | 100.30 | 104.42 | 4.33 |
| [121,126) | g | sieve | 48 | 101 | 85.5 | 1 | 1 | 43.8 | 0.080 | 0.0200 | 1088/754 | 100.37 | 104.76 | 4.27 |
| [126,131) | f | sieve | 48 | 96 | 96.0 | 2 | 0 | 29.4 | 0.084 | 0.0200 | 985/883 | 100.41 | 104.79 | 4.65 |
| [126,131) | g | sieve | 48 | 116 | 100.2 | 2 | 0 | 38.5 | 0.080 | 0.0200 | 1111/757 | 100.55 | 105.32 | 4.60 |
| [131,136) | f | sieve | 48 | 96 | 83.0 | 2 | 1 | 44.0 | 0.120 | 0.0200 | 1002/890 | 106.34 | 106.76 | 4.76 |
| [131,136) | g | sieve | 48 | 96 | 83.0 | 2 | 1 | 32.4 | 0.125 | 0.0200 | 1133/759 | 105.11 | 105.97 | 4.82 |
| [136,141) | f | sieve | 48 | 101 | 85.5 | 2 | 1 | 51.6 | 0.120 | 0.0200 | 1020/898 | 106.21 | 106.75 | 4.76 |
| [136,141) | g | sieve | 48 | 101 | 85.5 | 2 | 1 | 35.2 | 0.125 | 0.0200 | 1156/762 | 104.88 | 106.01 | 4.82 |
| [141,146) | f | sieve | 48 | 96 | 96.0 | 3 | 0 | 34.8 | 0.120 | 0.0200 | 1037/905 | 107.24 | 107.58 | 5.23 |
| [141,146) | g | sieve | 56 | 101 | 85.5 | 2 | 1 | 33.6 | 0.125 | 0.0200 | 1194/769 | 105.66 | 106.42 | 4.82 |
| [146,151) | f | sieve | 48 | 106 | 93.6 | 2 | 1 | 38.9 | 0.120 | 0.0200 | 1055/913 | 106.26 | 106.96 | 5.44 |
| [146,151) | g | sieve | 48 | 101 | 101.0 | 3 | 0 | 30.3 | 0.125 | 0.0200 | 1201/767 | 106.30 | 107.14 | 5.27 |
| [151,156) | f | sieve | 56 | 106 | 93.6 | 2 | 1 | 40.0 | 0.120 | 0.0200 | 1088/925 | 107.05 | 107.49 | 5.44 |
| [151,156) | g | sieve | 56 | 106 | 98.1 | 2 | 1 | 31.5 | 0.125 | 0.0200 | 1239/774 | 106.45 | 107.23 | 5.48 |
| [156,161) | f | sieve | 48 | 116 | 100.2 | 2 | 1 | 41.9 | 0.120 | 0.0200 | 1090/928 | 106.64 | 109.55 | 5.81 |
| [156,161) | g | sieve | 56 | 106 | 106.0 | 3 | 0 | 34.8 | 0.126 | 0.0200 | 1262/777 | 107.58 | 108.00 | 5.82 |
| [161,166) | f | sieve | 56 | 116 | 100.2 | 2 | 1 | 42.8 | 0.120 | 0.0200 | 1123/940 | 107.30 | 109.66 | 5.81 |
| [161,166) | g | sieve | 56 | 116 | 105.9 | 2 | 1 | 34.3 | 0.126 | 0.0200 | 1284/779 | 106.80 | 109.58 | 6.11 |
| [166,171) | f | sieve | 56 | 126 | 98.0 | 2 | 1 | 31.4 | 0.120 | 0.0200 | 1141/948 | 106.44 | 109.52 | 6.13 |
| [166,171) | g | sieve | 56 | 131 | 100.5 | 2 | 1 | 38.1 | 0.126 | 0.0200 | 1307/782 | 105.47 | 109.45 | 6.37 |
| [171,176) | f | sieve | 56 | 126 | 106.8 | 2 | 1 | 37.2 | 0.120 | 0.0200 | 1158/955 | 107.76 | 109.77 | 6.68 |
| [171,176) | g | sieve | 56 | 131 | 110.0 | 2 | 1 | 34.4 | 0.120 | 0.0200 | 1329/784 | 106.66 | 109.56 | 6.83 |
| [176,181) | f | sieve | 56 | 141 | 105.5 | 2 | 1 | 52.2 | 0.120 | 0.0200 | 1176/963 | 107.77 | 109.76 | 6.65 |
| [176,181) | g | sieve | 48 | 131 | 100.5 | 3 | 1 | 60.4 | 0.166 | 0.0200 | 1336/782 | 110.06 | 110.79 | 6.64 |
| [181,186) | f | sieve | 56 | 146 | 108.0 | 2 | 1 | 37.7 | 0.120 | 0.0200 | 1193/970 | 107.90 | 109.80 | 7.05 |
| [181,186) | g | sieve | 56 | 131 | 100.5 | 3 | 1 | 57.9 | 0.166 | 0.0200 | 1374/789 | 110.88 | 111.34 | 6.64 |
| [186,191) | f | sieve | 56 | 151 | 110.5 | 2 | 1 | 35.6 | 0.126 | 0.0200 | 1211/978 | 107.98 | 109.83 | 7.26 |
| [186,191) | g | sieve | 48 | 141 | 105.5 | 3 | 1 | 46.4 | 0.166 | 0.0200 | 1381/787 | 109.33 | 110.40 | 6.97 |
| [191,196) | f | sieve | 56 | 161 | 115.5 | 2 | 1 | 67.2 | 0.126 | 0.0200 | 1228/985 | 109.66 | 111.64 | 7.59 |
| [191,196) | g | sieve | 56 | 141 | 105.5 | 3 | 1 | 43.9 | 0.166 | 0.0200 | 1419/794 | 110.10 | 110.82 | 6.98 |
| [196,201) | f | sieve | 56 | 166 | 118.0 | 2 | 1 | 46.4 | 0.126 | 0.0200 | 1246/993 | 109.83 | 111.94 | 7.81 |
| [196,201) | g | sieve | 56 | 151 | 110.5 | 3 | 1 | 87.2 | 0.160 | 0.0200 | 1442/797 | 111.29 | 111.65 | 7.55 |
| [201,206) | f | sieve | 64 | 171 | 120.5 | 2 | 1 | 95.9 | 0.126 | 0.0200 | 1277/1005 | 111.79 | 112.68 | 8.00 |
| [201,206) | g | sieve | 48 | 146 | 119.9 | 3 | 1 | 47.4 | 0.166 | 0.0200 | 1448/794 | 110.60 | 112.01 | 7.68 |
| [206,211) | f | sieve | 56 | 181 | 125.5 | 2 | 1 | 56.3 | 0.126 | 0.0200 | 1281/1008 | 112.09 | 113.07 | 8.41 |
| [206,211) | g | sieve | 56 | 161 | 115.5 | 3 | 1 | 74.0 | 0.167 | 0.0200 | 1487/802 | 110.71 | 112.23 | 7.87 |
| [211,216) | f | sieve | 64 | 181 | 125.5 | 2 | 1 | 43.7 | 0.126 | 0.0200 | 1312/1020 | 112.28 | 113.17 | 8.41 |
| [211,216) | g | sieve | 56 | 166 | 118.0 | 3 | 1 | 56.9 | 0.167 | 0.0200 | 1509/804 | 110.40 | 112.34 | 8.11 |
| [216,221) | f | ecm | 16 | 171 | 136.3 | 5 | 1 | 17342.8 | 0.245 | 0.0200 | 1186/980 | 119.09 | 119.11 | 9.16 |
| [216,221) | g | sieve | 64 | 171 | 120.5 | 3 | 1 | 91.5 | 0.160 | 0.0200 | 1546/812 | 112.23 | 113.07 | 8.26 |
| [221,226) | f | ecm | 16 | 176 | 139.6 | 5 | 1 | 27263.1 | 0.246 | 0.0200 | 1203/987 | 119.93 | 119.94 | 9.36 |
| [221,226) | g | sieve | 56 | 181 | 125.5 | 3 | 1 | 96.4 | 0.167 | 0.0200 | 1554/809 | 111.49 | 112.96 | 8.66 |
| [226,231) | f | ecm | 16 | 171 | 148.9 | 5 | 1 | 77291.6 | 0.246 | 0.0200 | 1221/995 | 121.45 | 121.45 | 9.55 |
| [226,231) | g | sieve | 64 | 181 | 125.5 | 3 | 1 | 69.8 | 0.167 | 0.0200 | 1591/817 | 111.88 | 113.12 | 8.66 |
| [231,236) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 77422.5 | 0.280 | 0.0200 | 1238/1002 | 120.53 | 120.54 | 9.58 |
| [231,236) | g | sieve | 56 | 176 | 139.6 | 3 | 1 | 107.1 | 0.167 | 0.0200 | 1599/814 | 113.57 | 113.97 | 8.85 |
| [236,241) | f | ecm | 16 | 176 | 139.6 | 6 | 1 | 43092.0 | 0.286 | 0.0200 | 1256/1010 | 120.53 | 120.54 | 9.59 |
| [236,241) | g | sieve | 40 | 191 | 149.4 | 3 | 1 | 56.3 | 0.167 | 0.0200 | 1587/805 | 112.02 | 113.33 | 9.30 |
| [241,246) | f | ecm | 16 | 171 | 148.9 | 6 | 1 | 118958.9 | 0.280 | 0.0200 | 1273/1017 | 122.21 | 122.21 | 9.81 |
| [241,246) | g | sieve | 48 | 191 | 149.4 | 3 | 1 | 71.2 | 0.167 | 0.0200 | 1628/814 | 113.12 | 113.88 | 9.31 |
| [246,251) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 345204.5 | 0.280 | 0.0200 | 1291/1025 | 122.29 | 122.29 | 10.07 |
| [246,251) | g | sieve | 56 | 191 | 149.4 | 3 | 1 | 66.3 | 0.167 | 0.0200 | 1667/822 | 113.81 | 114.32 | 9.31 |
| [251,256) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 150292.0 | 0.280 | 0.0200 | 1308/1032 | 122.29 | 122.29 | 10.07 |
| [251,256) | g | sieve | 56 | 196 | 152.7 | 3 | 1 | 55.5 | 0.168 | 0.0200 | 1689/824 | 113.81 | 114.50 | 9.46 |
| [256,261) | f | ecm | 16 | 191 | 149.4 | 6 | 1 | 86349.0 | 0.286 | 0.0200 | 1326/1040 | 122.29 | 122.29 | 10.08 |
| [256,261) | g | sieve | 64 | 201 | 156.0 | 3 | 1 | 140.3 | 0.168 | 0.0200 | 1726/832 | 115.81 | 116.05 | 9.67 |
| [261,266) | f | ecm | 16 | 186 | 160.6 | 6 | 1 | 364423.0 | 0.280 | 0.0200 | 1343/1047 | 124.43 | 124.43 | 10.32 |
| [261,266) | g | sieve | 64 | 206 | 159.2 | 3 | 1 | 114.4 | 0.168 | 0.0200 | 1748/834 | 116.05 | 116.36 | 9.91 |

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work:** `2^105.614078` <= **`2^105.615`** core units (`2^105.597846` at `MESH = 384/384`), `-3.249` bits vs `c36` (`108.864`).
- **Metered cap:** `t = 2^117.5` core units; by Markov's inequality, `P(cap) <= 2^(105.615 - 117.5) = 2.643982e-04` (`2.908380e-04` at `+10%` headroom).
- **Node counter (`E[N]` live in `audit_cm.py`):** `E[N] = 2^14.2497`, hard abort `N_MAX = 2^25.8` -> `P(counter abort) <= 2^-11.5502`.
- **Repeat-ideal bound:** Smallest bin `[70, 76)` holds `>= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits; with `2^6` `1/p`-weighting slack, `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-10.06) <= 2^-10.05`.
