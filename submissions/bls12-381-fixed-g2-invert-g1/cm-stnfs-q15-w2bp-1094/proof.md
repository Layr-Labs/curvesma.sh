# recover-P: `(4,3)` `Q(zeta_15)` CM-STNFS + `Arch2-Prod-ZeroCopy-U32` 1-Owner Build + 1-Owner 4-Vec Arch B LA + Per-Environment `H5-Blind(a,b)` (`R_tree = 44`): expected work `2^109.448056241876`

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity `Q in G2`, `z in G_T`; output `P in G1` with `e(P, Q) = z` for exact supplied `Q` (`P = infinity` for `z = 1`). Unit: one dense `F_(p^12)` mul, `C = 20,000` baseline scalar ARMv8.0-A integer instructions (`simd = 0`, zero FP/NEON/SVE, zero `FEAT_LSE`/`DC ZVA`, zero external C measurements, zero black-box generator citations, zero amortization). Full `rho_dB` baked into Collection & Descent. All kernels (`KERNELS.md`, parsed & bit-exact-simulated baseline scalar ARMv8.0-A integer blocks), exact local-root Murphy-`alpha` & odd-integer `kap2_g` derivation (`p <= 10,000`), real-FB joint smoothness calibration (`800` coprime pairs), zero-column-loss `Z[zeta_12]` structured Gaussian filtering on `F_2core` (`2^16..2^24` x 6 seeds, `post_peel_cols_lost_max == 0`, `post_merge_cols_lost_max == 0`), 1-owner `< 2^64 B` 4-vector Arch B Block Wiedemann LA (`n = 112`, `m = 8192`, `te = 0`), per-environment stopped-process `FB`-blinded CRT descent (`lx = 71.0`, NO local per-scan cutoff, single deterministic tree cap `T_tree = 5.0*W_tree_1`, `R_tree = 44` conditionally `i.i.d.` restarts given `E in G_pre`), quantitative Law-of-Total-Variance & Gamma-Chernoff tail bridge (`h5_kappa_tail_bridge`), pure `Decimal(50)` `unitdomain.py`, and direct 200k-draw exact-`O_H`-coprime pure-product pool are self-contained (**no prior verdict or score is a premise**).

## Claim

`W` is total classical work on one target. Step 0: uniform 64b `u`; unless `u < K = 7195777280862110484` (least `K` with `K * (999785/1000000) >= 0.39 * 2^64`), abort (`< 1` u). Under declared exploratory STNFS model inequalities `H1/H2, H3, H-R, H-SM, H4, H4u, H-X, H5-Blind(a,b)` (§8; `RANK-BRIDGE.md`, `DESCENT.md`, `ADAPTIVE.md`), an opened run succeeds w.p. `>= 999785/1000000` (Boole, §6); overall success `>= 0.39` over uniform production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^110.8062 + 2^66) < 2^109.448056241876  dense-Fp12 units (C = 20,000).

`109.448056241876` is an **expected** total work bound `E[W]` (NOT every-run). The **ENTIRE opened run** runs under ONE hard aborting work meter `t_run = 2^123.0` u (block boundaries `<= 2^66` u): deterministic pointwise cap `W_actual <= 2^123.0 + 2^66` u, `E[min(W_open, t_run)] <= 2^110.8062`, Markov `p_meter <= 0.000214519364`. Analytic exploratory claim; all numbers match default `32/185` and refined `384/370` `audit_cm.py` (`checks_ok = true`, pure stdlib).

## 1. Reduction (fixed `Q in G2`, recover `P in G1`: two `G_T` logs) and field map

With pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `P0 in G1`:
0. **Run gate:** uniform 64b `u`; if `u >= K`, abort (`< 1` u).
1. **Validate:** `Q in G2 \ {infinity}` (`[r]Q = infinity`), `z^r = 1`; if `z = 1`, return `infinity`.
2. **Pairing & map:** `g = e(P0, Q)`, map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Randomize & solve:** independent uniform `e_1, e_2 in {0..r-1}`, core on `g^(e_1), z^(e_2)`, `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** `P = [a]P0`; return `P` iff `e(P, Q) = z`. Wrapper `< 2^66` u, failure `< 2^-48`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, Exact Murphy-`alpha`, & Order-12 Frobenius

BLS12-381: `u = -0xd201000000010000`, `r = u^4 - u^2 + 1 == 13 mod 15`, `p = (u-1)^2 r/3 + u == 7 mod 15`. Over PID `O_H = Z[zeta_15]` (`D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `O_H/P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` gives `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega+1)X^6 + X^3 - omega` (deg 72). `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`ord_r(p) = 12`) quotients into `Z[zeta_12]` 12-orbits (`algebra_bls12`). `check_algebra_bls12` verifies via exact polynomial GCD that `[2048, 8192]` contains `192` distinct unramified `O_H`-primes (`ram = 0`) admitting `f`-side degree-1 ideals for `B_blind` (`N_blind_pool = 192`). Pre-collection `W_f in F_r^(3x6), W_g in F_r^(1x2)` compress `P_chi Lambda` to `4` cols (`< 4/r` under `H-SM`).
- **Exact Local-Root Murphy-`alpha` & Odd-Integer Shift (`murphy_alpha_exact_live`):** `audit_cm.py` computes exact local valuation excess `sum_p (1/(p-1) - E[v_p(N)]) * log2(p)` across ramified `p in {2, 3, 5}` (`p=2`: `0` roots of `f` over `F_16`; `p=3`: exact `F_81` Hensel lifting giving `E[v_3(N_f)]=8/41, E[v_3(N_g)]=2/41`; `p=5`: `0` roots over `F_25`) and all `1,224` unramified primes `7 <= p <= 10,000`: **`alpha_f(<=10k) = 2.70436 <= 2.74`**, **`alpha_(g,>=3)(<=10k) = 1.10299 <= 1.28`** (`+0.2127b` conservative norm cushion, truncation yield bound `< 0.015b`), and exact odd-integer de Bruijn identity `1.0 + kappa_2 = 1.0525` (`max_err = 0.000102b`).
- **Live Real-FB Joint Smoothness Bridge (`real_fb_joint_coll_live`):** Direct trial-factorization smoothness check on `800` exact-`O_H`-coprime pairs (`B_fb=512, LP=8192`) confirms `obs_joint = 53` vs `mod_joint = 49.45` (`obs/mod = 1.0719`, `+0.1001b`).
- Live E2E (`toy_e2e_dlp`, `u_0=-5, p_0=7207, r_0=601`): `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps, `5/5` unseen DLPs.

## 3. Relation Collection: `Arch2-Prod-ZeroCopy-U32` (1-Owner Flat-Array Build) & Direct 200k Exact-`O_H`-Coprime Pool

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373`), excess `1.09` (`2^63.441701` rels).
- **Pure-product geometry & 2-adic `c(a)` (`239/239` exact-1):** `G = A_4 x A_2`, 10 Voronoi units (`f_Vor = 3/(5*sqrt(5)) = 0.268328`), pure product `q_a(a) <= s_{c(a)} R2_a`, `q_b(b) <= R2_b` (`A^2 = 317500.0`, `r_(ab) = 0.9815`, `s_c = (1.0, 1.5, 2.0)`, `STORAGE.md`). `N_a^can = 2^56.7697 < 0.90 * 2^57`, `L_max = 402.17 < W_tile = 416`, `V_kept = 2^120.2834`.
- **Direct 200k-Draw Exact-`O_H`-Coprime Pool (`P4_prod_2ad_prim_w2_317`, `n_draws = 200,000`):** Direct `gcd_O(a,b) == (1)` test on `200,000` draws (`f_prim = 0.95591`, `se_y_bits = 0.004684`). Retained yield: midpoint `2^63.704281`, upper-ret `2^63.698561` (**`margin_net = +0.247861`b**, `+19.49%` over `1.09*U'`), and **strict upper-endpoint smoothness + retention `2^63.558812` (`margin_strict_hi = +0.117111`b = `25.00 * SE`)**.
- **`Arch2-Prod-ZeroCopy-U32` (`W_tile=416, W_SB=32, P_m4=2^56.5, P_lo=P_f=2^59, P_g=2^64`):**
  - **1-owner simultaneous build (`<= 2^63.889 B < 2^64 B`):** `2^105.894017 <= 2^105.894517` u, total replay tables **`2^119.074 B < 2^120 B`**.
  - **Replay & Cofactorization:** Band-1a `rev_b1a_apply_u8` (`3.8`), Band-1b `rev_b1b_apply_u8` (`5.8`), Band-2 `rev_b2_prod_emit_1` + `u8` (`16.35`, `c_bd_chk=3.0`), `rev_scan_zero_u32` (`0.95`), fixed `g`-first cofactorization (`c_cof_2d_upper = 1.3784` ipc). Full-cell sieve (`ipc_tot = 13.8269`, `h2_allow = 1.08`) -> **`2^109.896079 <= 2^109.896579`** u.

## 4. Zero-Column-Loss `Z[zeta_12]` Filtering (`2^16..2^24`), 1-Owner 4-Vec Arch B BW LA, & Classic `PM-Basis`

6-bucket `Z[zeta_12] = Z[omega] (+) chi Z[omega]` (`filtsim_z12.py`, 6 seeds/scale, `B71.5_d04_e109`, `W = 200`, guarded unit-pivot + co-column-degree peel & exact `_dc` merge conservation):

| `U'` | Seeds | `N_l / U'` (SE) | `N_sw / U'` (SE) | `sd_l U'^(1/4)` | Depth | Unit % (Exp.) | `k5` (hex6) | `tab_r` | Lost Cols |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `2^16` | 6 | `0.003535` (`0.000303`) | `0.006900` (`0.000584`) | `0.0119` | `14.3` | `96.06%` (`89.83%`) | `10.8533` | `0.80693` | `0` |
| `2^17` | 6 | `0.003138` (`0.000114`) | `0.006177` (`0.000266`) | `0.0053` | `14.5` | `96.67%` (`91.09%`) | `10.6002` | `0.80579` | `0` |
| `2^18` | 6 | `0.002779` (`0.000081`) | `0.005581` (`0.000175`) | `0.0045` | `16.2` | `97.13%` (`91.73%`) | `10.4695` | `0.80307` | `0` |
| `2^19` | 6 | `0.002659` (`0.000073`) | `0.005323` (`0.000167`) | `0.0048` | `16.8` | `97.36%` (`92.27%`) | `10.3607` | `0.80024` | `0` |
| `2^20` | 6 | `0.002406` (`0.000048`) | `0.004887` (`0.000086`) | `0.0038` | `17.7` | `97.46%` (`92.65%`) | `10.2839` | `0.79896` | `0` |
| `2^21` | 6 | `0.002196` (`0.000040`) | `0.004448` (`0.000083`) | `0.0038` | `19.0` | `97.49%` (`92.83%`) | `10.2489` | `0.79695` | `0` |
| `2^22` | 6 | `0.002003` (`0.000024`) | `0.004083` (`0.000045`) | `0.0026` | `20.0` | `97.53%` (`93.01%`) | `10.2120` | `0.79592` | `0` |
| `2^23` | 6 | `0.001916` (`0.000006`) | `0.003888` (`0.000017`) | `0.0008` | `21.0` | `97.54%` (`93.24%`) | `10.1659` | `0.79510` | `0` |
| `2^24` | 6 | `0.001776` (`0.000006`) | `0.003603` (`0.000015`) | `0.0009` | `21.2` | `97.54%` (`93.41%`) | `10.1318` | `0.79480` | `0` |

- **Fixed-instance `H3` caps on `F_2core` (`cap_abort = 0`):** `N_l <= 0.0025 U'` (`1.41x` of `2^24`), `N_sw <= 0.00505 U'` (`1.40x` of `2^24`), `tab_ratio_charged = 1.00` (zero tabular compression discount), `k5_charged = 11.8` (`>= k5_live + 1.50 = 11.6364`, 6-bucket hexagonal `gen_run1` (`29.0` Ir/non-unit term)). Initial singleton purging determines the deterministic surviving 2-core factor base `F_known = F_2core subset F_full` (whose missing density `1 - peel_h` is explicitly charged in `eval_descent_live` via `miss_side`), and `H3`'s zero-column-loss guards conserve every column of `F_2core` through excess trimming and structured Gaussian merging. Live `filtsim_z12` verifies `repurge_backsub_errors == 0`, `merge_backsub_errors == 0`, bit-exact `2^12` parity (`N_sw=31, N_l=17, depth=7, exp_unit=0.712924, lost=0`), and `cols_2core_lost == 0` on algebraic `Q(zeta_15)` relations.
- **1-Owner `< 2^64 B` 4-Vec Arch B BW & `PM-Basis` (`RANK-BRIDGE.md`, `STORAGE.md`):** `n = 112 = 28 * 4`, `m = 8192`, `n_owners = 1`, `te = 0`, `128 B/slot`, double-buffered dense state `2 * N_sw * 128 B = 2^63.688 B < 2^64 B` (zero cross-owner transfers). **LA `= 2^107.869451 <= 2^107.869951`** u. Recursive `PM-Basis`: **Gen `= 2^92.847648 <= 2^92.848148`** u. Schwartz-Zippel `< 2^-66` (`| H-R, H-X`).

## 5. Individual Logarithms (`lx = 71.0`, Per-Environment Stopped-Process `FB`-Blinded CRT `L_{Q, Q'_t}`, Single Tree Cap `T_tree` + `R_tree = 44` Restarts)

Detailed in `DESCENT.md`, `ADAPTIVE.md`, & `KERNELS.md` §5 (`39` bins / `78` nodes, `lx = 71.0`, `boot = (265, 220.0, 16)`, `N_tree_max_1log = 3,376`, `min_mu = 32.15` [`384/370`: `32.76`], `p_bcap = 0.50`, `K_blind = 24`, `min_kappa_crit_pz = 0.9689`, `B_blind` has `192` exact unramified `O_H`-primes in `[2048, 8192]`, `N_blind_pool = 192`):
- **No Local Scan Cutoff:** Every blinding draw scans its full deterministic box `S_16(R)` with zero local per-scan truncation. The sole descent budget cap is `T_tree = 5.0 * W_tree_1`.
- **Adapted Stopping Time & Declared Per-Environment MODEL Inequality `H5-Blind(a)` (Supported, Not Proved):** On the good preprocessing event `G_pre` (`Pr[E not in G_pre] <= p_pre < 2^-50`, where deterministic collection + 2-core filtering determine `F_known = F_2core`), define the adapted stopping time `tau_exc := inf{j >= 1 : p_blind(Q_j; E) > 0.50}` (`1_{j < tau_exc}` is strictly `F_{j-1}`-measurable, determined before slot `j`'s fresh blinding coins!). `H5-Blind(a)` posits pointwise for **every `E in G_pre`** and jointly over `(e', blinding draws)`: `f_exc(E) := Pr_(e', blinding draws)[tau_exc < inf | E] <= eps_exc = 0.12`. Motivated quantitatively (not unconditionally proved at 381-bit scale) by the Law-of-Total-Variance decomposition `kappa(Q) = kappa_geom(Q) + kappa_arith(Q)` (`ADAPTIVE.md` §2; `h5_kappa_tail_bridge`): `prod_sq24_live` directly measures 381-bit `kappa_geom <= 0.1274` (`>= 42.04x` mean headroom -> `mu^-1 < 0.001`), while observed total `kappa_hat` on `(3,3)` (`2,160` balls, `mean = 0.0865, sd = 0.0563, max = 0.2010`, `15.69` SD below `0.9689`) and degree-`(72,24)` `(9,3)` (`432` balls, `mean = 0.2098, sd = 0.0489, max = 0.3163`, `15.53` SD below `0.9689`) give Gamma-Chernoff union bounds across all `3,376` worst-case knapsack slots `<= 3.44e-5 << eps_exc`.
- **Exact Stopped-Process Backward-Induction Closure (ZERO Conditioning on `not E_exc`, Pointwise in `E in G_pre`):** Fix any `E in G_pre`. Defining stopped pre-exceptional work `V_good(Q)` (`0` if `p_blind(Q; E) > 0.50`), backward induction on depth via the tower property (retaining exact unconditioned Newton-Maclaurin at every pre-exceptional slot `j < tau_exc`!) proves `E[V_good | E] <= W_tree_1`, `Pr[{tau_exc = inf} cap {any slot empty} | E] <= P_NM_tree = 0.0002012`, unconditional Markov `Pr[{tau_exc = inf} cap {W_tree > T_tree} | E] <= Pr[V_good > T_tree | E] <= 1/5.0`, and pointwise `W_att = min(W_tree, T_tree) <= V_good + T_tree * 1_{tau_exc < inf}` gives **`E[W_att | E] <= boot_own_exp + (1 + 5.0 * eps_exc) * W_tree_1`**. Hence **pointwise for every `E in G_pre`**, `f(E) <= P_fail_1 = 0.32020 < 0.45`, residual failure after `R_tree = 44` conditionally `i.i.d.` fresh `(e'_r, blinding)` restarts is `f(E)^44 <= P_fail_1^44 <= 1.0000e-18 < 2^-57` (in `fail_core`), and `E[W_desc | E] <= E[W_att | E] / (1 - P_fail_1)` gives expected two-log work `2^109.081096 <= 2^109.082100` u (`384/370`: `2^109.073111` u).

## 6. Success Probability & Unconditional Boole Union Bound

80-digit `Decimal` in `audit_cm.py`: BW SZ `< 2^-66` (`| H-R,H-X`), Schirokauer `< 2^-252` (`| H-SM`), `H3` caps (`0`), ECM Cheb `< 2^-70`, MR `<= 2^-80`, descent `R_tree = 44` residual `<= 1.0000e-18 < 2^-57`, unified meter `t_run = 2^123.0` Markov `p_meter <= 0.000214519364`, degen `<= 6/r`, wrapper `< 2^-48`. Floor: **`0.99978548 >= 999785/1000000`**. Gate `K = 7195777280862110484` satisfies `(K/2^64) * (999785/1000000) >= 0.39 > ((K-1)/2^64) * (999785/1000000)`.

## 7. Complete Expected-Work Ledger (Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Prod-ZeroCopy-U32`, `P4_prod_2ad_prim_w2_317`, `ipc=13.8269`, `c_cof=1.3784`, `m_strict_hi=+0.1171`b, `h2_allow=1.08`) | <= 109.896579 |
| replay `a`-tables build (1-owner 4-bank flat-array `<= 2^63.889 B < 2^64 B`) | <= 105.894517 |
| linear algebra (1-owner `< 2^64 B` 4-vec Arch B `Z[zeta_12]` BW, `n=112, m=8192, te=0`, `nl_cap=0.0025, sw_cap=0.00505, k5=11.8 (hex6), tab_r=1.00`) | <= 107.869951 |
| self-contained classic `PM-Basis` generator (`m=8192, n=112`) | <= 92.848148 |
| both individual logs (`lx=71.0`, `boot=(265,220,16)`, `N_blind_pool=192`, `p_bcap=0.50`, `K_blind=24`, NO local scan cutoff, per-`E` `eps_exc=0.12`, `T_tree=5.0*W_tree_1`, `R_tree=44`, `p_desc_loss<=1.0000e-18`; live `32/185`: `109.081096`, `384/370`: `109.073111`) | <= 109.082096 |
| band split (`<=100.735`), ECM (`<=96.2`), can `a` (`<=86`), 4-col SM (`<=84`), filt (`<=80`), `P_g/P_f/P_71` (`<=75.5`), FB (`<=75`), cert (`<=50`), ver (`<=11`) | <= 100.81 |
| wrapper (§1; `<= 2^80` Ir = `2^65.71` u) | <= 66 |

Core sum: `2^110.80599111 <= 2^110.8062` u. Meter `t_run = 2^123.0` (`p_meter <= 0.000214519364`) behind `K/2^64`: **`E[W] <= 1 + (K/2^64)(2^110.8062 + 2^66) < 2^109.448056241876`** dense-`F_(p^12)` u (`< 2^112.769984336763` G1 u).

## 8. Declared Exploratory Model Inequalities (Supported, Not Proved)

- **`H1/H2` (Collection Yield & Sieve Cost — Declared MODEL Inequality, Supported, Not Proved):** Two-sided `rho_dB` on exact-`O_H`-coprime pure-product pairs (`P4_prod_2ad_prim_w2_317`, `+0.117111`b strict-hi margin = `25.00 * SE`, `+19.49%` over `1.09*U'`), `h2_allow = 1.08`, exact local-root Murphy-`alpha` (`p <= 10,000`, `alpha_f = 2.704 <= 2.74`, `alpha_g_ge3 = 1.103 <= 1.28`, `kap2_g` err `0.0001b`) + live real-FB joint smoothness calibration (`800` pairs, `1.072x`), wrapped in `t_run = 2^123.0`.
- **`H3, H-SM, H-R, H-X` (Fixed-Instance Filtering on `F_2core` & Rank — Declared MODEL Inequality, Supported, Not Proved):** `r` unramified (proved), `P_chi Lambda` injective (`<4/r`), deterministic 2-core zero-column-loss filtering into fixed-instance caps (`nl_cap = 0.0025` = `1.41x`, `sw_cap = 0.00505` = `1.40x`, `k5_charged = 11.8`, `tab_ratio_charged = 1.00`, `cap_abort = 0`, deterministic surviving 2-core factor base `F_known = F_2core` on `G_pre` with missing density `1 - peel_h` explicitly charged in `miss_side`), fixed-instance full rank `H-R` and nonzero BW minors `H-X`.
- **`H4, H4u`:** Raw-opcode-verified & bit-exact-simulated baseline scalar ARMv8.0-A integer kernels (`KERNELS.md`, `simd = 0`, ZERO FP/LSE); `Z[zeta_12]` unit share `93.41%` at `2^24`.
- **`H5-Blind(a)` (Per-Environment Stopped-Process Special-`q` Blinding Descent — Declared MODEL Inequality, Supported, Not Proved):** Pairwise-`O_H`-coprime `FB`-blinded CRT `L_{Q, Q'_t}` (`p_sat = 0` proved, `192` exact unramified primes in `[2048, 8192]`), NO local scan cutoff, adapted stopping time `tau_exc := inf{j >= 1 : p_blind(Q_j; E) > 0.50}` with pointwise-in-`E in G_pre` bound `f_exc(E) := Pr_(e', blinding)[tau_exc < inf | E] <= eps_exc = 0.12`, closed via exact stopped-process backward induction on `V_good` (`E[W_att | E] <= boot_own_exp + (1 + 5.0*eps_exc)*W_tree_1`, ZERO conditioning on `not E_exc`), `T_tree = 5.0 * W_tree_1`, and `R_tree = 44` conditionally `i.i.d.` restarts (`p_desc_loss <= 1.0000e-18 < 2^-57`), motivated quantitatively by the Law-of-Total-Variance & Gamma-Chernoff bridge (`ADAPTIVE.md` §2; `h5_kappa_tail_bridge`).

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0 B` (`2^119.074 B` replay `a`-tables; 1-owner build `<= 2^63.889 B < 2^64 B`; 1-owner 4-vec Arch B LA `2^63.688 B < 2^64 B`; `P_71` across 128 owners `< 2^62 B`/owner; `STORAGE.md`). **Provenance:** Builds on Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, Boudot (2017/758). **No prior verdict or score is a premise.**
