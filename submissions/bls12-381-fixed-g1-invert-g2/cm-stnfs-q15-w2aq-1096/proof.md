# recover-Q: `(4,3)` `Q(zeta_15)` CM-STNFS + `Arch2-Prod-ZeroCopy-U16` 1-Owner Build + 1-Owner 4-Vec Arch B LA + Per-Environment `H5-Blind(a,b)` (`R_tree = 44`): expected work `2^109.610993760501`

Track `bls12-381-fixed-g1-invert-g2`. Input: nonidentity `P in G1`, `z in G_T`; output `Q in G2` with `e(P, Q) = z` for exact supplied `P` (`Q = infinity` for `z = 1`). Unit: one dense `F_(p^12)` mul, `C = 20,000` baseline scalar ARMv8.0-A integer instructions (`simd = 0`, zero FP/NEON/SVE, zero `FEAT_LSE`/`DC ZVA`, zero external C measurements, zero black-box generator citations, zero amortization). Full `rho_dB` baked into Collection & Descent. All kernels (`KERNELS.md`, parsed & bit-exact-simulated baseline scalar ARMv8.0-A integer blocks), exact local-root Murphy-`alpha` & odd-integer `kap2_g` derivation (`p <= 10,000`), real-FB joint smoothness calibration (`800` coprime pairs), zero-column-loss `Z[zeta_12]` structured Gaussian filtering (`2^16..2^24` x 6 seeds, `post_peel_cols_lost_max == 0`, `post_merge_cols_lost_max == 0`), 1-owner `< 2^64 B` 4-vector Arch B Block Wiedemann LA (`n = 112`, `m = 8192`, `te = 0`), per-environment pairwise-`O_H`-coprime `FB`-blinded CRT descent (`lx = 71.0`, NO local per-scan cutoff, single deterministic tree cap `T_tree = 4.5*W_tree_1`, `R_tree = 44` conditionally `i.i.d.` fresh-`e'` restarts given `E in G_pre`), live actual-`U_BLS` `(9,3)` CRT-blinded counting, live 100% per-parent `p_hat(Q)` audit on 8 depth-5 adaptive trees, pure `Decimal(50)` `unitdomain.py`, and direct 200k-draw exact-`O_H`-coprime pure-product pool are self-contained (**no prior verdict or score is a premise**).

## Claim

`W` is total classical work on one target. Step 0: uniform 64b `u`; unless `u < K = 7195964416171022347` (least `K` with `K * (999759/1000000) >= 0.39 * 2^64`), abort (`< 1` u). Under declared exploratory STNFS model inequalities `H1/H2, H3, H-R, H-SM, H4, H4u, H-X, H5-Blind(a,b)` (§8; `RANK-BRIDGE.md`, `DESCENT.md`, `ADAPTIVE.md`), an opened run succeeds w.p. `>= 999759/1000000` (Boole, §6); overall success `>= 0.39` over uniform production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^110.9691 + 2^66) < 2^109.610993760501  dense-Fp12 units (C = 20,000).

`109.610993760501` is an **expected** total work bound `E[W]` (NOT every-run). The **ENTIRE opened run** runs under ONE hard aborting work meter `t_run = 2^123.0` u (block boundaries `<= 2^66` u): deterministic pointwise cap `W_actual <= 2^123.0 + 2^66` u, `E[min(W_open, t_run)] <= 2^110.9691`, Markov `p_meter <= 0.000240161998`. Analytic exploratory claim; all numbers match default `32/185` and refined `384/370` `audit_cm.py` (`checks_ok = true`, pure stdlib).

## 1. Reduction (fixed `P in G1`, recover `Q in G2`: two `G_T` logs) and field map

With pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `Q0 in G2`:
0. **Run gate:** uniform 64b `u`; if `u >= K`, abort (`< 1` u).
1. **Validate:** `P in G1 \ {infinity}` (`[r]P = infinity`), `z^r = 1`; if `z = 1`, return `infinity`.
2. **Pairing & map:** `g = e(P, Q0)`, map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Randomize & solve:** independent uniform `e_1, e_2 in {0..r-1}`, core on `g^(e_1), z^(e_2)`, `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** `Q = [a]Q0`; return `Q` iff `e(P, Q) = z`. Wrapper `< 2^66` u, failure `< 2^-48`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, Exact Murphy-`alpha`, & Order-12 Frobenius

BLS12-381: `u = -0xd201000000010000`, `r = u^4 - u^2 + 1 == 13 mod 15`, `p = (u-1)^2 r/3 + u == 7 mod 15`. Over PID `O_H = Z[zeta_15]` (`D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `O_H/P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` gives `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega+1)X^6 + X^3 - omega` (deg 72). `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`ord_r(p) = 12`) quotients into `Z[zeta_12]` 12-orbits (`algebra_bls12`). `check_algebra_bls12` verifies via exact polynomial GCD that `[2048, 8192]` contains `192` distinct unramified `O_H`-primes (`ram = 0`) admitting `f`-side degree-1 ideals for `B_blind` (`N_blind_pool = 192`). Pre-collection `W_f in F_r^(3x6), W_g in F_r^(1x2)` compress `P_chi Lambda` to `4` cols (`< 4/r` under `H-SM`).
- **Exact Local-Root Murphy-`alpha` & Odd-Integer Shift (`murphy_alpha_exact_live`):** `audit_cm.py` computes exact local valuation excess `sum_p (1/(p-1) - E[v_p(N)]) * log2(p)` across ramified `p in {2, 3, 5}` (`p=2`: `0` roots of `f` over `F_16`; `p=3`: exact `F_81` Hensel lifting giving `E[v_3(N_f)]=8/41, E[v_3(N_g)]=2/41`; `p=5`: `0` roots over `F_25`) and all `1,224` unramified primes `7 <= p <= 10,000`: **`alpha_f(<=10k) = 2.70436 <= 2.74`**, **`alpha_(g,>=3)(<=10k) = 1.10299 <= 1.28`** (`+0.2127b` conservative norm cushion, truncation yield bound `< 0.015b`), and exact odd-integer de Bruijn identity `1.0 + kappa_2 = 1.0525` (`max_err = 0.000102b`).
- **Live Real-FB Joint Smoothness Bridge (`real_fb_joint_coll_live`):** Direct trial-factorization smoothness check on `800` exact-`O_H`-coprime pairs (`B_fb=512, LP=8192`) confirms `obs_joint = 53` vs `mod_joint = 49.45` (`obs/mod = 1.0719`, `+0.1001b`).
- Live E2E (`toy_e2e_dlp`, `u_0=-5, p_0=7207, r_0=601`): `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps, `5/5` unseen DLPs.

## 3. Relation Collection: `Arch2-Prod-ZeroCopy-U16` (1-Owner Flat-Array Build) & Direct 200k Exact-`O_H`-Coprime Pool

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373`), excess `1.09` (`2^63.441701` rels).
- **Pure-product geometry & 2-adic `c(a)` (`239/239` exact-1):** `G = A_4 x A_2`, 10 Voronoi units (`f_Vor = 3/(5*sqrt(5)) = 0.268328`), pure product `q_a(a) <= s_{c(a)} R2_a`, `q_b(b) <= R2_b` (`A^2 = 321600.0`, `r_(ab) = 0.9815`, `s_c = (1.0, 1.5, 2.0)`, `STORAGE.md`). `N_a^can = 2^56.8437 < 0.90 * 2^57`, `L_max = 404.76 < W_tile = 416`, `V_kept = 2^120.4314`.
- **Direct 200k-Draw Exact-`O_H`-Coprime Pool (`P4_prod_2ad_prim_v17`, `n_draws = 200,000`):** Direct `gcd_O(a,b) == (1)` test on `200,000` draws (`f_prim = 0.95642`, `se_y_bits = 0.005654`). Retained yield: midpoint `2^63.859129`, upper-ret `2^63.853894` (**`margin_net = +0.403193`b**, `+32.24%` over `1.09*U'`), and **strict upper-endpoint smoothness + retention `2^63.714125` (`margin_strict_hi = +0.272424`b = `48.18 * SE`)**.
- **`Arch2-Prod-ZeroCopy-U16` (`W_tile=416, W_SB=32, P_m4=2^56.5, P_lo=P_f=2^59, P_g=2^64`):**
  - **1-owner simultaneous build (`<= 2^63.924 B < 2^64 B`):** `2^105.958590 <= 2^105.95909` u, total replay tables **`2^119.121 B < 2^120 B`**.
  - **Replay & Cofactorization:** Band-1a `rev_b1a_apply_u4` (`4.25`), Band-1b `rev_b1b_apply_u4` (`6.25`), Band-2 `rev_b2_prod_emit_1` + `u4` (`16.75`, `c_bd_chk=3.0`), `rev_scan_zero_u16` (`1.1`), fixed `g`-first cofactorization (`c_cof_2d_upper = 2.0535` ipc). Full-cell sieve (`ipc_tot = 14.8197`, `h2_allow = 1.08`) -> **`2^110.144204 <= 2^110.144704`** u.

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

- **Fixed-instance `H3` caps (`cap_abort = 0`):** `N_l <= 0.00285 U'` (`1.60x` of `2^24`), `N_sw <= 0.0058 U'` (`1.61x` of `2^24`), `tab_ratio_charged = 1.00` (zero tabular compression discount), `k5_charged = 11.8` (`>= k5_live + 1.50 = 11.6364`, 6-bucket hexagonal `gen_run1` (`29.0` Ir/non-unit term)). Under zero-column-loss `H3`, `F_known = F_full` is the unique deterministic full factor base on `G_pre`. Live `filtsim_z12` verifies `repurge_backsub_errors == 0`, `merge_backsub_errors == 0`, bit-exact `2^12` parity (`N_sw=31, N_l=17, depth=7, exp_unit=0.712924, lost=0`), and `cols_2core_lost == 0` on algebraic `Q(zeta_15)` relations.
- **1-Owner `< 2^64 B` 4-Vec Arch B BW & `PM-Basis` (`RANK-BRIDGE.md`, `STORAGE.md`):** `n = 112 = 28 * 4`, `m = 8192`, `n_owners = 1`, `te = 0`, `128 B/slot`, double-buffered dense state `2 * N_sw * 128 B = 2^63.888 B < 2^64 B` (zero cross-owner transfers). **LA `= 2^108.258013 <= 2^108.258513`** u. Recursive `PM-Basis`: **Gen `= 2^93.066595 <= 2^93.067095`** u. Schwartz-Zippel `< 2^-66` (`| H-R, H-X`).

## 5. Individual Logarithms (`lx = 71.0`, Per-Environment `FB`-Blinded CRT `L_{Q, Q'_t}`, Single Tree Cap `T_tree` + `R_tree = 44` Restarts)

Detailed in `DESCENT.md`, `ADAPTIVE.md`, & `KERNELS.md` §5 (`39` bins / `78` nodes, `lx = 71.0`, `boot = (265, 220.0, 16)`, `N_tree_max_1log = 3,376`, `min_mu = 32.15` [`384/370`: `32.76`], `p_bcap = 0.50`, `K_blind = 24`, `min_kappa_crit_pz = 0.9689`, `B_blind` has `192` exact unramified `O_H`-primes in `[2048, 8192]`, `N_blind_pool = 192`):
- **No Local Scan Cutoff:** Every blinding draw scans its full deterministic box `S_16(R)` with zero local per-scan truncation. The sole descent budget cap is `T_tree = 4.5 * W_tree_1`.
- **Declared Per-Environment Conditional MODEL Inequality `H5-Blind(a, b)` (Supported, Not Proved):** On the good preprocessing event `G_pre` (`Pr[E not in G_pre] <= p_pre < 2^-50`, where `H3` yields the unique deterministic full factor base `E*`), `H5-Blind(a)` posits pointwise for **every `E in G_pre`** that `f_exc(E) := Pr_(e' ~ Unif(Z_r^*))[E_exc(h*g^(e'), E) | E] <= eps_exc = 0.12`; `H5-Blind(b)` posits `p_blind(Q; E) <= 0.50` on non-exceptional trees. Supported (not unconditionally proved at 381-bit scale) by: (i) live 100% per-parent `p_hat(Q)` audit across ALL `43/43` adaptively encountered parents on 8 live `q ~ 2^24` depth-5 trees (`0/8` trees exceeding `0.50` or `0.35`, `max_Q p_hat <= 0.1667`), (ii) `2,160` observed `(3,3)` CRT balls (`kappa_hat_max_Q = 0.2010`), (iii) **`32` live degree-`(72,24)` `(9,3)` CRT balls using ACTUAL production `U_BLS = -0xd201000000010000` (`0/32` empty, `kappa_hat_max_Q = 0.0000` sub-Poisson!)**, (iv) `432` high-mean `u=-5` `(9,3)` CRT balls (`0/432` empty, `kappa_hat_max_Q = 0.3163 < 0.9689`), and (v) 48 live 381-bit LLL lattice-geometry stability checks (`ADAPTIVE.md` §2).
- **Exact Pointwise-in-`E` Tree-Level Closure (`R_tree = 44` Conditionally `i.i.d.` Restarts Given `E`):** Fix any `E in G_pre`. Newton-Maclaurin proves `P(any slot empty | not E_exc, E) <= P_NM_tree = 0.0002012` AND `E[W_tree | not E_exc, E] <= W_tree_1`; Markov gives `P(W_tree > 4.5 * W_tree_1 | not E_exc, E) <= 1/4.5`. Hence **pointwise for every `E in G_pre`**, `f(E) <= P_fail_1 = 0.34242 < 0.45`, residual failure after `R_tree = 44` conditionally `i.i.d.` fresh-`e'` restarts is `f(E)^44 <= P_fail_1^44 <= 1.0000e-18 < 2^-57` (explicitly in `fail_core`), and `E[W_desc | E] <= E[W_attempt | E] / (1 - P_fail_1)` gives expected two-log work `2^108.970155 <= 2^108.971200` u (`384/370`: `2^108.961237` u), with ZERO `E[f(E)^R]` vs `(E[f(E)])^R` conditioning gap.

## 6. Success Probability & Unconditional Boole Union Bound

80-digit `Decimal` in `audit_cm.py`: BW SZ `< 2^-66` (`| H-R,H-X`), Schirokauer `< 2^-252` (`| H-SM`), `H3` caps (`0`), ECM Cheb `< 2^-70`, MR `<= 2^-80`, descent `R_tree = 44` residual `<= 1.0000e-18 < 2^-57`, unified meter `t_run = 2^123.0` Markov `p_meter <= 0.000240161998`, degen `<= 6/r`, wrapper `< 2^-48`. Floor: **`0.99975983 >= 999759/1000000`**. Gate `K = 7195964416171022347` satisfies `(K/2^64) * (999759/1000000) >= 0.39 > ((K-1)/2^64) * (999759/1000000)`.

## 7. Complete Expected-Work Ledger (Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Prod-ZeroCopy-U16`, `P4_prod_2ad_prim_v17`, `ipc=14.8197`, `c_cof=2.0535`, `m_strict_hi=+0.2724`b, `h2_allow=1.08`) | <= 110.144704 |
| replay `a`-tables build (1-owner 4-bank flat-array `<= 2^63.924 B < 2^64 B`) | <= 105.95909 |
| linear algebra (1-owner `< 2^64 B` 4-vec Arch B `Z[zeta_12]` BW, `n=112, m=8192, te=0`, `nl_cap=0.00285, sw_cap=0.0058, k5=11.8 (hex6), tab_r=1.00`) | <= 108.258513 |
| self-contained classic `PM-Basis` generator (`m=8192, n=112`) | <= 93.067095 |
| both individual logs (`lx=71.0`, `boot=(265,220,16)`, `N_blind_pool=192`, `p_bcap=0.50`, `K_blind=24`, NO local scan cutoff, per-`E` `eps_exc=0.12`, `T_tree=4.5*W_tree_1`, `R_tree=44`, `p_desc_loss<=1.0000e-18`; live `32/185`: `108.970155`, `384/370`: `108.961237`) | <= 108.971155 |
| band split (`<=100.8898`), ECM (`<=96.2`), can `a` (`<=86`), 4-col SM (`<=84`), filt (`<=80`), `P_g/P_f/P_71` (`<=75.5`), FB (`<=75`), cert (`<=50`), ver (`<=11`) | <= 100.97 |
| wrapper (§1; `<= 2^80` Ir = `2^65.71` u) | <= 66 |

Core sum: `2^110.96887058 <= 2^110.9691` u. Meter `t_run = 2^123.0` (`p_meter <= 0.000240161998`) behind `K/2^64`: **`E[W] <= 1 + (K/2^64)(2^110.9691 + 2^66) < 2^109.610993760501`** dense-`F_(p^12)` u (`< 2^112.932921855389` G1 u).

## 8. Declared Exploratory Model Inequalities (Supported, Not Proved)

- **`H1/H2` (Collection Yield & Sieve Cost — Declared MODEL Inequality, Supported, Not Proved):** Two-sided `rho_dB` on exact-`O_H`-coprime pure-product pairs (`P4_prod_2ad_prim_v17`, `+0.272424`b strict-hi margin = `48.18 * SE`, `+32.24%` over `1.09*U'`), `h2_allow = 1.08`, exact local-root Murphy-`alpha` (`p <= 10,000`, `alpha_f = 2.704 <= 2.74`, `alpha_g_ge3 = 1.103 <= 1.28`, `kap2_g` err `0.0001b`) + live real-FB joint smoothness calibration (`800` pairs, `1.072x`), wrapped in `t_run = 2^123.0`.
- **`H3, H-SM, H-R, H-X` (Fixed-Instance Filtering & Rank — Declared MODEL Inequality, Supported, Not Proved):** `r` unramified (proved), `P_chi Lambda` injective (`<4/r`), deterministic zero-column-loss filtering into fixed-instance caps (`nl_cap = 0.00285` = `1.60x`, `sw_cap = 0.0058` = `1.61x`, `k5_charged = 11.8`, `tab_ratio_charged = 1.00`, `cap_abort = 0`, deterministic full factor base `F_known = F_full` on `G_pre`), fixed-instance full rank `H-R` and nonzero BW minors `H-X`.
- **`H4, H4u`:** Raw-opcode-verified & bit-exact-simulated baseline scalar ARMv8.0-A integer kernels (`KERNELS.md`, `simd = 0`, ZERO FP/LSE); `Z[zeta_12]` unit share `93.41%` at `2^24`.
- **`H5-Blind(a, b)` (Per-Environment Conditional Special-`q` Blinding Descent — Declared MODEL Inequality, Supported, Not Proved):** Pairwise-`O_H`-coprime `FB`-blinded CRT `L_{Q, Q'_t}` (`p_sat = 0` proved, `192` exact unramified primes in `[2048, 8192]`), NO local scan cutoff, pointwise-in-`E in G_pre` `H5-Blind(a)` exceptional-tree fraction `f_exc(E) <= eps_exc = 0.12`, `H5-Blind(b)` non-exceptional slot availability `p_blind(Q; E) <= 0.50` (proving `E[W_tree | not E_exc, E] <= W_tree_1`), closed pointwise in `E` via `T_tree = 4.5 * W_tree_1` and `R_tree = 44` conditionally `i.i.d.` fresh-`e'` restarts (`p_desc_loss <= 1.0000e-18 < 2^-57`), supported by Layers 1–4 in `ADAPTIVE.md`.

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0 B` (`2^119.121 B` replay `a`-tables; 1-owner build `<= 2^63.924 B < 2^64 B`; 1-owner 4-vec Arch B LA `2^63.888 B < 2^64 B`; `P_71` across 128 owners `< 2^62 B`/owner; `STORAGE.md`). **Provenance:** Builds on Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, Boudot (2017/758). **No prior verdict or score is a premise.**
