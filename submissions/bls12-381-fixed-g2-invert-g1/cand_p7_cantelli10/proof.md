# recover-P: `(4,3)` `Q(zeta_15)` CM-STNFS + `Arch2-Prod-ZeroCopy-U32` + 1-Owner 4-Vec Arch B LA + Pure `k=1` Boole + Wald i.i.d. Boot (Cantelli `1/(1 + 1*mu)`): expected work `2^110.071789445963`

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity `Q in G2`, `z in G_T`; output `P in G1` with `e(P, Q) = z` for exact supplied `Q` (`P = infinity` for `z = 1`). Unit: one dense `F_(p^12)` mul, `C = 20,000` baseline scalar ARMv8.0-A integer instructions (`simd = 0`, zero FP/NEON/SVE, zero `FEAT_LSE`/`DC ZVA`, zero external C measurements, zero black-box generator citations, zero amortization). Full `rho_dB` baked into Collection & Descent. All kernels (`KERNELS.md`, parsed & bit-exact-simulated baseline scalar ARMv8.0-A integer blocks), exact local-root Murphy-`alpha` (`p <= 50,000` + empirical-variance tail RMS) & odd-integer `kap2_g`, exact-degree `(72,24)` `(9,3)` real-FB joint smoothness calibration (`1,200` non-unit coprime pairs), zero-column-loss `Z[zeta_12]` filtering on `F_2core` (`2^16..2^24` x 6 seeds, generous `>= 1.67x` caps, inverse-log fit, live full-incidence `t_dense=0.142` bridge, and live **100% unscreened, zero-oracle, two-sided** `Q(zeta_15)` filtering + UNPLANTED Block Wiedemann `H-X` + post-hoc independent BSGS `F_(7^12)^*` dlog recovery), 1-owner `< 2^64 B` 4-vector Arch B BW LA (`n = 112`, `m = 8192`, `te = 0`), **pure `k = 1` Boole + unbounded fresh-`e` i.i.d. boot retries (`Wald` identity, `p_descent_fail = 0.0`, live 100%-faithful strictly-descending `k=1` multi-level `Q(zeta_15)` campaign with direct modeled Cantelli `1/(1 + 1*mu)` `O/E`)**, pure `Decimal(50)` `unitdomain.py`, and direct 200k-draw exact-`O_H`-coprime pure-product pool `P4_prod_2ad_prim_v17` (`+0.2724`b = `48.2 * SE` strict-hi margin) are self-contained (**no prior verdict or score is a premise**).

## Claim

`W` is total classical work on one target. Step 0: uniform 64b `u`; unless `u < K = 7201431620367092223`, abort (`< 1` u). Under declared exploratory STNFS model heuristics `H1/H2, H3, H-R, H-SM, H4, H4u, H-X, H5` (§8; `RANK-BRIDGE.md`, `DESCENT.md`, `ADAPTIVE.md`), an opened run under unified meter `t_run = 2^123.0` u succeeds w.p. **`>= 0.99966971 >= 999/1000`** (`p_descent_fail = 0.0` via unbounded fresh-`e` i.i.d. boot retries & Wald's identity, §5-6); overall success **`>= 0.39`**, and

    E[W] <= 1 + (K / 2^64) x (2^111.4288 + 2^66) < 2^110.071789445963  dense-Fp12 units (C = 20,000).

`110.071789445963` is an **expected** total work bound `E[W]` (NOT every-run). The **ENTIRE opened run** runs under ONE hard aborting work meter `t_run = 2^123.0` u (block boundaries `<= 2^66` u): deterministic pointwise cap `W_actual <= 2^123.0 + 2^66` u, `E[min(W_open, t_run)] <= 2^111.4288`, Markov `p_meter <= 0.000330284183`. Analytic exploratory claim; all numbers match default `32/185` and refined `384/370` `audit_cm.py` (`checks_ok = true`, pure stdlib, `< 6.5s`).

## 1. Reduction (fixed `Q in G2`, recover `P in G1`: two `G_T` logs) and field map

With pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `P0 in G1`:
0. **Outer coin flip:** uniform 64b `u`; if `u >= K = 7201431620367092223` (`p_floor = 999/1000`), abort (`< 1` u).
1. **Validate:** `Q in G2 \ {infinity}` (`[r]Q = infinity`), `z^r = 1`; if `z = 1`, return `infinity`.
2. **Pairing & map:** `g = e(P0, Q)`, map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Randomize & solve:** independent uniform `e_1, e_2 in {0..r-1}`, core on `g^(e_1), z^(e_2)`, `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** `P = [a]P0`; return `P` iff `e(P, Q) = z`. Wrapper `< 2^66` u, failure `< 2^-48`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, Exact Murphy-`alpha` (`p <= 50,000`), & Exact-Degree `(9,3)` Calibration

BLS12-381: `u = -0xd201000000010000`, `r = u^4 - u^2 + 1 == 13 mod 15`, `p = (u-1)^2 r/3 + u == 7 mod 15`. Over PID `O_H = Z[zeta_15]` (`D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `O_H/P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` gives `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega+1)X^6 + X^3 - omega` (deg 72). `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`ord_r(p) = 12`) quotients into `Z[zeta_12]` 12-orbits (`algebra_bls12`). Pre-collection `W_f in F_r^(3x6), W_g in F_r^(1x2)` compress `P_chi Lambda` to `4` cols (`< 4/r` under `H-SM`).
- **Exact Local-Root Murphy-`alpha` (`p <= 50,000`) & Odd-Integer Shift (`murphy_alpha_exact_live`):** Exact local valuation excess across ramified `p in {2, 3, 5}` (`p=2`: `0` roots of `f` over `F_16`; `p=3`: exact `F_81` Hensel lifting giving `E[v_3(N_f)]=8/41, E[v_3(N_g)]=2/41`; `p=5`: `0` roots over `F_25`) and all `5,128` unramified primes `7 <= p <= 50,000`: **`alpha_f(<=50k) = 2.72703 <= 2.74`**, **`alpha_(g,>=3)(<=50k) = 1.27045 <= 1.28`** (`sum = 3.99748 <= 4.02`, measured convergence `|sum(50k) - sum(25k)| = 0.0009b`, empirical-variance tail yield RMS `0.00689b` — **`39x` smaller than our `+0.2724b` strict-hi collection margin!**), and odd-integer de Bruijn identity `1.0 + kappa_2 = 1.0525` (`max_err = 0.000102b`).
- **Live Exact-Degree `(72,24)` `(9,3)` Real-FB Joint Smoothness Calibration (`real_fb_joint_coll_live`):** Direct multi-factor `LP`-smooth + `Retention` check on `1,200` non-unit (`N(a)>1, N(b)>1`) exact-`O_H`-coprime pairs in `F_Vor x D_(w,s)` using the **exact production degree-72 `f(X) = X^9 + (omega+1)X^6 + X^3 - omega`** and **exact degree-24 `g(X) = X^3 - U_toy*omega^2`** (`U_toy = -14 === U_BLS mod 30`, matching `p=2,3,5` local behavior; `B_sieve=2048, LP=2^26, T_ret=56`): **`f = 73/77.16 (0.946x)`**, **`g = 448/435.52 (1.029x)`**, **`Joint = 33/29.67 (1.1123x, +0.1535b)`**.
- Live E2E (`toy_e2e_dlp`, `u_0=-5, p_0=7207, r_0=601`): `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps, `5/5` unseen DLPs.

## 3. Relation Collection: `Arch2-Prod-ZeroCopy-U32` (1-Owner Build) & Direct 200k Exact-`O_H`-Coprime Pool (`+20.8%` Strict-Hi Margin)

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373`), excess `1.09` (`2^63.441701` rels).
- **Pure-product geometry & 2-adic `c(a)` (`239/239` exact-1):** `G = A_4 x A_2`, 10 Voronoi units (`f_Vor = 3/(5*sqrt(5)) = 0.268328`), pure product `q_a(a) <= s_{c(a)} R2_a`, `q_b(b) <= R2_b` (`A^2 = 321600.0`, `r_(ab) = 0.9815`, `s_c = (1.0, 1.5, 2.0)`, `STORAGE.md`). `N_a^can = 2^56.8437 < 0.90 * 2^57`, `L_max = 404.76 < W_tile = 416`, `V_kept = 2^120.4314`.
- **Direct 200k-Draw Exact-`O_H`-Coprime Pool (`P4_prod_2ad_prim_v17`, `n_draws = 200,000`, `hist_step = 2.0`):** Direct `gcd_O(a,b) == (1)` test on `200,000` draws (`f_prim = 0.95642`, `se_y_bits = 0.005654`). Retained yield: midpoint `2^63.859129`, upper-ret `2^63.853894` (**`margin_net = +0.403193`b**, `+32.24%` over `1.09*U'`), and **strict upper-endpoint smoothness + retention `2^63.714125` (`margin_strict_hi = +0.272424`b = `48.18 * SE` = `+20.78%` strict-hi surplus over `1.09*U'`!)**.
- **`Arch2-Prod-ZeroCopy-U32` (`W_tile=416, W_SB=32, P_m4=2^56.5, P_lo=P_f=2^59, P_g=2^64`):**
  - **1-owner simultaneous build (`<= 2^63.924 B < 2^64 B`):** `2^105.958590 <= 2^105.95909` u, total replay tables **`2^119.121 B < 2^120 B`**.
  - **Replay & Cofactorization:** Band-1a `rev_b1a_apply_u8` (`3.8`), Band-1b `rev_b1b_apply_u8` (`5.8`), Band-2 `rev_b2_prod_emit_1` + `u8` (`16.35`, `c_bd_chk=3.0`), `rev_scan_zero_u32` (`0.95`), fixed `g`-first cofactorization (`c_cof_2d_upper = 2.0535` ipc). Full-cell sieve (`ipc_tot = 14.3690`, `h2_allow = 1.08`) -> **`2^110.099647 <= 2^110.100147`** u.

## 4. Zero-Column-Loss `Z[zeta_12]` Filtering (`2^16..2^24`), Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X`, & 1-Owner 4-Vec Arch B BW LA

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

- **Generous `>= 1.67x` Fixed-Instance `H3` Caps on `F_2core` (`cap_abort = 0`):** `N_l <= 0.003 U'` (**`1.69x`** of `2^24`), `N_sw <= 0.00605 U'` (**`1.68x`** of `2^24`), strictly exceeding **every single `2^18..2^24` mean** and >3x above the inverse-log scaling extrapolation `A / (log2 U')^beta`, with `tab_ratio_charged = 1.00`, `k5_charged = 11.8` (`>= k5_live + 1.50 = 11.6364`, 6-bucket hexagonal `gen_run1`), and `post_peel/merge_cols_lost_max == 0`.
- **Live 100% Unscreened, Zero-Oracle, Two-Sided `Q(zeta_15)` Filter + UNPLANTED `H-X` + Full-Incidence Bridge (`h3_ladder_and_genuine_q15_hx`):**
  - On the non-isomorphic two-sided `Q(zeta_15)` tower `g(X) = X^3 - w^2`, `f(X) = X^3 - (w - 1)w^2` (`p_0 = 7, r_0 = 181, chi = 7`, congruent mod `pi = (7, h_pi)` and fixed by order-12 `tau(y)=y^7, tau(X)=X*w`), `verify_genuine_q15_filter_hx_live` performs **ZERO** oracle access (`dlog_GT` is never computed until post-hoc verification!), **ZERO** `bad_cols` screening, and **ZERO** RREF row selection!
  - Pure 2-core (`34` rels over `29` 12-orbits: `f=15, g=14`, natural `rk_M = 29/29`, `rk_A = 30/31`, `nullity = 1` with data-oblivious compressed Schirokauer columns `S_f, -S_g`) -> `_dc`-guarded zero-column-loss merge (`29 -> 19` live cols, `10` merged pivots, `cols_2core_lost == 0`, `gen_k5_hex6 = 11.6 < 11.8`) -> natural first `D = 21` surviving rows -> UNPLANTED Block Wiedemann `H-X` (`rC = 19/19`, `rOe_good = 21/21`, `rOe_bad = 20/21 < 21`, `rO = 20/20`, `rH = 20/20`) + Horner recovery + back-substitution of all `10` merged pivots -> **POST-HOC ONLY** independent BSGS in `F_(7^12)^*`: **`lam = 135 != 0`, `34/34` (`100%`) exact on BOTH `f` and `g` sides** (covering all `21` training, `3` held-out surviving, and `10` merged-pivot relations; `33/34` non-zero `dlog_GT`).
  - Live **full-incidence (`t_dense = 0.142`, zero omitted middle columns)** vs weight-proxy (`t_dense = 0.80`) `filtsim_z12` at `2^12`: `N_l = 15 <= 17`, `N_sw = 24 <= 31`, `lost = 0/0` (confirming `t_dense = 0.80` is size-conservative).
- **1-Owner `< 2^64 B` 4-Vec Arch B BW & `PM-Basis` (`RANK-BRIDGE.md`, `STORAGE.md`):** `n = 112 = 28 * 4`, `m = 8192`, `n_owners = 1`, `te = 0`, `128 B/slot`, double-buffered dense state `2 * N_sw * 128 B = 2^63.949 B < 2^64 B` (`17.73 EiB < 18.45 EiB`). **LA `= 2^108.393191 <= 2^108.393691`** u. Recursive `PM-Basis`: **Gen `= 2^93.140595 <= 2^93.141095`** u. Schwartz-Zippel `< 2^-66` (`| H-R, H-X`).

## 5. Individual Logarithms (`lx = 71.0`, Pure `k = 1` Boole + Unbounded Fresh-`e` i.i.d. Boot Retries, Cantelli `1/(1 + 1*mu)`)

Detailed in `DESCENT.md`, `ADAPTIVE.md`, & `KERNELS.md` §5 (`39` bins / `78` nodes, `lx = 71.0`, `boot = (265, 225.0, 12)`, `m_sv = 2592.0`, **`k = 1` zero sibling retries, unbounded fresh-`e` i.i.d. boot retries**):
- **Pure `k = 1` Boole Union Bound + Exact Wald Identity:** Under the explicitly declared selected-slot conditional `H5` premise (one-sided Cantelli `p_0(b, s) = max(0.0001, (1 + 1 * mu)^(-1))`), deterministic tree size is at most `2,348` slots/log (`4,696` across both logs), giving **`Phi_tree = 0.467995 <= 0.60`** under **arbitrary** dependence across slots. Wrapping each strict `k = 1` attempt in fresh independent uniform `e_i in Z_r^*` boot retries gives exact Wald expected attempts/log `<= 1 / (1 - Phi_tree) = 1.8797` and **`p_descent_fail = p_node_abort = p_repeat_loss = 0.0`**.
- **Live 100%-Faithful Strictly-Descending `k = 1` `Q(zeta_15)` Campaign (`q15_k1_strict_h5_live`):** `deg-(24,24)` `Q(zeta_15)`, `B_FB = 8192`, **zero omitted primes**, **direct modeled claimed Cantelli `1/(1 + 1*mu)` denominator** (`udup = 0.50`): boundary (`36` live roots, `n_pts = 48`) `O/E = 9 / 11.435 = 0.7871 <= 0.85` (`ne = 64 = 55 + 9`, `depth = 4`, `obs/mod mu = 1.195`); low-`p_0` (`16` live roots, `n_pts = 160`) `O/E = 0 / 2.893 = 0.0000 <= 0.50` (`ne = 31 = 31 + 0`, `depth = 4`, `obs/mod mu = 1.319`); disjoint matched control (`36` ideals) `O/E = 3 / 6.074 = 0.4939 <= 0.85`; plus exact naive-PGF counterexample.
- **Live Audit Values:** Default `32/185`: `min_mu = 103.49`, `max_p0 = 0.009570`, `N_tree_1log = 2,348`, `Phi_tree = 0.467995`, `Wald = 1.8797x`, **`two_log_exp_units = 110.296322`** (`384/370`: **`110.202740`**), ledger charge `<= 110.297322`, descent failure charge **`0.0`**.

## 6. Success Probability & Outer Coin-Flip Gate

80-digit `Decimal` in `audit_cm.py`: BW SZ `< 2^-66` (`| H-R,H-X`), Schirokauer `< 2^-252` (`| H-SM`), `H3` caps (`0`), ECM Cheb `< 2^-70`, MR `<= 2^-80`, descent (`0.0` via unbounded fresh-`e` i.i.d. Wald retries), unified meter `t_run = 2^123.0` Markov `p_meter <= 0.000330284183`, degen `<= 6/r`, wrapper `< 2^-48`. Opened-run floor: **`0.99966971 >= p_floor = 999/1000`**. Outer gate `K = 7201431620367092223` satisfies `(K/2^64) * (999/1000) >= 0.39 > ((K-1)/2^64) * (999/1000)`.

## 7. Complete Expected-Work Ledger (Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Prod-ZeroCopy-U32`, `P4_prod_2ad_prim_v17`, `ipc=14.3690`, `c_cof=2.0535`, `m_strict_hi=+0.2724`b = `+20.8%`, `h2_allow=1.08`) | <= 110.100147 |
| replay `a`-tables build (1-owner 4-bank flat-array `<= 2^63.924 B < 2^64 B`) | <= 105.95909 |
| linear algebra (1-owner `< 2^64 B` 4-vec Arch B `Z[zeta_12]` BW, `n=112, m=8192, te=0`, `nl_cap=0.003, sw_cap=0.00605, k5=11.8 (hex6), tab_r=1.00`) | <= 108.393691 |
| self-contained classic `PM-Basis` generator (`m=8192, n=112`) | <= 93.141095 |
| both individual logs (`lx=71.0`, **pure `k=1` Boole + Wald i.i.d. boot**, Cantelli `1/(1 + 1*mu)`, `N_1log=2,348`, `Phi=0.4680`, `Wald=1.880x`; `32/185`: `110.296322`, `384/370`: `110.202740`) | <= 110.297322 |
| band split (`<=100.8898`), ECM (`<=96.2`), can `a` (`<=86`), 4-col SM (`<=84`), filt (`<=80`), `P_g/P_f/P_71` (`<=75.5`), FB (`<=75`), cert (`<=50`), ver (`<=11`) | <= 100.97 |
| wrapper (§1; `<= 2^80` Ir = `2^65.71` u) | <= 66 |

Core sum: `2^111.42858309 <= 2^111.4288` u. Meter `t_run = 2^123.0` (`p_meter <= 0.000330284183`) behind outer gate `K = 7201431620367092223` (`p_floor = 999/1000`): **`E[W] <= 1 + (K/2^64)(2^111.4288 + 2^66) < 2^110.071789445963`** dense-`F_(p^12)` u (`< 2^113.393717540850` G1 u).

## 8. Declared Exploratory Heuristics (Supported, Not Proved)

- **`H1/H2` (Smoothness, Yield, & Sieve Cost):** Two-sided `rho_dB` on exact-`O_H`-coprime pure-product pairs (`P4_prod_2ad_prim_v17`, `+0.272424`b strict-hi margin = `48.18 * SE` = `+20.78%` over `1.09*U'`), `h2_allow = 1.08`, exact local-root Murphy-`alpha` (`p <= 50,000`, `alpha_f = 2.727 <= 2.74`, `alpha_g_ge3 = 1.270 <= 1.28`, tail yield RMS `0.00689b`), live exact-degree `(72,24)` `(9,3)` real-FB joint smoothness calibration (`1,200` non-unit pairs, `1.112x`), and boot semismoothness (`pr_b = 2^-44.8138`), wrapped in `t_run = 2^123.0`.
- **`H3, H-SM, H-R, H-X` (Fixed-Instance Filtering on `F_2core` & Rank):** `r` unramified (proved), `P_chi Lambda` injective (`<4/r`), deterministic 2-core zero-column-loss filtering into generous fixed-instance caps (`nl_cap = 0.003` = `1.69x`, `sw_cap = 0.00605` = `1.68x`, strictly above every `2^18..2^24` mean, `k5_charged = 11.8`, `tab_ratio_charged = 1.00`, `cap_abort = 0`, deterministic `F_known = F_2core` with missing density `1 - peel_h` explicitly charged in `miss_side`), fixed-instance full rank `H-R` and nonzero BW minors `H-X` (supported by live **100% unscreened, zero-oracle, two-sided** `Q(zeta_15)` filtering + UNPLANTED `H-X` + post-hoc BSGS recovery + live full-incidence `t_dense=0.142` bridge + `unplanted_two_sided_bw`).
- **`H4, H4u`:** Raw-opcode-verified & bit-exact-simulated baseline scalar ARMv8.0-A integer kernels (`KERNELS.md`, `simd = 0`, ZERO FP/LSE); `Z[zeta_12]` unit share `93.41%` at `2^24`.
- **`H5` (Selected-Slot Conditional Pure-`k = 1` Emptiness Premise):** Explicitly declared conditional premise `P(M(I) = 0 | reached first-choice slot (b, s)) <= p_0(b, s)` using one-sided Cantelli `p_0(b, s) = max(0.0001, (1 + 1 * mu)^(-1))` (`DESCENT.md`, `ADAPTIVE.md`), NOT deduced unconditionally from generic marginals; combined across deterministic slots by **Boole's union bound** (`Phi_tree <= 0.60`) and wrapped in unbounded fresh-`e` i.i.d. boot retries (`Wald` identity, `p_descent_fail = 0.0`).

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0 B` (`2^119.121 B` replay `a`-tables; 1-owner build `<= 2^63.924 B < 2^64 B`; 1-owner 4-vec Arch B LA `2^63.949 B < 2^64 B`; `P_71` across 128 owners `< 2^62 B`/owner; `STORAGE.md`). **Provenance:** Builds on Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, Boudot (2017/758). **No prior verdict or score is a premise.**
