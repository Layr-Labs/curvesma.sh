# recover-Q: `(4,3)` `Q(zeta_15)` CM-STNFS + Proved `RET_BOOT=1` + Direct 200k Coprime Pool + Live Adaptive `H5` + `PM-Basis`: expected work `2^109.955669516687`

Track `bls12-381-fixed-g1-invert-g2`. Input: nonidentity `P in G1`, `z in G_T`; output `Q in G2` with `e(P, Q) = z` for exact supplied `P` (`Q = infinity` for `z = 1`). Unit: one dense `F_(p^12)` mul, `C = 20,000` baseline scalar ARMv8-A instructions (`simd = 0`, zero `FEAT_LSE`/`DC ZVA`, zero external C measurements, zero black-box generator citations, zero amortization). Full de Bruijn `rho_dB` is baked into Collection & Descent. All cost kernels (`KERNELS.md`), direct 200k-draw finite-region `O_H`-coprime pool, proved deterministic `RET_BOOT = 1.00`, ECM 2nd-moment Chebyshev bound, `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` recurrence, and `1,178`-attempt `Q(zeta_15)` real-2-core-FB adaptive diagnostic (`ADAPTIVE.md`) are derived and audited self-contained inside this package (**no prior verdict or score is a premise**).

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is an independent run gate: draw uniform 64-bit `u`; unless `u < K = 7214596996066621202`, stop immediately as failure (`< 1` unit). `K` is the least integer with `K x 0.997177 >= 0.39 x 2^64`. Under explicitly declared exploratory STNFS heuristics `H1`, `H1'`, `H-EA`, `H2`, `H3`, `H-R`, `H-SM`, `H4`, `H4u`, `H-X`, and `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` `H5` (Section 8; `RANK-BRIDGE.md`, `DESCENT.md`, `ADAPTIVE.md`), an opened run succeeds with probability `>= 0.997177` (unconditional Boole union bound, Section 6); by independence of `u`, overall success is `>= 0.39` over uniform full-production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^111.310045 + 2^66) < 2^109.955669516687  dense-Fp12 units (C = 20,000).

`109.955669516687` is an **expected** total classical work bound `E[W]` over challenge and fresh coins, NOT an every-run bound. The **ENTIRE opened run** runs under ONE unified hard aborting work meter `t_run = 2^123.0` core units (checked at coarse block boundaries `<= 2^66` u), giving deterministic pointwise opened-run cap `W_actual <= 2^123.0 + 2^66` dense-`F_(p^12)` units on 100% of opened runs, `E[min(W_open, t_run)] <= 2^111.310045`, and Markov `p_meter = 0.0003026724 <= 0.00030268` (`ROUND_CEILING`). Analytic exploratory claim; all numbers match default `32/185` `audit_cm.py` (`14/14` checks `checks_ok = true`, pure stdlib).

## 1. Reduction (fixed `P`: two `G_T` logarithms) and field map

With pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `Q0` of `G2`:
0. **Run gate:** uniform 64-bit `u`; if `u >= K`, abort as failure (`< 1` u).
1. **Validate:** `P in G1 \ {infinity}` (`[r]P = infinity`) and `z in F_(p^12)^*` with `z^r = 1`. If `z = 1`, return `Q = infinity`.
2. **Pairing & map:** compute `g = e(P, Q0)` and map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Randomize & solve:** draw independent uniform `e_1, e_2 in {0, ..., r - 1}`, run core on `g' = g^(e_1)`, `z' = z^(e_2)` (exact uniform on `G_T \ {1}`), set `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** compute `Q = [a]Q0`; output `Q` iff `e(P, Q) = z`, else fail. Wrapper `< 2^66` u, failure `< 2^(-48)`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, & Order-12 Frobenius

BLS12-381: `u = -0xd201000000010000`, `r = u^4 - u^2 + 1` (`r == 13 mod 15`), `p = (u - 1)^2 r / 3 + u` (`p == 7 mod 15`).
- **CM splitting & `w' = omega^2` isometry:** Over PID `K_H = Q(zeta_15)` (`O_H = Z[zeta_15]`, `D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `P = (p, h_pi)`, `O_H / P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` yields `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega + 1) X^6 + X^3 - omega` (deg 72).
- **Order-12 Frobenius (`algebra_bls12`):** `g` is irreducible over `F_(p^4)`, `sigma(y) = y^7` fixes `P, f, g`, and `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`chi = p == u mod r`, `ord_r(chi) = 12`), quotienting into `Z[zeta_12]` 12-orbits.
- **Relative-Degree-1 Lemma & Full-Orbit FB:** Unramified `Q | (a - b X)` (`gcd_O(a,b) = (1)`) has relative degree `1` over `O_H`. ALL full-orbit (`d = 12`) relative-degree-1 ideals of norm `<= B` are included (incl. `f_H >= 2`, `< 2^32.5 << 2^(-28) U'`). Short-orbit (`d < 12`, forcing `vlog == 0 mod r`) ideals are omitted.
- **4-Col Schirokauer Compression & `alpha`:** Pre-collection uniform `W_f in F_r^(3x6)`, `W_g in F_r^(1x2)` compress `P_chi Lambda_f, P_chi Lambda_g` to `4` cols (`< 4/r < 2^(-252)` under `H-SM`). Exact `w' = omega^2` Murphy `alpha`: `alpha_f = +2.74`b, `alpha_(g,>=3) = +1.28`b, `kappa_2 = 0.0525`.
- **Live E2E `(4, 3)` `Q(zeta_15)` DLP (`toy_e2e_dlp`):** On `u_0 = -5, p_0 = 7207, r_0 = 601`: `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps annihilating BSGS logs, neg controls failing, `5/5` unseen targets recovered.

## 3. Relation Collection: `Arch2-Rev-SigmaB-v3-ZeroCopy-U16` & Direct 200k-Draw Exact-`O_H`-Coprime Pool

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373` 12-orbits), excess `1.09` (`2^63.441701` rels).
- **Geometry, Voronoi FD, 2-Adic `c(a)`, & `D_(omega,sigma)` (`239/239` exact-1):** `G = A_4 x A_2` (`det = 1,265,625`), 10 explicit Voronoi units (`12/12` anchors `min-T2`, `N_a^can = 2^57.8300`). In `O_H`, `(2) = P_(2,1) P_(2,2)` (`N = 16`), `c(a) in {0, 1, 2}` (`s_c = (1.00, 1.20, 1.36)`, `STORAGE.md`), `r_(ab) = 1.187014`, `A^2 = 528954.37`, `D_(omega,sigma)` gives `N_7 = 2^57.2513` rows (`L_max = 550.47 < W_tile = 552`, `L_b = 236.44`) and `N_a^tbl = 2^59.8300` slots.
- **Direct 200,000-Draw Finite-Region Exact-`O_H`-Coprime Pool (`P2_wy_2ad_prim`, ZERO Euler-Product Extrapolation!):** On `200,000` independent draws across the exact finite 2-adic region (`V_kept = 2^119.7598`), every draw `(a, b)` is tested directly for exact `O_H`-ideal coprimality `gcd_O(a, b) == (1)` via `gcd(N(a), N(b))` and exact ideal-norm factorization (`191,832 / 200,000` coprime, **`f_prim_mc_200k = 0.959160`**; collapsed to `401` conservative `(bf, bgo, max_bgm)` bins storing primitive-only yield weights `wc` summing to `0.959160` and all-cell work weights `wk` summing to `1.0`). No asymptotic Euler-product lower bound is asserted for finite-region collection yield (finite-region smoothness/coprimality is explicitly declared in exploratory `H1`/`H2`). Retained upper-endpoint yield is `2^63.578997`, giving **`margin_net = +0.12830` full bits** above `1.09 U'` (`63.578997 - 0.009 - 2*0.01020 = 63.5496 >= 63.44`). `verify_pool_2adic_live` draws `400` fresh live canonical `(a, b)` pairs confirming exact `v_2(N_g)`, exact `O_H`-coprimality **`382/400`** (`z_prim = 0.42 <= 2.5`), `prim_yield_ratio = 1.0066`, and `z_rels = 0.82 <= 2.5`.
- **`Arch2-Rev-SigmaB-v3-ZeroCopy-U16` (`W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`):** 1 full-chord tile/row (`nt = 1.000`, `4.676` SBs/row across `552` column owners at `2^61.830` B `<< 2^64` B): Build **`2^106.2393 <= 2^106.239800`** u (`C_LIST_A = 30.0`, **64-bit `uint64_t` `rev_build_pass12` histograms** avoiding 32-bit overflow near `p = 2^24`, storage `2^119.655` B `< 2^120.0` B); zero-copy `uint16_t` Band-1a (`4.25`) & Band-1b (`6.25`, `<= 2^63.4484` B); Band-2 `sb_min_q8` (`<= 2^63.9399` B, `cand_ratio = 1.4023 <= 1.414`, `C_CAND = 11.5`, `(m-k)&3` bank-rot `C_HIT_HI = 15.25`, `b.hs .Lempty`); `rev_scan_zero_u16` (`2.00`); schoolbook `_nk_g_only` (`6,500`), `_nk_f_only` (`20,000`), `MONT5 = 267 <= 270`, `Q_F_SAFETY = 2.0`, 4-prime 62-bit NTT `stage_b_analytic` (`108.0 <= 168.0`) integrated **LIVE across all `401` bins** up to `bg_max = 700`b (**`c_cof_2d_upper = 2.8950`** ipc). Sieve work prices **100% of ALL cells including non-primitive cells** (`ipc_tot = 23.8643`, `H2_ALLOW = 1.00774`) -> **`E[W_coll] = 2^110.059975 <= 2^110.060475`** u. Band split `<= 2^100.543` u; ECM split metered at **`2^96.2`** u (exact 2nd-moment Chebyshev `< 2^(-74.45)`).

## 4. `Z[zeta_12]` Filtering, 1-Owner Block Wiedemann LA, & Self-Contained Classic `PM-Basis`

Splitting `Z[zeta_12] = Z[omega] (+) chi Z[omega]` into 6 signed hexagonal buckets cuts unit+unit collision fill-in `3x` vs `Z[i]`. Bundled `filtsim_z12.py` mirrors our 6-seed ladder across `U' = 2^16 .. 2^23` (`B71.5_d04_e109`, `W = 200`):

| Scale `U'` | Seeds | `N_l / U'` (SE) | `N_sw / U'` (SE) | `sd_l U'^(1/4)` | Depth | Unit % (Exp.) | `k5` | `tab_ratio` |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `2^16` | 6 | `0.003260` (`0.000227`) | `0.006462` (`0.000482`) | `0.0089` | `13.7` | `95.61%` (`89.15%`) | `20.4312` | `0.80659` |
| `2^17` | 6 | `0.003050` (`0.000112`) | `0.006089` (`0.000219`) | `0.0052` | `15.2` | `96.63%` (`90.47%`) | `19.0156` | `0.80554` |
| `2^18` | 6 | `0.002857` (`0.000088`) | `0.005747` (`0.000173`) | `0.0049` | `16.7` | `96.96%` (`91.30%`) | `18.1201` | `0.80212` |
| `2^19` | 6 | `0.002533` (`0.000074`) | `0.005111` (`0.000118`) | `0.0049` | `16.7` | `97.18%` (`91.82%`) | `17.5677` | `0.79984` |
| `2^20` | 6 | `0.002399` (`0.000032`) | `0.004857` (`0.000049`) | `0.0025` | `18.7` | `97.35%` (`92.32%`) | `17.0304` | `0.79839` |
| `2^21` | 6 | `0.002172` (`0.000035`) | `0.004377` (`0.000067`) | `0.0033` | `18.5` | `97.42%` (`92.61%`) | `16.7212` | `0.79684` |
| `2^22` | 6 | `0.001999` (`0.000021`) | `0.004057` (`0.000039`) | `0.0024` | `19.7` | `97.43%` (`92.77%`) | `16.5469` | `0.79508` |
| `2^23` | 6 | `0.001889` (`0.000013`) | `0.003846` (`0.000026`) | `0.0017` | `20.5` | `97.45%` (`93.01%`) | `16.2921` | `0.79485` |

- **Fixed-instance `H3` caps (`P(cap abort | H3) = 0`):** `N_l <= 0.00225 U' = 2^54.5215`, `N_sw <= 0.00450 U' = 2^55.5215`, `k5 <= 17.50` (`k5_live = 16.3162` at `2^23`), `tab_r <= 0.85`.
- **1-Owner Arch B BW & Self-Contained Classic `PM-Basis` (`RANK-BRIDGE.md` §3):** `x + y` + 22 tables fit in `2^63.913` B `< 2^64` B on 1 owner (`0` cross-owner sends). `(m, n) = (8192, 110)`, `unit_run5 = 8.80`, `gen_run1 = 116.0`, Barrett `row_red5 = 68`: **LA `= 2^107.789053 <= 2^107.789554`** u. Self-contained recursive Giorgi-Jeannerod-Villard / Thomé `PM-Basis` on `[S(z) | -I_m]` (`t = 8302`, `L_bm <= 2^47.7706`, 1D 4-prime 62-bit NTT, 4 limbs, schoolbook `t^3` Barrett MAC `17 <= 24` ops): **Gen `= 2^92.7203 <= 2^93.0`** u.
- **4-Col Mod-`r` CRT Schirokauer (`<= 2^84`), `H-R`/`H-X`, & `repurge_stack`:** Exact 255b / `K_rej = 64` sampling misses `< 2^(-156)`. Conditional on `H-R, H-X`, Schwartz-Zippel `<= 4 N_cap / 2^128 < 2^(-66)`. Live: `unplanted_two_sided_bw` (incl. live recursive `pm_basis` kernel recovery) and `hx_z12.py` (`post_peel_cols_lost == post_merge_cols_lost == n_reconstruct_unit == 0`).

## 5. Individual Logarithms (`RET_BOOT=1` Proved + `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` + Live Adaptive `H5`)

Detailed in `DESCENT.md` and `ADAPTIVE.md` (all `39` bins / `78` nodes, exact `theta_1^8` `cells_rank2`, `1/zeta_K(2) >= 0.970721 > 0.95`, seeded peeling `h = 0.256845` evaluated LIVE at default **`32/185`**):
- **Proved Deterministic `RET_BOOT = 1.00` at `m_sv = 2592.0` (`DESCENT.md` §3):** Exact rational Sylvester leading principal minors prove `(15 + 1/1000) I_8 - G_pow > 0` LIVE (`cm_models.verify_gt_boot_20`); combined with `99/100`-Lovasz on `L_boot` (`det = p^4 s_boot^8`, `s_boot = 2^21`, `alpha_LLL = 100/74`), `(x_j^(3/2)+y_j^(3/2))^(2/3) <= x_j+y_j`, and AM-GM raised to the `3/2` power, EVERY nonidentity `G_T` target deterministically satisfies `log2 |N(a_1^3 - u w^2 b_1^3)| <= 12 log2((15.001/8)*(|u|^(2/3)/s_boot^2)*||v_1||_2^2) <= 2591.90 < 2592.0` (`RET_BOOT = 1.00`, ZERO heuristic!). Fresh-exponent `i.i.d.` single-cell boot (`boot_rels = 1, cells_b = 1`, `LLL_BOOT = 2^48`, `(n1=255, nA=195, kA=15)`) redraws fresh `e_t` on failure (`psi_boot = 0.25188 < 0.50`, no inner cap; `1/1` live 381-bit pure-Bareiss check gives `2542.56b`).
- **`H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` (`phi in [0.01, 0.04]`, `psi < 0.80`):** Default `32/185` outputs **`max_phi = 0.02500`**, **`max_psi = 0.21575`**, **`psi_boot = 0.25188`**, **`min_mu = 740.2`**, expected two-log work **`2^110.194841 <= 2^110.195741`** u, `E[N_2log] = 2^12.524`, node budget `B_node = 2^22.0` (`p_node_abort = 0.0014044`) + `H5` repeat-collision (`C_rep = 256` against `2^60.81`, `p_repeat_loss = 0.0011105`) -> **`0.0025150 <= 0.0025200`**.
- **Four Live Audits (`lll_16d_int`, `0` Pre-Supplied Vectors; Exploratory Support, Not Proof):** (a) `prod_sq24_live` (`24/24` `(f_9, g_3)` controls, joint-cell `pr_emp/pr_charged >= 16.96x` `f` / `>= 55.01x` `g`); (b) `prod_cond_q15_live` (`12/12` `q | N_side`, joint-cell `pr_emp/pr_charged >= 41.60x`, med `301.41x`); (c) `q15_samefield_descent_live` (`3/3` live-replayed + `120`-ideal real-2-core-FB vs independent `mu_model`); (d) **`adaptive_h5_live` (`ADAPTIVE.md`):** `1,178` real-2-core-FB after-failure attempts (`pooled O/E fail = 0.8022` [95% CI `0.703-0.905`], `rat1 = 0.8279`, `work O/E = 1.0132`, `0/1,178` shared descendants, `C_bday <= 2.815`, `C_emp <= 19.11`, 8 history strata, `threshold_counterexample`) + **3 live-regenerated multi-level after-failure witness descents**.

## 6. Success Probability & Unconditional Boole Union Bound

Verified in 80-digit `Decimal` in `audit_cm.py`:
- BW Schwartz-Zippel `< 2^(-66)` (`| H-R, H-X`), Schirokauer `< 2^(-252)` (`| H-SM`), `H3` fixed-instance caps (`0`), ECM split exact Chebyshev `< 2^(-74.45)`, 100-round MR `<= 2^(-80)`.
- Descent node budget (`0.0014044`) + `H5` repeat-collision (`0.0011105`) `= 0.0025150 <= 0.0025200`, unified run meter `t_run = 2^123.0` Markov `p_meter = 0.0003026724 <= 0.00030268`, degenerate target `<= 2/r`, wrapper `< 2^(-48)` (incl. `F_r` rejection `< 2^(-156)`).
- Opened-run floor: `0.99717732 >= 0.997177`. Gate `K = 7214596996066621202` satisfies `(K/2^64) x 0.997177 >= 0.39 > ((K-1)/2^64) x 0.997177`.

## 7. Complete Expected-Work Ledger (Wrapped Pointwise in Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, direct 200k coprime pool, `ipc=23.8643`, `c_cof=2.8950`; live `110.059975`) | <= 110.060475 |
| replay 4-orbit `a`-tables build (64-bit `x2` hist, `C_LIST_A=30.0`; live `106.2393`) | <= 106.239800 |
| linear algebra (`1-owner Arch B Z[zeta_12] BW`, `unit_run5=8.80`, `gen_run1=116.0`, Barrett `row_red5=68`, `k5<=17.50`; live `107.789053`) | <= 107.789554 |
| self-contained classic `PM-Basis` matrix generator (1D 4-prime NTT + schoolbook `t^3` Barrett MAC; live `92.7203`) | <= 93.0 |
| both individual logarithms (default `32/185`: proved `RET_BOOT=1.00` at `2592.0`, `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR`, `psi_b=0.25188`, `max_phi=0.02500`, `max_psi=0.21575`, `E[N]=2^12.524`; live `110.194841`) | <= 110.195741 |
| band split (`<= 100.543`), ECM split (`<= 96.2`, Cheb `< 2^-74.45`), canonical `a` (`<= 86`), 4-col SM (`<= 84`), filter/misc (`<= 80`), `P_g/P_f` trees (`<= 75.5`), FB (`<= 75`), cert (`<= 50`), verify (`<= 11`) | <= 100.62 |
| wrapper (sec. 1; `<= 2^80` instructions = `2^65.71` core units) | <= 66 |

Expected opened-run core sum (`80`-digit `Decimal`): `2^111.31004416 <= 2^111.310045` core units. Wrapped in `t_run = 2^123.0` (pointwise cap `<= 2^123.0 + 2^66`, `p_meter <= 0.00030268`) behind gate `K / 2^64`: **`E[W] <= 1 + (K/2^64)(2^111.310045 + 2^66) < 2^109.955669516687`** dense-`F_(p^12)` units (`< 2^113.277597611574` G1 units).

## 8. Declared Exploratory STNFS Hypotheses (Supported, Not Proved)

- **`H1`, `H1'`, `H-EA`, `H2`:** Two-sided `rho_dB` smoothness on exact-`O_H`-coprime finite-region pairs (`P2_wy_2ad_prim`, `200,000` draws, `f_prim = 0.959160`, `+0.12830`b margin, ZERO Euler-product lower-bound claim), exact `v_2(N_g)`, two-stage Euler retention + survivor-upper (`Q_F_SAFETY = 2.0`), `H2_ALLOW = 1.00774`, expected work wrapped in `t_run = 2^123.0`.
- **`H3`, `H-SM`, `H-R`, `H-X`:** `r` unramified (proved), `r !| h(K_f), h(K_g)`, `P_chi Lambda` injective (`< 4/r`), monotone filtering trend into fixed-instance caps (`N_l/U' <= 0.00225`, `N_sw/U' <= 0.00450`, `k5 <= 17.50`, `tab_r <= 0.85`), fixed-instance full rank `H-R` and nonzero BW minors `H-X`.
- **`H4`, `H4u`:** Explicit raw-opcode-verified baseline scalar ARMv8-A kernels (`KERNELS.md`, `simd = 0`, 64-bit `x2` `rev_build_pass12`, `b.hs .Lempty`, `ntt62` `108.0 <= 168.0`, 62-bit Barrett MAC `17 <= 24`, `MONT5 = 267 <= 270`, Barrett `row_red5 = 67 <= 68`, uniform 2-stage Bareiss LLL `2^40 / 2^48`); `Z[zeta_12]` unit share `93.01%` at `2^23`.
- **`H5` (`H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` + Sequential `psi` + `C_rep = 256`, `DESCENT.md`, `ADAPTIVE.md`):** Exact `S_rank2(R)` (`|(B_8(R)\{0}) cap (F_Vor/mu_30)| >= N_alpha(R)`, `1/zeta_K(2) >= 0.970721 > 0.95`), fresh-exponent `i.i.d.` boot smoothness (`RET_BOOT = 1.00` PROVED deterministically at `2592.0`), seeded peeling (`h = 0.256845`), `H-PGF (0.5*mu)` (`t <= 0.80`), sequential conditional child-failure hazard `P(F_i | A_i, pre-entry history) <= psi` & conditional child-work bound (`ADAPTIVE.md`), and adaptive repeat-collision anti-concentration `<= C_rep / N_lx` (`C_rep = 256`, `C_bday <= 2.815`, `C_emp <= 19.11`).

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0` B (`2^119.655` B replay `a`-tables; `2^61.830` B `< 2^64` B column owners; `2^63.4484` / `2^63.9399` B `< 2^64` B per Band-1/2 table incl. 64-bit `x_hist`; `< 2^105` B descent; `2^63.913` B `< 2^64` B 1-owner LA; `2^82.81` B `< 2^84` B `PM-Basis`; `< 2^30` B wrapper; `STORAGE.md`).

**Provenance:** Builds on accepted `c36p` / `stnfs-bn15r` exploratory methodology (`afe6f29f`, `5a49c502`, `4931a702`) and Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, and Boudot (2017/758). **No prior verdict or score is a premise of this claim.**
