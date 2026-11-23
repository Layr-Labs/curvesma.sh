# `recover-P` (`bls12-381-fixed-g2-invert-g1`) through two `G_T` logarithms with `(4, 3)` `Q(zeta_15)` CM-STNFS behind a run gate: expected work `2^109.733089445958` units

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity `Q in G_2` and `T in G_T`; output `P in G_1` with `e(P, Q) = T` for the exact supplied `Q` (`P = O` for `T = 1`). The score unit is one dense `F_{p^12}` multiplication, `C = 20,000` pure baseline scalar ARMv8.0-A integer instruction equivalents (`simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64/x18; AAPCS64 `sub/add sp` `MONT5 = 279 <= 280`; `mont2_cios = 47 <= 48 < 58.84` Ir for `L = 2`; zero undefined ASM labels and per-kernel disjoint GPR maps across all 30 hot-loop blocks). No advice, smaller targets, or amortization are used. **No prior verdict or score is a premise.**

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is a run gate: draw a fresh 64-bit word `u` from an independent entropy stream **drawn before and strictly isolated from the challenge `(Q, T)` and all internal solver PRNG seeds** (blinding exponents `e`, Schirokauer/BW matrices, boot cell choices); unless `u < K = 7201431620367092223`, stop as a failure (`< 1` u). `K` is the least integer with `K * 0.9990 >= 0.39 * 2^64`. Under exploratory premises `H1`–`H5`, `H-R`, `H-SM`, `H-X`, and `H-FFT` (section 9), the algorithm returns the exact `P` with probability `>= 0.9990` once past the gate (section 7, retaining `+7.14e-4` headroom under `+10%` meter sensitivity + `2^-48` wrapper, and `+6.30e-4` modeled-success headroom under simultaneous `1.5x` `H5` variance-floor stress); `u` is independent of the challenge and all solver coins, so cold-run success is `>= 0.39` and

    E[W] <= 1 + (K / 2^64) * (2^111.0901 + 2^66) < 2^109.733089445958  dense-Fp12 units (C = 20,000).

`2^111.0901` is the sum of section 8's rows (every-run bounds for collection and linear algebra plus expected work for the two Pure `k = 1` Cantelli `K1_CANTELLI_P0005_MU50` individual logarithms `2^109.9092` and ECM relation splitting — every major row carries positive slack above its live recomputed bound, `>= +0.0094`b on collection and `>= +0.0160`b on all other major rows). The score `109.733089445958` is an **expected** work bound, not a worst-case bound on every run: every run is bounded pointwise behind the hard aborting meter `T_run = 2^123` u plus wrapper:

    W <= 2^123 + 2^66 < 2^(123 + 2^-56)  units.

Analytic, exploratory-v1 claim. Under native `python3 -I -B` (pure Python 3.11+ stdlib), `audit_cm.py` recomputes live from scratch all `16/16` checks: algebraic/CRT/NTT identities, Murphy-alpha, real-FB joint collisions + non-vacuous `(f, g)` survival coupling, toy E2E DLP, `sigma/B` full-slab Band-2 exactness, 400-draw live `R+7`-envelope 2-adic cross-check of the 1,000,000-draw `A^2 = 321,500` histogram, Voronoi + 16D LLL + `4x` 381-bit `G_T` `F_{p^6}`-subfield elimination, `srt_toy`, `0.25`-bit 2D `Retention`/`c_cof_2d_upper` integration with static 3-class `T_g[c] = (299, 287, 275)` + constant `T_f_const = 118`b pre-gates + full `0.02`-bit-grid `pool_hist_2d` collection-cut `ret_mult` conditioning (live-asserted `>= 0.90 / 0.82`) + deterministic relation-split cap `n_rels_upper = 2^63.530` & Cauchy-Schwarz correlated random-sum ECM Chebyshev `< 2^-70`, `2^12` filtering (Stage-B row-rejection + expanded `W_max = 4,096` guard) + H-X, `2^14` `[6b]` Stage-B-capped collection-cut-conditioned `missing_recon`, unscreened two-sided `Q(zeta_15)` H-X (`R0 = 181` BSGS + algebraic `gamma_0`-blinding ratio + full-matrix `exp_unit_hex6`), production-aligned structured `PM-Basis` (`prod_eq_schur = True`, leading `[z^{s_out[j]}] q_j(z)[:s_sm] != 0` 1-pass selector + zero-SM-lead negative control), raw-opcode ARMv8.0-A parsing (`30` hot-loop blocks incl. `rev_b2_row_k0_step4` and `47`-op `mont2_cios`, `0` undefined labels, disjoint physical register maps) + 64-bit simulations (incl. `mont2_cios` exact 128-bit odd-modulus verification, `rev_b2_row_k0_step4` at `p = 2^64 - 59` & `Q(zeta_15)` orbits, `64 KiB`-buffered Band-2 full & row-end drains, `47`-op scan, AAPCS64 `MONT5 = 279`, 1-vec Horner `lsr #2`, 4D Fincke-Pohst DFS + ambient sign canonicalization, 55b `SQ_LUT`, half-ball rank `0..17,071`, `Decimal` `+6` `bmul_6` witness), exact 384-mesh 76-node `K1_CANTELLI_P0005_MU50` raw-`L_q` descent evaluation + live `1.5x` variance-floor stress + live regeneration of `4` production-size `71..240`-bit raw-`L_q` `deg-(72, 24)` `Q_CM`-LLL reductions with `32` exact `B_16(16)` norm & `>= 2x` yield-headroom checks + live regeneration of `4` raw-`L_q` strict-descent `deg-(32, 24)` parents and `4` reached children + `+-16`-ULP Dickman enclosure (15 points) + live integer `BOUND_CELL (+6)` + live rejection of the `(135, g)` counterexample + live `L_min = 2` branch-tree DP / 12-bit ladder / Suyama-Newton guard + live `>= 4` distinct encountered children + live 384-entry Part C recomputation + live multi-bin exact-DP `UDUP = 0.40` recomputation + 31k strict-descent raw-`L_q` O/E & Fano + 15.2k production raw-`L_q` archive assertions + `q15_retain16`, and the 80-digit `Decimal` ledger/gate + manifest verification.

## 1. Reduction (Fixed `Q`: Calibrated `gamma_0`-Blinded Two-Logarithm Inversion) and Field Map

Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md` §5) recovers a nonzero homogeneous virtual-logarithm functional `V_g(alpha) = lambda * log_g0(phi_iso^(-1)(alpha)^{e_cof}) mod r` (`e_cof := (p^12 - 1)/r`, `gcd(r, e_cof) == 1` asserted live). Pick ONE Stage-1 relation `(a_0, b_0)` with `V_g(a_0 - b_0 X) != 0 mod r`, and define `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^{e_cof} in G_T \ {1}` (`<= 2^11` u, ZERO third descent) with known `V_0 := e_cof * V_g(a_0 - b_0 X) mod r != 0`. Given `(Q, T)` with `Q = [a] Q_0 in G_2 \ {O}` and `T = e(P, Q) = g_0^{a b} in G_T`:
  0. Run gate (`Claim`): draw fresh independent 64-bit `u`; if `u >= K_gate`, stop as failure (`< 1` u).
  1. Check `Q` (on twist, `[r]Q = O`, `Q != O`) and `T` (`T^r = 1`). If `T == 1`, return `P = O`.
  2. Compute `h := e(P_0, Q) = g_0^a in G_T` (`<= 2^11` u) and map `h, T` into the `Q(zeta_15)` tower (`field-map.md`, wrapper `<= 2^66` u, failure `< 2^-48`).
  3. For each target `z in {h, T}`: run unbounded `i.i.d.` fresh-`e` Wald boot renewals (`z' := z * gamma_0^e`, Guillevic `F_{p^6}`-subfield elimination with `W_1 != 0` unconditionally, `1` uniform half-cell of `L_boot(z')` in `B_16(2)/{+-1}`, Pure `k = 1` zero-retry raw-`L_q` tree descent) until a complete tree yields `V(z') = (lambda * e_cof) * log_g0(z) + e * V_0 mod r`, then unblind `V(z) := (V(z') - e * V_0) mod r`.
  4. Compute `b := V(T) * V(h)^(-1) mod r` (cancelling `lambda * e_cof`) and output `P := [b] P_0 in G_1` if `e(P, Q) == T`, otherwise fail.

## 2. Tower, Polynomials, and `Aut = C_12` Symmetry

Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, lcm(p^4-1, p^6-1)) == 1` and `gcd(r, (p^12-1)/r) == 1` asserted live). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`, monogenic maximal order `O_H = Z[zeta_15]`, `h_K = 1`), `g(x) = x^3 - u omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = <sigma, omega> = Z/12Z`. Smoothness bounds: **`B_f = 2^70.79`, `B_g = 2^71.19`, `log2(U') = 62.827711`**, target excess `exc = 1.09`.

## 3. Relation Collection (`P5_prod_2ad_prim_v5_1m_r83`, `Arch2-Prod-ZeroSlice-Buf64K-B2-U16-Streamed`)

- **Region Volume (`A^2 = 321,500, r_ab = 0.83`):** `V_a = 2^55.8744 < 0.92 * 2^56` (`N_a^tbl = 4 V_a = 2^57.8744`), `V_b = 2^64.5535` (`L_max = 440.08 < W_tile = 448`, `L_b = 189.02`, `N_7 = 2^54.9910`, `V_kept = N_a^tbl * N_7 * L_b = 2^120.4279`, `f_prim = 0.955571`).
- **1,000,000-Draw 2D Pool & Static Pre-Gates:** 100% deterministic `c_gbatch` (`P_{71.19} + P_{70.79}` + retained `g`-batch candidate product trees across `2,048` owners at `3.091 EiB << 16.00 EiB` each, ALL THREE `L-11` subtrees `<= 1.053 EiB < 2^62 B`, `live = 71.39 <= 75.5` u) -> `det_credited = 0.99`. Static 2-adic class thresholds `T_g[c] = (299, 287, 275)` and constant `T_f_const = 118`b give `margin_net_rethi = +0.10531`b, `margin_strict_hi = +0.03474`b (`15.1` SE), `log2_rels_2se_live = 62.982177 >= 62.980 = LOG2_RELS`, conditioned odd-norm means `(547.280, 671.016) >= (546.0, 670.0)`, `c_cof_2d_upper_live = 1.0557` Ir/cell.
- **Per-Band `ipc_tot` Decomposition (`P_min = 32, p_m4_s = 52.0, p_m4 = 57.5, p_lo = p_f = 58.0, p_g = 63.5`, `H2_ALLOW = 1.01`):**

| Component | `g`-side (Ir/cell) | `f`-side (Ir/cell) | Notes (`KERNELS.md`, `STORAGE.md`) |
|---|---:|---:|---|
| Band-1a-long `[2^32, 2^52.0)` hits | `2.3848` | `2.3848` | `rev_b1_stream_u16` `1.125` + `rev_b1a_apply_u16` `3.625` <= `4.78` |
| Band-1a-short `[2^52.0, 2^57.5)` hits | `1.0332` | `1.0332` | `rev_b1_stream_tail1` `4.0` + `rev_b1a_apply_1` `6.0` = `10.0` |
| Band-1b `[2^57.5, 2^58.0)` hits | `0.1068` | `0.1068` | `rev_b1_stream_tail1` `4.0` + `rev_b1b_apply_1` `8.0` = `12.0` |
| Band-1 4-bank lookups + guards + row init | `1.0720` | `1.0720` | `4 * 96 = 384` (1a-long), `4 * 10 = 40` (1a-short), `12` (1b) |
| Band-2 `[2^58.0, 2^63.5)` row `k_0` + probe/SB hdr | `0.9767` | `0.0000` | `8.25 <= 10.0` `rev_b2_row_k0_step4` + `19 <= 20` probe + `cbz`-guarded `10` SB hdr |
| Sequential `64 KiB`-buffered Band-2 `(x4, w_lp)` | `2.2215` | `0.0000` | `15` emit `[x7],#16` + `6` recv `[x_rx],#16` + `0.004` flush = `21.004 <= 22.0` + `28` row-end drains |
| 64-cell `ORR` scan + blunt `1,000`-Ir cold | `0.8400` | *(combined)* | Live `sga_mean`-tied `0.8231 <= 0.84` (`47`-op `cbz`/`.Lscan_done: ret` + `1,000` Ir cold) |
| 2D Cofactorization (`c_cof_2d_upper_live`) | `1.0557` | *(combined)* | Static 3-class `T_g=(299,287,275)` + `T_f=118`b, `MONT5=280`, `ntt3_F` |
| **Total `ipc_tot` (exact sum, before `1.01` `H2_ALLOW`)** | **`14.2874`** | *(total)* | `V_kept + log2(ipc_tot * 1.01 / 20000)` = **`109.991165`** <= **`110.0006`** (`+0.0094`b slack) |

  `replay_a_tables <= 2^105.3423` u (`live = 105.3263`, `+0.0160`b slack), 1-owner sieve RAM `2^63.087 B`, 1-prime-per-owner Band-2 build peak `<= 14.981 EiB` & resident `<= 7.491 EiB`, total tables `2^118.349 B < 2^120 B`.

## 4. Filtering and Linear Algebra (`ArchB-1Dense-16Proj-16Shard-4Vec-ZeroCopy`, Stage-B Row-Rejected Ladder, Conservative Caps `nl_cap = 0.00278`, `sw_cap = 0.00568`, `nonunit_ub = 0.25`)

- **Stage-B Row-Rejected Conditioned Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0`, conditioned norms `[547.0, 671.0]`, `sbf = 375`, `sbg = 315`, `W_max = 4,096`):** Raw smooth rows (`ret = 1.0`, zero prime relocation between rows) whose residual exceeds `row_bf > 375` or `row_bg > 315` are discarded outright (`stageb_row_keep_frac = 0.6483`, automatically inducing `ret_f_realized = 0.8246`, `ret_g_realized = 0.7770` at `U = 2^25`). At `U = 2^25`, `nl_25 = 0.0004961` vs conservative `nl_cap = 0.00278` (`5.60x` headroom), `sw_25 = 0.0008920` vs conservative `sw_cap = 0.00568` (`6.37x` headroom), `eu_25 = 0.918213`, exact expanded `avg_w_25 = 200.0182 < 201.0`, `max_row_w_25 = 1,596.05 <= 4,096`, `merge_stack_mean_w = 73.29 < 201.0`. Note on unit-fraction accounting: the earlier `0.7659` figure was the 4-tuple ring-entry fraction whereas the Krylov kernel executes expanded 6-term hexagonal entries (`exp_unit_hex6 = 0.93443` at `U = 1024`, climbing to `0.9618` at `U = 4096`), and charging **`nonunit_ub = 0.25`** (`k5_charged = 12.00` for 4-vec Krylov vs `k5_live = 8.8284`; `k5_horner_charged = 14.00` for 1-vec Horner vs `k5_horner_live = 10.7825`) conservatively covers even `1 - 0.7659 = 0.2341`.
- **Complete `< 16.00 EiB` Owner Map (`m = 32,768, n = 128 = 32 * 4, n_sm = 4`, ALL `< 2^64 B` asserted live):** Dense owner `2^63.368 B = 10.323 EiB` (`2.581 EiB` Horner reuse); `16` projection owners holding `Z + Y=[S,BZ] + X` at `2^63.359 B = 10.263 EiB`; `2,048` purge owners at `2^62.552 B = 5.866 EiB` (`1,040 B/rel`); `8,192` mutable filtering pools at `2^62.459 B = 5.498 EiB` (`nonunit_ub = 0.25`); `2,048` Krylov sequence owners at `2^63.354 B = 10.224 EiB`; `2^22` PM-Basis owners at `2^63.108 B = 8.623 EiB` (final `F(z)` `0.639 EiB`); `16` read-only matrix shard owners at `2^62.124 B = 4.358 EiB`; aggregate streaming `c_stream_mat = 73.0`, `c_stream_mat_1vec = 291.0`, `c_xproj = 70.0`.
- **Explicit LA Total Work (`audit_cm.py`):** With `N = 0.00278 * U' = 2^54.3370`, `N_sw = 0.00568 * U' = 2^55.3678`, `k5_charged = 12.00`, `k5_horner_charged = 14.00`, `1`-pass Horner justified by `PM-Basis` `[z^{s_out[j]}] q_j(z)[:s_sm] != 0`: **`linear_algebra = 107.232346 <= 107.2484` u (`+0.0160`b slack)**, **`matrix_generator = 98.2557 <= 98.2717` u (`+0.0160`b slack)**.

## 5. Main-Phase Success & Missing-Log Reconstruction (`[6b]`)

Collection ECM relation splitting deterministically truncates at `n_rels_upper = 2^63.530` (`Pr[N_rels > n_rels_upper] == 0`), using prime-exponent Newton perfect-power precheck and Cauchy-Schwarz correlated random-sum variance under post-`2^48` collection `p_r = 1e-3`, giving `ew_l2 = 94.5038 <= 95.3025` u and Chebyshev `< 2^-70.62 < 2^-70` (asserted live). Block Wiedemann fails with probability `<= 2^-66`, misc `< 2^-80`, Schirokauer injectivity `<= 2/r` (`+ 4/r` descent = `6/r`). `[6b]` `missing_recon_6b_live` (synthetic `2^14`-column Stage-B-capped peeling diagnostic) confirms `meas_miss <= mod_miss` across all occurrence bands.

## 6. Individual Logarithms (`K1_CANTELLI_P0005_MU50`, Pure `k = 1` Zero-Retry Boole on Raw `L_q` + Exact One-Sided Cantelli + Unbounded `i.i.d.` Fresh-`e` Wald Boot)

`DESCENT.md` (`verify_q15_descent.py` + `q15_ecm16.py` + `q15_bound_cell.py`): `76` 100%-eager uncapped `dsieve16` nodes (`BOUND_CELL = 5,000` Ir charged on both sieve and direct-ECM branches, `LB = LX = 68.5`, `bins [68.5, 75), [75, 80)..[255, 260)`), **single-scan (`K = 1`) raw special-`q` lattice `L_q`** (`bp = (0.0, 0.0)`), exact one-sided Cantelli `p_0(b) = max(0.0005, (1 + 0.5 * mu(b))^{-1})` under explicitly declared exploratory STNFS variance premise `H5` (`Var(M | slot) <= max(2 * mu, 0.0005 * mu^2)`, `min_mu = 85.64 >= 50.0`, `p0_max = 0.022820`), linear missing-leaf retention `max(0.0, 1.0 - miss)`, branch-tree-capped `L_min = 2` `47`-op `mont2_cios` + `12`-bit geometric ECM ladder (`SEPARATION = 0.5`, explicit failed-replay retry `p_r = 0.0043`, `2^26` Newton precheck) + unbounded `i.i.d.` fresh-`e` Wald boot (`BOOT16 = (230, 312, 15, 2)`, `BOOT_KEEP_FRAC = 0.90 < CP95 LB = 0.9395`, `BALL_PAIR = 0.0`, `N_slots = 682 <= 1,000`, `Phi_tree = 0.445595 <= 0.50`, `E[runs/log] = 1.8037 <= 2.0`, **`boot_phi = 0.0`**, **`p_descent_fail = 0.0`**).

Expected two-log work `E[D] <= 2^109.9092` u (`live = 2^109.893179` u, `+0.0160`b slack). The unified hard meter caps the entire opened run at `T_run = 2^123` u: `E[min(W_run, T_run)] <= 2^111.0901` u, and by Markov's inequality **`P(meter) <= E[W_run] / T_run <= 0.00025988 <= 0.00027`** (`<= 0.00028587` at `+10%`, `<= 0.00037023` under simultaneous `1.5x` variance-floor stress + `+10%`).

## 7. Success

Core opened-run success (80-digit `Decimal`, **`p_descent_fail = 0.0`** unconditionally under unbounded Wald renewal):

    P(success | open) >= (1 - 2^-66 - 2^-70 - 2^-80 - 2/r) * (1 - 0.00025988 - 4/r)
                      >= 0.9997401   (exact unrounded: 0.99974012)

with `+10%` meter sensitivity (`<= 0.00028587`) AND wrapper failure (`2^-48`), **`P(success | open) >= 0.99971413 >= 0.9990`** (`+7.14e-4` headroom). As a separate modeled opened-run success / meter headroom diagnostic (explicitly NOT a substitute for `H5` in the charged work ledger, since raising the variance ceiling increases `two_log_units` from `109.8932` to `110.6323`), under simultaneous `1.5x` variance-floor stress (`P0_FLOOR = 0.00075`, `FANO_CAP = 3.0`, `Phi_tree = 0.667853 < 0.70`, `E[runs/log] = 3.0107`) + `+10%` stressed meter + `2^-48` wrapper, **`succ_head_stress = 0.99962977 >= 0.9990 = p_floor`** (`+6.30e-4` headroom, asserted live in `audit_cm.py`). Past the gate `P(success | open) >= 0.9990`; cold-run success `>= (K / 2^64) * 0.9990 >= 0.39`, and `E[W] <= 1 + (K / 2^64) * (2^111.0901 + 2^66) < 2^109.733089445958`.

## 8. Ledger (`C = 20,000` Ir/u)

| Component | Live `log2(u)` | Charged `log2(u)` | Slack |
|---|---:|---:|---:|
| `relation_collection` | `109.991165` | `110.0006` | `+0.0094`b |
| `replay_a_tables` | `105.3263` | `105.3423` | `+0.0160`b |
| `linear_algebra` | `107.232346` | `107.2484` | `+0.0160`b |
| `matrix_generator` | `98.2557` | `98.2717` | `+0.0160`b |
| `both_individual_logs_exp` | `109.893179` | `109.9092` | `+0.0160`b |
| `band_split` | `100.0057` | `100.0217` | `+0.0160`b |
| `ecm_split_expected` | `94.5038` | `95.3025` | `+0.7987`b |
| `canonical_a_and_sort` | `<= 86.0` | `86.0` | — |
| `schirokauer_4col_crt` | `<= 84.0` | `84.0` | — |
| `filtering_and_misc` | `<= 73.5` | `80.0` | `+6.5`b |
| `pg_pf_product_trees` | `71.3876` | `75.5` | `+4.1`b |
| `factor_base_roots_ideals` | `<= 75.0` | `75.0` | — |
| `numerical_certificate` | `<= 50.0` | `50.0` | — |
| `verification` | `<= 11.0` | `11.0` | — |
| **Opened-Run Expected Sum (`E[W_run]`)** | **`111.090071`** | **`111.0901`** | — |
| **Hard Pointwise Meter (`T_run + wrapper`)** | — | **`< 123 + 2^-56`** | — |
| **Gated Cold-Run Expected Claim (`E[W]`, `p_cold >= 0.39`)** | **`109.733089445957`** | **`109.733089445958`** | — |

## 9. Heuristics (Explicit Exploratory-v1 Premises)

- **`H1` (Smoothness & Split Density):** Dickman-de Bruijn / Murphy-alpha smoothness (`alpha_f = 2.74, alpha_g_total = 0.77`, exact 2-adic `v_2(N_g)` class law, 1,000,000-draw `A^2 = 321,500, r_ab = 0.83` norm pool cross-checked online), within-bin conditional `(f, g)` survival independence (accompanied at toy scale by `real_fb_joint_coll_live`), exploratory `0.99 / 1.01` `li(x)` ideal-density brackets, and split-prime fraction `<= KH = 0.08`; 1-cell-per-fresh-blind boot has `BALL_PAIR = 0.0` exact.
- **`H2` (Collection & Sieve Cost):** Uniform 1,000,000-draw `r_ab = 0.83` pool with `H2_ALLOW = 1.01`, `det_credited = 0.99`, static 3-class `T_g[c] = (299, 287, 275)` + constant `T_f_const = 118`b pre-gates, explicit Band-1a/1b/2 streaming (`rev_b2_row_k0_step4` `8.25 <= 10.0`, `64 KiB`-buffered Band-2 `21.004 <= 22.0`), `47`-op scan + blunt `1,000`-Ir cold bound (`<= 0.84`), generous scalar build bounds (`32.0 / 10.0`), and `SIEVE_HIT = 2,600` Ir/hit.
- **`H3` (Matrix & Filtering):** `r` unramified and prime to class numbers; conservative zero-column-loss filtering caps `nl_cap = 0.00278`, `sw_cap = 0.00568` (`> 5.6x` above the Stage-B row-rejected `2^25` ladder; extrapolating from toy `Z[zeta_12]` scales `U = 2^16..2^25` to production `U' = 2^62.83`), deterministic `W_max = 4,096` guard, `nonunit_ub = 0.25`, `k5_charged = 12.00` (4-vec Krylov), `k5_horner_charged = 14.00` (1-vec Horner), `c_stream_mat = 73.0`, `c_stream_mat_1vec = 291.0`, `c_xproj = 70.0`.
- **`H-R` & `H-SM` (`RANK-BRIDGE.md` §1):** Filtered `A = [M | S]` has `rank(M) = N`, `rank(A) = D - 1`, and data-oblivious Schirokauer compression loses injectivity with probability `<= 4/r` (`6/r` total).
- **`H4` (Kernel Cost):** Pure baseline scalar ARMv8.0-A integer instruction counts (`KERNELS.md`, `STORAGE.md`, `30` hot-loop blocks incl. `rev_b2_row_k0_step4` and `47`-op `mont2_cios`; abstract 64-bit owner address model).
- **`H-X` (Block Wiedemann):** Full Krylov/Hankel ranks for block-diagonal `X` and `4` Schirokauer columns at charged `(d, e)` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx`, `unplanted_two_sided_bw`).
- **`H-FFT` (Local SSA):** Odd-`k`-exact `ntt3_F(N)` / `fft_F(N)` bounds local and 4-step distributed large-integer multiplication past `2^26` limbs.
- **`H5` (Pure `k = 1` Raw-`L_q` Cantelli Descent — Explicitly Declared Exploratory STNFS Variance Premise):** `384`-mesh lower sums (`+-16`-ULP `rho_lo < rho < rho_hi`), `100%` eager cofactorization, `BOUND_CELL = 5,000` Ir on both sieve and direct-ECM branches, linear union-bound missing-leaf retention `max(0, 1 - miss)` at `LX = 68.5`, `LOG2_RELS = 62.980`, `RET_MULT = (0.90, 0.82)`, `q15_retain16` omitted-mass laws, `UDUP = 0.40`, `BOOT_KEEP_FRAC = 0.90 < CP95 LB = 0.9395`, branch-tree-capped `L_min = 2` `12`-bit geometric ECM ladder (`sep = 0.5`, `p_r = 0.0043`, `2^26` Newton precheck), and the explicitly declared exploratory STNFS variance premise `Var(M | slot) <= max(2 * mu(b), 0.0005 * mu(b)^2)` on raw-`L_q` reached-slot usable-relation counts (`min_mu = 85.64 >= 50.0`), implying via exact one-sided Cantelli `P(M = 0 | slot) <= p_0(b) = max(0.0005, (1 + 0.5 * mu(b))^{-1})` across `N_slots = 682` tree slots (ZERO sibling retry, ZERO cross-slot independence; declared as an exploratory STNFS variance premise rather than deduced unconditionally from mean yield, toy Clopper-Pearson bounds, or lattice geometry alone, and empirically accompanied by 31k raw-`L_q` strict-descent `deg-(32, 24)` O/E & within-bin Fano `1.11..1.25 << 2.0`, 15.2k production `q >= 2^70` raw-`L_q` `deg-(72, 24)` geometry/norm checks, live regeneration inside `verify_q15_descent.py` with `>= 2x` modeled yield headroom, and production Part C `384/384`).

## 10. Memory & Reproducibility

Peak storage `< 2^120` B (`2^118.349` B replay tables during collection; all per-owner spaces `< 2^64` B asserted live in `audit_cm.py`). Reproduction: `python3 -I -B audit_cm.py <<< '{}'`.
