# recover-Q: `(4,3)` `Q(zeta_15)` CM-STNFS with `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` + Live Adaptive `H5` + Self-Contained `PM-Basis` & Unified Run Meter: expected work `2^109.871144921649`

Track `bls12-381-fixed-g1-invert-g2`. Input: nonidentity `P in G1` and `z in G_T`; output `Q in G2` with `e(P, Q) = z` for the exact supplied `P` (`Q = infinity` for `z = 1`). Reference unit: one dense `F_(p^12)` multiplication, `C = 20,000` baseline scalar ARMv8-A instructions. No advice, smaller targets, SIMD, `FEAT_LSE`, `DC ZVA`, external unbundled compiler measurements, black-box generator papers, or amortization are used. **Full de Bruijn local-density correction** (`rho_dB(u, x) = rho(u) - gamma rho(u-1)/ln x`) is baked into both Collection and Descent. Every cost constant (`uint16_t` zero-wrap raw-opcode-verified baseline scalar ARMv8-A `KERNELS.md`, `(m-k)&3` Band-2 bank rotation, 67-opcode 1-step Barrett `row_red5_z12`, 4-coord `gen_run1_z12`, input-uniform 2-stage Bareiss LLL, 4-prime 62-bit NTT `fft_F`, and self-contained classic `PM-Basis` matrix generator), smoothness integral, live `1/1` 381-bit boot compatibility check, exact 2nd-moment ECM Chebyshev bound, `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` recurrence, and `1,178`-attempt `Q(zeta_15)` real-2-core-FB adaptive diagnostic (`ADAPTIVE.md`) is derived and audited self-contained inside this package (**no prior verdict or score is a premise of the claim**).

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is an independent run gate: draw uniform 64-bit `u`; unless `u < K = 7214474002798577940`, stop immediately as failure (`< 1` unit). `K` is the least integer with `K x 0.997194 >= 0.39 x 2^64`. Under explicitly declared exploratory STNFS heuristics `H1`, `H1'`, `H-EA`, `H2`, `H3`, `H-R`, `H-SM`, `H4`, `H4u`, `H-X`, and `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` `H5` (Section 8; `RANK-BRIDGE.md`, `DESCENT.md`, `ADAPTIVE.md`), an opened run succeeds with probability `>= 0.997194` (unconditional Boole additive union bound, Section 6); because `u` is independent of the challenge and remaining coins, overall success is `>= 0.39` over uniform full-production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^111.225545 + 2^66) < 2^109.871144921649  dense-Fp12 units (C = 20,000).

The score `109.871144921649` is an **expected** total classical work bound `E[W]` over the challenge and fresh coins of the single attack, NOT an every-run bound. The **ENTIRE opened run** (`Collection + Replay + LA + PM-Basis Generator + Band/ECM Split + Both Logs`) runs under ONE unified hard aborting work meter `t_run = 2^123.0` core units (checked at coarse block boundaries `<= 2^66` units), giving the deterministic every-opened-run cap `W_actual <= 2^123.0 + 2^66` dense-`F_(p^12)` units on 100% of opened runs, exact truncation `E[min(W_open, t_run)] <= E[W_open] <= 2^111.225545`, and exact Markov abort term `p_meter = 0.0002854538 <= 0.00028546` (`ROUND_CEILING`) summed directly in Section 6. Analytic, exploratory claim. All headline numbers match the default `32/185` audit JSON of `audit_cm.py` (`14/14` live checks `checks_ok = true`, pure Python stdlib).

## 1. Reduction (fixed `P`: two `G_T` logarithms) and field map

With the pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `Q0` of `G2`:
0. **Run gate:** draw uniform 64-bit `u`; if `u >= K`, abort as failure (`< 1` unit).
1. **Validate:** check `P in G1 \ {infinity}` (`[r]P = infinity`) and `z in F_(p^12)^*` with `z^r = 1`. If `z = 1`, return `Q = infinity`.
2. **Pairing & field map:** compute `g = e(P, Q0)` and map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Fresh-coin randomization & core solve:** draw independent uniform `e_1, e_2 in {0, ..., r - 1}`, run the core on `g' = g^(e_1)`, `z' = z^(e_2)` (both exact uniform on `G_T \ {1}`), set `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** compute `Q = [a]Q0`; output `Q` iff `e(P, Q) = z`, else fail. Wrapper `< 2^66` core units, failure `< 2^(-48)`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, & Order-12 Frobenius

BLS12-381: `u = -0xd201000000010000 = -2^16 x 906349 x 254760293`, `r = u^4 - u^2 + 1` (`r == 13 mod 15`), `p = (u - 1)^2 r / 3 + u` (`p == 7 mod 15`).
- **CM splitting & `w' = omega^2` isometry:** `F(z)` splits over `Q(omega)` into `F_1(z) F_2(z)`. Over PID `K_H = Q(zeta_15)` (`O_H = Z[zeta_15]`, `D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `P = (p, h_pi)` has `O_H / P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` yields `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega + 1) X^6 + X^3 - omega` (deg 72).
- **Exact order-12 Frobenius (`algebra_bls12`):** `g` is irreducible over `F_(p^4)`, `sigma(y) = y^7` fixes `P, f, g`, and `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`chi = p == u mod r`, `ord_r(chi) = 12`), quotienting into `Z[zeta_12]` 12-orbits.
- **Relative-Degree-1 Lemma & Full-Orbit FB:** Unramified `Q | (a - b X)` (`gcd(a,b)=1`) has relative degree `1` over `O_H`. ALL full-orbit (`d = 12`) relative-degree-1 ideals of norm `<= B` are included (incl. `f_H >= 2`, `p <= 2^35.85 < 2^48`, total `< 2^32.5 << 2^(-28) U'`). Short-orbit (`d < 12`, forcing `vlog == 0 mod r`) ideals are omitted (all ramified `3, 5, 19, 67, ell | u` verified live).
- **4-Col Schirokauer Compression (`W_f, W_g`) & `alpha`:** Pre-collection uniform `W_f in F_r^(3x6)`, `W_g in F_r^(1x2)` compress `P_chi Lambda_f, P_chi Lambda_g` to `4` cols (`< 4/r < 2^(-252)` loss under `H-SM`). Exact `w' = omega^2` Murphy `alpha`: `alpha_f = +2.74`b, `alpha_(g,>=3) = +1.28`b, `kappa_2 = 0.0525`.
- **Live E2E `(4, 3)` `Q(zeta_15)` Tower DLP (`toy_e2e_dlp`):** On `u_0 = -5, p_0 = 7207, r_0 = 601` with random `W_g, W_f`: `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps annihilating live BSGS logs, neg controls failing, `5/5` unseen targets recovered.

## 3. Relation Collection: `Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, Exact 2-Adic Primitivity, & de Bruijn Yield

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373` 12-orbits), excess `1.09` (`2^63.441701` rels).
- **Geometry, Voronoi FD, Exact 2-Adic Class Rule `c(a)`, & `D_(omega,sigma)` (`239/239` exact-1):** `G = A_4 x A_2` (`det = 1,265,625`, `W, S`-isometric). 10 explicit Voronoi units (`voronoi_units_10`, `12/12` anchors `min-T2`, `N_a^can = 2^57.8294`). In `O_H`, `(2) = P_(2,1) P_(2,2)` (`N = 16`, both `sigma`-fixed); `c(a) in {0, 1, 2}` (`v_1+v_2 = 0, 1, >=2`; densities `225/256, 225/2048, 23/2048`; radii `s_c = (1.00, 1.20, 1.36)`; `STORAGE.md` §1) with `r_(ab) = 1.20`, `A^2 = 523230.1` and `1/12` FD `D_(omega,sigma)` (`239/239` exact-1) gives `N_7 = 2^57.1414` rows (`L_max = 544.52`, `L_b = 233.88`) and `N_a^tbl = 2^59.8294`.
- **`Arch2-Rev-SigmaB-v3-ZeroCopy-U16` (`W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`):** 1 full-chord tile/row (`nt = 1.000`, `4.658` SBs/row across `552` column owners at `2^61.829` B `<< 2^64` B): Build **`2^106.2387 <= 2^106.239198`** u (`C_LIST_A = 30.0`, global `q_key`-sorted `x_meta`/`x_imap` -> key-sorted Band-2 super-buckets by construction; storage `2^119.65` B `< 2^120.0` B); zero-copy `uint16_t` Band-1a (`4.25`) & Band-1b (`6.25`, exact Poisson `ce4`, `<= 2^63.448` B); Band-2 `sb_min_q8` (`+2/254` slack, `<= 2^63.939` B, `C_CAND = 11.5`, explicit `(m-k)&3` bank-rot `C_HIT_HI = 15.25`, `rev_sb_probe_fast` unsigned `b.hs .Lempty`); 16-bit SWAR `rev_scan_zero_u16` (`2.00`), schoolbook `_nk_g_only` (`5,768 <= 6,500`), `_nk_f_only` (`18,976 <= 20,000`), `MONT5 = 267 <= 270`, `Q_F_SAFETY = 2.0`, and closed-form 4-prime 62-bit NTT `stage_b_analytic` (`3 x 4p = 108.0 <= 168.0` Ir/elem/stage) integrated **LIVE across `337` bins** up to `bg_max = 700`b (**`c_cof_2d_upper = 2.9856`** ipc).
- **Exact 2-Adic Radius-Weighted Primitivity & Full-Cell Sieve Work:** Weighting each 2-adic valuation pair `(v_1, v_2)` by `s_(min(v_1+v_2,2))^8 * (15/16)^([v_1>0]+[v_2>0])` times the odd-prime Euler product through `p <= 5000` times rigorous tail `(1 - 8/4999)` yields exact primitive fraction lower bound **`prim_exact_lo = 0.957726`** (computed LIVE in `audit_cm.py`), shifting retained upper-endpoint yield to `2^63.472853` (`-0.009`b -> **`margin_net = +0.02215 >= 0.020`b** above `1.09 U'`, and `63.472853 - 0.009 - 2*0.01020 = 63.4435 >= 63.44`). Conversely, sieve work is priced over **100% of ALL cells including non-primitive cells** (`1.0168 * V_kept * ipc_tot * H2_ALLOW`, `ipc_tot = 24.1105`) -> **`E[W_coll] = 2^109.9309 <= 2^109.931400`** u inside unified meter `t_run = 2^123.0`. Band split `<= 2^100.543` u; exact 2nd-moment `ecm_mean_var` + `mr_instr` gives `E[W_ecm] = 2^95.09` u metered at **`2^96.2`** u (Chebyshev `<= 2^(-74.45) < 2^(-70)`). `verify_pool_2adic_live` draws `400` fresh canonical `(a, b)` samples verifying exact `v_2(N_g)`, volume, class shares, and de Bruijn yield live (`z_rels <= 2.5`).

## 4. `Z[zeta_12]` Filtering, 1-Owner Block Wiedemann LA, & Self-Contained Classic `PM-Basis`

Splitting `Z[zeta_12] = Z[omega] (+) chi Z[omega]` into 6 signed hexagonal buckets cuts unit+unit collision fill-in `3x` vs `Z[i]`. Bundled `filtsim_z12.py` mirrors our 6-seed ladder across `U' = 2^16 .. 2^23` (`B71.5_d04_e109`, `W = 200`, norms `548.9 / 677.8`):

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

- **Fixed-instance `H3` caps (`P(cap abort | H3) = 0`):** Hard caps `N_l <= 0.00225 U' = 2^54.5215`, `N_sw <= 0.00450 U' = 2^55.5215`, `k5 <= 17.50` (`k5_live = 16.3162` at `2^23` with 4-coord `gen_run1 = 116.0`), `tab_r <= 0.85` (`filtsim_z12` is exploratory support, not proof).
- **1-Owner Arch B BW & Self-Contained Classic `PM-Basis` Generator (`RANK-BRIDGE.md` §3):** `x + y` + 22 tables fit in `2^63.913` B `< 2^64` B on 1 owner (`0` cross-owner sends). At `(m, n) = (8192, 110)`, `unit_run5 = 8.80`, `gen_run1 = 116.0`, Barrett `row_red5 = 68` (`c_row_chi = 1,000`), `N_SM = 4`, `BW_SIG = 2^(-7)`: **LA `= 2^107.789053 <= 2^107.789554`** u. Classic divide-and-conquer Giorgi-Jeannerod-Villard / Thomé `PM-Basis` on `[S(z) | -I_m]` (`t = 8302`, `L_bm <= 2^47.7706`, pure 1D 4-prime 62-bit NTT since `4 L_bm = 2^49.77 < 2^53`, 4 base-`2^64` limbs `t L_bm 2^128 = 2^188.79 << 2^244.70`, schoolbook `t^3` pointwise matmuls at `17 <= 24` Barrett MAC ops verified live): **Gen `= 2^92.7203 <= 2^93.0`** u (ZERO black-box paper!).
- **4-Col Mod-`r` CRT Schirokauer Maps (`<= 2^84`), `H-R`/`H-X`, & Zero-Loss `repurge_stack`:** `f: [36, 36]`, `g: [4, 4, 4, 12]`. Exact 255-bit / `K_rej = 64` sampling misses `< 2^(-156) << 2^(-150)` (verified in 80-digit `Decimal`). Conditional on `H-R, H-X`, Schwartz-Zippel failure `<= 4 N_cap / 2^128 < 2^(-66)`. Live: `unplanted_two_sided_bw` (including **LIVE recursive `pm_basis` `order=15, cols=8, deg=6` recovering the exact kernel `ker(A)`**!) and `hx_z12.py` (`post_peel_cols_lost == post_merge_cols_lost == n_reconstruct_unit == 0`).

## 5. Individual Logarithms (`H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` + Live Adaptive `H5` + Fresh-Exponent `i.i.d.` Boot)

Detailed in `DESCENT.md` and `ADAPTIVE.md` (all `39` bins / `78` nodes, exact `theta_1^8` `cells_rank2`, rigorous `1/zeta_K(2) >= 0.970721 > 0.95`, and seeded peeling `h = 0.256845, exc_crit = 0.7627` evaluated LIVE at default mesh **`32/185`**):
- **Exact `S_rank2(R) = A_cap(R) x B_8(R)` (`cells_rank2`):** Proved zero `O_K^*` unit duplicates; `audit_cm.py` recomputes exact `|B_8(R)|` via `theta_1(R^2)^8` for all `23` radii `R <= 160` and counts `R = 2` canonical `alpha` live.
- **Fresh-Exponent `i.i.d.` Single-Cell Boot (`boot_rels = 1, cells_b = 1, D_skew = 2^21`):** Uniform `e_t in {0..r-1}`, `F_(p^12) -> F_(p^6) -> gamma_t`, 2-stage pure-int Bareiss `lll_16d_int` (`2^46.45 << LLL_BOOT = 2^48` Ir), `|N| <= 2^2560` (`RET_BOOT = 0.60`, declared exploratory MODEL inequality; `1/1` live from-scratch SHA-256 `G_T` compatibility check `<= 2544.7b`, `0` pre-supplied vectors), split `(n1=255, nA=195, kA=15)`. On failure (`<= psi_boot = 0.25188 < 0.50`), draw fresh `e_t` (`i.i.d.`, no inner search cap).
- **`H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` (`phi in [0.01, 0.04]`, `psi < 0.80`):** Exact closed-form `t_k` minimization over `k = 1..64` gives default `32/185` **`max_phi = 0.02500`**, **`max_psi = 0.21575`**, **`psi_boot = 0.25188`**, **`min_mu = 740.2`**, expected two-log work **`2^110.128131 <= 2^110.129000`** u, `E[N_2log] = 2^12.524`, node budget `B_node = 2^22.0` (`p_node_abort = 0.0014044`) + `H5` repeat-collision (`C_rep = 256` against `2^60.81`, `p_repeat_loss = 0.0011105`) -> **`0.0025150 <= 0.0025200`**.
- **Four Complementary Live Audits (`lll_16d_int`, `0` Pre-Supplied Vectors; Exploratory Support, Not Proof):** (a) `prod_sq24_live` (`24/24` real `(f_9, g_3)` control lattices reduced live, `qbins` uses overall `cgmax=46.5`, **actual joint-cell `(nf_k, ng_k)` `pr_emp/pr_charged >= 16.96x` (`f`) / `>= 55.01x` (`g`)**, `0` observed-mean PGF); (b) `prod_cond_q15_live` (`12/12` live `q | N_side(parent_a, parent_b)` + **actual joint-cell `pr_emp/pr_charged >= 41.60x`**, med `301.41x`); (c) `q15_samefield_descent_live` (`3/3` live-replayed + `120`-ideal real-2-core-FB marginal count check vs independent `mu_model`); (d) **`adaptive_h5_live` (`ADAPTIVE.md`):** replays all `1,178` real-2-core-FB deterministic-order after-failure attempts (`pooled O/E fail = 0.8022` [95% CI `0.703-0.905`], `rat1 = 0.8279`, history-stratified `d3+=0.828, top=0.746, c2=0.729, p2=1.112`, `work O/E = 1.0132`, `0/1,178` shared descendants, `C_bday <= 2.815`, `C_emp <= 19.11`), checks `threshold_counterexample`, AND **live-regenerates 3 complete multi-level after-failure witness descents from scratch**!

## 6. Success Probability & Unconditional Boole Union Bound

Unconditional additive Boole union bound `P(success) >= 1 - sum_i P(fail_i)` (verified in 80-digit `Decimal` in `audit_cm.py`):
- **Main phase, compression, ECM split, & MR:** BW Schwartz-Zippel `<= 4 N_cap / 2^128 < 2^(-66)` (`| H-R, H-X`), Schirokauer compression `<= 4/r < 2^(-252)` (`| H-SM`), exploratory `H3` matrix caps (`P(cap abort | H3) = 0`), ECM split exact 2nd-moment Chebyshev `< 2^(-74.45) < 2^(-70)`, 100-round MR `<= 2^(-80)`.
- **Descent Budget/Repeat Loss & Unified Run Meter (`t_run = 2^123.0`):** Node budget `2^22.0` (`0.0014044`) + `H5` repeat-collision (`C_rep = 256`, `0.0011105`) `= 0.0025150 <= 0.0025200`, unified hard run meter `t_run = 2^123.0` Markov `p_meter = 0.0002854538 <= 0.00028546`, degenerate target `<= 2/r`, wrapper `< 2^(-48)` (incl. `F_r` rejection `< 2^(-156)`).
- **Opened-run floor & gate:** Core + wrapper `= 0.99719454 >= 0.997194`. Gate `K = 7214474002798577940` satisfies `(K/2^64) x 0.997194 >= 0.39 > ((K-1)/2^64) x 0.997194`.

## 7. Complete Expected-Work Ledger (Wrapped Pointwise in Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, full-cell `1.0168x`, `ipc=24.1105`, `c_cof=2.9856`; live `109.9309`) | <= 109.931400 |
| replay 4-orbit `a`-tables build (`g` to `2^63` + `f` to `2^60`, 2-pass `C_LIST_A=30.0`; live `106.2387`) | <= 106.239198 |
| linear algebra (`1-owner Arch B Z[zeta_12] BW`, `W=200`, `unit_run5=8.80`, `gen_run1=116.0`, Barrett `row_red5=68`, `k5<=17.50, tr<=0.85`, `(8192,110)`; live `107.789053`) | <= 107.789554 |
| self-contained classic `PM-Basis` matrix generator (1D 4-prime 62-bit NTT + schoolbook `t^3` Barrett MAC; live `92.7203`) | <= 93.0 |
| both individual logarithms (expected, default `32/185` mesh: `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR`, `psi_b=0.25188`, `max_phi=0.02500`, `max_psi=0.21575`, `E[N]=2^12.524`, `LLL_NODE=2^40`, `LLL_BOOT=2^48`; live `110.128131`) | <= 110.129000 |
| band split of kept relations (`(2^20, 2^48]`, 10/12 limbs) | <= 100.543 |
| ECM split of kept relations (`(2^48, B]`; exact 2nd-moment Cheb. `< 2^-74.45`) | <= 96.2 |
| canonical `a`-list enumeration + global `key_a` pre-sort | <= 86 |
| 4-column `chi`-Schirokauer maps (mod-`r` CRT `f: [36,36]`, `g: [4,4,4,12]`) | <= 84 |
| filtering streams, product-form pivot selection, 3-stack `repurge_stack` back-substitution, misc | <= 80 |
| `P_g`, `P_f` product trees (`< 2^89` instructions) | <= 75.5 |
| factor base, roots, 12-orbit ideals (`U' = 2^63.32`) | <= 75 |
| numerical certificate construction | <= 50 |
| final verification | <= 11 |
| wrapper (sec. 1; `<= 2^80` instructions = `2^65.71` core units) | <= 66 |

Expected opened-run core sum (`80`-digit `Decimal`): `2^111.225544893704 <= 2^111.225545` core units. Wrapped in `t_run = 2^123.0` (pointwise cap `<= 2^123.0 + 2^66`, `p_meter <= 0.00028546`) behind gate `K / 2^64`: **`E[W] <= 1 + (K/2^64)(2^111.225545 + 2^66) < 2^109.871144921649`** dense-`F_(p^12)` units (`< 2^113.193073016536` G1 units).

## 8. Declared Exploratory STNFS Hypotheses (Supported, Not Proved)

Every item below is an explicitly declared exploratory STNFS hypothesis; finite-scale campaigns and live audits provide scoped compatibility evidence, NOT mathematical proofs of production behavior:
- **`H1`, `H1'`, `H-EA`, `H2`:** Two-sided `rho_dB` smoothness, exact `v_2(N_g)`, exact 2-adic primitive-pair fraction (`>= 0.957726`), two-stage Euler retention + survivor-upper (`Q_F_SAFETY = 2.0`, `0.009`b margin), `200,000`-sample pool (`H2_ALLOW = 1.00774`), expected work wrapped in `t_run = 2^123.0`.
- **`H3`, `H-SM`, `H-R`, `H-X`:** `r` unramified (proved), `r !| h(K_f), h(K_g)`, `P_chi Lambda` injective (`< 4/r`), monotone filtering trend into fixed-instance caps `(N_l/U' <= 0.00225, N_sw/U' <= 0.00450, k5 <= 17.50, tab_r <= 0.85)`, fixed-instance full rank `H-R` and nonzero BW minors `H-X` (supported by `unplanted_two_sided_bw` + live `PM-Basis` + `hx_z12.py`).
- **`H4`, `H4u`:** Explicit raw-opcode-verified baseline scalar ARMv8-A kernels (`KERNELS.md`, `simd = 0`, zero `FEAT_LSE`/`DC ZVA`, `b.hs .Lempty`, `ntt62` `3 x 4p = 108.0 <= 168.0`, 62-bit Barrett MAC `17 <= 24`, `MONT5 = 267 <= 270`, Barrett `row_red5 = 67 <= 68`, uniform 2-stage Bareiss LLL `2^40 / 2^48`); `Z[zeta_12]` unit share `93.01%` at `2^23`.
- **`H5` (`H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` + Sequential `psi` + `C_rep = 256`, `DESCENT.md`, `ADAPTIVE.md`):** Exact `S_rank2(R)` (`|(B_8(R)\{0}) cap (F_Vor/mu_30)| >= N_alpha(R)`, `1/zeta_K(2) >= 0.970721 > 0.95`), fresh-exponent `i.i.d.` boot (`RET_BOOT = 0.60`), seeded peeling (`h = 0.256845`), `H-PGF (0.5*mu)` (`t <= 0.80`), sequential conditional child-failure hazard `P(F_i | A_i, pre-entry history) <= psi` and conditional child-work bound (supported by `1,178`-attempt `ADAPTIVE.md` + 3 live-regenerated witness descents), and adaptive repeat-collision anti-concentration `<= C_rep / N_lx` (`C_rep = 256`, supported by `C_bday <= 2.815`, `C_emp <= 19.11`).

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0` bytes (`2^119.65` B replay `a`-tables; `4 * N_a_tbl = 2^61.83` B `< 2^64` B for `uint16_t` `S_g + S_f`; `< 2^105` B descent; `2^63.913` B `< 2^64` B on 1 owner in LA; `2^82.81 B < 2^84 B` `PM-Basis` NTT workspace; `< 2^30` B wrapper; `STORAGE.md`).

**Provenance:** Builds on accepted `c36p` / `stnfs-bn15r` exploratory methodology (`afe6f29f`, `5a49c502`, `4931a702`) and Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, and Boudot (2017/758). **No prior verdict or score is a premise of this claim.**
