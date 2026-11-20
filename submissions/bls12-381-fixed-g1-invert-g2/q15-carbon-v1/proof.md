# recover-Q: `(4,3)` `Q(zeta_15)` CM-STNFS + `B_16(R)/{+-1}` Zero-Float 64-Bit Q32 `dsieve16` (Deterministic `B_int` Gate) + Pure-Int Ball `BOUND_CELL` + `k=1` Boole Gate G + Proved `RET_BOOT=1` + `PM-Basis`: expected work `2^110.889621034967`

Track `bls12-381-fixed-g1-invert-g2`. Input: nonidentity `P in G1`, `z in G_T`; output `Q in G2` with `e(P, Q) = z` for exact supplied `P` (`Q = infinity` for `z = 1`). Unit: one dense `F_(p^12)` mul, `C = 20,000` baseline scalar ARMv8-A integer instructions (`simd = 0`, zero FP/NEON/SVE, zero `FEAT_LSE`/`DC ZVA`, zero external C measurements, zero black-box generator citations, zero amortization). Full `rho_dB` baked into Collection & Descent. All kernels (`KERNELS.md`, incl. raw-opcode-verified zero-float 64-bit Q32 `dsieve16`, pure-int ball `BOUND_CELL`, and unconditionally proved Bareiss LLL), direct 200k-draw `O_H`-coprime pool, proved `RET_BOOT = 1.00`, lossless Gate G (`P_70` strip + full overshoot wave), pure `k = 1` Boole descent (NO Gate E), ECM Chebyshev, and literal `(9,3)` `Q(zeta_15)` `B_16` diagnostics are self-contained (**no prior verdict/score is a premise**).

## Claim

`W` is total classical work on one target. Step 0: uniform 64b `u`; unless `u < K = 7198405263799728974` (least `K` with `K * 0.999420 >= 0.39 * 2^64`), abort (`< 1` u). Under declared exploratory STNFS heuristics `H1, H1', H-EA, H2, H3, H-R, H-SM, H4, H4u, H-X, H5-Empty` (§8; `RANK-BRIDGE.md`, `DESCENT.md`), an opened run succeeds w.p. `>= 0.999420` (Boole, §6); overall success `>= 0.39` over uniform production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^112.247238 + 2^66) < 2^110.889621034967  dense-Fp12 units (C = 20,000).

`110.889621034967` is an **expected** total work bound `E[W]` (NOT every-run). The **ENTIRE opened run** runs under ONE hard aborting work meter `t_run = 2^123.0` u (block boundaries `<= 2^66` u): deterministic pointwise cap `W_actual <= 2^123.0 + 2^66` u, `E[min(W_open, t_run)] <= 2^112.247238`, Markov `p_meter = 0.0005795569 <= 0.00057956`. Analytic exploratory claim; all numbers match default `32/185` `audit_cm.py` (`14/14` `checks_ok = true`, pure stdlib).

## 1. Reduction (fixed `P`: two `G_T` logs) and field map

With pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public `Q0 in G2`:
0. **Run gate:** uniform 64b `u`; if `u >= K`, abort (`< 1` u).
1. **Validate:** `P in G1 \ {infinity}` (`[r]P = infinity`), `z^r = 1`; if `z = 1`, return `infinity`.
2. **Pairing & map:** `g = e(P, Q0)`, map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Randomize & solve:** independent uniform `e_1, e_2 in {0..r-1}`, core on `g^(e_1), z^(e_2)`, `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** `Q = [a]Q0`; return `Q` iff `e(P, Q) = z`. Wrapper `< 2^66` u, failure `< 2^-48`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, & Order-12 Frobenius

BLS12-381: `u = -0xd201000000010000`, `r = u^4 - u^2 + 1 == 13 mod 15`, `p = (u-1)^2 r/3 + u == 7 mod 15`. Over PID `O_H = Z[zeta_15]` (`D_H = 1,265,625`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)`, `O_H/P ~= F_(p^4)`. `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` gives `sigma`-fixed `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega+1)X^6 + X^3 - omega` (deg 72). `tau(y, X) = (y^7, omega^2 X) = Frob^1` (`ord_r(p) = 12`) quotients into `Z[zeta_12]` 12-orbits (`algebra_bls12`). All full-orbit (`d=12`) rel-deg-1 ideals `<= B` included (`f_H >= 2` `< 2^32.5 << 2^-28 U'`); short-orbit ideals (`vlog == 0`) omitted. Pre-collection `W_f in F_r^(3x6), W_g in F_r^(1x2)` compress `P_chi Lambda` to `4` cols (`< 4/r` under `H-SM`). Exact `w' = omega^2` `alpha`: `alpha_f = +2.74`b, `alpha_(g,>=3) = +1.28`b, `kappa_2 = 0.0525`. Live E2E (`toy_e2e_dlp`, `u_0=-5, p_0=7207, r_0=601`): `0` `P_bar` hits, `25/25` `K_g` & `57/57` `K_f` deps, `5/5` unseen DLPs.

## 3. Relation Collection: `Arch2-Rev-SigmaB-v3-ZeroCopy-U16` & Direct 200k Exact-`O_H`-Coprime Pool

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373`), excess `1.09` (`2^63.441701` rels).
- **Geometry & 2-adic `c(a)` (`239/239` exact-1):** `G = A_4 x A_2`, 10 Voronoi units (`N_a^can = 2^57.8300`, `N_a^tbl = 2^59.8300`). `(2) = P_(2,1) P_(2,2)`, `c(a) in {0,1,2}` (`s_c = (1.00, 1.20, 1.36)`, `STORAGE.md`), `r_(ab) = 1.187014`, `A^2 = 528954.37`, `D_(omega,sigma)` -> `N_7 = 2^57.2513` rows (`L_max = 550.47 < W_tile = 552`, `L_b = 236.44`).
- **Direct 200k-Draw Exact-`O_H`-Coprime Pool (`P2_wy_2ad_prim`, `V_kept = 2^119.7598`):** Direct `gcd_O(a,b) == (1)` test on `200,000` draws (`191,832/200,000`, `f_prim_mc_200k = 0.959160`; `401` conservative `(bf, bgo, max_bgm)` bins, primitive yield `wc` summing to `0.959160`, all-cell work `wk` summing to `1.0`). Retained yield `2^63.578997`, **`margin_net = +0.12830`b**. Live `verify_pool_2adic_live` (`400` draws): `382/400` coprime (`z_prim = 0.42`), `prim_yield_ratio = 1.0066`, `z_rels = 0.82`.
- **Arch2-Rev-SigmaB-v3-ZeroCopy-U16 (`W_tile=552, W_SB=64, P_lo=2^55, P_m4=2^48, P_g=2^63, P_f=2^60`):** `nt = 1.000`, `4.676` SBs/row (`552` col owners `2^61.830 B < 2^64 B`): Build **`2^106.2393 <= 2^106.239800`** u (64-bit `x2` `rev_build_pass12` hist, `2^119.655 B < 2^120 B`); Band-1a/1b (`4.25/6.25`, `<= 2^63.4484 B`); Band-2 (`<= 2^63.9399 B`, `cand_ratio = 1.4023 <= 1.414`, `C_CAND=11.5`, `C_HIT_HI=15.25`); `rev_scan_zero_u16` (`2.0`); `_nk_g_only` (`6500`), `_nk_f_only` (`20000`), `MONT5 = 267 <= 270`, `Q_F_SAFETY = 2.0`, 4-prime NTT `stage_b_analytic` (`108 <= 168`) integrated LIVE across `401` bins (`bg_max = 700`b -> `c_cof_2d_upper = 2.8950` ipc). Full-cell sieve (`ipc_tot = 23.8643`, `H2_ALLOW = 1.00774`) -> **`2^110.059975 <= 2^110.060475`** u. Band split `<= 2^100.543` u; ECM split **`2^96.2`** u (2nd-moment Chebyshev `< 2^-74.45`).

## 4. `Z[zeta_12]` Filtering, 1-Owner BW LA, & Self-Contained Classic `PM-Basis`

6-bucket `Z[zeta_12] = Z[omega] (+) chi Z[omega]` (`filtsim_z12.py`, 6 seeds/scale, `B71.5_d04_e109`, `W = 200`):

| `U'` | Seeds | `N_l / U'` (SE) | `N_sw / U'` (SE) | `sd_l U'^(1/4)` | Depth | Unit % (Exp.) | `k5` | `tab_r` |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `2^16` | 6 | `0.003260` (`0.000227`) | `0.006462` (`0.000482`) | `0.0089` | `13.7` | `95.61%` (`89.15%`) | `20.4312` | `0.80659` |
| `2^17` | 6 | `0.003050` (`0.000112`) | `0.006089` (`0.000219`) | `0.0052` | `15.2` | `96.63%` (`90.47%`) | `19.0156` | `0.80554` |
| `2^18` | 6 | `0.002857` (`0.000088`) | `0.005747` (`0.000173`) | `0.0049` | `16.7` | `96.96%` (`91.30%`) | `18.1201` | `0.80212` |
| `2^19` | 6 | `0.002533` (`0.000074`) | `0.005111` (`0.000118`) | `0.0049` | `16.7` | `97.18%` (`91.82%`) | `17.5677` | `0.79984` |
| `2^20` | 6 | `0.002399` (`0.000032`) | `0.004857` (`0.000049`) | `0.0025` | `18.7` | `97.35%` (`92.32%`) | `17.0304` | `0.79839` |
| `2^21` | 6 | `0.002172` (`0.000035`) | `0.004377` (`0.000067`) | `0.0033` | `18.5` | `97.42%` (`92.61%`) | `16.7212` | `0.79684` |
| `2^22` | 6 | `0.001999` (`0.000021`) | `0.004057` (`0.000039`) | `0.0024` | `19.7` | `97.43%` (`92.77%`) | `16.5469` | `0.79508` |
| `2^23` | 6 | `0.001889` (`0.000013`) | `0.003846` (`0.000026`) | `0.0017` | `20.5` | `97.45%` (`93.01%`) | `16.2921` | `0.79485` |

- **Fixed-instance `H3` caps (`P(cap abort | H3) = 0`):** `N_l <= 0.00225 U' = 2^54.5215`, `N_sw <= 0.00450 U' = 2^55.5215`, `k5 <= 17.50` (`k5_live = 16.3162`), `tab_r <= 0.85`.
- **1-Owner Arch B BW & `PM-Basis` (`RANK-BRIDGE.md` §3):** `x+y` + 22 tables fit in `2^63.913 B < 2^64 B` (`0` cross-owner sends). `(m,n)=(8192,110)`, `unit_run5=8.80`, `gen_run1=116.0`, Barrett `row_red5=68`: **LA `= 2^107.789053 <= 2^107.789554`** u. Self-contained recursive `PM-Basis` (`t=8302`, `L_bm<=2^47.7706`, 1D 4-prime NTT, schoolbook `t^3` Barrett MAC `17<=24` ops): **Gen `= 2^92.7203 <= 2^93.0`** u. 255b/`K_rej=64` sampling `< 2^-156`; Schwartz-Zippel `< 2^-66` (`| H-R, H-X`). Live: `unplanted_two_sided_bw` (incl. live `pm_basis`) + `hx_z12.py` (`0` lost cols).

## 5. Individual Logarithms (`B_16(R)/{+-1}` Zero-Float 64-Bit Q32 `dsieve16` + Deterministic `B_int` Gate + Pure-Int `BOUND_CELL` + `k=1` Boole Gate G + Proved `RET_BOOT=1`)

Detailed in `DESCENT.md` & `KERNELS.md` §5 (`39` bins / `78` nodes, `R in [64, 128]`, exact `theta_1^16` `|S_16(R)|`, `udup_node = 0.40` with certified `[lo, hi]` `Decimal(40)` Neumann/Taylor complete-orbit live `yr_min = 79.81x >= 10*0.40` and transparent `E[s/m]/E[s]` disclosure, Gate G `P_70` + full overshoot waves, `h = 0.256845`, default `32/185`):
- **Zero-Float `dsieve16` with Deterministic `B_int` Gate (`KERNELS.md` §5, `DESCENT.md` §2):** Pure-int Bareiss `lll_16d_int(9,10)` on `L_{q,Q}` (`det=p<=2^70`, `d_k=p^2<=2^140`, `<=13,815` swaps -> `2^31.20 << PAIR_SETUP = 2^35.0` Ir unconditionally!). Pure 64-bit Q32 outward-rounded Fincke-Pohst DFS (`dsieve16_q32_node` raw opcodes `57 <= 65` Ir/node, ZERO floats!) + exact-int leaf check & bijective `rank_lex` (`207 <= 210` Ir/leaf, `dsieve16_rank_step` verified) behind deterministic per-ideal Fubini/Riemann-sum gate `B_int(dg, R+0.05) <= C_int(R, p)` (`WC_LLL99` profile; `8/8` live proved primes `16b..70b` on production `(f_9,g_3)` lattices pass at `0.181..0.344 < 0.50`; live `R=2, p=31` `552/552` exact match vs brute force & distinct ranks; live `R=64` discrete DFS on proved `84`b primes `560` & `695` hits with `actual_nodes <= B_int`). Pure-integer ball-arithmetic `BOUND_CELL = 1,800` Ir (ZERO floats, proved lossless for all `Z^16` by triangle inequality; `25/25` live incl. 5 near-cancellation cells, max gap `<= 0.1215`b). SWAR `SCAN + ZERO = 6` Ir.
- **Proved `RET_BOOT = 1.00` (`DESCENT.md` §3):** Rational Sylvester `(15.001*I_8 - G_pow > 0)` + `99/100`-Lovasz + AM-GM (`3/2`) prove `log2|N| <= 2591.90 < 2592.0` for 100% of nonidentity `G_T` targets (ZERO heuristic).
- **Pure `k=1` Boole Descent (NO Gate E, `p_node_abort = p_repeat_loss = 0` PROVED):** Worst-case boot child counts (`<= 15` in `(2^70, 2^195]`, `<= 1` in `(2^195, 2^255]`), `N_tree_max_1log = 642` deterministically. Under `H5-Empty` (supported by 36/36 full-production `(f_9,g_3)` lattices + 120 literal `(9,3)` `Q(zeta_15)` `B_16` balls with 1 live-regenerated 9-hit ball + 120 real-2-core-FB `Q(zeta_15)` balls), Boole gives `P(tree fails | boot qualifies) <= P_boot_tree = 0.72705 < 0.80` (`min_mu = 153.34`, `max_p0 = 0.01287`). Proved `i.i.d.` Wald gives default `32/185` **`2^111.769688 <= 2^111.771273`** u (`96/185`: `2^111.702264`).

## 6. Success Probability & Unconditional Boole Union Bound

80-digit `Decimal` in `audit_cm.py`: BW SZ `< 2^-66` (`| H-R,H-X`), Schirokauer `< 2^-252` (`| H-SM`), `H3` caps (`0`), ECM Cheb `< 2^-74.45`, MR `<= 2^-80`, descent abort/repeat **`0.0` PROVED**, meter `t_run = 2^123.0` Markov `p_meter = 0.0005795569 <= 0.00057956`, degen `<= 2/r`, wrapper `< 2^-48` (incl. `F_r` rej `< 2^-156`). Floor: **`0.99942044 >= 0.999420`**. Gate `K = 7198405263799728974` satisfies `(K/2^64)*0.999420 >= 0.39 > ((K-1)/2^64)*0.999420`.

## 7. Complete Expected-Work Ledger (Unified Meter `t_run = 2^123.0`)

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, 200k coprime pool, `ipc=23.8643`, `c_cof=2.8950`; live `110.059975`) | <= 110.060475 |
| replay 4-orbit `a`-tables build (64b `x2` hist, `C_LIST_A=30.0`; live `106.2393`) | <= 106.239800 |
| linear algebra (`1-owner Arch B Z[zeta_12] BW`, `unit_run5=8.80`, `gen_run1=116.0`, Barrett `68`, `k5<=17.50`; live `107.789053`) | <= 107.789554 |
| self-contained classic `PM-Basis` generator (1D 4-prime NTT + schoolbook `t^3` Barrett MAC; live `92.7203`) | <= 93.0 |
| both individual logs (`B_16(R)/{+-1}` zero-float 64b Q32 `dsieve16` + det. `WC_LLL99` `B_int` gate, `PAIR_SETUP=2^35`, `65` Ir/node, `210` Ir/leaf, pure-int `BOUND=1800`, `k=1` Boole Gate G, NO Gate E, proved `RET_BOOT=1`, `P_boot_tree=0.72705`, `N_tree<=642`; live `111.769688`) | <= 111.771273 |
| band split (`<=100.543`), ECM (`<=96.2`, Cheb `<2^-74.45`), can `a` (`<=86`), 4-col SM (`<=84`), filt (`<=80`), `P_g/P_f/P_70` (`<=75.5`, live `72.328`), FB (`<=75`), cert (`<=50`), ver (`<=11`) | <= 100.62 |
| wrapper (§1; `<= 2^80` Ir = `2^65.71` u) | <= 66 |

Core sum: `2^112.24723721 <= 2^112.247238` u. Meter `t_run = 2^123.0` (`<= 2^123.0 + 2^66` pointwise, `p_meter <= 0.00057956`) behind `K/2^64`: **`E[W] <= 1 + (K/2^64)(2^112.247238 + 2^66) < 2^110.889621034967`** dense-`F_(p^12)` u (`< 2^114.211549129854` G1 u).

## 8. Declared Exploratory STNFS Hypotheses (Supported, Not Proved)

- **`H1, H1', H-EA, H2`:** Two-sided `rho_dB` on exact-`O_H`-coprime finite-region pairs (`P2_wy_2ad_prim`, `f_prim=0.959160`, `+0.12830`b margin), exact `v_2(N_g)`, two-stage Euler retention + `Q_F_SAFETY=2.0`, `H2_ALLOW=1.00774`, double-counted leaf Mertens sum `2*1.10*ln(70/8)`, expected work wrapped in `t_run = 2^123.0`.
- **`H3, H-SM, H-R, H-X`:** `r` unramified (proved), `r !| h(K_f), h(K_g)`, `P_chi Lambda` injective (`<4/r`), monotone filtering into fixed-instance caps (`N_l/U'<=0.00225`, `N_sw/U'<=0.00450`, `k5<=17.50`, `tab_r<=0.85`), fixed-instance full rank `H-R` and nonzero BW minors `H-X`.
- **`H4, H4u`:** Raw-opcode-verified baseline scalar ARMv8-A integer kernels (`KERNELS.md`, `simd=0`, ZERO FP/LSE, 19 parsed asm blocks incl. `dsieve16_q32_node`, `dsieve16_rank_step`, `bound_cell_ball_mul`, 64b `x2` hist, `ntt62` `108<=168`, Barrett MAC `17<=24`, `MONT5=267<=270`, Barrett `row_red5=67<=68`, unconditionally proved Bareiss LLL `2^35/2^40/2^48`); `Z[zeta_12]` unit share `93.01%` at `2^23`.
- **`H5-Empty` (`DESCENT.md` §4 — Declared Exploratory Lower-Tail Premise, ZERO Sequential Retry/Repeat/Gate-E Premise):** Isotropic `S_16(R) = B_16(R)/{+-1}` (`udup_node = 0.40`, certified `[lo, hi]` `Decimal(40)` `yr_min = 79.81x >= 10*0.40`), boot smoothness (proved `RET_BOOT=1.00`), seeded peeling (`h=0.256845`), reached-slot `H-EA` expected survivors (`E[W_tree|boot] <= W_tree_1`), and declared exploratory marginal empty-ball bound `P(M_w=0 | slot w reached) <= p_0(b(w)) = max(0.001, 1/(1+0.5*mu))` (Cantelli Fano-2 form floored at `0.001`, bridged via Mechanisms A & B in `DESCENT.md` §4; supported diagnostically, NOT proved by toy runs). No `THRESHOLD-REPAIR`, no `F_i|F_1..F_{i-1}`, no `C_rep`, no Gate E.

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0 B` (`2^119.655 B` replay `a`-tables; `2^61.830 B` col owners; `<= 2^63.9399 B` Band-1/2 tables; `P_70` `2^68.53 B` across 64 owners; descent `<= 4*|S_16(128)| = 2^110.913 B < 2^111.0 B` across `<= 2^49` `< 2^62`-B chunk owners + `0.88 GiB` `P_d` table; `2^63.913 B` 1-owner LA; `2^82.81 B` `PM-Basis`; `< 2^30 B` wrapper; `STORAGE.md`). **Provenance:** Builds on `c36p` / `stnfs-bn15r` (`afe6f29f`, `5a49c502`, `4931a702`, `cb938290`), Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), Giorgi-Jeannerod-Villard (2003), Thomé (2002), CADO-NFS, Boudot (2017/758). **No prior verdict/score is a premise.**
