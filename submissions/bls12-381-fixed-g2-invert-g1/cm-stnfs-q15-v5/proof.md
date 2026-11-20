# recover-P: (4,3) CM-STNFS over Q(zeta_15) with Arch2-Rev-SigmaB-v3 behind a run gate: expected work `2^109.249832659662`

Track `bls12-381-fixed-g2-invert-g1`. Input: a nonidentity Q in G2 and z in GT; output P in G1 with e(P,Q) = z for the exact supplied Q (P = infinity for z = 1). Reference unit: one dense F_p12 multiplication, C = 20,000. No advice, smaller targets, SIMD, `FEAT_LSE`, hand-written assembly, or amortization are used. **Full de Bruijn local-density correction** (`rho_dB(u, x) = rho(u) - gamma rho(u-1)/ln x`) is baked into both Collection and Descent. Every cost constant (including the closed-form non-recursive 1D/2D 4-prime 62-bit radix-2 NTT `fft_F`), smoothness integral, and descent recurrence is derived and verified self-contained inside this package (no prior verdict or external constant is a premise).

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is a run gate (Section 9): draw uniform 64-bit `u`; unless `u < K = 7205708882997339893`, stop immediately as failure (`< 1` unit). `K` is the least integer with `K x 0.998407 >= 0.39 x 2^64`. Under heuristics `H1`, `H1'`, `H-EA`, `H2`, `H3`, `H-R`, `H-SM`, `H4`, `H4u`, `H-X`, and `H5` (Section 8; `RANK-BRIDGE.md`, `THRESHOLD-REPAIR.md`), an opened run succeeds with probability `>= 0.998407` (strict additive union bound, Section 6); because `u` is independent of the challenge and remaining coins, overall success is `>= 0.39` over uniform full-production challenges and fresh coins, and

    E[W] <= 1 + (K / 2^64) x (2^110.6059865849 + 2^66) < 2^109.249832659662  dense-Fp12 units (C = 20,000).

Here `2^110.6059865849` core units (`C = 20,000` scalar AArch64 instructions) sums Section 7's every-run bounds with expected two-log descent (`2^107.449`, Section 5). The score `109.249832659662` is an expectation over the challenge and fresh coins of the single attack, not a bound on every opened run; every opened run is strictly metered at:

    W <= 2^119.5026900659 + 2^66 < 2^119.502690065901  dense-Fp12 units (C = 20,000).

Analytic, exploratory claim. All declared numbers round conservatively and are audited live by `audit_cm.py` (`checks_ok = true`, pure Python stdlib). Neither finite-scale simulator outcomes nor statistical confidence intervals are substituted for the algorithmic success probability; `H3` and `H5` are declared exploratory heuristics supported by the monotone `2^16..2^23` campaign and live audits.

## 1. Reduction (fixed `Q`: two `G_T` logarithms) and field map

With the pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generator `P0` of `G1`:
0. **Run gate:** draw uniform 64-bit `u`; if `u >= K`, abort as failure (`< 1` unit).
1. **Validate:** check `Q in G2 \ {infinity}` (`[r]Q = infinity`) and `z in F_(p^12)^*` with `z^r = 1`. If `z = 1`, return `P = infinity`.
2. **Pairing & field map:** compute `g = e(P0, Q)` and map `(g, z)` into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
3. **Fresh-coin randomization & core solve:** draw independent uniform `e_1, e_2 in {1, ..., r - 1}`, run the core on `g' = g^(e_1)`, `z' = z^(e_2)` (both exact uniform on `G_T \ {1}`), set `a = (e_1 log(z')) (e_2 log(g'))^(-1) mod r`.
4. **Output & check:** compute `P = [a]P0`; output `P` iff `e(P, Q) = z`, else fail. Wrapper `< 2^66` core units, failure `< 2^(-48)`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, & Order-12 Frobenius

For BLS12-381, `u = -0xd201000000010000 = -2^16 x 906349 x 254760293`, `r = u^4 - u^2 + 1` (255-bit prime, `r == 13 mod 15`), `p = (u - 1)^2 r / 3 + u` (381-bit prime, `p == 7 mod 15`).

- **CM splitting & `w' = omega^2` isometry:** BLS12-381's `F(z) = (z - 1)^2 (z^4 - z^2 + 1) + 3 z` splits over `Q(omega)` into `F_1(z) F_2(z)` (`F_1(z) = z^3 - z^2 - (omega + 1) z - omega`). Over `K_H = Q(zeta_15)` (`O_H = Z[zeta_15]`, `D_H = 3^4 x 5^6 = 1,265,625`, PID `h = 1`, `omega = y^5`, `zeta_5 = y^3`), `Phi_15 == h_pi h_pibar (mod p)` (`deg 4` each), `P = (p, h_pi)` has `O_H / P ~= F_(p^4)`. The unimodular `T_2`-isometry `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` turns `w = y` into **100% `sigma`-fixed** `g(X) = X^3 - omega^2 u` (deg 24) and `f(X) = X^9 + (omega + 1) X^6 + X^3 - omega` (deg 72).
- **Exact order-12 Frobenius (`audit_cm.py`: `algebra_bls12`):** `(u omega^2)^((p^4-1)/3) != 1` (so `g` is irreducible over `F_(p^4)` and `F_(p^12) ~= O_H[X]/(P, g)`), `sigma(y) = y^7` fixes `P, f, g` and equals `y^p mod h_pi`, and `tau(y, X) = (y^7, omega^2 X)` equals exact `Frob^1` (`chi = p == u mod r`, `ord_r(chi) = 12`, `Phi_12(chi) == 0 mod r`), quotienting relations and factor base into `Z[zeta_12]` 12-orbits.
- **Split-base localization & Pigeonhole Zero-Log Theorem:** `Psi_f, Psi_g` localize only at `P_f, P_g` above `P` (fixed by `tau`); conjugate `P_bar` ideals have norm `>= p^4 ~= 2^1524`. Every ramified/non-maximal prime `ell in {2, 3, 5, 19, 67, 906349, 254760293}` (and every `ell !== 1, 11, 14 mod 15`) has `tau`-orbit length `d < 12`, where `gcd(chi^d - 1, r) = 1` forces `vlog(Q) == 0 (mod r)` identically.
- **4 Schirokauer columns & Murphy `alpha`:** `<tau> ~= C_12` acts freely on the 36 complex places of `K_f` and 12 of `K_g` (`dim((O_H^x/r)^chi) = 0`), giving `3 + 1 = 4` Schirokauer columns. Exact `w' = omega^2` Murphy `alpha`: `alpha_f = +2.718`b (charged `+2.74`b), `alpha_g = +0.762`b (`alpha_(g,>=3) = +1.267`b, charged `+1.28`b; `kappa_2 = 0.0525`).
- **Live E2E `(4, 3)` `Q(zeta_15)` Tower DLP (`audit_cm.py`: `toy_e2e_dlp`):** On toy BLS12 `u_0 = -5, p_0 = 7207 == 7 mod 15, r_0 = 601` (`P_bar` present): `0` `P_bar` hits, `27` omitted short-orbit factors (`0` errors), `25/25` `K_g` deps annihilating live BSGS logs, neg controls failing (`20/26` drop `SM_g`; `40/40` wrong `chi`), `5/5` unseen target recoveries, and `57/57` `K_f` (deg-72, `3` `SM_f` cols) deps annihilating live BSGS logs (`62/62` fit).

## 3. Relation Collection: `Arch2-Rev-SigmaB-v3` & Exact de Bruijn Yield

Operating point **`B71.5_d04_e109`**: `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = 2^63.317373` 12-orbits), filtering excess `1.09`, collection target `1.148 U' = 2^63.516495` relations.

### Geometry, 10-Facet Voronoi FD, & `100%` `b_0`-Invariant `1/12` FD `D_(omega,sigma)`

In `B_(1,omega) = {1, omega, zeta_5, omega zeta_5, zeta_5^2, omega zeta_5^2, zeta_5^3, omega zeta_5^3}`, the `T_2` Gram matrix is `G = A_4 x A_2` (`det G = 1,265,625`, `W^T G W = S^T G S = G`, verified live).
- **Canonical `a` (10 explicit Voronoi units, verified LIVE in `audit_cm.py`):** `U_10 = {(y^3+1)^(+-1), (y-1)^(+-1), (y^2-1)^(+-1), ((y^3+1)(y-1))^(+-1), ((y^3+1)^(-1)(y^2-1))^(+-1)}` (`cm_data.json`: `voronoi_units_10`); all 10 have norm `+-1` and all `12/12` anchors satisfy `T_2(a) <= T_2(u a)` (`p_can = 0.268857`, `N_a^can = 2^57.8294`).
- **2-adic regions & `D_(omega,sigma)` (`239/239` exact-1):** Because `v_2(N_g)` depends solely on `a` (and is `sigma`-invariant), class radii `(1.00, 1.20, 1.36)` with `r_(ab) = 1.20`, `A^2 = 523230.1` yield unified key `key_a(a) = q_a(a) - s_c(a) <= -q_b(b)`. Moving `<sigma>` onto `b` via the `b_0/b_1`-invariant `1/12` FD `D_(omega,sigma)` (verified exhaustively on `Q <= 30`: `2,940` pts, `239/239` exact-1, `0` dups, `0` misses) gives `N_7 = 2^57.1414` rows (`L_max = 544.52`, mean `L_b = 233.88`) and `N_a^tbl = 4 N_a^can = 2^59.8294` slots.

### `Arch2-Rev-SigmaB-v3` (`W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

One 4-orbit `a`-table built at `s_0` over `a_slot = (m << 58) | oid` serves all 4 `sigma`-conjugates simultaneously (`m' = (m - k) mod 4`). Since `L_max = 544.52 < W_tile = 552`, every `b_0`-row fits in **one single full-chord tile** (`nt = 1.000`, `4.658` SBs/row across `552` column owners at `2^61.829` B `<< 2^64` B each):
- **Build (`C_LIST_A = 35.0`, `C_LBKT = 4.0`):** **`2^106.4594`** units, peak table storage **`2^119.65`** B (`< 2^120.0` B).
- **Band 1a (`[2^24, 2^48)`, `m`-partitioned) & Band 1b (`[2^48, 2^55)`, `rev_b1_early4_find`):** `m`-partitioning eliminates conjugate rotation (`C_HIT_LO_M4 = 6.25`); `rev_b1_early4_find` exits empty Band-1b buckets in `11` Ir (`C_HIT_LO_U = 8.0`, `<= 2^63.448` B/table `< 2^64` B).
- **Band 2 (`[2^55, P_hi)`):** `W_SB = 64` with per-SB tight `q_min_j` and **1-byte `sb_min_q8` gate** (`+1/254` quant slack, `9` B/SB, `<= 2^63.939` B/table `< 2^64` B) skipping **`85.4%`** of spanned SBs (`C_SB_SKIP = 6.0`, `C_SB_ENTER = 20.0`; verified live with `0` mismatches in `audit_cm.py`). Exact 16D ellipsoid integration gives `C/H = 1.4071` (charged **`1.422`**, `C_CAND = 13.0`), with 8B `log_p`-grouped Band-2 hit buffers (`C_HIT_HI = 14.5`).
- **Zero 9D walk & live 2D-histogram closed-form NTT Stage B:** Active-prefix `rev_scan_zero_u4` (`C_SCAN_ZERO = 2.95`), `_nk_g_only` (`4,500`), `_nk_f_only` (`6,500`), and closed-form 4-prime 62-bit NTT Stage B (`stage_b_analytic`, `KERNELS.md` §3, `(s_a_g, s_b_g) = (345, 295)`, `(s_a_f, s_b_f) = (405, 345)`) are integrated **LIVE across all `1,090` bins** up to `bg_max = 700`b inside `audit_cm.py` (`r_eff = 0.68001`, **`c_cof_2d_upper = 3.2887`** ipc).

### Pool Yield, Band/ECM Split + Chebyshev, & Collection Cost

Evaluated on the pre-registered `200,000`-sample duplicate-free pool `P2_wy_2ad` (seed `900003`; 12 anchors verified live to `< 5e-11`b) with `rho_dB` and per-row `Retention` integrated jointly:
- `V = 2^119.5919` (`SE_c = 0.00290`, `SE_w = 0.00303`, `H2_ALLOW = 1.00774 >= 1 + 2 max(SE_c, SE_w)`), mean norms `(550.13, 678.22)`b (relation-weighted `(544.38, 673.75)`b).
- Retained yield `2^63.516495` (`SE_y = 0.01020`b); strict upper-endpoint retention `2^63.513963` minus `0.009`b gives net margin **`0.06326`b (`>= 0.060`b)** above `1.09 U'`, and its 2-SE lowered count exceeds descent's `2^63.46`.
- **Per-cell cost:** `g` `12.98` + `f` `9.65` (Band 1 `12.79`, gate+probe `2.29`, Band 2 `7.55`) + scan/zero `2.95` + 2D-upper cofactorization `3.29` = **`28.87` ipc** -> **`2^110.166694`** core units.
- **Band split & ECM split:** Across `<= 2^64.066` kept relations: band split `<= 2^100.543` units; Suyama-Montgomery ECM + limb-scaled `mr_instr` gives `E[W_ecm] = 2^95.6179` units, metered at **`2^96.5`** units with Chebyshev tail `<= 2^-77.14 < 2^(-70)`.

## 4. `Z[zeta_12]` Filtering & 1-Owner Block Wiedemann Linear Algebra

Splitting `Z[zeta_12] = Z[omega] (+) chi Z[omega]` into 6 signed hexagonal buckets cuts unit+unit collision fill-in `3x` vs `Z[i]`. Bundled `filtsim_z12.py` (hosted experiment) mirrors our compiled C `filtsim_z12` campaign across 6 seeds at every doubling `U' = 2^16 .. 2^23` (`B71.5_d04_e109`, `W = 200`, conservative cell-mean norms `548.9 / 677.8`):

| Scale `U'` | Seeds | `N_l / U'` (SE) | `N_sw / U'` (SE) | `sd_l U'^(1/4)` | Depth | Unit % (Exp.) | `k5` | `tab_ratio` |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `2^16` | 6 | `0.003260` (`0.000227`) | `0.006462` (`0.000482`) | `0.0089` | `13.7` | `95.61%` (`89.15%`) | `10.5342` | `0.80659` |
| `2^17` | 6 | `0.003050` (`0.000112`) | `0.006089` (`0.000219`) | `0.0052` | `15.2` | `96.63%` (`90.47%`) | `10.2017` | `0.80554` |
| `2^18` | 6 | `0.002857` (`0.000088`) | `0.005747` (`0.000173`) | `0.0049` | `16.7` | `96.96%` (`91.30%`) | `9.9914` | `0.80212` |
| `2^19` | 6 | `0.002533` (`0.000074`) | `0.005111` (`0.000118`) | `0.0049` | `16.7` | `97.18%` (`91.82%`) | `9.8616` | `0.79984` |
| `2^20` | 6 | `0.002399` (`0.000032`) | `0.004857` (`0.000049`) | `0.0025` | `18.7` | `97.35%` (`92.32%`) | `9.7354` | `0.79839` |
| `2^21` | 6 | `0.002172` (`0.000035`) | `0.004377` (`0.000067`) | `0.0033` | `18.5` | `97.42%` (`92.61%`) | `9.6628` | `0.79684` |
| `2^22` | 6 | `0.001999` (`0.000021`) | `0.004057` (`0.000039`) | `0.0024` | `19.7` | `97.43%` (`92.77%`) | `9.6218` | `0.79508` |
| `2^23` | 6 | `0.001889` (`0.000013`) | `0.003846` (`0.000026`) | `0.0017` | `20.5` | `97.45%` (`93.01%`) | `9.5620` | `0.79485` |

- **Fixed-instance `H3` structural premise & deterministic hard caps (`P(cap abort | H3) = 0`):** Collection's region, polynomials, and ECM curves are fixed by a public seed, so the finite relation set and filtered matrix are **deterministic for the fixed production instance**. Both `N_l/U'` and `N_sw/U'` fall at every doubling (`N(250)/U' = 0.003552` at `2^23`), explicit unit share rises monotonically to `93.01%` (`k5_exact = 9.5620`), `tab_ratio` slope `+ 2 SE = -0.00160 < 0`, and lowering `w_mid` from `1.0 -> 0.5 -> 0.0` in `filtsim_z12.py` (`log2U = 14..15`, 4 seeds) lowers `N_l/U'` (`-5.5%` at `2^14`) and `N_sw/U'` (`-20.3%` at `2^14`, `-17.1%` at `2^15`), putting production (`w_mid = 0`) at or below the finite-scale trend. The algorithm enforces hard abort caps on all four parameters: **`N_l <= 0.00225 U' = 2^54.5215`**, **`N_sw <= 0.00450 U' = 2^55.5215`**, **`k5 <= 10.50`**, and **`tab_r <= 0.85`** (any run above a cap aborts). **`H3`** is declared explicitly as a fixed-instance structural premise asserting that the fixed production matrix satisfies all four caps (`P(cap abort | H3) = 0`; no synthetic filtering-seed probability is substituted for challenge/coin randomness). Note: `filtsim_z12` samples column incidences independently at the production occurrence profile, so the scaling ladder is exploratory support for `H3`, not a proof for jointly smooth retained relations.
- **1-Owner Arch B BW & Closed-Form NTT Generator:** Input `x` + output `y` + all 22 dense tables fit in **`2^63.913` B `< 2^64` B** on ONE owner (`0` cross-owner sends/combines). At `(m, n) = (8192, 110)`, `C_ROW + C_CHI = 785`, `N_SM = 4`, `BW_SIG = 2^(-7)`: **LA `<= 2^107.174671`**, **Generator `<= 2^98.584`** (using closed-form 1D/2D 4-prime 62-bit NTT `fft_F`, `KERNELS.md` §3).
- **4-Column Mod-`r` CRT Schirokauer Maps (`<= 2^84`, verified live):** `Phi_15 == h_(r,1) h_(r,2) mod r` (`[4, 4]`); `K_f` stays degree-9 irreducible over both components (`[36, 36]`, `eps_f = r^36 - 1`); `K_g` stays irreducible on `w_1 = -u^2` and splits into 3 linear factors on `w_2 = u^4` (`u w_2^2 = (u^3)^3`), giving **`[4, 4, 4, 12]`** (`eps_g = r^12 - 1`).
- **Finite Relation Rank (`H-R`), Block Wiedemann (`H-X`), & 1D-Kernel Horner Diagnostic (`RANK-BRIDGE.md`, `hx_z12.py`):** `RANK-BRIDGE.md` separates fixed-instance `H-R` (`rank(M) = N`, `rank(A) = D - 1`, 1D virtual-log kernel `(w, c)` with `c != 0`), finite `H-X` (nonzero `N x N` control and `(D-1) x (D-1)` observation minors), and the proved conditional Schwartz-Zippel bound `<= 4 N_cap / 2^128 < 2^(-66)` given `H-R` and `H-X`. `hx_z12.py` runs `filtsim_z12` live unscreened (ALL `N` live cols, first `D = N + 4` live rows), plants an exact 1D kernel `(w_true, c_true)` (testing `c_0 = 0` on even seeds), evaluates the actual 2-step recorded-pivot product operator `op(x) = (M_live - E_live L_rec^(-1) B_rec) x`, solves the projected Hankel kernel `H f = 0`, and verifies exact Horner recovery `A(w_rec, c_rec) = 0` proportional to `(w_true, c_true)` plus 1-block pad negative control failure (`6/6` on `log2U in 12..14`). This tests the product operator and 1D-kernel BW recovery on a constrained synthetic system; `H-R` and production `H-X` remain declared heuristic premises.

## 5. Individual Logarithms

Detailed in `DESCENT.md`, `THRESHOLD-REPAIR.md`, and `ADAPTIVE.md` (all `39` bins / `78` nodes computed **100% LIVE** from `plan_struct` by `cm_models.eval_descent_live`). Using constructive 16D integer LLL -> rank-2 submodule `(v_1 = b_1, v_2 = b_k)` (`a - b X`, exact relative norm `a^3 - omega^2 u b^3`, `RET_BOOT = 0.50` decile filter at `m_max = 2600`b, `BOOT_RELS = 48`), half-scale NB-PGF **`T_k + psi^k`** union bound (`r = 2, alpha = 0.5, DISC = 0`, NO `PGF(psi)` substitution; `THRESHOLD-REPAIR.md`), **100% full survivor work** (`1.0x`), **full worst-case multi-factor ECM extraction loops** (`0.999^K` detection), and **limb-scaled `mr_instr(nb)`** (`5.14e9` Ir at `2600`b):
- **Live lattice, counterexample, & reduced-tower verification (`audit_cm.py`, `ADAPTIVE.md`):** Constructive 16D LLL on `p_0 = 7207` (`k = 2, index = 1`); 4 bundled 381-bit BLS12-381 `w'=omega^2` boot bases (`det = p^4`, `68/68` norms `< 2600`b, `4/4` full 10-decile passes); live `M=1/M=100` counterexample refuting `PGF(psi)` (`0.1 > 1/28`, repaired `0.132258 >= 0.1`); and a 1,600-root uncensored reduced-`Q(zeta_15)` tower campaign (`ADAPTIVE.md`, `4,292` nodes, `731` after-fail attempts including fixed 5%-missing shared-DB runs: all `PGF_2 <= 0.950 < 1.0`, every row-specific `fail O/E` 95% CI bracketing `1.0`, every `work O/E < 0.85`; 4-root slice replayed live in `audit_cm.py`).
- **Expected two-log work:** `2^107.448555` <= **`2^107.449`** core units (`2^107.383549` at `{"descent_mesh": 384}`), metered at `t = 2^119.5` (`P(cap) <= 2.356609e-04`), `phi_boot = 1.34e-21 <= 1.0e-12`, `E[N] = 2^14.4757`.

## 6. Success Probability & Strict Additive Union Bound

Strict additive union bound `P(success) >= 1 - sum_i P(fail_i)` (verified in 80-digit `Decimal` inside `audit_cm.py`):
- **Main phase & ECM meter:** BW Schwartz-Zippel `<= 4 N_cap / 2^128 < 2^(-66)` (proved conditional on fixed-instance `H-R`, `H-X`; `RANK-BRIDGE.md`), fixed-instance `H3` caps (`P(cap abort | H3) = 0`), ECM Chebyshev `< 2^(-77.1) < 2^(-70)`.
- **Descent & target (2 logs, `H5`):** Two boots `<= 2 x 1.0e-12`, metered cap `<= 2.356609e-04` (`2.592270e-04` at `+10%`), node counter abort (`N_MAX = 2^25.8`) `<= 2^-11.3242`, `tau`-orbit repeats `<= 2^-10.05`, degenerate target `<= 2/r`.
- **Wrapper:** `< 2^(-50) < 2^(-48)` (`field-map.md`).
- **Opened-run floor & gate:** Core `>= 0.998431`; with `+10%` cap headroom & wrapper, `>= 0.998407 >= 0.998407`. Gate `K = 7205708882997339893` satisfies `(K/2^64) x 0.998407 >= 0.39 > ((K-1)/2^64) x 0.998407`.

## 7. Complete Ledger

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (`Arch2-Rev-SigmaB-v3`, sec. 3, `H-EA`; est. `110.155570` + `0.774%` H2 allow.) | <= 110.166694 |
| replay 4-orbit `a`-tables build (`g` to `2^63` + `f` to `2^60`, 3-pass `a_0`-chord; sec. 3) | <= 106.4594 |
| linear algebra (`1-owner Arch B Z[zeta_12] BW`, `W=200`, hard caps `k5<=10.5, tr<=0.85`, `(8192,110)`; sec. 4) | <= 107.174671 |
| matrix Berlekamp-Massey generator (closed-form 4-prime 62-bit NTT; sec. 4) | <= 98.584 |
| both individual logarithms (metered cap; expected `107.449`; sec. 5) | <= 119.5 |
| band split of kept relations (`(2^20, 2^48]`, 10/12 limbs; sec. 3) | <= 100.543 |
| ECM split of kept relations (`(2^48, B]`; metered cap; exp. `95.62`, Cheb. `< 2^-77`; sec. 3) | <= 96.5 |
| canonical `a`-list enumeration + global `key_a` pre-sort | <= 86 |
| 4-column `chi`-Schirokauer maps (mod-`r` CRT `f: [36,36]`, `g: [4,4,4,12]`; meas. `<= 82.8`) | <= 84 |
| filtering streams, product-form pivot selection, back-substitution, misc | <= 80 |
| `P_g`, `P_f` product trees (`< 2^89` instructions) | <= 75.5 |
| factor base, roots, 12-orbit ideals (`U' = 2^63.32`) | <= 75 |
| numerical certificate construction | <= 50 |
| final verification | <= 11 |
| wrapper (sec. 1; `<= 2^80` instructions = `2^65.71` core units) | <= 66 |

Core cap sum (`80`-digit `Decimal`): `2^119.502690065885 <= 2^119.5026900659` core units (`< 2^119.502690065901` P/Q = `< 2^122.824618160788` G1). Replacing `2^119.5` by expected descent `2^107.449` gives opened-run expected core work `2^110.605986584879 <= 2^110.6059865849` core units, and behind the run gate `K / 2^64 = 0.3906222613`, `< 2^109.249832659662` dense-`F_(p^12)` units (`C = 20,000`) and `< 2^112.571760754550` G1 units (`C = 2,000`).

## 8. Heuristics

- **`H1` (Smoothness & Local Laws):** Two-sided smoothness with explicit `rho_dB(u, x)` on Collection and Descent, Murphy `alpha = (+2.74, +1.28)`, exact `v_2(N_g)` and `kappa_2`, boot pair ratio `kappa <= 2^24` (`[11]`: `0.7076`).
- **`H1'` & `H-EA` (Retention & Survivor Work):** Two-stage Euler-product retention and survivor-upper laws (`cm_models.Retention` + closed-form 4-prime NTT `stage_b_analytic`, `0.009`b band-miss margin; `f/g` dependence `q_f = min(1, 10 rho(u_f))`).
- **`H2` (Pool Estimation):** Pre-registered `200,000`-sample duplicate-free pool `P2_wy_2ad` (`H2_ALLOW = 1.00774 >= 1 + 2 max(SE_c, SE_w)`); exact 16D ellipsoid concentration (`+1%` discrete safety on `C/H`).
- **`H3` & `H-SM` (Fixed-Instance Matrix Caps & Class Numbers):** `r` unramified (proved), `r !| h(K_f), h(K_g)`, `P_chi Lambda` injective on the 4-dim `chi`-unit space (`toy_e2e_dlp`); **`H3`** is a fixed-instance structural premise asserting that the fixed production filtered matrix satisfies all four hard caps `(N_l/U' <= 0.00225, N_sw/U' <= 0.00450, k5 <= 10.50, tab_r <= 0.85)` (`P(cap abort | H3) = 0`), supported by the monotone `2^16..2^23` independent-incidence scaling ladder and `w_mid -> 0` trend.
- **`H-R` & `H-X` (Finite Relation Rank & Block Wiedemann, `RANK-BRIDGE.md`):** Fixed-instance `H-R` (`rank(M) = N`, `rank(A) = D - 1`, 1D virtual-log kernel with `c != 0`) and finite `H-X` (nonzero `N x N` control and `(D-1) x (D-1)` observation minors on sparse block-diagonal `X` placing the 4 pad coords in distinct blocks). Neither follows from smooth counts, class numbers, or synthetic replays; conditional on `H-R` and `H-X`, Schwartz-Zippel gives failure `<= 4 N_cap / 2^128 < 2^(-66)` (`hx_z12.py` tests 1D-kernel recorded-pivot Hankel+Horner recovery live).
- **`H4` & `H4u` (Kernels & Unit Share):** Baseline scalar AArch64 counts (`KERNELS.md`, zero SIMD, zero `FEAT_LSE`); `Z[zeta_12]` 6-bucket unit share (`93.01%` explicit at `2^23`, `k5 <= 10.50`, `tab_r <= 0.85`).
- **`H5` (Descent Model, `DESCENT.md`, `THRESHOLD-REPAIR.md`, `ADAPTIVE.md`):** Lower-sum semismoothness with `rho_dB`, deterministic 50% boot-ball decile pass rate (`RET_BOOT = 0.50` at `m_max = 2600`b), half-scale NB-PGF count envelope `G(t) = (1 + 0.25 mu (1 - t))^(-2)` (`r = 2, alpha = 0.5, DISC = 0`) paired with separate disclosed conditional child-failure hazard `P(F_i | A_i, prehistory) <= psi` (`phi_node = max(0.02, min_k(T_k + psi^k + 1e-10))`, `phi_boot = psi^48`) and conditional child-work envelope, full `1.0x` survivor work, independent subtrees of distinct ideals (`[7e]`), uniform orbit collision rate, Poisson missing-ideal model, random-integer ECM orders, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`.

## 9. Memory & Provenance

**Peak memory:** `< 2^120.0` bytes (`2^119.65` B replay `a`-tables; `< 2^105` B descent; `2^63.913` B `< 2^64` B on 1 owner in LA; `< 2^30` B wrapper; `STORAGE.md`).

**Provenance:** CM-split cyclotomic towers over `Q(zeta_15)` build on Kim-Barbulescu (2016), Guillevic-Singh (2020), Zhu-Liu (2025), and APT26; filtering follows CADO-NFS and Boudot (2017/758). Prior CurveSmash lineage (`c24..c36`, `stnfs-bn15r`; gbotrel, Claude Opus 5.5 / Fable 5.2) and credited community ideas: run gate, expected work, high descent cap (kaankolcu `c613ae01`, `3359b9e2`), Shoup `chi`-mul (soispoke `cf131637`, and `P(M<k)+psi^k` union bound `cfdcfdc4`), carry-chain combine (yukonautominer `9a7df12c`), and `(m, n)` sizing (mpjunior92 `ea928ab4`). No prior verdict is a premise.
