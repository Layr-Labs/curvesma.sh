# bn254-pairing-inversion: Pure-Integer `C6-m184_w216-rec` `(6,2)` STNFS — Expected Work `2^86.620236211` Units

Track `bn254-pairing-inversion`. Input: `z` uniform in `GT` (`g = e(P0, Q0)`); output `(P, Q)` in `G1 x G2` with `e(P, Q) = z` (`P = P0, Q = [x]Q0`, G2 scalar mult `output <= 2^11.761903` u; if `z = 1`, `(P0, inf)`). Score unit: `C = 20,000` pure baseline scalar ARMv8.0-A integer instructions (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, `FEAT_LSE`, `DC ZVA`, `fp64`, external C; `p_num = 0`). 100% self-contained across `proof.md`, `RANK-BRIDGE.md`, `STORAGE.md`, `KERNELS.md`, `ADAPTIVE.md`, and live `audit_bn_cm.py` (zero imported prior verdicts). No advice, smaller targets, or amortization. Analytic, exploratory claim.

## Claim

Under declared exploratory STNFS premises `H1`, `H2`, `H3`, `H5`, `H-EA`, `H-SM`, `H-R`, `H-X`, `H-SEP` (§6, `RANK-BRIDGE.md`, `ADAPTIVE.md`; zero `H-FFT`, Band-3 `H1'`, Poisson Band-2 occupancy, `THRESHOLD-REPAIR` `psi^k`), using exact 3-class `rho_dB` on pure-product `s = 0.80, kap = 0.0` (`200k` exact-HNF cells), pure-int Band-1/2 direct-apply sieve (`W_SB = 128`, 1-bit `occ_bm[]`, deterministic `p_occ = min(1.0, mu_sb)`, 64-cell `ORR <= 0.96`, zero-sort 2-pass lists, 48-op Radix-4 + 46-op (`<= 48`) Radix-3 `< 2^62` Shoup NTT, deterministic zero-GCD `P_14` strip `C_SB = 56,000`, `7`-limb `f`-SRT), fixed-instance `H3` 4-vec `Z[omega]` LA (`2^12..2^24` all-4-monotonic ladder, caps `0.00120/0.00078/13.80/1.00`, live `NU_H = 0.9837`, ZERO TABM), 5-SM shifted Beckermann–Labahn $\sigma$-basis extraction, analytically input-uniform 2-stage Bareiss 12D LLL (`<= 2^39.81 < 2^40` Boot, `<= 2^32.89 << 2^36` Special-`q <= 2^220`, `<= 2^25.98 < 2^26.2 = PAIR_SETUP` for ALL `p <= 2^59`), and `62`-row `LB = 59` **pure `k = 1` zero-retry Boole + Wald `i.i.d.` fresh-`e` boot** (`p_0 = max(0.001, (1+0.5*mu)^{-1})`, `Phi_tree = 0.515420 <= 0.60`, `E[runs/log] <= 2.063642`, `p_node_abort = p_repeat_loss = p_descent_fail = 0`, 6-bit ECM ladder `266` DAG / `15` Boot rungs, unified meter `t_run = 2^112`), the attack succeeds w.p. `>= 0.39` (opened-run `>= 0.999996` via `FIX-U`), and

    E[W] <= 1 + (K / 2^64) * 2^87.978684411 <= 2^86.620236211 units (C = 20,000 Ir),
    K = 7194258965782588261.

Open iff fresh 64-bit coin `u < K` (least integer with `(K/2^64)*0.999996 >= 0.39`; closed run fails at `1` unit). `86.620236211` bounds strictly **expected** work `E[W]`. Unified hard aborting meter `t_run = 2^112` encloses the opened run pointwise (`Markov = 5.873e-8`; split `2^112` `8.306e-10`). Work rounds up; `success` rounds down; finite prototypes are scoped evidence, not proofs.

## 1. Tower, Reciprocal Pair, Exact `sigma(w) == w`, & Mod-`r` CRT

BN254: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.
- **`K_h` & `w`:** Totally imaginary sextic over `Q(sqrt(-2))` (unit rank `2`, order-3 `sigma`): `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -2^9 * 19^2`, `h_K = 1`, `w_K = 2`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`, `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3`, `w(y) = -y^5 + y^4 - y^2 - 2y` (`N(w) = 216`), `(u*w)^((p^6-1)/2) = -1` in `F_{p^6}`, **`w(sigma(y)) == w(y)`** (live).
- **Reciprocal pair & `tau`:** `f_rec(x) = (-2y^5+2y^4-2y^2-4y+3)x^8 + (3y^5-3y^4+3y^2+6y)x^6 - 12x^4 + (-3y^5+3y^4-3y^2-6y+12)x^2 + (2y^5-2y^4+2y^2+4y-5)`, `g_rec(x) = u*w*x^2 - 1` (`N(lead)=N(const)=729`, `f_rec(x^2->(u*w)^-1) = 0 mod (p,h)`), `tau(y,x) = (sigma(y), -x)` (`Frob^10`, `chi = p^10 mod r`), mod-`r` DDF `h:[1^6], g_rec:[2], f_rec:[2^4]`. Pure-int `G_sym`/`D_fi` Voronoi FD certified for ALL radii via Dedekind, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), and `Decimal(55)` index-1 & hexagon theorem (`KERNELS.md` §0).
- **Exact 3-Class 2-Adic Split & Map:** `2 = P_2^2` (`N=8, v_{P_2}(w)=1`) partitions primitive `(a,b)` into `(v_2(N_f), v_2(N_g)) = (0,0)` (`7/9`), `(0,3)` (`1/9`), `(6,0)` (`1/9`); `alpha_f_odd = ALPHA_F_RAW + 2/3`, `alpha_g_odd = ALPHA_G_RAW + 1/3`. `phi_iso: F_{p^{12}} -> F_p[y,x]/(h, uwx^2-1)` (`wrapper <= 2^53.712288`); if `z = 1` output `(P0, inf)`; else uniform `e in [0,r-1]` -> `T = phi_iso(z*g^e) = c_0 + c_1 x`, `R = c_1/c_0 mod (p,h)`, `x = (vlog(T)-e*vlog(g))*vlog(g)^-1 mod r`, output `P = P0, Q = [x]Q0` (`output <= 2^11.761903` u), checked by `e(P, Q) == z` (`verify <= 2^11.000001` u).

## 2. Pure-Product (`kap = 0.0`) Deterministic-Occ Band-1/2 Direct-Apply Collection (`12.668` Ir/Cell)

- **200k exact-HNF pool (`s = 0.80, kap = 0.0`):** `B_f = 2^61`, `B_g = 2^60` (`lx2 = P_mid = 2^54.5`, `sl = 160, sb = 165`), `U' = 2^53.641526`, target `0.82*U'` (`+0.0666`b). Exact 3-class `rho_dB` + class-matched `Retention` across `136,563` prim + `63,437` nonprim == `200,000` cells (live `1,000`-sample HNF `|z| <= 1.15`): `V = 2^97.8362` (`SE = 0.0046`b), `V_b = 2^52.3194`, `N_a = 2^45.5168`, `N_5 = 2^43.6888`, `L_b = 396.33`, `L_max = 807.40`, `N_f = 402.69`b (**AM-GM `sup log2|N_f| <= 441.66 < 448`b = 7 limbs**), `N_g = 471.47`b, `r_g = 0.5966`.
- **Single-Arena Band-1/2 Sieve (`P_min=14, P_cut=45.0, P_mid=54.5, W_SB=128`, ZERO Band 3; `KERNELS.md` §1):** Band 1 (`rep_lo = 5.024`): 1-owner arena `2^{56.18} B`, `APPLY8 = 3.75`, `C_B1_ROW = 24.0`. Band 2 (`probe=0.417, occ_sb=1.208, hi_ap=3.124`): 1-bit `occ_bm[]` (`12/14/7/7`), 64-cell `ORR` `59/64 <= 0.96` (`0/128k` errs). **Deterministic `p_occ = min(1.0, mu_sb)` (ZERO Poisson):** `#(non-empty SBs) <= min(S_p, N_a)` deterministically (`det_occ_bound_ok = True`), giving entered-SB rate `<= min(1.0, mu_sb)` (mean `0.6067`) under `H2`.
- **Zero-Sort 2-Pass Lists (`lists <= 2^84.688287`):** `C_LIST_B1 = 20, C_LIST_B2 = 24, C_LBUCKET = 6` (`2^52.32 B`, `zero_sort_exact_ok = True`).
- **Stage B (`cof = 1.935` Ir/cell):** `C_NORMS = 30k`, deterministic zero-GCD `417`-chunk `P_14` strip (`55,949 <= 56,000`, `0/20` errs incl. prime cubes), `4` disjoint LP windows (`0` fn), 48-op Radix-4 + 46-op (`<= 48`) Radix-3 `< 2^62` lazy `[0, 2p)` Shoup NTT (`60*n*k + 220*n` / `60*n_3*(k-2) + 364*n_3`, live `n=16` & `n=24` convolutions), `7`-limb `f`-SRT (`Q_F_SAFETY = 2.0`), nonprim `0.0396`. Total: **`12.668` Ir/cell** `=> collection <= 2^87.216306`, global RAM **`2^98.394 B < 2^100 B`**.

## 3. Fixed-Instance `H3` `Z[omega]` Filtering (`2^12..2^24` Ladder + `NU_H`), 4-Vec LA, & Shifted $\sigma$-Basis

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `k5_exp` (`8.00/58`) | `tab_r` | `w_mid` |
|---|---|---:|---:|---:|---:|---:|
| `2^12` | **LIVE Python** (`W<=161, kern=0, lost=0`) | `0.014160` | `0.009277` | `20.23` | `0.7760` | `14.07` |
| `2^13` | **LIVE Python** (`W<=161, lost=0`) | `0.006836` | `0.004028` | `19.51` | `0.7806` | `13.71` |
| `2^20` | C `filtsim` (`n=4`) | `0.001271` | `0.000751` | `13.85` | `0.7482` | `11.54` |
| `2^21` | C `filtsim` (`n=8`) | `0.001031` | `0.000630` | `13.62` | `0.7436` | `11.25` |
| `2^22` | C `filtsim` (`n=4`) | `0.000923` | `0.000570` | `13.35` | `0.7411` | `10.96` |
| `2^23` | C `filtsim` (`n=4`, max `0.000808/0.000510/13.20`) | `0.000761` | `0.000489` | `13.15` | `0.7390` | `10.67` |
| `2^24` | C `filtsim` (`n=4`, max `0.000695/0.000452/12.92`) | `0.000673` | `0.000441` | `12.90` | `0.7363` | `10.38` |
| **Fixed-Inst. `H3`** | **Deterministic Abort (`W<=161, P(abort|H3)=0`)** | **`0.001200`** | **`0.000780`** | **`13.80`** | **`1.0000`** | **`~0.0`** |

- **Fixed-Instance `H3` Support (`verify_h3_scaling_and_nu_h` LIVE):** Collection and filtering are deterministic on the fixed BN254 instance, so the four core dimensions are fixed constants declared under exploratory `H3` (`P(cap exceeded | H3) = 0`, ZERO synthetic Chebyshev tail claim). Diagnostic support: (1) all 4 metrics decrease strictly across `2^20..2^24` (`2^13..2^24`), caps sit **`1.782x`** (`N_sw/U'`) and **`1.768x`** (`N_l/U'`) above `2^24` means; (2) empirical inverse-log fit `39.597/(log2 U')^{3.458}` (`max_rel_err = 0.028`) projects `4.10e-5` at `2^{53.64}` (`29.3x` below cap); (3) Chebotarev `C6` `<tau>`-orbit splitting up to `2^{14}` (`419` `d=3` + `102` `d=6` = `623` orbits -> **`NU_H = 0.9837`**, `0.9992` by `2^{16}`) + live `185`-rel genuine `C6` filter (`0` cols lost, `98.1%` units).
- **4-Vec `Z[omega]` LA & Shifted $\sigma$-Basis (`RANK-BRIDGE.md`, `hx_z6.py`):** `unit4 = 8.00`, `gen1 = 56 <= 58`, `tail = 722 <= 740` (ZERO TABM), `bw = (1024, 35)`, `nslots = 7`: **`la <= 2^84.626125`**, **`gen <= 2^77.417545`**, 1-owner RAM **`2^59.649 B << 2^64 B`** (`lazy = 0`). `S_full` (`sm <= 2^73.827705`), Schwartz-Zippel `3.656e-6`; live unplanted `F_1187^*` `18x18` shifted $\sigma$-basis (`4` SM-nonzero, `4/4` `Hf=0`, `[407]` across `68/68` rels) + Schur.

## 4. Pure `k = 1` Zero-Retry Boole + Wald `i.i.d.` Boot (`LB = 59`) + 6-Bit ECM Ladder + `ADAPTIVE.md`

- **Pure `k = 1` Boole + Wald `i.i.d.` Fresh-`e` Boot (`ADAPTIVE.md` §1):** All `62` rows (`plan59`) run pure `k = 1` zero-retry descent (`M = 0` aborts root attempt; ZERO internal retries, `THRESHOLD-REPAIR`, or `B_max`/`T_boot_max`, `p_node_abort = p_repeat_loss = p_descent_fail = 0`). Each slot charges `p_0(b) = max(0.001, (1+0.5*mu(b))^{-1})` (`mu_min = 713.91, p0_max = 2.794e-3`); unconditional Boole gives **`Phi_tree = 0.515420 <= 0.60`** over **`c_det_one = 411`** slots (`boot = [163, 88, 12, 4]`, `ea = [5, 1250]`, `x_EA = 64.4167`, `2*P_x = 2^62.945 B`). Fresh uniform `e in [0, r-1]` + input-uniform 2-stage Bareiss 12D LLL (`<= 2^{39.81} < 2^{40}` for EVERY `R`, `<= 2^{32.89} << 2^{36}` for EVERY `q <= 2^{220}`, `<= 2^{25.98} < 2^{26.2} = PAIR_SETUP` for EVERY `p <= 2^{59}`; `KERNELS.md` §§0,4) gives **`E[runs/log] <= 1/(1-Phi_tree) = 2.063642`**, **`E[trials/log] <= 1.1552e17`**, **`E[two_logs] <= 2^85.623497`** (`two_dag <= 2^84.6544`, `two_boot_exp <= 2^84.5919`), enclosed in unified meter `t_run = 2^112` (`5.873e-8`).
- **Four-Pillar Live Support (`26/26`):** (1) Real-2-core-FB `C6` `O_h` `ADAPTIVE.md` (`u_t=2, p_t=13`, `B_FB=2^{14}`, `9,870` roots, `1,856` after-fail: **`6/6` `O/E fail <= 0.871` (pooled `0.7935`)**, `Work O/E = 0.9985`, `10/10` bins `C_bday <= 2.620`, **3 live witnesses**); (2) Above-cutoff (`B_FB=1,024`) deg-8 `f_rec` mini-BN254 (`100/100` non-empty, Fano `0.807/1.106`, + 2 depth-2 chains); (3) `24/24` prod yield floor (`min 2.43x` cons / `7.4x` actual Murphy `alpha`); (4) `224`-ideal depth-2 PGF (`0/224` empty).
- **6-Bit ECM Ladder (`KERNELS.md` §4):** Upper-endpoint `[lo, hi_end]` enclosure across **266 DAG** (`0.6174`) & **15 Boot** (`0.6081 > 0.50`) rungs (`63/63` Suyama-12 live).

## 5. Itemized 16-Row Ledger, `FIX-U`, & Exact Rational Gate

| Row | Comp. (`log2` u) | Decl. (`log2` u) | Row | Comp. (`log2` u) | Decl. (`log2` u) |
|---|---:|---:|---|---:|---:|
| `collection` | `87.216306` | `87.216306` | `fb` | `65.038073` | `65.038073` |
| `chi_setup` | `65.712288` | `65.712288` | `split` | `81.834890` | `81.834891` |
| `la` | `84.626124` | `84.626125` | `gtree` | `59.185169` | `59.185169` |
| `gen` | `77.417544` | `77.417545` | `cert` | `50.000001` | `50.000001` |
| `lists` | `84.688286` | `84.688287` | `verify` | `11.000001` | `11.000001` |
| `two_logs` | `85.623496` | `85.623497` | `wrapper` | `53.712288` | `53.712288` |
| `canon` | `74.029127` | `74.029127` | `output` | `11.761902` | `11.761903` |
| `sm` | `73.827705` | `73.827705` | `misc` | `69.544407` | `69.544407` |
| **`exp_run`** | **`87.978684410`** | **`87.978684411`** | **`E[W]` (claim)** | **`86.620236211`** | **`86.620236211`** |

**`FIX-U`:** Schwartz-Zippel `3.656e-6` + Fixed-instance `H3` given `H3` `0.0` + Split `2^112` `8.306e-10` + Descent Wald/node/repeat `0.0` + Unified meter `t_run = 2^112` `5.873e-8` + MR `1.735e-18` + Group/boot/SM `1.056e-59` + Numerical `0.0` `=>` **`success = 0.99999628 >= 0.999996`**.

## 6. Declared Exploratory Heuristics

- **H1 / H2:** Full `rho_dB` with odd-part Murphy `alpha`, 3-class 2-adic split, all-`200k`-cell Euler retention/survivors (`lx2=2^54.5`, zero Band 3); uniform Voronoi FD sampling (`SE=0.0046`b) & `b`-row residue equidistribution mod `p` (paired with deterministic `p_occ <= min(1.0, mu_sb)`).
- **H3 (Fixed-Instance Deterministic Filtering):** On fixed BN254 `C6-m184_w216-rec` at `U'=2^53.64` (`exc=0.82, W<=161.0`), deterministic purge+merge satisfies the 4 fixed-instance caps (`0.00120/0.00078/13.80/1.00`, ZERO TABM; `P(cap exceeded | H3)=0`), supported diagnostically by the all-4-monotonic `2^12..2^24` ladder (`1.78x` headroom), inverse-log fit (`29.3x` below cap), live `NU_H=0.9837`, & live `185`-rel `C6` filter; LA/gen/RAM priced at cap (`lazy=0`).
- **H-SM / H-R / H-X (`RANK-BRIDGE.md`):** Non-degenerate 5-SM $\chi$-projections; `rank M_prod = N_l`, 1D kernel of `[M_prod|S_prod]`; controllability/observability of padded `(B,Y,X)` yielding Schwartz-Zippel `3.656e-6` & deterministic shifted $\sigma$-basis extraction.
- **H5 / H-EA / H-SEP:** `H-EA`: `QBINS`/`BOOT_BINS` lower-bound smoothness and `QBINS_LO`/`BOOT_LO_OFF` upper-bound survivors; `H5`: fresh-`e` `i.i.d.` boot trial succeeds w.p. `>= prob_b` and at every reached pure-`k=1` DAG slot `P(M=0|slot) <= max(0.001,(1+0.5*mu)^{-1})`, propagated via unconditional Boole (`Phi_tree=0.515420, c_det_one=411`) & Wald `i.i.d.` boot retries (`E[runs/log]<=2.063642`, meter `t_run=2^112`); `H-SEP`: 6-bit ECM separation `>= 0.50` (`min_dag=0.6174, min_boot=0.6081, 63/63` Suyama-12).
