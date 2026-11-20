# bn254-gt-dlp: 100% Pure-Integer `C6-m184_w216-rec` `(6,2)` STNFS (Band-0/1/2 Flat-Array Sieve, Shoup `< 2^62` NTT, `7`-Limb `f`-SRT, `2^23` `Z[omega]` LA, `H-PGF(0.5) + THRESHOLD-REPAIR` Descent) — Expected Work `2^86.992667169` Units

Track `bn254-gt-dlp`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `x`. Score unit: `C = 20,000` pure baseline scalar ARMv8.0-A integer instructions (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, `p_num = 0`). Every component is derived analytically here and in `RANK-BRIDGE.md`, `STORAGE.md`, `KERNELS.md`, and verified live in `audit_bn_cm.py` (self-contained; no prior verdict or score is a premise). No advice, smaller targets, or amortization are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared exploratory STNFS premises `H1`, `H2`, `H3`, `H5` (`H-PGF(0.5)+THRESHOLD-REPAIR`), `H-EA`, `H-SM`, `H-R`, `H-X`, and `H-SEP` (§6, `RANK-BRIDGE.md`; zero `H-FFT`, zero Band-3 `H1'`), using exact 3-class 2-adic `rho_dB` yield at `s = 0.80, kap = 0.322` (`89,200` exact-HNF cells), 100% pure-integer Band-0+1+2 sieve (`W_SB = 128`, flat-array zero-dict 3-pass + exact 200-bin Poisson heapsort + post-sort `sb_u7` list build, Shoup `< 2^62` 3-prime NTT, `7`-limb `f`-SRT), fixed-instance `H3` 1-owner SoA `Z[omega]` plain-sparse LA (ZERO TABM, `2^21..2^23` caps), 5-SM, and 62-row `LB = 59` `H-PGF(0.5)+THRESHOLD-REPAIR` + 6-bit ECM ladder + fresh-`e` `i.i.d.` boot (`K = 16`, `NODE_BUDGET = 2e6`, `t = 2^112`), the attack succeeds with probability `>= 0.39` over uniform challenges and fresh coins (opened-run success `>= 0.998933` via `FIX-U`), and

    E[W] <= 1 + (K / 2^64) * 2^88.349580962 <= 2^86.992667169 units (C = 20,000 instructions),
    K = 7201914631658704969.

- **Run gate & expected-work reading:** Open iff a fresh 64-bit coin `u < K = 7201914631658704969` (least integer with `(K/2^64)*0.998933 >= 0.39`; closed run fails at cost `1` unit). `86.992667169` is strictly an **expected** work bound `E[W]` (not every-run). Combined two-log work runs under hard aborting meter `t = 2^112` (`E[min(D, t)] <= E[D] <= 2^84.801522`, Markov `6.493e-9`), explicit `NODE_BUDGET = 2,000,000` (`7.250e-4`), and repeated-ideal abort (`1.383e-4`); split has a `2^112` meter. Work rounds up; `success` rounds down. Finite prototypes (`audit_bn_cm.py`) provide scoped compatibility evidence, not production proofs.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Exact `sigma(w) == w`, & Live Mod-`r` CRT

BN254 parameters: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h`:** Totally imaginary sextic over `Q(sqrt(-2))` (signature `(0,3)`, unit rank `2`, order-3 `sigma`): `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -2^9 * 19^2`, `h_K = 1`, `w_K = 2`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`.
- **Exact `sigma`-invariant twist:** `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3`, `w(y) = -y^5 + y^4 - y^2 - 2y`, `N(w) = 216`, `|w_j| = sqrt(6)` (`iso_err < 1e-10`), `(u*w)^((p^6 - 1)/2) = -1` in `F_{p^6}`, and **`w(sigma(y)) == w(y)` identically in `Q[y]/(h(y))`** (verified live).
- **Reciprocal pair:** `f_rec(x) = (-2y^5+2y^4-2y^2-4y+3)x^8 + (3y^5-3y^4+3y^2+6y)x^6 - 12x^4 + (-3y^5+3y^4-3y^2-6y+12)x^2 + (2y^5-2y^4+2y^2+4y-5)`, `g_rec(x) = u*w*x^2 - 1` (`N(lead)=N(const)=729`, `f_rec(x^2 -> (u*w)^-1) = 0 mod (p, h)`).
- **Order-6 Frobenius, Mod-`r` CRT, & Pure-Integer FD:** `tau(y, x) = (sigma(y), -x)` (`Frob^10`, `chi = p^10 mod r`, `chi^2 - chi + 1 = 0 mod r`); live DDF mod `r` gives `h:[1^6]`, `g_rec:[2]`, `f_rec:[2^4]`. Pure-integer `G_sym`/`D_fi` Voronoi FD certified for ALL radii via Dedekind `O_lat == O_h`, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), and outward-rounded `Decimal(55)` Rouché/Cholesky index-1 & convex-hexagon theorem (`KERNELS.md` §0).
- **Exact 3-Class 2-Adic Split & Field Map:** `2 = P_2^2` (`N(P_2) = 8, v_{P_2}(w) = 1`) partitions primitive `(a, b)` into `(v_2(N_f), v_2(N_g)) = (0,0)` (`7/9`), `(0,3)` (`1/9`), `(6,0)` (`1/9`); odd-part `alpha_f_odd = ALPHA_F_RAW + 2/3`, `alpha_g_odd = ALPHA_G_RAW + 1/3`. Canonical `phi_iso: F_{p^{12}} -> F_p[y,x]/(h, uwx^2 - 1)` (`wrapper <= 2^53.712288`); fresh uniform `e in [0, r-1]` maps `T = phi_iso(z*g^e) = c_0 + c_1 x`, `R = c_1/c_0 mod (p, h)` (`c_0 != 0` since `gcd(r, p^6-1) = 1`; `c_1 = 0` w.p. `1/r` in `FIX-U`), giving `vlog(a + bx) = vlog(z) + e*vlog(g) mod r` and scalar recovery `x = (vlog(T) - e*vlog(g))*vlog(g)^-1 mod r`, certified by `g^x == z` (`verify <= 2^11.000001`).

---

## 2. 100% Pure-Integer Band-0+1+2 Collection (`20.579` Ir/Cell, Flat-Array Lists, Shoup NTT, `7`-Limb `f`-SRT)

- **Asymmetric bounds & `s = 0.80` 200k-pool exact 3-class 2-adic yield:** `B_f = 2^61`, `B_g = 2^60` (sieved to `lx2 = P_mid = 2^54.75`, Stage-B `sl = 162, sb = 167`), `U' = 2^53.641526` `tau`-orbits, target `0.82 * U'` (`+0.0666`b margin). Exact 3-class `rho_dB` + class-matched `Retention(ng - v_{2g}, 162, 167)` across `777` exact-HNF Voronoi bins (`60,606 + 28,594 == 89,200`, live `1,500`-sample HNF audit `|z| <= 1.22`): `V = 2^97.5887` (`SE = 0.0063`b), `V_b = 2^52.7781`, `N_a = 2^45.9755`, `N_5 = 2^44.0710`, `L_b = 417.90`, `L_max = 851.34`, `N_f = 400.84`b (**analytic AM-GM ceiling `sup_{(a,b) in R} log2|N_f(a,b)| <= 445.33 < 448 bits = 7 limbs`**!), `N_g = 471.61`b, `r_g = 0.6307`.
- **Band-0+1+2 Sieve (`P_min = 14, P_0 = 30.0, P_b1_dense = 41.0, P_cut = 44.5, P_mid = 54.75, W_SB = 128`, ZERO Band 3):**
  - *Band 0 `[2^14, 2^30)` (`rep_b0 = 3.300` Ir/cell):* `B_a = 2014` owners of `2^35` slots (`2^62.72 B << 2^64 B`), 8x-unrolled `APPLY8 = 3.75` Ir/hit + `24.0` tail.
  - *Band 1 `[2^30, 2^44.5)` (`rep_lo = 7.444` total Ir/cell):* Dense `[2^30, 2^41)` uses 4x-unrolled `emit4` (`2.50`) + explicit `+24.0` tail per `(p,b)` row + `COPY64` (`1.25`) + `APPLY8` (`3.75`) = `7.50` Ir/hit; sparse `[2^41, 2^44.5)` uses direct 1-hit emit (`7.0`) + `COPY64` + `APPLY8` = `12.00` Ir/hit (zero tail); both use exact `min(bsearch, c_e4)`.
  - *8B-Packed Band 2 `[2^44.5, 2^54.75)` (`probe = 0.764, occ_sb = 1.291, hi_ap = 5.860` Ir/cell):* `W_SB = 128` (`max_sb = 8 <= 8`, `cr2 = 1.356`), 10-op `b2_swar_probe <= 12.0`, `b2_sb_enter = 16.0` (`11` setup + `<= 5` loop-exit), 12-op 8B-packed `b2_cand = 12.0` (enforcing BOTH `x_p < x_end` AND `rec < x_msb`, exact `M_t[t]` + `u64` wrap, `0` errs live), `C_HIT_HI = 12.00` (`7.0 + 1.25 + 3.75`), 32-cell `ORR` scan+zero `1.78125 <= 1.95` Ir/cell (`0/16,000` errs live).
- **Explicit Flat-Array Zero-Dict List Build (`lists <= 2^86.057995`):** Band 0/1 (`< 0.2%` of primes) use 2-pass `Q(a)`-order split-`2^{28}` 6-coord residue + prefix-sum + scatter (`C_LIST_A_LO = 52.0` Ir/entry, `0/64` errs live). Band 2 (`99.8%` of primes) uses 3-pass flat arrays (`b2_count_line` 14-op `<= 15.5` packing `tmp_sb` & `tmp_rec` + `b2_prefix_sum2` `4.5` Ir/SB zeroing `cnt[]` & `sb_u7[]` + `b2_scatter` `8 <= 8.5` = `24.0 == C_LIST_B2_BASE`, plus Pass-3 dispatch `5.0` Ir/SB -> `9.5 <= C_LBUCKET = 10.0` Ir/SB, exact 200-bin Poisson heapsort + `8.0/K` post-sort `sb_u7[sb] = qa_u7_tbl[entries[off[sb]] >> 8]`, avg `39.36` Ir/entry; builder RAM `2^52.86 B << 2^64 B`; ZERO dicts, live `flat_sbu7_exact_ok = True` & heapsort tests `0` errs!).
- **Stage B & Proved `7`-Limb `f`-SRT (`cof_b = 3.270` Ir/cell):** `C_NORMS = 30,000` (with early `gcd_{O_h}(a,b) == (1)`), `417`-chunk `P_14` strip (`~18.3k <= C_SB = 24,000`), `4` disjoint LP windows (gaps `[841, 757, 673]`, `0` fn live), Shoup `< 2^62` 3-prime NTT (`84*n*k + 260*n + 512*n`, 18-op `ntt62_butterfly` + 3-prime Garner CRT live) `7`-limb `f`-SRT (`sup_nf = 445.33 < 448`b, `gbatch_cost(61, leaf=7) ~ 8.34M` Ir), `Q_F_SAFETY = 4.0`, exact-HNF nonprimitive `0.0571` Ir/cell (`6` live witnesses). Totals: **`20.579` Ir/cell** `=> collection <= 2^87.668881`, global RAM **`2^99.078 B < 2^100 B`**.

---

## 3. Production `Z[omega]` Filtering (`2^12..2^23` 6-Scale Ladder), Fixed-Instance `H3` 4-Cap LA, & Unplanted 5-SM Rank Bridge

- **6-Scale Production Ladder (`2^12 -> 2^13` LIVE Python -> `2^20 -> 2^21 -> 2^22 -> 2^23` C, `W <= 161.0`, `(61, 60, exc = 0.82)`):** Every metric decreases strictly monotonically across all 6 scales (`2^21..2^23` C runs all strictly below hard caps; `w_mid -> 0.0` at `U' = 2^53.64`):

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `k5_exp` (`9.60/66`) | `tab_r` | `w_mid` |
|---|---|---:|---:|---:|---:|---:|
| `2^12` | **LIVE Python** (`W <= 161`, `kernel_err = 0`, `cols_lost = 0`, `backsub_err = 0`) | `0.014160` | `0.009277` | `23.16` | `0.7760` | `14.07` |
| `2^13` | **LIVE Python** (`W <= 161`, `cols_lost = 0`) | `0.006836` | `0.004028` | `22.38` | `0.7806` | `13.71` |
| `2^20` | C `filtsim` (`n = 4` seeds) | `0.001271` | `0.000751` | `16.20` | `0.7482` | `11.54` |
| `2^21` | C `filtsim` (`n = 8` seeds, max `0.001128 / 0.000679 / 16.29`) | `0.001031` | `0.000630` | `15.94` | `0.7436` | `11.25` |
| `2^22` | C `filtsim` (`n = 4` seeds, max `0.000965 / 0.000591 / 15.66`) | `0.000923` | `0.000570` | `15.64` | `0.7411` | `10.96` |
| `2^23` | C `filtsim` (`n = 4` seeds, max `0.000735 / 0.000473 / 15.59`) | `0.000690` | `0.000452` | `15.51` | `0.7390` | `10.67` |
| **Hard Caps (`H3`)** | **Deterministic Abort (`W <= 161.0`, `lazy = 0`, ZERO TABM)** | **`0.001150`** | **`0.000700`** | **`16.50`** | **`1.0000`** | **`~0.0`** |

- **Fixed-Instance `H3` (Explicit Exploratory Premise), Register-Resident `Z[omega]` LA, & Zero-Loss Excess Trim:** The monotone decay beyond `2^23` to `U' = 2^53.64` (under the tower `tau`-orbit splitting factor `6`) is a named exploratory STNFS premise (`H3`, §6 — NOT claimed proved from `2^23` alone). Every single run in `2^21..2^23` lies strictly below the 4 hard caps (`N_sw <= 0.00115*U'`, `N_l <= 0.00070*U'`, `k5 <= 16.50`, `tab_r <= 1.00`, `1.55x–1.67x` headroom above `2^23` mean!). Plain-sparse SoA register-resident `unit_run5_reg_z6 = 9.40 <= 9.60`, `gen_run1_z6 = 66.0`, `row_tail_3bucket_z6 = 722 <= 740.0` (1-step Barrett, ZERO TABM), `bw = (1024, 35)`, `nslots = 7`: **`la <= 2^84.601544`**, **`gen <= 2^78.214178`**, 1-owner RAM **`2^59.65 B << 2^64 B`** (`lazy = 0`). Degree-preserving excess trim + LIFO `repurge_stack` back-sub (`< 2^55` ops in `misc`) loses `0` post-peel columns (verified live at `2^13` and on `185` genuine `C6` rels).
- **5-SM & Unplanted 5-SM `F_1187^*` BW + DLP (`RANK-BRIDGE.md`, `hx_z6.py`):** $S_{\mathrm{full}} \in \mathbb{F}_r^{R \times 5}$ (`sm <= 2^73.831685`, rank loss `<= 5/r` under `H-SM`), padded `B = [M_prod | 0]`, 5-block-separated `X`, conditional Schwartz-Zippel `4 N_sw_cap / 2^64 = 3.504e-6`; verified live on unplanted `F_1187^*` (`r = 593`) + `Z[omega]` Schur identity.

---

## 4. Live 62-Row `H-PGF(0.5)` + `THRESHOLD-REPAIR` (`phi >= 0.01`) + 6-Bit ECM Ladder + `i.i.d.` Boot Descent (`LB = 59`)

- **Reached-Ideal `H-PGF(0.5)` Calibrated Envelope (ZERO Exact NB Claim!):** Let `M` be the number of distinct qualifying relations found at a reached node before child descents, and let `mu` be `evaluate_descent`'s Dickman/Murphy, retention, detection, and norm-envelope scale (`mu_min = 711.97`). `mu` is **NOT** asserted to equal the measured mean of `M`, and no exact negative-binomial distribution is claimed. For `t in [0, 0.8]`, exploratory STNFS hypothesis `H-PGF` bounds the reached-ideal ensemble PGF by the calibrated conservative envelope (`alpha = 0.5`):
  $$E[t^M] \le G(t) = \left(1 + \frac{0.5\,\mu\,(1 - t)}{r}\right)^{-r}, \qquad r = \begin{cases} 0.5 & \text{if effective window }\le 1\text{ prime},\\ 1.0 & \text{otherwise.}\end{cases}$$
- **`THRESHOLD-REPAIR` (ZERO Count/Solvability Independence Assumed!):** Let `psi = kA*(pm("f", nA) + pm("g", nA)) + kB*(pm("f", n) + pm("g", n))` be the union bound of child failure envelopes for a relation. Under the sequential hazard hypothesis `P(F_i | M >= i, F_1, ..., F_{i-1}, pre-entry history) <= psi`, for **every** integer `k >= 1`,
  $$P(\text{node fails}) \le P(M < k) + \psi^k \le \frac{G(t)}{t^{k-1}} + \psi^k \quad (\forall\,t \in (0, 0.8],\; \text{and } G(0) + \psi \text{ at } k = 1).$$
  Minimizing over the fixed finite certificate grid `k = 1..64, t = 0.01..0.80` and flooring at `PHI_FLOOR = 0.01` gives `phi = max(0.01, min_{k,t}(G(t)/t^{k-1} + psi^k))` (neither substituting `G(psi)` nor assuming count/solvability independence!). Across all 62 rows (`LB = 59`), **`0.01 <= phi <= 0.01731 <= 0.04`** and **`max psi = 0.19839 < 0.80`**. Expected attempted relations per node are at most `1 / (1 - psi)`, giving full-ball cofactor work (ZERO lazy cofactor discount!) `w_exp = own + child_work / (1 - psi)` and expected invocations `c_det = 1 + child_inv / (1 - psi)`.
- **6-Bit Geometric ECM Ladder & B-Window Quotient (`KERNELS.md` §4):** Every node and boot trial explicitly charges `extraction_ladder(lo, nA, nbits, 0.999, step=6) + (kA + kB)*PRIME`, stripping all `<= kA` window-A primes rung-by-rung with exact `_ecm_best` separation **`min_dag = 0.6184 > 0.50`** across all 266 DAG rungs (`T0` & `T`) and **`min_boot = 0.6389 > 0.50`** across all 16 Boot rungs up to `156.17`b (`tau = 1250`, verified live!) and certifying the `<= 1` remaining B-window quotient prime directly via Miller-Rabin (`PRIME`).
- **Fresh-Exponent `i.i.d.` Boot (`boot = [165, 88, 12, 4]`, `ea = [5, 1250]`, `K_boot = 16`), Explicit Node Budget (`2,000,000`), & Repeated-Ideal Loss:** Build the 1-owner boot prime-product tree (`x_EA = 64.52`, `2*P_x = 2^63.05 B < 2^64 B`) once (`tree_b`), then draw fresh uniform `e in [0, r-1]` whenever a radius-4 ball is exhausted and try up to `K_boot = 16` distinct boot relations sequentially (`psi_boot = 0.19136`, `E[tries] <= 1.2367`, `P(all 16 fail) <= psi_boot^16 = 3.234e-12`, plus per-log boot-search cap `1e-4` -> `phi_log_two = 2.000e-4`). Combined expected work across both logs: **`E[D] <= two_logs = 2^84.801522` units** (`two_dag = 2^84.035`, `two_boot = 2^83.523`), metered by hard cap **`t = 2^112`** (`Markov 6.493e-9`). Expected total node invocations across both logs: **`E[N_two] = 1,450.0`**. Aborting at explicit **`NODE_BUDGET = 2,000,000`** adds Markov tail **`7.250e-4`**, and aborting on any repeated ideal over bins `>= 2^59` (`>= 2^53.5` ideals/bin) adds birthday loss `<= (2e6)^2 / (2 * 2^53.5) = 1.383e-4` — both charged explicitly in `FIX-U`!
- **Scoped Empirical Support (24+8 Production Lattices & 224-Ideal Depth-2 `C6` Campaign):** `audit_bn_cm.py` verifies all `24` prod special-`q` bases (live 2-stage `11/40->9/10` Bareiss LLL, `24/24` paired `>= 1.74x`, `vol/2 hiR_min = 1.44x`) + `8` `G_T` boot certs (`1.12x`), and the 224-ideal exact-ideal-primitive depth-2 `C6` campaign (`100` d1 + `24` d2 + `100` ctrl, `8` replayed live incl. `2` 2-hop depth-2 chains: `0/224` empty, `M/mu = 1.2395`, `emp_pgf(t) <= G_{0.5}(t)` across `t in {0, 0.2, 0.5, 0.8}` for both `r = 1.0` and `r = 0.5` across all 5 cohorts, `S_obs/sv_up = 0.1961 < 1`). Neither geometry nor toy data proves production transfer; `H-PGF` extrapolates across degree, prime size, and adaptive depth as an explicit exploratory STNFS hypothesis.

---

## 5. Itemized 15-Row Ledger, Strict Additive Union Bound (`FIX-U`), & Exact Rational Run Gate

| Ledger Row | Computed (`log2` u) | Declared Upper Bound (`log2` u) |
|---|---:|---:|
| `collection` | `87.668880` | `87.668881` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `84.601544` | `84.601544` |
| `gen` | `78.214177` | `78.214178` |
| `lists` | `86.057995` | `86.057995` |
| `two_logs` | `84.801521` | `84.801522` |
| `canon` | `74.487831` | `74.487831` |
| `sm` | `73.831684` | `73.831685` |
| `misc` | `69.548386` | `69.548387` |
| `fb` | `65.038073` | `65.038073` |
| `split` | `81.838870` | `81.838870` |
| `gtree` | `59.594331` | `59.594332` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| **Expected opened run (`exp_run`)** | **`88.349580962`** | **`88.349580962`** |
| **Gated claim (`E[W]`)** | **`86.992667169`** | **`86.992667169`** |

**Strict additive union bound on opened-run success (`FIX-U`, under declared heuristics):**
- Main-phase conditional Schwartz-Zippel (`4 N_sw_cap / 2^64`): `3.504e-06`
- Filtering 4-cap abort given fixed-instance `H3`: `0.0`
- Split `2^112` budget abort: `8.329e-10`
- Fresh-exponent `i.i.d.` boot union (`2 * (psi_boot^16 + 1e-4)`): `2.000e-04`
- Explicit descent node-budget abort (`E[N_two] / 2,000,000`): `7.250e-04`
- Explicit descent repeated-ideal abort (`(2e6)^2 / (2 * 2^53.5)`): `1.383e-04`
- Combined two-log hard meter `t = 2^112` Markov tail (`E[D] / 2^112`): `6.493e-09`
- 64-round Miller-Rabin (`2^69 * 2^-128`): `1.735e-18`
- Group-order, boot `c_1 = 0`, & random SM compression (`8/r`): `3.655e-76`
- Numerical / floating-point error (`simd = 0`, zero `fp64`, zero Band 3): **`0.0`**
- **Opened-run success:** `1 - sum(failures) = 0.99893316 >= 0.998933`.

---

## 6. Declared Exploratory Heuristics (Named Production-Scale STNFS Model Premises)

- **H1:** Standard NFS smoothness under full de Bruijn `rho_dB` with odd-part Murphy `alpha`, exact 3-class disjoint 2-adic split, and live all-`89,200`-active-cell norm-histogram-weighted Euler-product retention/survivor laws (Band-0+1+2 up to `lx2 = P_mid = 2^54.75`, zero Band 3).
- **H2:** Uniform sampling over the certified 6-facet Voronoi FD at `s = 0.80, kap = 0.322` (`200k` pool, `89,200` active = `60,606` prim + `28,594` nonprim exact-HNF cells, `SE = 0.0063`b).
- **H3 (Exploratory Filtering Monotonicity / Fixed-Instance 4-Cap Premise):** Extending the 6-scale monotone trend (`2^12 -> 2^13 -> 2^20 -> 2^21 -> 2^22 -> 2^23`, tower `tau`-orbit factor `6`) to production `U' = 2^53.64` (`exc = 0.82`, `W <= 161.0`, zero-loss excess trim) is an explicit exploratory STNFS premise: the production `Z[omega]` matrix conforms to the 4 deterministic hard caps (`N_sw <= 0.00115*U'`, `N_l <= 0.00070*U'`, `k5 <= 16.50`, `tab_r <= 1.00` ZERO TABM). LA (`bw = (1024, 35)`), generator, and 1-owner RAM (`2^59.65 B`) are priced deterministically at the 4-cap (`lazy = 0`).
- **H-SM / H-R / H-X (`RANK-BRIDGE.md`):** `H-SM`: non-degenerate 5-SM $\chi$-eigenspace projections mod $r$; `H-R`: `rank M_prod = N_l`, 1D kernel of `[M_prod | S_prod]`; `H-X`: controllability/observability of padded `(B, Y, X)`, giving conditional Schwartz-Zippel `4 N_sw_cap / 2^64`.
- **H5 (`H-PGF(0.5)` + `THRESHOLD-REPAIR`) / H-EA / H-SEP:** At `LB = 59`, `(0)` `QBINS`/`BOOT_BINS` lower-bound paired smoothness & `QBINS_LO`/`BOOT_LO_OFF` upper-bound expected survivors (`H-EA`); `(i) H-PGF(0.5)`: reached-ideal ensemble PGF `E[t^M] <= (1 + 0.5*mu*(1-t)/r)^(-r)` on `t in [0, 0.8]` (`r = 0.5` for `<= 1` window prime, `r = 1` otherwise; `0.5` is a calibrated envelope parameter, no exact NB claimed), sequential after-failure child hazard `P(F_i | M >= i, F_1..F_{i-1}, history) <= psi` via `THRESHOLD-REPAIR` `P(fail) <= P(M < k) + psi^k` (`phi >= 0.01`, zero count/solvability independence assumed), adaptive repeat-collision bound `<= C_rep / N_ideals(bin)` (`N_eff >= li(2^59) = 2^53.5` in `[59, 65)` -> `repeat_prob <= 1.383e-4`), `NODE_BUDGET = 2e6`, and `t = 2^112` meter; `(ii) H-SEP`: 6-bit geometric ECM ladder isolates every window-A prime with per-curve cond separation `>= 0.50` (`min_dag = 0.6184` across all 266 DAG rungs, `min_boot = 0.6389` across all 16 Boot rungs to `156.17`b under exact `_ecm_best`, `63/63 = 100%` live Suyama-12 isolation; 2-addend `MONT5 = 367.76`, `KERNELS.md` §4). All finite-scale prototypes give scoped support; production scale transfer rests on these named premises.
