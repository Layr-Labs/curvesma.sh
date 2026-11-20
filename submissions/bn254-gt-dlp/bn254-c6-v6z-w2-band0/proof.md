# bn254-gt-dlp: 100% Pure-Integer `C6-m184_w216-rec` `(6,2)` STNFS (Band-0 Direct-Apply + Band-1/2 SWAR, Zero Band 3) — Expected Work `2^87.626011643` Units

Track `bn254-gt-dlp`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `x`. Score unit: `C = 20,000` pure baseline scalar AArch64 integer instructions (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, **ZERO floating-point instructions anywhere in the entire algorithm**, `p_num = 0` unconditionally). Every component is derived analytically here and in `RANK-BRIDGE.md`, `STORAGE.md`, `KERNELS.md`, and verified live in `audit_bn_cm.py` with zero external citations. No advice, smaller targets, or amortization across targets are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared exploratory heuristics `H1`, `H2`, `H3`, `H5`, `H-EA`, `H-SM`, `H-R`, `H-X`, and `H-SEP` (§6, `RANK-BRIDGE.md`; zero `H-FFT`, zero Band-3 `H1'`, closed-form `< 2^62` 3-prime NTT + Karatsuba, `KERNELS.md` §2), using exact 3-class 2-adic yield (`61,051` exact-HNF prim + `28,100` exact-HNF nonprim = `89,151` cells, `rho_dB`), **100% pure-integer Band-0 local-owner direct-apply (`[2^14, 2^30)`, `B_a = 567`, `N_a_blk = 2^35`, `2^62.72 B << 2^64 B`, `C_HIT_B0 = 4.00 + 24.0` tail, `0.4375` Ir/cell 64B `score_b0` stream) + buffered Band-1 `[2^30, 2^43)` & Band-2 `[2^43, 2^55.5)` SWAR `sb_u7` (ZERO Band 3)** with strict `< 2^64 B` per-owner buffered hits (`C_HIT_LO = 13.25`, `C_HIT_HI = 14.5`) & `P_14` Stage B (`sl = 164, sb = 169`), fixed-instance `H3` 1-owner `Z[omega]` plain-sparse LA (`bw = (1024, 35)`, ZERO TABM, `W <= 161.0`, `lazy = 0`), 5-SM, and 62-row `LB = 59` `k = 1` DAG + `B_max = 10` i.i.d. root-only boot retries (`NB(r = 1)`, `t = 2^112`), the attack succeeds with probability `>= 0.39` over uniform challenges and fresh coins (opened-run success `>= 0.999993` via `FIX-U`), and

    E[W] <= 1 + (K / 2^64) * 2^88.984455515 <= 2^87.626011643 units (C = 20,000 instructions),
    K = 7194280548710566105.

- **Run gate:** Open the run iff a fresh 64-bit coin `u < K = 7194280548710566105`, the least integer with `(K / 2^64) * 0.999993 >= 0.39` (`audit_bn_cm.py` checks in exact rationals that `K` passes and `K - 1` fails). A closed run fails at cost `1` unit.
- **Expected-work reading & hard phase meters:** `87.626011643` is strictly an **expected** work bound `E[W]` — **not** an every-run bound (modeled `Ns`, `gated`, `surv_b` are expected bounds under `H-EA`/`H1`). Each root run has deterministic `C_det = 663` and `T_boot_max`, at most `B_max = 10` i.i.d. root runs per log, and the combined two-log phase runs under hard aborting meter `t = 2^112` (`E[min(D, t)] <= E[D] <= 2^84.785146`, Markov tail `6.420e-9 < 1e-8` in `FIX-U`); split has a `2^112` aborting meter. Work bounds round upward; `success` rounds downward.
- **Scope of prototypes:** Finite-scale prototypes (`filtsim_z6`, genuine-`C6` filter, unplanted 5-SM `F_1187^*` BW+DLP, pure-int `G_sym`/`D_fi` Voronoi FD `382/382` + Dedekind maximality `O_lat == O_h` + exact `Fraction` `L D L^T` enum + outward-rounded `[lo, hi]` `Decimal(55)` Rouché/Cholesky index-1 & global hexagon proof, Band-2 exact `M_t[t]` `sb_u7` + Stage-B, unified inductive `delta = 9/10` Bareiss 12D LLL live on all `24` prod special-`q` + `q = 2^{220}-1` + `8` boot + exact-`(nA<=1,nB<=1)` ideal-primitive `224`-ideal campaign, live 32/32 `MONT5` 2-addend CIOS `335 <= 345 <= 367.76` + live 5-limb `Z[omega]` LA & 1-step Barrett `row_red5_z6` + per-curve `92/93 = 98.92%` Suyama-12 ECM isolation, live `1,500`-sample HNF pool audit) provide **compatibility evidence (supported, not proved — not certification)** for §6; production success follows from §6 via `FIX-U`.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Exact `sigma(w) == w`, & Live Mod-`r` CRT

BN254 parameters: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h`:** Totally imaginary sextic over `K_2 = Q(sqrt(-2))` (signature `(0, 3)`, unit rank `2`, order-3 `sigma`): `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -184,832 = -2^9 * 19^2`, `h_K = 1`, `w_K = 2`, `Reg = 2.374167`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`.
- **Exact `sigma`-invariant twist:** `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3` satisfies `sigma(y) = y^(p^4) = y^(p^10) mod (p, h)`. The twist `w(y) = -y^5 + y^4 - y^2 - 2y` has `N(w) = 216`, constant modulus `|w_j| = sqrt(6)` across all three complex places (`sigma`-isometry of `Q_a, Q_b`, `iso_err < 1e-10`), `(u*w)^((p^6 - 1)/2) = -1` in `F_{p^6}`, and **`w(sigma(y)) == w(y)` identically in `Q[y]/(h(y))`** (verified live in exact rationals).
- **Content-reduced reciprocal pair:** Dividing `w^4 x^8 F_BN(1/(w x^2))` by `c_0 = 4(1 + sqrt(-2))^2` (`N(c_0) = 2^12 * 3^6`, `21.509775`b) gives:
  `f_rec(x) = (-2y^5 + 2y^4 - 2y^2 - 4y + 3)x^8 + (3y^5 - 3y^4 + 3y^2 + 6y)x^6 - 12x^4 + (-3y^5 + 3y^4 - 3y^2 - 6y + 12)x^2 + (2y^5 - 2y^4 + 2y^2 + 4y - 5)`, `g_rec(x) = u * w * x^2 - 1`.
  Leading/constant coeffs of `f_rec` have norm `3^6 = 729`, `f_rec(x^2 -> (u*w)^-1) = 0 mod (p, h)`, `N(disc(f_rec)) = 2^96 * 3^96`.
- **Order-6 Frobenius, Mod-`r` CRT, & Pure-Integer `G_sym`/`D_fi` Voronoi FD:** `tau(y, x) = (sigma(y), -x)` acts as `Frob^10` (`chi = p^10 mod r`, `chi^2 - chi + 1 = 0 mod r`). Live DDF mod `r` confirms `h:[1^6]`, `g_rec:[2]`, `f_rec:[2^4]`, so `eps = r^2 - 1` annihilates all 6 residue-field unit groups. All geometry uses exact integer `sigma`-invariant `G_sym` (`G_SC = 196,608`, rel diff `7.17e-6`) and `D_fi = U_fi^T G_sym U_fi - G_sym` (`330` Ir/raw `a` inside `canon`, `KERNELS.md` §0; Dedekind `O_lat == O_h` + exact `Fraction` `L D L^T` enum `5,670` pts / `382/382` orbits + outward-rounded `[lo, hi]` `Decimal(55)` Rouché/Cholesky index-1 & global hexagon proof!).
- **Exact 3-Class 2-Adic Split & Odd-Part Alpha:** Since `2 = P_2^2` (`N(P_2) = 8, v_{P_2}(w) = 1`), every primitive `(a, b)` lies in one of 3 disjoint classes: `(v_2(N_f), v_2(N_g)) = (0,0)` (`7/9`), `(0,3)` (`b in P_2`, `1/9`), `(6,0)` (`a+b in P_2`, `1/9`). With `+0.08`b pad: `alpha_f = +1.580`b, `alpha_g = -0.039`b; odd-part `alpha_f_odd = ALPHA_F_RAW + 2/3`, `alpha_g_odd = ALPHA_G_RAW + 1/3`, summed explicitly over all 3 classes.
- **Explicit Field Map, Challenge-Dependent Boot, & Scalar Reconstruction:**
  1. Compute canonical `phi_iso: F_{p^{12}} -> F_p[y, x]/(h(y), u w x^2 - 1)` (`< 2^50` ops inside `wrapper <= 2^53.712288`).
  2. Draw fresh uniform `e in [0, r-1]` per trial and map `T = phi_iso(z * g^e) = c_0(y) + c_1(y) x`. Because `gcd(r, p^6 - 1) = 1` (verified live), `c_0 = 0` is impossible in `G_T`, and `c_1 = 0` iff `T = 1` (exact prob `1/r`, in `FIX-U`). Set `R(y) = c_1/c_0 mod (p, h)`; every `(a, b)` in `b == R a mod (p, h)` satisfies `a + b x = (a/c_0) T` in `F_{p^{12}}^*` (verified live), so `vlog(a + b x) = vlog(z) + e * vlog(g) mod r`.
  3. *Scalar reconstruction & verification:* Descend `a + b x` to recover `vlog(T)` and independently descend a randomized `g^{e'}` to recover `vlog(g)`, compute `x = (vlog(T) - e * vlog(g)) * vlog(g)^{-1} mod r`, and certify `g^x == z` (`verify <= 2^11.000001`).

---

## 2. 100% Pure-Integer Band-0+1+2 Score-Only Collection (`40.368` Instr/Cell, ZERO Band 3, Exact-HNF Nonprimitive Model, `< 2^62` 3-Prime NTT)

- **Asymmetric bounds, live 200k-pool exact 3-class 2-adic yield, & live `1,500`-sample HNF audit:** `B_f = 2^61` (batch-tested), `B_g = 2^60` (sieved up to `lx2 = P_mid = 2^55.5`, `(2^55.5, 2^60]` recovered in Stage-B `g`-factorization via `sl = 164, sb = 169`), `U' = 2^53.641526` `tau`-orbits, target `0.82 * U'` (`+0.0666`b margin). Exact 3-class 2-adic summation under `rho_dB` & class-matched `Retention(ng - v_{2g}, 164, 169)` across `901` Voronoi bins partitioned by `12x6` `O_h` HNF (`61,051 + 28,100 == 89,151`, `frac_reg = 0.445755` vs `0.44535`, `0.99816` guard): `V = 2^97.5260` (`SE = 0.0075`b, radii within `2.2%`), `V_b = 2^54.5448`, `N_a = 2^44.1469`, `N_5 = 2^45.5433`, `L_b = 512.52`, `L_max = 1044.10`, `N_f = 404.80`b, `N_g = 468.43`b, `r_g = 0.7281`. Live `1,500`-sample embedded Voronoi-FD + `12x6` `O_h` HNF audit (`641` active = `441` prim + `200` nonprim) matches `c6_data.json` within `|z| <= 1.50`.
- **Exact Integer Active-Region Cutoffs `M_t[t]`, Single-Owner Score Arena (`2^55.18 B << 2^64 B`), & Streamed 2-Arena SWAR Scan+Zero (`2.65` Ir/cell):** Canonical `a`-points are sorted ONCE by exact integer `Q(a) = a^T G_sym a` inside `ins_pre` (`lists <= 2^84.92`), giving exact contiguous active prefixes `a_slot < M_t[t] = bisect_right(Q_a_sorted, qa_max_t[t]) <=> Q(a) <= qa_max_t[t]` (`2^48.36 << canon = 2^72.66` units, `<= 8.4 KiB` L1 table/row, `KERNELS.md` §0). Row-major `uint16_t` storage + `bias_table` + `Q_a_sorted` occupies `2^55.1846 B << 2^64 B` (`STORAGE.md`); `+48` (`P_2`) + 8x-unrolled SWAR scan+zero (`0.4375 + 2.15625 = 2.594 <= 2.65` Ir/cell (explicit 64B `score_b0` stream+zero + `score_b12` scan+zero), `0/10,000` errors live) touches strictly `V` active cells.
- **Band-0+1+2 100% pure-integer score-only sieve (`P_min = 14`, `P_0 = 30.0`, `P_cut = 43.0`, `P_mid = lx2 = 55.5`, `slack = 164`, `slack_b = 169`; ZERO Band 3, ZERO `fp64`, zero prime recording):**
  - **Band 0 `[2^14, 2^30)` (`rep_b0 = 3.492` Ir/cell) + Band 1 `[2^30, 2^43)` (`rep_lo = 9.184` total Ir/cell):** `B_a = 567` independent `< 2^64 B` block owners (`2^62.72 B << 2^64 B`, `C_HIT_B0 = 4.00 + 24.0` tail, 2-block test `0` mismatches); Band 1 Parallel `8B keys[] + 8B a_slots[]` (`16 B * N_a`, `STORAGE.md`), per-row `x_row_base = t * N_a`, 7-op flush-checked 64-bit emit + `2.0` copy + `4.25` unrolled-4 apply = `C_HIT_LO = 13.25` Ir/hit + exact `min(bsearch, 11 + 14*(1 - exp(-mu*fa)) + mu*fa)` (`KERNELS.md` §1).
  - **Band 2 `[2^43, 2^55.5)` (`probe = 5.508`, `occ_sb = 2.014`, `hi_ap = 8.932` Ir/cell):** 8-lane SWAR `sb_u7` (`W_SB = 256`, `max_sb = 6 <= 8`, `C_SWAR_PROBE = 18.0`, `cr2 = 1.577`), exact 16B `(res_u64, a_slot_u64)` records tested via `a_slot < M_t[t]` in `C_CAND = 13.0` Ir (`1523/1523`, `0` mismatches live), buffered `C_HIT_HI = 14.5 >= 13.25` Ir (`KERNELS.md` §1). ZERO Band 3 (`p_num = 0`).
- **Score-Only `P_14` Stage B, Direct Empirical Exact-HNF Nonprimitive Model (`cof_b = 12.081` instr/cell under `H1`/`H2`), & Closed-Form `< 2^62` 3-Prime NTT:**
  - Early ideal primitivity reject (`gcd_{O_h}(a,b) == (1)` inside `C_NORMS = 30,000`), `417`-chunk `2^14`-smooth strip (`~18.3k << C_SB = 80,000`), `4` disjoint large-prime windows (`k = 0..3`, `0` fn live), `< 2^62` `ntt3_F`-priced SRT/factorization (zero `H-FFT`, zero Montgomery carry overflow), 64-round MR (`2^-128`). Nonprimitive `(a, b)` (`28,100` exact-HNF cells, 68 floored half-bit bins, 6 live HNF witnesses) contribute `0.09724` Ir/cell (`cof_b = 12.0806`).
- **Totals:** `40.368` instr/cell (`* 1.0033`) `=> collection <= 2^88.578201`; `lists <= 2^84.920671`; global RAM `2^99.36 B < 2^100 B` (`swar_probe = 22` sensitivity `2^87.659`, `STORAGE.md`).

---

## 3. Production `Z[omega]` Filtering, Fixed-Instance `H3` 4-Cap LA, & Unplanted 5-SM Rank Bridge

- **5-Scale Production Ladder (`2^13 -> 2^14` LIVE Python -> `2^20 -> 2^21 -> 2^22` C, `W <= 161.0`, `(61, 60, exc = 0.82)`):** `verify_bn254_v6.py` runs hexagonal `Z[omega]` filtering live at `2^13` (`kernel_errors = 0`) and `2^14`, while `2^20 -> 2^21 -> 2^22` C runs confirm `ALL4` monotonicity and `compute_w_mid_at_scale` verifies monotone decay of `w_mid -> 0.0` at `U' = 2^53.64`:

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `k5_exp` | `tab_r` | `w_mid` |
|---|---|---:|---:|---:|---:|---:|
| `2^13` | **LIVE Python** (`W <= 161`, `kernel_err = 0`, `cols_lost = 0`, `backsub_err = 0`) | `0.007568` | `0.004517` | `25.06` | `0.7783` | `13.71` |
| `2^14` | **LIVE Python** (`W <= 161`, `cols_lost = 0`) | `0.003967` | `0.002441` | `26.21` | `0.7678` | `13.37` |
| `2^20` | C `filtsim` (`n = 4` seeds) | `0.001271` | `0.000751` | `19.90` | `0.7482` | `11.54` |
| `2^21` | C `filtsim` (`n = 8` seeds) | `0.001031` | `0.000630` | `19.67` | `0.7436` | `11.25` |
| `2^22` | C `filtsim` (`n = 4` seeds) | `0.000923` | `0.000570` | `19.39` | `0.7411` | `10.96` |
| **Hard Caps (`H3`)** | **Deterministic Abort (`W <= 161.0`, `lazy = 0`, ZERO TABM)** | **`0.001600`** | **`0.000980`** | **`22.00`** | **`1.0000`** | **`~0.0`** |

- **Live Genuine-`C6`-Relation Filtering (`audit_bn_cm.py [live_genuine_c6_rel_filter]`):** Generates `220` genuine smooth `C6` relations from `g(x) = w x^2 - 1` and `f(x) = x^2 - (4 + w)` over BN254 `K_h, w` (`50` post-peel cols, `25` order-6 `tau`-orbits, `490` `Z[omega]` unit hits) and runs singleton peel + degree-preserving excess trim (`22` re-purged cols) + unit-pivot Markowitz merge live: `post_peel_cols_lost = 0`, `repurge_backsub_errors = 0`, `unit_entry_frac = 0.9593`, `kernel_errors = 0` mod `1,000,033`.
- **Fixed-Instance `H3`, Zero-TABM Pricing, & Zero-Loss Excess Trim:** Guard `totw_post <= 161.0*(N_alive - 1)` + deterministic abort if `N_sw > 0.00160*U'`, `N_l > 0.00098*U'`, `k5 > 22.00` (plain sparse `C_UNIT_5 = 13.80, C_GEN_1 = 66.0, c_tail_z6 = 950.0`, ZERO TABM `tab_r_cap = 1.00`). Degree-preserving excess trim + LIFO `repurge_stack` reverse back-sub (`< 2^55` ops inside `misc`) guarantees `0` post-peel `<= 2^59` cols lost (verified live at `U' = 2^13` on `2,890` re-purged cols and on `220` genuine `C6` rels). Pricing `la <= 2^85.961759`, `gen <= 2^79.052015` (`bw = (1024, 35)`, `nslots = 7`), 1-owner RAM `2^59.68 B << 2^64 B` (`lazy = 0`), `P(cap abort | H3) = 0`.
- **5-SM & Unplanted 5-SM `F_1187^*` BW + DLP (`RANK-BRIDGE.md`, `hx_z6.py`):** $S_{\mathrm{full}} = [W_f(P_\chi \Lambda_f) \mid W_g(P_\chi \Lambda_g)] \in \mathbb{F}_r^{R \times 5}$ (`sm <= 2^73.834205`, rank loss `<= 5/r` under `H-SM`). `RANK-BRIDGE.md` separates `H-R`, `H-X` (padded `B = [M_prod | 0]`, `Y = [S_1..S_5, B Z]`, 5-block-separated `X`), and Recovery Theorem (`4 N_sw_cap / 2^64 = 4.875e-6`). `hx_z6.py` verifies live on an unplanted 5-SM `F_1187^*` (`r = 593`) instance (`91` rels, `N = 59`, `D = 64`, `4` held-out): `rank(M)=59`, `rank(A)=63`, 4-SM control `ker_dim=0`, `rank(C_d)=59`, `rank(O_e A)=63`, 1-block-pad `60 < 63`, single ratio `[450]` across all `68` purged rels, DLP `275^468 == 679 mod 1187`, and `Z[omega]` Schur `prod_eq_schur = True`.

---

## 4. Live 62-Row Zero-Retry-DAG (`k = 1`) + `B_max = 10` Strictly I.I.D. Root-Only Boot Retries (`NB(r = 1)`, `t = 2^112`, `LB = 59`)

- **Zero-Retry-DAG (`k = 1`, `boot_rels = 1`) & Per-Run Cache (`P(repeat) = 0`):** Every node (62 rows, `q > 2^59`) evaluates `B(R_v)` ONCE and takes the FIRST USABLE relation (`k = 1`, zero internal retries); missing FB leaves `<= 2^59` discarded via `max(0, 1 - miss)`. Acyclic (`n_i < b_lo`, `C_det = 663`/run); per-run cache gives `P(repeat) = 0`, cleared between root runs `1..10` for exact i.i.d. independence.
- **Ultra-Weak Standard Geometric `NB(r = 1.0, alpha = 1.0)` `H5` (`p0_max = 1.404e-3`) + `B_max = 10` Strictly I.I.D. Root Runs + `t = 2^112` Meter:**
  - Under reached-ideal `H5` (`P(boot usable) >= prob_b` and `P(M_usable(q_sigma) = 0 | reached) <= (1 + mu_v)^(-1)`, `mu_min = 711.57`, `p0_max = 1.404e-3`), unconditional Boole gives `Phi_dag_1 = 0.245124`; fresh uniform `e in [0, r-1]` gives `E[min(N_boot, T_boot_max)] <= 1/prob_b` and boot-cap tail `<= 1e-4`, so **`Phi_1 = 0.245224 <= 0.35`**.
  - Build the single-owner boot prime-product tree (`x_EA = 64.42`, `2*P_x = 2^62.95 B < 2^64 B`) **ONCE** (`tree_b`), then execute up to `B_max = 10` strictly i.i.d. fresh-`e_t` root runs per log:
    - **Failure amplification:** `P(all 10 root runs fail) <= Phi_1^10 = 7.863e-7`/log `=>` **`2 * Phi_1^10 = 1.573e-6`** across both logs.
    - **Expected work (exact i.i.d. linearity, zero conditioning bias):** Because `I_t = 1{runs 1..t-1 failed}` is strictly independent of run `t`'s work `W_t` (`E[W_t] <= W_1 = trial_b/prob_b + over_b + W_dag_one` under `H-EA`):
      $$E[W_{\mathrm{log}}] \le \mathrm{tree}_b + E[W_1] \sum_{t=0}^9 \Phi_1^t = \mathrm{tree}_b + \frac{1 - \Phi_1^{10}}{1 - \Phi_1} E[W_1] \le \mathrm{tree}_b + 1.3249 \cdot E[W_1].$$
    - Across both logs: **`E[D] <= two_logs = 2^84.785146` units** (`two_dag = 2^84.322`). Wrapping the entire combined two-log phase in hard aborting meter **`t = 2^112` units** preserves `E[min(D, t)] <= E[D]` with Markov tail **`E[D]/2^112 = 6.420e-9 < 1e-8`** in `FIX-U`.
- **Live 12D Lattice Audit (24 Special-`q` + 8 `G_T` Boot) & 224-Ideal Usable-Relation + Depth-2 `C6` Campaign:**
  1. `audit_bn_cm.py [descent_12d_lattices]` verifies all `24` production special-`q` bases (`12f + 12g`, `2^66..2^213`, reduced live by unified `bareiss_lll_12d`) and `8` genuine `G_T` boot witnesses (`(1 + R_k x)^r in F_{p^6}`, `|det(B)| == p^6`, Bareiss-certified) over `B_12(4)/{+-1}` (`13,315,520` cells): `24/24` paired ratios `>= 1.74x >= 1.0`, `hi_radius` (exact `vol/2` mod `+-v`, `64` samples) min `1.44x >= 1.0`, `6/6` `cf`/`cg` sextiles `<= QBINS`, `6/6` `QBINS_LO`, `5/5` boot quintiles `<= BOOT_BINS` (`1.12x`).
  2. **224-Ideal Exact-Ideal-Primitive Depth-2 `C6` Campaign (`audit_bn_cm.py [c6_pgf_campaign]`):** On the `tau`-invariant `(6, 2)` tower over exact BN254 `K_h, w` (`u_t = 13, p_t = 79`), `B_12(sqrt(3))/{+-1}` (`1,024` pts), `LB = 13, nA = 15, nB = 16, T_surv = 18`, enforcing BOTH exact `nA <= 1 and nB <= 1` (`kA = kB = 1` matching `rows_2w`) AND exact `O_h` ideal primitivity (`a O_h + b O_h == O_h`, `model_env` scaled by `BASE_PRIM = 0.683542`), missing-FB discard (`29` primes) & per-factor ECM thinning (`p_ecm = 0.996`) across `100` depth-1 (`50g + 50f`) + `24` depth-2 (`12g + 12f`) + `100` controls (`16` balls incl. `4` 2-hop chains replayed live, `exact = True`): `0 / 224` empty (`mu_u in [27.64..57.67]`, `M_u in [26..95]`, `M/mu = 1.2395 >= 1.0`, `emp_p0 = 0 <= nb2 = 0.00875 <= nb1 = 0.02445`); BOTH `NB(r=1)` & `NB(r=2)` PGFs hold at `t in {0, 0.2, 0.5, 0.8}` across all 5 cohorts; `d1/ctrl = 1.0757`, `d2/ctrl = 1.0074`; `S_obs / sv_up = 0.1961` (`max = 0.3555 < 1.0`). Scoped compatibility evidence (supported, not proved); production transfer rests on `H5` and `H-EA`.

---

## 5. Itemized 15-Row Ledger, Strict Additive Union Bound (`FIX-U`), & Exact Rational Run Gate

| Ledger Row | Computed (`log2` u) | Declared Upper Bound (`log2` u) |
|---|---:|---:|
| `collection` | `88.578200` | `88.578201` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `85.961759` | `85.961759` |
| `gen` | `79.052015` | `79.052015` |
| `lists` | `84.920670` | `84.920671` |
| `two_logs` | `84.785145` | `84.785146` |
| `canon` | `72.659219` | `72.659220` |
| `sm` | `73.834204` | `73.834205` |
| `misc` | `69.550907` | `69.550908` |
| `fb` | `65.038073` | `65.038073` |
| `split` | `81.841390` | `81.841391` |
| `gtree` | `60.372479` | `60.372479` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| **Expected opened run (`exp_run`)** | **`88.984455515`** | **`88.984455515`** |
| **Gated claim (`E[W]`)** | **`87.626011643`** | **`87.626011643`** |

**Strict additive union bound on opened-run success (`FIX-U`, under declared heuristics `H1`, `H2`, `H3`, `H5`, `H-SM`, `H-R`, `H-X`):**
- Main-phase conditional Schwartz-Zippel (`4 N_sw_cap / 2^64`): `4.875e-06`
- Filtering 4-parameter hard-cap abort under fixed-instance `H3`: `0.0`
- Split `2^112` budget abort: `8.344e-10`
- Zero-retry-DAG + `B_max = 10` root-only i.i.d. boot retries (`2 * Phi_1^10`, `H5` `NB(r = 1.0, alpha = 1.0)` + per-run `1e-4` boot cap): `1.573e-06`
- Combined two-log phase hard meter `t = 2^112` Markov tail (`E[D] / 2^112`): `6.420e-09`
- Descent repeat-ideal abort (acyclic DAG + per-run cache): `0.0`
- 64-round Miller-Rabin (`2^69 * 2^-128`): `1.735e-18`
- Group-order, boot `c_1 = 0`, & random SM compression (`8/r`): `3.655e-76`
- Numerical / floating-point error (100% `-mgeneral-regs-only` integer, zero Band 3, zero `fp64`): **`0.0`**
- **Opened-run success:** `1 - sum(failures) = 0.99999354 >= 0.999993`.

---

## 6. Declared Exploratory Heuristics (Named Production-Scale Model Inequalities)

- **H1:** Standard NFS smoothness under full de Bruijn `rho_dB` with odd-part Murphy `alpha`, exact 3-class disjoint 2-adic split (`7/9:(0,0), 1/9:(0,3), 1/9:(6,0)`), and live all-`89,151`-active-cell norm-histogram-weighted Euler-product score-only retention/survivor laws (`61,051` exact-HNF primitive + `28,100` exact-HNF nonprimitive cells, Band-0+1+2 sieve up to `lx2 = P_mid = 2^55.5`, zero Band 3).
- **H2:** Uniform sampling over the certified 6-facet Voronoi FD (`200k`-sample pool, `89,151` active cells = `61,051` exact-HNF primitive + `28,100` exact-HNF nonprimitive, `SE = 0.0075`b; sample radii within `2.2%` of production).
- **H3 (Fixed-Instance Deterministic 4-Cap Conformance):** On the fixed public-seed production matrix (`U' = 2^53.64`, `exc = 0.82`), after filtering halts under deterministic `W_avg <= 161.0` (with zero-column-loss excess-row trim + reverse `repurge_stack` back-sub), the filtered `Z[omega]` matrix satisfies all four deterministic hard caps (`N_sw <= 0.00160*U'`, `N_l <= 0.00098*U'`, `k5 <= 22.00`, `tab_r <= 1.00` (ZERO TABM)). LA (`bw = (1024, 35)`), generator, and 1-owner RAM (`2^59.68 B`) are priced deterministically at the 4-cap (`lazy = 0`).
- **H-SM / H-R / H-X (`RANK-BRIDGE.md`):** `H-SM`: $\chi$-eigenspace SM projections $P_\chi \Lambda_f$ (rank 4), $P_\chi \Lambda_g$ (rank 1) non-degenerate mod $r$ (random $W_f, W_g$ rank-loss `<= 5/r`); `H-R`: `rank M_prod = N_l`, `A = [M_prod | S_prod]` has 1D kernel `(w_true, c_true)`, `c_true != 0`; `H-X`: finite controllability `rank(C_d) = N_l` and observability `O_e` injective on `Im(A)` on padded `(B, Y, X)`, giving conditional Schwartz-Zippel `4 N_sw_cap / 2^64` (supported, not proved, by live unplanted 5-SM `F_1187^*` DLP in `hx_z6.py`).
- **H5 / H-EA / H-SEP (Model Inequalities — Supported, Not Proved):** At `LB = 59`, `(0)` `QBINS`/`BOOT_BINS` lower-bound paired smoothness and `QBINS_LO`/`BOOT_LO_OFF` upper-bound two-sided residual survivors; `(i) H5`: (a) for each fresh-exponent i.i.d. uniform boot target `T = z * g^e in G_T`, `P(boot trial produces a usable relation) >= prob_b` (`E[min(N_boot, T_boot_max)] <= 1/prob_b`, tail `<= 1e-4`), and (b) for every DAG slot `sigma` reached along the `k = 1` first-usable-relation rule within a root run, `P(M_usable(q_sigma) = 0 | slot sigma reached) <= (1 + mu_model)^(-1)` (`NB(r = 1.0, alpha = 1.0)`, `p0_max = 1.404e-3`), combined via unconditional Boole's inequality (`Phi_1 = 0.245224 <= 0.35`, zero internal retries, zero intra-DAG independence) and amplified across `B_max = 10` strictly i.i.d. fresh-`e_t` root-only boot runs (`2 * Phi_1^10 = 1.573e-6`, `E[runs] <= 1.3249x`); `(ii) H-EA / H-SEP`: slot-conditional expected two-sided survivor counts (`Ns`, `gated`, `surv_b`) are upper-bounded by `QBINS_LO`/`BOOT_LO_OFF`, and Suyama-12 stored-prefix Stage-1/Stage-2 ECM isolates a smooth cofactor per curve with conditional probability `>= SEPARATION = 0.5` (`ceil(log1p(-0.999)/log1p(-0.5*p))`: stored-prefix binary search isolates any two smooth factors deterministically whenever their largest prime-power steps differ `s_1 != s_2`, prob `> 0.90 >> 0.50` across Stage-1 `pi(B_1)` steps and `> 0.99` across Stage-2 giant steps; supported live by `92/93 = 98.92%` per-curve single-prime isolation across `96` curves on `12` real 256-bit 3-factor composites incl. `55/56` multi-factor splits, and 2-addend scalar CIOS/SOS `335 <= 345 <= MONT5 = 367.76` Ir, `KERNELS.md` §4; `ntt3_F`, zero `H-FFT`). The 24+8 lattices and exact-`(nA<=1,nB<=1)` ideal-primitive 224-ideal depth-2 `C6` campaign provide scoped compatibility evidence (not certification); production transfer rests on `H5`, `H-EA`, and `H-SEP`.
