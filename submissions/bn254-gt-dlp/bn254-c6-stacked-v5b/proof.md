# bn254-gt-dlp: Self-Contained Pure-Scalar Content-Reduced Reciprocal CM-Sextic `(6,2)` STNFS (Full de Bruijn + 4-Cap `Z[omega]` LA) — Expected Work `2^87.535215439` Units

Track `bn254-gt-dlp`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `x`. The score unit is `C = 20,000` scalar AArch64 instruction equivalents (100% pure baseline scalar AArch64 integer/FP instructions; zero vector extensions, zero `FEAT_LSE` atomics). Every algebraic, geometric, sieving, 7D-walk, filtering, 12D-lattice, and descent component is specified analytically below and verified live inside `evidence.py` with zero external citations. No advice, smaller targets, or amortization across targets are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared exploratory heuristics H1-H5, H-EA, H-FFT, H-R, H-X, and H-SEP (Section 6), with **full de Bruijn local-density correction `rho_dB`**, **4-parameter deterministic hard-cap linear-algebra pricing `(N_sw_cap, N_l_cap, k5_cap, tab_r_cap)`** (`W_avg <= 201.0`, `lazy = 0`), **parameter-free analytical 7D Gaussian-volume-ratio Band-3 Fincke-Pohst pricing**, **full un-discounted Schirokauer maps**, **pure coarse `QBINS`/`BOOT_BINS` 12D descent using `min_k (P(M < k) + psi^k)` with an explicit conditional candidate-failure hazard**, and **100% pure baseline scalar AArch64 kernels**, the attack succeeds with probability at least `0.39` over uniform challenges and fresh coins (an opened run succeeds with probability at least `0.995968` by the strict additive union bound `FIX-U`; the run gate opens at most a share `0.391579` of runs), and

    E[W] <= 1 + (K / 2^64) * 2^88.887840704 <= 2^87.535215439 units (C = 20,000 instructions),
    K = 7223354755119366416.

- **Run gate:** Read one fresh 64-bit coin `u` and open the run iff `u < K = 7223354755119366416`, the least integer with `(K / 2^64) * 0.995968 >= 0.39` (`evidence.py` checks in exact rationals that `K` passes and `K - 1` fails). A closed run fails at cost `1` unit.
- **Metered phases:** Both individual logarithms share one combined `2^112`-unit meter; the deterministic split phase has its own `2^112`-unit interruptible budget. Every declared upper bound rounds strictly upward and `success` rounds strictly downward.
- **Scope of prototypes & heuristics:** Finite-scale simulations (`filtsim_z6`, 7D walk regimes, and 6-base/2-boot 12D lattice diagnostics) are exploratory diagnostic checks of algorithm mechanics and parameter conservatism. **Neither finite-scale simulator outcomes nor sample confidence intervals are substituted for the attack's production success probability**, which is derived explicitly from the declared exploratory premises in Section 6 via the strict additive union bound in Section 5.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Order-6 Frobenius, & Live Mod-`r` CRT

BN254 parameters: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h`:** Rank-2 CM sextic over `K_2 = Q(sqrt(-2))`:
  `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -184,832 = -2^9 * 19^2`, signature `(0, 3)`, unit rank `2`, class number `h_K = 1`, `w_K = 2`, `Reg = 2.374167`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`.
- **Order-3 automorphism & twist:** `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3` fixes `K_2` and satisfies `sigma(y) = y^(p^4) = y^(p^10) mod (p, h)`. The `sigma`-invariant twist `w = sqrt(-2)*(1 + sqrt(-2)) = -y^5 + y^4 - y^2 - 2y` has norm `N(w) = 216`, constant modulus `|w_j| = sqrt(6)` across all three complex places (giving exact `sigma`-isometry of `Q_a` and `Q_b`), and `(u*w)^((p^6 - 1)/2) = -1` in `F_{p^6}`.
- **Content-reduced reciprocal pair:** Dividing `w^4 * x^8 * F_BN(1 / (w * x^2))` by `c_0 = 4(1 + sqrt(-2))^2` (`N(c_0) = 2^12 * 3^6`, `21.509775` bits) yields the primitive balanced polynomial:
  `f_rec(x) = (-2y^5 + 2y^4 - 2y^2 - 4y + 3)x^8 + (3y^5 - 3y^4 + 3y^2 + 6y)x^6 - 12x^4 + (-3y^5 + 3y^4 - 3y^2 - 6y + 12)x^2 + (2y^5 - 2y^4 + 2y^2 + 4y - 5)`,
  `g_rec(x) = u * w * x^2 - 1`.
  Both leading and constant coefficients of `f_rec` have norm `729 = 3^6`, `f_rec(x^2 -> (u*w)^-1) = 0 mod (p, h)`, and `N(disc(f_rec)) = 2^96 * 3^96`.
- **Order-6 Frobenius:** `tau(y, x) = (sigma(y), -x)` acts on `F_{p^12}` as `Frob^10` (`(u*w)^((p^10 - 1)/2) = -1`), with eigenvalue `chi = p^10 mod r = 36u^3 + 18u^2 + 6u + 2` (`192` bits, `chi^2 - chi + 1 = 0 mod r`, `chi^3 = -1 mod r`). Every ramified prime of `O_{K_g}` lies above `{2, 3, 19}` and has `tau`-orbit length `< 6`.
- **Live Mod-`r` CRT Certification (`evidence.py [algebra]`):** Distinct-degree factorization modulo `r` executes live over **all 6 roots** of `h(y) mod r`: `h(y) mod r` has degree pattern `[1, 1, 1, 1, 1, 1]`, and at **every root**, `g_rec(x) mod r` has degree `[2]` and `f_rec(x) mod r` splits as `[2, 2, 2, 2]`. Thus `eps = r^2 - 1` (`508` bits) annihilates every residue-field unit group modulo `r` identically across all 6 components.
- **Exhaustive Small-Radius Fundamental Domain Bijection (`evidence.py [exhaustive_fd]`):** On `T_2(v) <= 50.0`, Fincke-Pohst enumeration finds all `1552` nonzero lattice points in `O_h`, grouping into `107` exact `<+-1, O_h^x, sigma>` orbits. With the 6-facet Voronoi cone + exact integer unit-matrix boundary tie-break (`UMATS[0..5]`), **`107/107`** orbits have exactly 1 canonical representative (`missed_orbits = 0`, `duplicate_orbits = 0`).
- **Small primes & Murphy alpha:** `2 = P_2^2` with `N(P_2) = 8` gives exact 2-adic classes `(v_2(N_f), v_2(N_g))` of `(0,0)` (`7/9`), `(0,3)` (`b in P_2`, `1/9`), and `(6,0)` (`a in P_2`, `1/9`). Murphy alpha plus `+0.08`b uniform safety margin: `alpha(f_rec) = +1.580`b, `alpha(g_rec) = -0.039`b.

---

## 2. Self-Contained Pure-Scalar Collection (`Arch B-Hybrid` + Analytical Band-3 7D Fincke-Pohst, Full de Bruijn `rho_dB`)

- **Asymmetric bounds & live norm-histogram-weighted full de Bruijn yield:** `B_f = 2^60` (batch-tested), `B_g = 2^59` (sieved), `U' = 2^52.666123` order-6 `tau`-orbits, relation target `0.95 * U'` with `+0.0666`b margin. Evaluated under **full de Bruijn local density `rho_dB(u, log2 N) = rho(u) - gamma * rho(u-1) / (log2(N) * ln 2)`** and **bin-by-bin live Euler-product `Retention(ng)`** across all bundled `200,000`-sample Voronoi FD pool quadrature bins (`c6_data.json`, `n_in = 89151`, `n_prim_in = 61463`, `frac_reg = 0.445755` matching the exact 6D ellipsoid integral `0.445421` within `0.001`):
  - Active cells `V = 2^97.5979` (`SE = 0.0077`b), `V_b = 2^54.5807`, canonical `a`-points `N_a = 2^44.1829`, `b`-rows `N_5 = 2^45.5733` (exact 5D Schur complement of `G_b`, average chord `L_b = 514.7`, exact max chord `L_max = 2 * R_b / sqrt(G_b[0][0]) = 1048.45`), relation-weighted norms `N_f = 404.77`b, `N_g = 468.50`b, effective `r_g = 0.9277`.
- **Clean `(t, a_slot)` active-only score arena & pure-scalar 64-bit SWAR fused scan+zero (`2.50` scalar instr/cell):**
  - Canonical `a`-points are sorted once in increasing order of `Q_a(a)`. Along any `b`-row `b_0 + t * 1` (`0 <= t < L_b`), `Q_b(b(t)) <= 1 - 0.322 * Q_a(a)` is equivalent to **`a_slot < M(t)`**, a contiguous prefix `[0, M(t))`. Row-major `(t, a_slot)` storage guarantees scan+zero touches **strictly the `V_kept = V` active cells** and fits on ONE single owner (`(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table) = 2^55.22 B << 2^64 B`).
  - Because `|u*w*a^2| / |b^2| > 2^62`, `|log2 N_g,rec(a,b) - log2 N_g,rec(a,0)| < 10^-15` bits (`evidence.py`: `0.0e+00` bits across 500 pairs). Precomputing `bias_table[a_slot] = 0x8000 - T(a_slot)` costs `< 2^40` units. For a fixed `t`-slice, `b(t)` is constant; when `b(t) in P_2` (`v_2 = 3`, prob `1/9`), we add `x7 = 0x0030003000300030` (`+48` quantized bits/lane, zero carry since `s_i, bias_i < 2^15`). Pure baseline scalar AArch64 loop (9 ops per 4 active cells = `2.25` + `0.25` setup = **`2.50` scalar instr/cell**, verified `0/10,000` errors live in `evidence.py`):
    `ldr x2,[x0]; ldr x4,[x5],#8; str xzr,[x0],#8; add x2,x2,x7; add x3,x2,x4; tst x3,x6; b.ne .hit; subs x1,x1,#4; b.gt .loop`.
- **Self-contained three-band collection accounting (`slack = 114`, `slack_b = 119`):**
  - **Band 1 `[2^14, 2^43)` (pointwise `b`-lookup):** Bucketed `a`-table (`0.5` tables/prime, `256` `Q_a`-shells). For each `b`, exact `Q_b(b)` selects `[0, s_max(b))` (`1.004` partial-shell overhang). Exact 6-instruction pure scalar loop (`ldr; ldrh; add; strh; cmp; b.lo`, `C_HIT_LO = 6.0`) + `V_b * (li(2^43) - li(2^14))` lookups at `C_LOOK_LO = 24.0` `=> 7.65` instr/cell.
  - **Band 2 `[2^43, 2^53.5)` (`Q_a`-sorted super-bucket `a`-table, `W_sb = 128`):** Every `b`-row chord (`L_max = 1048.45 < 64 * 128 = 8192`) spans `avg_sb = 5.01` super-buckets and at most 2 adjacent 64-bit bitmap words (exact expectation `bmw = 1.0627` words/probe at `C_PROBE = 18.0` `=> 1.48` instr/cell). Live Poisson integration of `1 - exp(-N_a * W_sb / p)` charges `C_OCC_SB = 13.0` scalar instr per occupied super-bucket (`2.22` instr/cell). Exact 6D ellipsoid integration gives `E[cands]/E[hits] = 1.8066` at `C_CAND = 13.0` (`qa <= qa_max`, `t < L_b`, `a_slot < M[t]`) + `C_HIT_HI = 9.0`/active hit `=> 7.10` instr/cell. Verified `0` mismatches over `1523` hits live in `evidence.py`.
  - **Band 3 `[2^53.5, 2^58.0)` (Self-Contained 7D Fincke-Pohst Coset Walk with Exact Analytical Gaussian-Volume-Ratio Scaling):**
    - For each Band-3 ideal `(p, s, +-rho)`, one 7D LLL basis `Bred` is computed for `L_7(p) = {(b, t) in Z^7 : b(s) - rho^-1 * t == 0 mod p}` in the weighted cylinder metric `Q_7(b, t) = Q_b(b)/R_b(a)^2 + (w_t * (2t/L_a - 1))^2 <= 1 + w_t^2` with `w_t = 0.40` (analytic minimum `1/sqrt(6) = 0.40825`, budget `C_LLL7 = 2^20` scalar ops/ideal `=> bases <= 2^57.350707` units).
    - **Analytical Scale-Invariant Volume-Ratio Derivation:** The enclosing 7D ellipsoid `B_7(R_7)` (`R_7 = sqrt(1 + w_t^2)`) has exact volume ratio over the target cylinder `Cyl_7 = B_6(1) x [-w_t, w_t]`:
      `VR_7 = Vol(B_7(R_7)) / Vol(Cyl_7) = ((16*pi^3/105) * (1 + w_t^2)^(7/2)) / ((pi^3/6) * 2*w_t) = 1.921299`.
      For any ideal prime `p` and run length `L_a`, the expected hit count per coset is the dimensionless ratio `h(p, L_a) = V_b_max * L_a / p` (spanning `h in [4.6, 103]` across `[2^53.5, 2^58.0)`), and `det(L_7) = Vol(Cyl_7) / h = Vol(B_7(R_7)) / (VR_7 * h)`. Under the standard 7D Gaussian volume heuristic for LLL-reduced bases (`||b_j*|| ~ det^(1/7)`), `E[leaves]/E[hits] = VR_7 = 1.9213` (bounded conservatively by `FP_LEAVES_PER_HIT = 2.15`, `+12%` pad) and the expected internal Fincke-Pohst nodes per hit depends only on `h`:
      `nu_node^Gauss(h) = (1/h) * (1 + sum_{d=1}^6 (V_d / V_7^(d/7)) * (VR_7 * h)^(d/7))`.
      Integrating `nu_node^Gauss(h(p, L_a))` live across `[2^53.5, 2^58.0)` with hit-density weight `dp / (p ln p)` and the `(fge16=0.989, 1-fge16)` run-length mixture (`2^5.52` vs `2^2.50`) yields **`fp_nodes_int = 4.7904`** internal nodes/hit.
    - **Live 3-Regime Verification Spanning Band-3 Hit Densities (`evidence.py [walk7_prototype]`):** Across `2510` hits (`0` missed), empirical `leaves/hit <= 1.95 <= 2.15` and empirical `nodes/hit` sits at `0.68x-0.78x` of `nu_node^Gauss(h)`:

| Band-3 Density Regime | `h_emp` (hits/coset) | Hits | Missed | `leaves/hit` | `nodes/hit` (Emp.) | `nu_node^Gauss(h)` (Pricing) |
|---|---:|---:|---:|---:|---:|---:|
| `sparse_band3_hi` | `6.57` | `197` | `0` | `1.8985` (`<= 2.15`) | `5.1878` | `6.8682` (`0.76x`) |
| `mid_band3` | `24.11` | `434` | `0` | `1.9954` (`<= 2.15`) | `2.8986` | `4.2546` (`0.68x`) |
| `dense_band3_lo` | `104.39` | `1879` | `0` | `1.9260` (`<= 2.15`) | `2.0000` | `2.6964` (`0.74x`) |

    - Per internal node `C_FP_NODE = 20.0`, per leaf `C_FP_LEAF = 20.0`, `1.05x` run overhang `=> fp_walk = 11.77` instr/cell; per-(run, ideal) 7D Babai/GS projection `C_RUN_INIT = 500.0 => babai_init = 6.86` instr/cell; bucket write (`12`) + owner-to-owner transfer (`8`) + scalar apply (`5`) = `25.0`/hit `=> bkt_top = 2.02` instr/cell.
  - **Live norm-histogram-weighted Stage-B / Cofactorization (`cof = 0.59` instr/cell, analytical derivations in Section 2.1):** Integrated across all bundled pool norm bins (`sf = 7.549e-07`, `s0 = 4.354e-08`) with `C_NORMS = 27,000`, `C_SB = 80,000` (score-only full 8-limb `P_14` strip, no prime recording in Band 1–3), `gbatch_cost(60)`, `C_TDIV`, `C_BATCH`, `gbatch_cost(59)`.
  - **Totals:** `42.184` scalar instr/cell (times `1.0033` work allowance) `=> collection <= 2^88.713573` units; 3-pass pure-C `a`-table build (`C_LIST_A = 35.0`) `=> lists <= 2^82.568935` units; peak RAM `2^97.84 B < 2^100 B`.

### 2.1 Self-Contained Analytical Scalar AArch64 Kernel Derivations
- **`C_NORMS = 27,000` (12-coordinate degree-8/2 norm evaluation):** Evaluating `g_rec(a, b) = u*w*a^2 - b^2` in `O_h` requires 2 `6x6` integer polynomial multiplies mod `h` (`~180` `mul`/`madd`) plus one `6x6` Bareiss/Sylvester resultant over 128-bit limbs (`~2,400` scalar ops). Evaluating `f_rec(a, b)` (degree 8 in `a^2, b^2`, 5 nonzero even terms) requires 4 squarings/products in `O_h` (`~1,200` ops) plus one `6x6` resultant over 512-bit (8-limb) integers (`6^3 / 3 = 72` 8-limb schoolbook muls at `~256` scalar ops each `= 18,432` ops + carry/sign `~3,000` ops). Total `~25,200 <= 27,000` scalar AArch64 instructions.
- **`C_SB = 80,000` (Score-Only Full 8-Limb Stage-B `P_14` Strip + Disjoint Large-Prime Window Test, ZERO Prime Recording in Band 1–3):** Band 1's 6-instruction hit loop ONLY accumulates the 16-bit quantized sieve score `s = sum_{p in [2^14, 2^58), p | N_g} floor(16 * log2 p)` (no sieved primes are recorded in Band 1, Band 2, or Band 3). Every first-stage survivor (`floor(16 * log2 N_g) - s <= 16 * 114`) enters Stage B with its full un-reduced `g`-norm `N_g` (`n_g = log2 N_g ~ 469 <= 512` bits `= 8` 64-bit limbs; 8-limb Montgomery multiplication `mulmod_8 = 360` scalar ops). Stage B strips ONLY primes `p < 2^14` using the fixed precomputed product `P_14 = prod_{p < 2^14} p` (`1,900` primes, `~23,440` bits `= 367` 64-bit words): (a) remainder `y = P_14 mod N_g` via `367` 8-limb Barrett/schoolbook steps (`367 * 120 = 44,040 <= 45,000` scalar ops), (b) nine 8-limb modular squarings `y_9 = y^(2^9) mod N_g` (`9 * 360 = 3,240` ops, covering prime powers up to `p^512`), (c) one 8-limb binary GCD `m = gcd(N_g, y_9)` (`<= 15,000` ops), and (d) **exact score-only disjoint large-prime window test** on `O = floor(16 * log2 N_g) - s - floor(16 * log2 m)`:
  - *Quantization rounding lemma:* An `n_g`-bit integer has at most `c = floor(ceil(n_g) / 14)` prime factors `>= 2^14` (at `n_g = 469`, `c = 33`), so `0 <= 16 * log2 S - s < c`, `0 <= 16 * log2 m - floor(16 * log2 m) < 1`, and `-1 < floor(16 * log2 N_g) - 16 * log2 N_g <= 0`. Hence ANY squarefree `2^59`-smooth norm `N_g = m * S * Q` with `k in {0, 1, 2}` large primes in `[2^58, 2^59)` MUST satisfy `O <= 16 * 119` AND the **disjoint large-prime window condition**:
    `O <= c + 1` (`k = 0`), **OR** `16 * 58 * k - 1 <= O <= min(16 * 59 * k + c + 1, 16 * 119)` for `k in {1, 2}`
    (at `n_g = 469, c = 33`, the three strictly disjoint intervals are `O <= 34`, `927 <= O <= 978`, and `1855 <= O <= 1904`; checking `O` against these 3 intervals takes `<= 6` scalar `cmp`/`b.hs` ops inside the `16,760`-op reserve for log/normalization/control).
  - Total `45,000 + 3,240 + 15,000 + 16,760 = 80,000` scalar AArch64 instructions at full 8-limb entry width, with zero sieved-prime division and zero `119`-bit pre-Stage-B cofactor assumption. Because `Retention.survivor_upper` (`model_bn254_c6.py`) evaluates the exact enclosing disjoint-window tolerance `(c + 2)/16 + 0.25` around `k * [58, 59]`, Stage B rejects all gap composites and lowers the survivor rate from `sf = 7.549e-07` to `s0 = 4.354e-08` with **zero false negatives** on squarefree smooth norms (verified live in `evidence.py [stageB_disjoint_window]`: `stageB_window_disjoint = True`, `retention_encloses_stageB = True`, `stageB_window_false_neg = 0/150` across `k in {0, 1, 2}`, `gap_composites_rejected = 100/100`). Only Stage-B survivors advance to `gbatch_cost(60)`, `C_TDIV`, `C_BATCH`, and `gbatch_cost(59)`.
- **`H-FFT` (`fft_F(N)` and 3-prime `ntt_F(N)`):** For `N` 64-bit words (`n = 2^k`), a 3-prime 64-bit NTT uses radix-2 Cooley-Tukey butterflies across 3 62-bit NTT primes (`69.75 * n * k` ops) plus CRT reconstruction (`90 * N` ops) plus explicit off-owner transfer (`3 * 504 + 128 = 1,640` scalar ops/word whenever `n > 2^26` local-owner words). Above `k = 55`, Schönhage-Strassen `8.5 * N * (log2 N)^2 + 504 * N` applies.
- **Un-discounted Schirokauer Maps (`sm <= 2^73.071624`):** Each relation computes `alpha^eps mod r^2` (`eps = r^2 - 1`, `508` bits) in `O_{K_f} / r^2 O_{K_f}` (degree 48) and `O_{K_g} / r^2 O_{K_g}` (degree 12) with **zero discount factor**.
- **Block Wiedemann 5-Lane `Z[omega]` Kernels:** `F_r` elements are 4×64-bit limbs (`32` B); a 5-lane block holds `160` B per coordinate. Unit entry `C_UNIT_5 = 7.8`; general non-unit entry `C_GEN_1 = 33.0`; projection `C_PROJ = 28.0`; dense Schirokauer column `C_FULL = 97.0`; 64-bit bucket gather `C_Z64 = 40.0`; 9-limb Barrett reduction `C_RED9 = 270.0`.

---

## 3. Production `Z[omega]` Filtering, Fixed-Instance Premise `H3`, 4-Cap LA, & Rank Bridge (`H-R`, `H-X`)

- **Complete 5-Scale Monotonic Production Ladder (`2^13 -> 2^14` LIVE Python -> `2^20 -> 2^21 -> 2^22` C, `W <= 201.0`, `(60, 59, exc=0.95)`):**
  `verify_bn254_c6.py` executes the **exact production `Z[omega]` filtering pipeline live inside `evidence.py`** (Dickman-`kappa` occurrence profile at `Lf=60, Lg=59, Nf=405.5, Ng=468.7, exc=0.95`, singleton + union-find clique purge, Markowitz merge `kmax=64, cmax=255` with deterministic pre-merge stop guard `W_avg <= 201.0`, Boudot order-labeled recorded pivots, and dense-band TABM subset-sum grouping; also callable on demand via `part="filtsim"` for any `log2U in [10, 15]` and `seed`), plus exact `Z[omega]` Markowitz-merge kernel preservation (`kernel_errors = 0`) and a **live `w_mid = 0` conservatism check**:

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `exp_unit` | `k5_exp` | `tab_r` |
|---|---|---:|---:|---:|---:|---:|
| `2^13` (`w_mid > 0`) | **LIVE Python** (`W = 196.01 <= 201`, `kernel_errors = 0`) | `0.00549` | `0.00354` | `71.27%` | `15.04` | `0.7752` |
| `2^13` (`w_mid = 0`) | **LIVE Python** (`W = 31.46 <= 201`, zero unsimulated mid-band) | `0.00317` | `0.00256` | `43.28%` | `22.09` | `0.8474` |
| `2^14` (`w_mid > 0`) | **LIVE Python** (`W = 199.33 <= 201`) | `0.00488` | `0.00287` | `79.09%` | `13.07` | `0.7735` |
| `2^20` (`w_mid > 0`) | C `filtsim` (`n = 4` seeds) | `0.00211` | `0.00118` | `90.00%` | `10.32` | `0.7500` |
| `2^21` (`w_mid > 0`) | C `filtsim` (`n = 8` seeds) | `0.00187` | `0.00107` | `90.41%` | `10.22` | `0.7466` |
| `2^22` (`w_mid > 0`) | C `filtsim` (`n = 4` seeds) | `0.00163` | `0.00094` | `90.80%` | `10.12` | `0.7441` |
| **Hard Caps** | **Deterministic Abort (`W <= 201.0`, `lazy = 0`)** | **`0.00310`** | **`0.00175`** | — | **`11.50`** | **`0.8500`** |

- **Deterministic Pre-Merge Weight Guard & 4-Parameter Post-Filtering Hard Caps (Zero Jensen / Composition Gap):**
  - Before committing any Markowitz merge step, the filter checks the exact post-merge average row weight `totw_post <= W_ROW * (N_alive - 1)` (`W_ROW = 201.0`) and halts before any merge that would exceed `W_ROW`, guaranteeing **`W_avg <= 201.0` deterministically on EVERY run** (verified live: `avg_w = 196.01` at `2^13`, `199.33` at `2^14`).
  - After filtering stops with `W_avg <= 201.0`, the algorithm **deterministically aborts** if ANY of `N_sw > 0.00310 * U' = 2^44.3326` (`+91%` above `2^22` mean), `N_l > 0.00175 * U' = 2^43.5077` (`+87%` above `2^22` mean), `k5 > 11.50` (`+40%` above `k5_ev = 8.19`), or `tab_r > 0.85` (`+14.2%` above `0.744`) occurs.
  - With `lazy = 0`, `W_LA`, `W_gen`, and 1-owner RAM are **strictly monotonically increasing in `(N_sw, N_l, k5, tab_r, W_ROW)`**. Pricing `la <= 2^84.915528`, `gen <= 2^78.513159`, and 1-owner RAM `<= 2^59.69 B < 2^64 B` **at the 4-parameter hard cap** holds deterministically on EVERY non-aborted run.
- **Fixed-Instance Structural Premise `H3` (`P(cap abort | H3) = 0`) & Scope of `filtsim_z6`:**
  - Because relation collection and filtering run deterministically on the fixed BN254 instance, **H3** is declared cleanly as a **fixed-instance structural premise**: the fixed production filtered matrix satisfies all four caps (`N_sw <= 0.00310*U'`, `N_l <= 0.00175*U'`, `k5 <= 11.50`, `tab_r <= 0.85`, `W_avg <= 201.0`), so `P(cap abort | H3) = 0` (no synthetic probability is assigned to a deterministic matrix, while the algorithm still enforces all four deterministic abort caps and prices LA/gen/RAM at the 4-cap worst case!).
  - At finite scales `2^13..2^22`, unsimulated mid-band primes contribute positive `w_mid > 0` that vanishes identically at production `U' = 2^52.67` (`lsmall > t_dense * L_g`); the live `u13_wmid0` check proves `N_sw/U'(w_mid=0) = 0.00317 < 0.00549` and `N_l/U'(w_mid=0) = 0.00256 < 0.00354`, confirming the finite-scale ladder is conservative in `w_mid`. Note that `filtsim_z6` generates column incidences independently at the production Dickman-`kappa` occurrence profile (providing exploratory support for fixed-instance premise `H3`, not a proof for jointly smooth retained relations).
- **Rank Bridge (`H-R`, `H-X`, & Conditional Schwartz-Zippel, Live `evidence.py [hx_z6_rank_bridge]`):**
  To separate deterministic matrix properties from fresh-coin Block Wiedemann failure probability cleanly:
  1. **Fixed-instance relation saturation (`H-R`):** On the filtered `D = N + s` (`s = 5` Schirokauer columns) system `A = [M | S]` modulo `r`, `rank(M) = N` and `rank(A) = D - 1`, so `ker(A) = span((w_log, c_log))` is strictly 1-dimensional with nonzero Schirokauer component `c_log != 0`.
  2. **Fixed-instance Krylov controllability & observability (`H-X`):** On the sparse bucketed starting support `Y = A Z` (`n` lanes) and sparse projection `X` (`m` lanes), controllability `rank(C_d) = D - 1` (full `N`-rank on sparse coordinates) and observability `O_e` is injective on `im(A)` (`rank(O_e A) = D - 1`).
  3. **Conditional fresh-coin failure bound:** Conditioned on `(H-R, H-X)`, fresh random Krylov starting vectors and projections fail to isolate the unique 1D kernel with probability at most the Schwartz-Zippel bound `4 * N_sw_cap / 2^64` (charged explicitly in Section 5 as `main_matrix`).
  - **Live 1D-Kernel Prototype (`verify_hx_z6`):** Over `Z[omega] mod r_toy` (`r_toy = 1000003`, `N = 35, s = 5, D = 40, m = 8, n = 8`), `evidence.py` builds an exact 1D-kernel `Z[omega]` matrix `A = [M | S]` and verifies live `rank(M) = 35`, `rank(A) = 39`, `rank(C_d) = 39` (`rank(C_d[:N]) = 35`), `rank(O_e A) = 39`, an explicit deficient-observability failure witness `rank(O_e_bad A) = 5 < 39`, and exact projected Block Wiedemann Hankel recovery of `(w_true, c_true)` (`rec_ok = True`).

---

## 4. Complete Live 62-Row 12D Descent (`LB = 59`, Pure `QBINS`/`BOOT_BINS`, `ADAPTIVE` Campaign, & `THRESHOLD-REPAIR`)

- **Finite-Sample Paired 12D Lattice Diagnostic Across `2^66..2^213` & Full `B_12(4)` (`evidence.py [descent_12d_lattices]`):**
  `c6_data.json` bundles **`6`** exact `12x12` LLL-reduced integer bases for representative `C6-m184_w216-rec` special-`q` ideals spanning `[60, 215]` bits across **BOTH sides** (`3` `g`-side + `3` `f`-side: `g:2^66, f:2^86, g:2^111, f:2^141, g:2^176, f:2^213`) and **`2`** exact `12x12` boot bases `mod p`. `evidence.py` verifies live:
  1. Exact Bareiss `12x12` integer determinants (`|det(B)| == q` for all 6 special-`q` bases; `|det(B)| == p^6` for both boot bases) and exact big-int congruences (`a(s) == rho * b(s) mod q`, `b == R * a mod (p, h)`) on all 12 columns of all 8 bases;
  2. **Exact `6x6` big-int algebraic resultant norms (`det_int` of multiplication matrices)** cross-checked against complex embeddings across all 6 special-`q` bases (`max_err = 3.90e-11` bits `< 10^-8`) AND exact integer divisibility **`q | N_f(a,b)`** on all 3 `f`-side bases and **`q | N_g(a,b)`** on all 3 `g`-side bases (`exact_side_divisibility_ok = True`);
  3. **Exact uniform sampling over the FULL discrete 12D ball `B_12(4) = {c in Z^12 \ {0} : ||c||_2^2 <= 16}`** via exact square partitions (`tot_pts / 2 == CELLS[4] = 13,315,520`, including coordinates `|c_i| in {3, 4}`, `46/400` sampled boot vectors have `|c_i| >= 3`);
  4. **Direct paired joint smoothness diagnostic across BOTH `f` and `g` sides (NO separately-sorted marginal substitution, NO upward sub-bin inflation!):**
     - `descent_bn254_c6.py` evaluates node smoothness **directly from pure 6-bin `QBINS`** and boot smoothness **directly from pure 5-quintile `BOOT_BINS`**.
     - Directly on the sampled `(a_i, b_i)` pairs (`n = 120` FULL `B_12(4)` vectors per special-`q` base = `720` pairs, subtracting `log2(q)` on the **`f`-side** for the 3 `f`-side bases and on the **`g`-side** for the 3 `g`-side bases), the empirical **joint** paired smoothness exceeds the `QBINS` paired model on every tested base (`g:2^66` (`4.19x`), `f:2^86` (`5.88x`), `g:2^111` (`6.87x`), `f:2^141` (`6.25x`), `g:2^176` (`7.32x`), `f:2^213` (`5.91x`)), and across large radii `R in {16, 64, 128}` (`n = 270` pairs, min paired ratio `2.95x >= 1.0x`), while all 6 `QBINS_LO` survivor-cost lower bounds lie below the empirical lower quantiles.
     - On `n = 400` exact uniform FULL `B_12(4)` boot vectors, all 5 empirical quintiles lie strictly below `BOOT_BINS` (`min = -196.91 >= BOOT_LO_OFF = -200.0`, `med = -181.86 <= -179.5`, `max = -174.45 <= -171.0`), and empirical boot smoothness exceeds `BOOT_BINS` (`1.15x >= 1.0x` at `L = 110`).
- **Node Failure Bound `min_{k >= 1} (P(M < k) + psi^k)` & `THRESHOLD-REPAIR` Counterexample (`evidence.py [h5_adaptive_and_threshold]`):**
  - Let `M` be the number of usable descent candidates generated at a node and let `E_j` be the event that candidate `j` fails to descend. Under declared heuristic `H5(ii)` (explicit conditional candidate-failure hazard):
    `P(E_j | E_1, ..., E_{j-1}, M >= j) <= psi` for every `j >= 1`.
  - For EVERY integer `k >= 1`, the elementary union bound gives (with NO independence assumption on `M` and NO invalid `PGF(psi)` substitution):
    `P(node fails) = P(all M candidates fail) <= P(M < k) + P(M >= k, E_1, ..., E_k) <= P(M < k) + psi^k`.
  - **Live `THRESHOLD-REPAIR` Counterexample Witness:** Taking `M = 1` with prob `0.10` (sole candidate fails) and `M = 100` with prob `0.90` (immediate success) gives true failure `0.10` while satisfying the `mu = 60, shape = 1` PGF envelope `E[t^M] <= (1 + 30(1 - t))^-1` on `[0, 0.8]`. The invalid `PGF(0.10) = 0.035714` strictly understates true failure `0.10`, whereas `descent_bn254_c6.failure_bound(60, 1, 0.10)` returns **`0.132258 >= 0.10`** (at `k = 1`), verified live in `evidence.py` (`threshold_counterexample.ok = True`).
  - Under declared heuristic `H5(i)`, the usable-candidate count `M` has lower tail bounded by the overdispersed Gamma-Poisson / Negative Binomial law with mean `mu / 2` and shape `r_sh in {0.5, 1.0}` (`P(M < k) <= inf_{0 < t <= 0.8} t^(-(k-1)) (1 + (mu / (2 r_sh))(1 - t))^(-r_sh)`). Minimizing over `k in {1, ..., 64}` in `descent_bn254_c6.py` gives `phi_node = max(0.01, min_k (P(M < k) + psi^k)) <= 0.04` on all 62 rows (`phi_max = 0.03077`, `psi_max = 0.39318 < 0.8`, `boot_psi = 0.44227`).
- **Real Reduced-Tower `C6` Recursive Descent Campaign (`ADAPTIVE`, `c6_adaptive.json` replayed LIVE in `evidence.py`):**
  To test `H5(i)–(iii)` directly on genuine recursive 12D descent trees of our tower `C6-m184_w216-rec` (reduced pair `f_4(x) = x^4 + w x^2 + w^2`, `g_2(x) = 4 w x^2 - 1` over `O_h`, smoothness bound `B_0 = 32768`, exact `O_h`-ideal primitivity via `rank_mod_p_6x12 == 6`, exact degree-1 root matching on every child, exact norm deduplication, **ZERO censoring**), we ran **`1360` roots across 4 settings (`f` and `g` sides)** yielding **`959` after-failure recursive attempts** across **`519` independent root clusters**, computing child-size-class row-specific `O/E` (calibrated strictly on first attempts):

| Setting | Roots OK | `mean(M)` | `P_emp(M<1) <= T_1` | After Fail | `OE_fail` | `ratio_1` | `OE_work` | Shared / Colls |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `f19_d1` (`f`, `19b`, `d=1`) | `339/400` | `3.73` | `0.0375 <= 0.3489` | `113/217` | `0.729` | `0.843` | `0.998` | `0` / `0` |
| `f20_d1` (`f`, `20b`, `d=1`) | `286/320` | `4.94` | `0.0031 <= 0.2882` | `169/270` | `0.847` | `1.093` | `1.004` | `0` / `0` |
| `f21_d1` (`f`, `21b`, `d=1`) | `301/320` | `6.53` | `0.0000 <= 0.2345` | `238/389` | `0.823` | `0.922` | `1.067` | `1` / `3` |
| `g19_d1` (`g`, `19b`, `d=1`) | `223/320` | `1.97` | `0.1469 <= 0.5035` | `40/83` | `0.832` | `1.008` | `0.954` | `0` / `0` |
| **POOLED (`959` after)** | **`1149/1360`** | — | **`12/12` `M<1,2,3` OK** | — | **`0.8095`** `[0.7685, 0.8434]` | **`0.9607`** | **`1.0307`** `[0.9970, 1.0632]` | **`1` / `3`** |

  `evidence.py` replays all `959` raw after-failure records from `c6_adaptive.json` live: all `12/12` empirical `M`-tails (`M < 1, 2, 3`) stay strictly inside `T_k(mean(M), alpha=0.5, r=1)`; pooled post-failure `OE_fail = 0.8095` (95% root-cluster bootstrap CI upper bound `0.8434 < 1.00`) and `ratio_1 = 0.9607 <= 1.00` show zero post-failure degradation; pooled child-work `OE_work = 1.0307` (95% CI upper bound `1.0632 < CHILD_X = 1.25`); and shared child ideals across retries occur in only `1/959` (`0.104%`) with `3` total tree collisions across all `1360` roots. (As with all finite-scale campaigns, these reduced-tower diagnostics support the plausibility of declared heuristic `H5`; they are NOT substituted for the production success probability.)
- **12D `dsieve`, Full Survivor Work, & Multi-Factor Stored-Prefix ECM Pricing (`descent_bn254_c6.py`, `T0` strictly at `LB = 59`):**
  - Per 12D cell: `HITS_PER_CELL = 2 * ln(59/8) = 4.00` hits x `SIEVE_HIT = 1,100` + `SCAN = 10` + `ZERO = 2` + `BOUND_CELL = 2,400`, plus `PAIR_SETUP = 2^25.1` and `LLL_NODE = 2^36`. Full survivor work is charged using `QBINS_LO` and `BOOT_LO_OFF = -200.0`.
  - `extraction_price(fbits, nbits, 0.999)` explicitly prices stored-prefix binary-search GCD isolation for ALL up to `F_ch = ceil(nbits / 20)` factors simultaneously (`SEPARATION = 0.5`), and `det` / `det_b` apply the exact Bernoulli union bound `1 - 0.001 * k_total <= 0.999^(k_total)` across all small and large-window factors (`T0` strictly at `LB = 59`, `T0 = 136` on `[75, 80)`).
- **Declared `H5(iii)` History-Conditional Orbit-Collision Bound:**
  - Under declared heuristic `H5(iii)`, conditioned on the prior descent history `F_{t-1}`, the `t`-th visited ideal above `2^59` lands in any specific previously visited `tau`-orbit with probability `<= C_peak / M_59`, where `M_59 = li(2^59) / 6 = 2^51.0977` and `C_peak = 4.0` (`+2.0` bits safety pad over uniform). Summing `sum_{t=1}^{N_max} (t - 1) * (C_peak / M_59)` gives `P(repeat) <= (N_max^2 / 2) * (C_peak / M_59) = 2^-9.5417` and Markov node-budget abort `P(N > N_max) <= E[N] / N_max = 2^-8.5410` (`E[N] = 2^11.7370`, `N_max = 2^20.278`, `two_logs = 2^84.018234 <= 2^84.018234` units).

<details>
<summary>Complete 62-Row Production Descent Plan Table (`LB = 59`, executed live in `descent_bn254_c6.py`)</summary>

| Bin | Side | `n` | `iA` | `kA` | `kB` | `R` | `thin` |
|---|---|---:|---:|---:|---:|---:|---:|
| `[59,65)` | `f` | `59` | `0` | `0` | `0` | `56` | `1` |
| `[59,65)` | `g` | `59` | `0` | `0` | `0` | `56` | `1` |
| `[65,70)` | `f` | `59` | `0` | `0` | `0` | `64` | `1` |
| `[65,70)` | `g` | `59` | `0` | `0` | `0` | `64` | `1` |
| `[70,75)` | `f` | `59` | `0` | `0` | `0` | `80` | `1` |
| `[70,75)` | `g` | `59` | `0` | `0` | `0` | `80` | `1` |
| `[75,80)` | `f` | `71` | `69` | `2` | `0` | `32` | `1` |
| `[75,80)` | `g` | `67` | `17` | `1` | `1` | `40` | `1` |
| `[80,85)` | `f` | `74` | `35` | `1` | `1` | `32` | `1` |
| `[80,85)` | `g` | `69` | `44` | `2` | `0` | `48` | `0` |
| `[85,90)` | `f` | `72` | `57` | `2` | `0` | `48` | `0` |
| `[85,90)` | `g` | `69` | `49` | `1` | `1` | `48` | `0` |
| `[90,95)` | `f` | `71` | `31` | `1` | `1` | `48` | `0` |
| `[90,95)` | `g` | `69` | `31` | `1` | `1` | `56` | `0` |
| `[95,100)` | `f` | `73` | `26` | `1` | `1` | `56` | `1` |
| `[95,100)` | `g` | `78` | `45` | `2` | `0` | `64` | `0` |
| `[100,105)` | `f` | `73` | `22` | `1` | `1` | `64` | `0` |
| `[100,105)` | `g` | `70` | `30` | `1` | `1` | `80` | `0` |
| `[105,110)` | `f` | `71` | `27` | `1` | `1` | `80` | `0` |
| `[105,110)` | `g` | `76` | `22` | `1` | `1` | `80` | `0` |
| `[110,115)` | `f` | `80` | `48` | `2` | `0` | `96` | `0` |
| `[110,115)` | `g` | `82` | `21` | `1` | `1` | `80` | `0` |
| `[115,120)` | `f` | `84` | `24` | `1` | `1` | `80` | `0` |
| `[115,120)` | `g` | `98` | `21` | `1` | `1` | `64` | `0` |
| `[120,125)` | `f` | `98` | `22` | `1` | `1` | `64` | `0` |
| `[120,125)` | `g` | `99` | `25` | `1` | `1` | `80` | `1` |
| `[125,130)` | `f` | `99` | `25` | `1` | `1` | `80` | `0` |
| `[125,130)` | `g` | `97` | `23` | `1` | `1` | `96` | `0` |
| `[130,135)` | `f` | `108` | `20` | `1` | `1` | `80` | `0` |
| `[130,135)` | `g` | `102` | `31` | `1` | `1` | `96` | `0` |
| `[135,140)` | `f` | `104` | `33` | `1` | `1` | `80` | `0` |
| `[135,140)` | `g` | `96` | `27` | `2` | `1` | `96` | `0` |
| `[140,145)` | `f` | `104` | `43` | `1` | `1` | `96` | `0` |
| `[140,145)` | `g` | `103` | `27` | `2` | `1` | `96` | `0` |
| `[145,150)` | `f` | `111` | `24` | `1` | `1` | `112` | `0` |
| `[145,150)` | `g` | `99` | `46` | `2` | `1` | `96` | `0` |
| `[150,155)` | `f` | `145` | `56` | `2` | `0` | `80` | `0` |
| `[150,155)` | `g` | `104` | `47` | `2` | `1` | `96` | `0` |
| `[155,160)` | `f` | `124` | `56` | `1` | `1` | `80` | `0` |
| `[155,160)` | `g` | `108` | `66` | `2` | `1` | `80` | `0` |
| `[160,165)` | `f` | `129` | `62` | `1` | `1` | `80` | `0` |
| `[160,165)` | `g` | `109` | `64` | `2` | `1` | `96` | `0` |
| `[165,170)` | `f` | `129` | `62` | `1` | `1` | `96` | `0` |
| `[165,170)` | `g` | `134` | `24` | `2` | `1` | `112` | `0` |
| `[170,175)` | `f` | `119` | `72` | `5` | `1` | `14` | `0` |
| `[170,175)` | `g` | `124` | `64` | `2` | `1` | `80` | `0` |
| `[175,180)` | `f` | `118` | `72` | `5` | `1` | `16` | `0` |
| `[175,180)` | `g` | `127` | `64` | `2` | `1` | `96` | `0` |
| `[180,185)` | `f` | `119` | `72` | `5` | `1` | `16` | `0` |
| `[180,185)` | `g` | `142` | `48` | `2` | `1` | `96` | `0` |
| `[185,190)` | `f` | `126` | `72` | `5` | `1` | `14` | `0` |
| `[185,190)` | `g` | `144` | `42` | `2` | `1` | `128` | `0` |
| `[190,195)` | `f` | `125` | `72` | `5` | `1` | `16` | `0` |
| `[190,195)` | `g` | `153` | `48` | `2` | `1` | `96` | `0` |
| `[195,200)` | `f` | `126` | `72` | `6` | `1` | `16` | `0` |
| `[195,200)` | `g` | `160` | `48` | `2` | `1` | `96` | `0` |
| `[200,205)` | `f` | `125` | `72` | `6` | `1` | `16` | `0` |
| `[200,205)` | `g` | `163` | `48` | `2` | `1` | `112` | `0` |
| `[205,210)` | `f` | `132` | `72` | `5` | `1` | `16` | `0` |
| `[205,210)` | `g` | `170` | `48` | `2` | `1` | `112` | `0` |
| `[210,215)` | `f` | `139` | `72` | `5` | `1` | `14` | `0` |
| `[210,215)` | `g` | `151` | `72` | `6` | `1` | `14` | `0` |

</details>

---

## 5. Itemized Ledger, Strict Additive Union Bound (`FIX-U`), & Run Gate

| Ledger Row | Computed (`log2` u) | Declared Upper Bound (`log2` u) |
|---|---:|---:|
| `collection` | `88.713572` | `88.713573` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `84.915527` | `84.915528` |
| `gen` | `78.513159` | `78.513159` |
| `lists` | `82.568934` | `82.568935` |
| `two_logs` | `84.018234` | `84.018234` |
| `canon` | `72.695158` | `72.695159` |
| `sm` | `73.071623` | `73.071624` |
| `misc` | `68.788325` | `68.788326` |
| `fb` | `64.061920` | `64.061920` |
| `bases` | `57.350707` | `57.350707` |
| `split` | `80.275459` | `80.275460` |
| `gtree` | `60.174753` | `60.174754` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| **Expected opened run (`exp_run`)** | **`88.887840704`** | **`88.887840704`** |
| **Gated claim (`E[W]`)** | **`87.535215438`** | **`87.535215439`** |

**Strict additive union bound on opened-run success (`FIX-U`, zero survival-factor multiplication; under declared heuristics H1-H5, H-R, H-X):**
- Main-phase conditional Schwartz-Zippel (`4 N_sw_cap / 2^64`, given `H-R, H-X`): `4.804e-06`
- Filtering 4-parameter hard-cap abort (under fixed-instance premise `H3`): `0.000000` (`== 0.000000`)
- Split `2^112` budget abort: `2.818e-10`
- Log-phase `2^112` meter Markov bound: `3.773e-09`
- Descent node-budget abort (`H5`, `2^-8.5410`): `2.685e-03`
- Descent history-conditional repeated-ideal collision (`H5(iii)`, `C_peak = 4.0`, `2^-9.5417`): `1.342e-03`
- Group-order & numerical (`2/r + 2^-48`): `3.553e-15`
- **Opened-run success:** `1 - sum(failures) = 0.99596882 >= 0.995968`.

---

## 6. Declared Exploratory Heuristics

- **H1 / H1':** Standard NFS smoothness under full de Bruijn local density `rho_dB` with Murphy `alpha`, exact 2-adic split, live norm-histogram-weighted Euler-product retention/survivor laws on `N_f`, `N_g`, and the 7D Gaussian volume heuristic for Band-3 LLL-reduced coset bases (verified diagnostically across 3 hit-density regimes `h in [6, 105]` where empirical `nodes/hit` sits at `0.68x-0.78x` of the analytical Gaussian bound).
- **H2:** Uniform sampling over the certified 6-facet Voronoi fundamental domain (`200,000`-sample bundled pool quadrature, `SE = 0.0077`b).
- **H3 (Fixed-Instance Production Filtering 4-Cap Conformance):** On the fixed deterministic BN254 production instance (`U' = 2^52.67`), after filtering halts under the deterministic pre-merge weight guard `W_avg <= 201.0`, the filtered hexagonal `Z[omega]` matrix satisfies all four deterministic hard caps simultaneously — `N_sw <= 0.00310 * U'`, `N_l <= 0.00175 * U'`, `k5 <= 11.50`, and `tab_r <= 0.85` (so `P(cap abort | H3) = 0`; the algorithm still enforces all four deterministic abort caps and prices LA, generator, and 1-owner RAM deterministically at the 4-cap worst case with `lazy = 0`). Supported qualitatively by the monotonic 5-scale ladder `2^13 -> 2^14 -> 2^20 -> 2^21 -> 2^22`, live `w_mid = 0` conservatism check, `+87%..+91%` dimension / `+14%..+40%` per-row cap margins, and super-critical peeling `e_crit = 0.718 << 0.95`.
- **H4 / H-R / H-X (Rank Bridge):** Fixed-instance relation saturation `H-R` (`rank(M) = N`, `rank(A) = D - 1`, 1D kernel with `c_log != 0`) and fixed-instance Krylov controllability/observability `H-X` (`rank(C_d) = D - 1`, `O_e` injective on `im(A)` on sparse `X, Y` support; illustrated live in `verify_hx_z6`), reducing Block Wiedemann failure probability cleanly to the conditional Schwartz-Zippel bound `4 N_sw_cap / 2^64`.
- **H5 / H-EA / H-FFT / H-SEP (STNFS Descent Norm, Tail, Conditional Hazard, & Collision Premise):** At `B_g = 2^59`, (0) across production special-`q` (`60..215` bits) and boot lattices, paired candidate smoothness is lower-bounded by pure `QBINS` / `BOOT_BINS` and survivor rates are upper-bounded by `QBINS_LO` / `BOOT_LO_OFF = -200.0` (supported diagnostically on 6 special-`q` bases `2^66..2^213`, radii `R = 4..128`, and 2 boot bases on exact FULL `B_12(4)`); (i) usable candidate counts `M` satisfy the overdispersed Negative Binomial lower-tail bound (`mean = mu/2`, `r_sh in {0.5, 1.0}`); (ii) successive candidates at a node satisfy the explicit conditional failure hazard `P(E_j | E_1, ..., E_{j-1}, M >= j) <= psi` (yielding `phi <= min_k (P(M < k) + psi^k)`, never `PGF(psi)`); and (iii) conditioned on prior descent history, any next child ideal hits a specific previously visited orbit above `2^59` with probability `<= 4.0 / (li(2^59)/6)` — all three components tested directly on `1360` reduced-tower `C6` recursive descent roots (`959` after-failure attempts in `c6_adaptive.json`), together with scaled remainder tree FFT pricing (`H-FFT`) and stored-prefix multi-factor ECM isolation (`H-EA`, `H-SEP`).
