# bn254-pairing-inversion: Pure-Scalar `C6-m184_w216-rec` `(6,2)` STNFS (`sigma`-Folded 12D Band-3 + 8-Lane SWAR `sb_u7` + Score-Only Stage B + Fixed-Instance `H3` 4-Cap `Z[omega]` LA) — Expected Work `2^86.548991570` Units

Track `bn254-pairing-inversion`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `(P0, [x]Q0)` in `G1 x G2` satisfying `e(P0, [x]Q0) = e(P0, Q0)^x = z`. After recovering `x` via the GT-DLP algorithm below, one fixed-base `G2` scalar multiplication `[x]Q0` (`<= 69,456,896` scalar AArch64 instructions, `output <= 2^11.761903` units) produces the valid preimage `(P0, [x]Q0)`. The score unit is `C = 20,000` pure baseline scalar AArch64 instruction equivalents (`-mgeneral-regs-only`, `simd = 0`, zero vector extensions, zero `FEAT_LSE` atomics). Every algebraic, geometric, sieving, 12D-walk, filtering, rank-bridge, 12D-lattice, and descent component is specified analytically below (and in companion files `RANK-BRIDGE.md`, `STORAGE.md`, `KERNELS.md`) and verified live inside `audit_bn_cm.py` with zero external citations. No advice, smaller targets, or amortization across targets are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared exploratory heuristics `H1-H5`, `H-EA`, `H-FFT`, `H-R`, `H-X`, and `H-SEP` (Section 6 and `RANK-BRIDGE.md`), with **full de Bruijn local-density correction `rho_dB`**, **100% score-only 3-band `g`-side sieving & `P_14` Stage B (zero sieve prime recording)**, **`sigma`-orbit-folded 12D Ellipsoid Fincke-Pohst Band-3 enumeration**, **fixed-instance `H3` 4-parameter deterministic hard-cap linear-algebra pricing (`W_avg <= 161.0`, `lazy = 0`)**, **full un-discounted Schirokauer maps**, **pure coarse `QBINS`/`BOOT_BINS` 12D descent using `min_k (P(M < k) + psi^k)`**, and **100% pure baseline scalar AArch64 kernels**, the attack succeeds with probability at least `0.39` over uniform challenges and fresh coins (an opened run succeeds with probability at least `0.997079` by the strict additive union bound `FIX-U`; the run gate opens at most a share `0.391143` of runs), and

    E[W] <= 1 + (K / 2^64) * 2^87.903225261 <= 2^86.548991570 units (C = 20,000 instructions),
    K = 7215306097858570014.

- **Run gate:** Read one fresh 64-bit coin `u` and open the run iff `u < K = 7215306097858570014`, the least integer with `(K / 2^64) * 0.997079 >= 0.39` (`audit_bn_cm.py` checks in exact rationals that `K` passes and `K - 1` fails). A closed run fails at cost `1` unit.
- **Metered phases:** Both individual logarithms share one combined `2^112`-unit meter; the deterministic split phase has its own `2^112`-unit interruptible budget. Every declared upper bound rounds strictly upward and `success` rounds strictly downward.
- **Scope of prototypes & heuristics:** Finite-scale prototypes (`filtsim_z6`, `hx_z6`, `sb_u7`, `sigma`-folded 12D walk, and 6-base/2-boot 12D lattice diagnostics) are exploratory diagnostic checks of algorithm mechanics and parameter conservatism. **Neither finite-scale simulator outcomes nor sample confidence intervals are substituted for the attack's production success probability**, which is derived explicitly from the declared exploratory premises in Section 6 via the strict additive union bound in Section 5.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Exact `sigma(w) == w`, & Live Mod-`r` CRT

BN254 parameters: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h`:** Rank-2 CM sextic over `K_2 = Q(sqrt(-2))`:
  `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -184,832 = -2^9 * 19^2`, signature `(0, 3)`, unit rank `2`, class number `h_K = 1`, `w_K = 2`, `Reg = 2.374167`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`.
- **Order-3 automorphism & exact `sigma`-invariant twist:** `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3` fixes `K_2` and satisfies `sigma(y) = y^(p^4) = y^(p^10) mod (p, h)`. The content-reduced twist `w(y) = -y^5 + y^4 - y^2 - 2y` has norm `N(w) = 216`, constant modulus `|w_j| = sqrt(6)` across all three complex places (giving exact `sigma`-isometry of `Q_a` and `Q_b`, `iso_err < 1e-10`), `(u*w)^((p^6 - 1)/2) = -1` in `F_{p^6}`, and satisfies **`w(sigma(y)) == w(y)` identically in `Q[y]/(h(y))`** (verified live in exact rational arithmetic in `audit_bn_cm.py [algebra]`: `sigma_fixes_w_exact = True`).
- **Content-reduced reciprocal pair:** Dividing `w^4 * x^8 * F_BN(1 / (w * x^2))` by `c_0 = 4(1 + sqrt(-2))^2` (`N(c_0) = 2^12 * 3^6`, `21.509775` bits) yields:
  `f_rec(x) = (-2y^5 + 2y^4 - 2y^2 - 4y + 3)x^8 + (3y^5 - 3y^4 + 3y^2 + 6y)x^6 - 12x^4 + (-3y^5 + 3y^4 - 3y^2 - 6y + 12)x^2 + (2y^5 - 2y^4 + 2y^2 + 4y - 5)`,
  `g_rec(x) = u * w * x^2 - 1`.
  Both leading and constant coefficients of `f_rec` have norm `729 = 3^6`, `f_rec(x^2 -> (u*w)^-1) = 0 mod (p, h)`, and `N(disc(f_rec)) = 2^96 * 3^96`.
- **Order-6 Frobenius:** `tau(y, x) = (sigma(y), -x)` acts on `F_{p^12}` as `Frob^10`, with eigenvalue `chi = p^10 mod r = 36u^3 + 18u^2 + 6u + 2` (`192` bits, `chi^2 - chi + 1 = 0 mod r`, `chi^3 = -1 mod r`).
- **Live Mod-`r` CRT & Exhaustive Voronoi FD (`audit_bn_cm.py`):** Distinct-degree factorization mod `r` confirms `h` splits `[1,1,1,1,1,1]` and at every root `g_rec` has degree `[2]` and `f_rec` splits `[2,2,2,2]`, so `eps = r^2 - 1` (`508` bits) annihilates all 6 residue-field unit groups identically. On `T_2(v) <= 50.0`, Fincke-Pohst verifies all `1,552` nonzero lattice points group into `107/107` exact `<+-1, O_h^x, sigma>` orbits (`missed = 0`, `dups = 0`).
- **2-adic classes & Murphy alpha:** `2 = P_2^2` (`N(P_2) = 8`) gives `(v_2(N_f), v_2(N_g))` of `(0,0)` (`7/9`), `(0,3)` (`b in P_2`, `1/9`), `(6,0)` (`a in P_2`, `1/9`). Murphy alpha with `+0.08`b safety pad: `alpha(f_rec) = +1.580`b, `alpha(g_rec) = -0.039`b.

---

## 2. Pure-Scalar Score-Only Collection (`29.855` Instr/Cell, Full de Bruijn `rho_dB`)

- **Asymmetric bounds & live 200k-pool full de Bruijn yield:** `B_f = 2^61` (batch-tested), `B_g = 2^60` (sieved), `U' = 2^53.641511` order-6 `tau`-orbits, relation target `0.82 * U'` with `+0.0666`b margin. Evaluated under **full de Bruijn `rho_dB(u, log2 N) = rho(u) - gamma * rho(u-1) / (log2(N) * ln 2)`** and live bin-by-bin Euler-product `Retention(ng)` across all bundled `200,000`-sample Voronoi FD pool bins (`c6_data.json`, `n_in = 89,151`, `frac_reg = 0.445755` vs exact 6D integral `0.4454`):
  - Active cells `V = 2^96.9622` (`SE = 0.0073`b), `V_b = 2^54.2629`, canonical `a`-points `N_a = 2^43.8650`, `b`-rows `N_5 = 2^45.3084`, average chord `L_b = 496.10`, max chord `L_max = 1010.65`, relation-weighted norms `N_f = 402.58`b, `N_g = 467.90`b, effective `r_g = 0.9325`.
- **Single-owner `(t, a_slot)` `uint16_t` score arena (`2^54.85 B << 2^64 B`) & 64-bit SWAR scan+zero (`2.15` instr/cell):**
  - Canonical `a`-points are sorted once in increasing `Q_a(a)` order; active cells along any `b`-row form the contiguous prefix `a_slot < M(t)`. Row-major `(t, a_slot)` `uint16_t` storage occupies `2^54.8504 B << 2^64 B` (`STORAGE.md`). Because `|log2 N_g,rec(a,b) - log2 N_g,rec(a,0)| < 10^-15` (`0.0` bits across 500 live pairs), `bias_table[a_slot] = 0x8000 - T(a_slot)` is precomputed once, with `+0x0030003000300030` added when `b(t) in P_2`. Pure baseline scalar AArch64 8x-unrolled 64-bit SWAR scan+zero (`KERNELS.md`, `1.91 <= 2.15` Ir/cell, verified `0/10,000` errors live in `audit_bn_cm.py`) touches strictly the `V` active cells.
- **Self-contained three-band score-only sieve (`P_min = 14`, `slack = 114`, `slack_b = 119`; ZERO prime recording!):**
  - **Band 1 `[2^14, 2^44)` (pointwise `b`-lookup + `early4_find` + unrolled-4 `uint16_t` apply):** Unrolled-4 pure-scalar `bn254_apply_8b_u16_u4` (`C_HIT_LO = 4.25` Ir/hit, `KERNELS.md`) + `min(bsearch, bn254_b1_early4_find)` integrated live over `[2^14, 2^44)` `=> rep_lo = 6.635` instr/cell.
  - **Band 2 `[2^44, 2^53)` (8-lane 64-bit SWAR `sb_u7` gate + per-SB tight `q_min_j`, `W_SB = 128`):** Each super-bucket stores a 1-byte quantized `sb_u7[q] in [0..127]` (`0` = empty, `1..127` = `127 - floor(126 * min_qa)`). An 8-phase per-row `bias_tbl[0..7]` tests up to 8 consecutive SBs in ONE 64-bit pure-scalar SWAR load (`C_SWAR_PROBE = 18.0` Ir, `probe = 1.280` instr/cell), skipping `45.2%` of occupied SBs (`P_occ = 0.4636 -> P_ent = 0.2541`, `C_SB_ENTER = 14.0 => occ_sb = 1.234` instr/cell). Exact 6D ellipsoid integration gives per-SB tight `cr2 = 1.301` (vs `1.807` row-wide) at `C_CAND = 13.0`, `C_HIT_HI = 9.0 => hi_ap = 4.822` instr/cell. Verified live (`1523/1523` hits, `0` mismatches) in `audit_bn_cm.py`.
  - **Band 3 `[2^53, 2^59)` (`sigma`-Orbit-Folded 12D Ellipsoid Fincke-Pohst Walk):**
    - Because **`sigma(w) == w` identically in `K_h`** and `sigma` is an exact isometry of `(Q_a, Q_b)`, every prime ideal `P = (p, s, rho)` of `g_rec` satisfies `(a, b) in L_P <=> (sigma^k(a), sigma^k(b)) in L_{sigma^k(P)}`.
    - Thus we compute ONE 12D LLL basis per `sigma`-orbit (`1/3` as many ideals; generous `C_LLL12 = 2^22` ops/orbit `=> bases <= 2^58.786037` units) and execute ONE 12D Fincke-Pohst enumeration on the exact 12D sieve ellipsoid `Q_b(b)/R_B^2 + kap * Q_a(a)/R_A^2 <= 1` with half-space symmetry `v <-> -v` at the root (`VR_12 = 1.000000`, `babai_init = 0.000000`, zero `a`-run overhang).
    - At each leaf, we check the `sigma`-invariant Voronoi unit-FD (`f_FD_units = 0.3891`, so `1/0.3891 = 2.5700` leaves per canonical hit, bounded conservatively by `MULT_FD_12D = 2.62`) and map passing leaves via `sigma^k` to their unique canonical `(a, b)` hit (emitting all 3 `sigma`-images on the `< 10^-9` measure-zero `sigma`-fixed subspace).
    - Under the parameter-free 12D Gaussian volume law at dimensionless density `h_12(p) = 3 * V / (0.3891 * p) in [2^40.9, 2^46.9]`, live integration across `[2^53, 2^59)` gives `fp12_nodes_int = 0.129` internal nodes/leaf (`C_FP12_NODE = 26.0`, `C_FP12_LEAF_FD = 36.0 => fp12_walk = 11.057` instr/cell; 8B `log_p`-grouped bucket `C_BKT_B3 = 16.0 => bkt_top = 1.716` instr/cell).
    - **Live Verification (`audit_bn_cm.py [walk12_sigma_folded]`):** Audited across 4 real BN254 `C6` prime `sigma`-orbits (`p in {379, 419, 569, 619}`): **`602 / 602` canonical hits recovered, `0` missed, `0` duplicates, `leaves/can_hit = 2.3688 <= 2.62`**.
- **Score-Only Stage B Window Lemma & Batch Cofactorization (`cof = 0.961` instr/cell; ZERO sieve prime recording!):**
  - As derived in `KERNELS.md` §2 (matching accepted `bn15r` `c5-gt`), Bands 1–3 only accumulate `1/16`-bit quantized `uint16_t` scores `s = sum floor(16 log2 p)` (zero sieve prime recording). On each Stage-A (`sf`) survivor (`L - s <= 16 * 114`, where `L = floor(16 log2 |N_g|)`), Stage B computes exact norms (`C_NORMS = 27,000`) and strips the **exact `2^14`-smooth part `m` of `N_g`** via `417` 64-bit prime-product chunks covering all `1,900` primes `< 2^14` (`~18,344 << C_SB = 80,000` scalar ops; verified `0` errors live in `audit_bn_cm.py`).
  - **Window Lemma:** Let `L_m = floor(16 log2 m)` and `O = L - s - L_m`. Since at most `c = floor(ceil(n_g) / 14)` scored ideals in `[2^14, 2^59)` each round `16 log2 p` down by `< 1`, `O` lies within `(-1, c + 1]` of `16` times the true outside bit-mass (`log2 |N_g| - log2 m - sum_{p in [2^14,2^59)} log2 p`). Because `S_A = 114 < S_B = 119`, `O <= 16 * 119` alone would pass all Stage-A survivors; the **`~18x` Stage-B rejection (`sf = 1.11e-6 -> s0 = 6.23e-8`)** comes from enforcing **BOTH**:
    1. `O <= 16 * S_B = 1904`, **AND**
    2. `O` lies in one of the **disjoint large-prime windows** for `k in {0, 1, 2}` unsieved primes in `(2^59, 2^60]`:
       `O <= c + 1` (`k = 0`), OR `16 * 59 * k - 1 <= O <= min(16 * 60 * k + c + 1, 16 * 119)` (`k = 1, 2`) — matching `Retention.survivor_upper` identically (verified live in `audit_bn_cm.py`).
  - Only passing `s0` survivors enter `f`-side exact SRT (`gbatch_cost(61)`) and (if `f` passes, `<= q_f`) full `g`-side factorization (`C_TDIV + C_BATCH + gbatch_cost(60)`).
- **Totals:** `29.855` scalar instr/cell (times `1.0033` allowance) `=> collection <= 2^87.579177` units; `lists <= 2^81.911151` units; global RAM `2^97.40 B < 2^100 B` (and `< 2^64 B` per-owner buffered sensitivity `2^86.6856`, `STORAGE.md`).

---

## 3. Production `Z[omega]` Filtering, Fixed-Instance `H3` 4-Cap LA, & Rank Bridge

- **Complete 5-Scale Production Ladder (`2^13 -> 2^14` LIVE Python -> `2^20 -> 2^21 -> 2^22` C, `W <= 161.0`, `(61, 60, exc = 0.82)`):**
  `verify_bn254_v6.py` runs the exact hexagonal `Z[omega]` filtering pipeline live at `2^13` (with exact kernel check, `kernel_errors = 0`) and `2^14`, AND verifies monotonic decay of the aggregate mid-band weight `w_mid` (`13.71 -> 13.37 -> 11.54 -> 11.25 -> 10.96 -> 0.19`, reaching `0.0` per-side at `U' = 2^53.64`):

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `k5_exp` | `tab_r` | `w_mid` |
|---|---|---:|---:|---:|---:|---:|
| `2^13` | **LIVE Python** (`W <= 161`, `kernel_errors = 0`) | `0.004883` | `0.002686` | `14.12` | `0.7729` | `13.71` |
| `2^14` | **LIVE Python** (`W <= 161`) | `0.002319` | `0.001343` | `15.17` | `0.7743` | `13.37` |
| `2^20` | C `filtsim` (`n = 4` seeds) | `0.001271` | `0.000751` | `10.75` | `0.7482` | `11.54` |
| `2^21` | C `filtsim` (`n = 8` seeds) | `0.001031` | `0.000630` | `10.63` | `0.7436` | `11.25` |
| `2^22` | C `filtsim` (`n = 4` seeds) | `0.000923` | `0.000570` | `10.50` | `0.7411` | `10.96` |
| **Hard Caps (`H3`)** | **Deterministic Abort (`W <= 161.0`, `lazy = 0`)** | **`0.001600`** | **`0.000980`** | **`12.00`** | **`0.8500`** | **`~0.0`** |

- **Fixed-Instance Premise `H3` & Deterministic 4-Parameter Hard Caps:**
  - Before any Markowitz merge step, the filter enforces the deterministic pre-merge guard `totw_post <= 161.0 * (N_alive - 1)`. After filtering stops with `W_avg <= 161.0`, the algorithm **deterministically aborts** if ANY of `N_sw > 0.00160 * U' = 2^44.3538`, `N_l > 0.00098 * U' = 2^43.6466`, `k5 > 12.00`, or `tab_r > 0.85` occurs.
  - Because `LA`, `gen`, and 1-owner RAM (`2^59.68 B << 2^64 B`) with `lazy = 0` are strictly monotonically increasing in `(N_sw, N_l, k5, tab_r, W_ROW)`, pricing `la <= 2^84.872876` and `gen <= 2^78.666335` (`bw = (768, 35)`) **at the 4-parameter hard cap** holds deterministically on every non-aborted run.
  - Because the full-production BN254 relation matrix is a **fixed deterministic instance** (given the fixed public sieve region and public seed), `H3` is stated explicitly as a **fixed-instance exploratory premise** (supported qualitatively — not proved — by the strictly monotonic `ALL4` 3-scale C ladder `2^20 -> 2^21 -> 2^22` with `+72%..+73%` dimension and `+14%` per-row margins, `w_mid -> 0`, and super-critical peeling `e_crit = 0.676 << 0.82`, `h = 0.481`): under `H3`, the fixed production filtered matrix satisfies all four caps, so `P(cap abort | H3) = 0`.
- **Exact 1D-Kernel Rank Bridge (`RANK-BRIDGE.md` & `hx_z6.py`):**
  `RANK-BRIDGE.md` separates `(a) H-R` (1D kernel), `(b) H-X` (finite controllability & observability at charged `d, e` on sparse block-diagonal `X` — non-pad coordinates partitioned across `m` blocks, Schirokauer pad coordinates in distinct blocks `0..s_sm-1` — and dense random `Z` generating `Y = [S, B Z]`, matching production), and `(c)` conditional Schwartz-Zippel (`4 N_sw_cap / 2^64 = 4.875e-6`). `audit_bn_cm.py` executes `hx_z6.py` live on an unscreened filtered `Z[omega]` recorded-pivot system (`N = 16, D = 20, rank_A = 19`), verifying `rank(C_d) = 16`, `rank(O_e A) = 19`, one-block control `17 < 19`, and exact recovery of the planted 1D kernel `(w_true, c_true)` via projected Block-Wiedemann Hankel `H f = 0`. (Production-scale `H-R` and `H-X` remain explicitly declared fixed-instance exploratory premises supported diagnostically, not proved, by `hx_z6.py`.)

---

## 4. Live 62-Row 12D Descent (`LB = 59 <= min(Lf, Lg)`, Certified Pure `QBINS`/`BOOT_BINS`)

- **Certified 100%-Dominated C6 Envelopes (`descent_bn254_v6.py`, ZERO `PRCAP`):**
  - `QBINS = [(1/6, -1.23, -1.00), (1/6, 8.10, 2.72), (1/6, 14.96, 5.42), (1/6, 20.79, 7.70), (1/6, 27.67, 10.30), (1/6, 55.33, 18.25)]` (pure 6-bin, zero sub-bin `PRCAP`!), conservative `QBINS_LO = [(-95,-55), (-10,-55), (0,-55), (8,-55), (14,-55), (20,-55)]`, pure 5-quintile `BOOT_BINS = [(0.2,-184.45), (0.2,-181.82), (0.2,-180.07), (0.2,-178.14), (0.2,-170.68)]`, `BOOT_LO_OFF = -200.0`.
  - Audited live in `audit_bn_cm.py [descent_12d_lattices]` on 6 exact special-`q` bases (`3 f + 3 g`, `2^66..2^213`) and 2 boot bases using exact uniform FULL `B_12(4)` sampling (`13,315,520` cells including `|c_i| in {3, 4}`), exact Bareiss dets, exact `6x6` big-int norms (`max_err = 3.9e-11`b), and exact `q | N_f` & `q | N_g` divisibility: paired joint smoothness beats `QBINS` on all 6 bases (`[2.01x, 2.78x, 3.03x, 2.70x, 2.98x, 2.35x]`) and across `R in {16, 64, 128}` (`min = 1.12x >= 1.0`), 6/6 `QBINS_LO` pass, and all 5 boot quintiles pass (`boot_smooth_ratio_L110 = 1.13x >= 1.0`).
- **Node Failure Bound `min_k (P(M < k) + psi^k)`, `THRESHOLD-REPAIR` Counterexample, & Multi-Factor ECM:**
  - Uses `phi_node = max(0.01, min_k (P(M < k) + psi^k)) <= 0.04` (explicit conditional hazard `H5(ii)`, never `PGF(psi)`), multi-factor stored-prefix ECM isolation (`SEPARATION = 0.5`), and history-conditional orbit collision `C_peak = 4.0` (`M_59 = li(2^59)/6`), yielding `two_logs <= 2^83.472242` units, `node_abort = 1.945e-3`, `repeat_ideal = 9.714e-4` (since `B_g = 2^60 > 2^59`, `LB = 59` is a strict conservative lower bound).
  - **Live `THRESHOLD-REPAIR` Counterexample Witness (`audit_bn_cm.py [h5_adaptive_and_threshold]`):** For `M = 1` w.p. `0.10` (sole candidate fails) and `M = 100` w.p. `0.90` (immediate success), true failure is `0.10` while `E[t^M] <= (1 + 30(1-t))^-1` holds on `[0, 0.8]`. Invalid `PGF(0.10) = 0.035714 < 0.10` understates failure, whereas `descent_bn254_v6.failure_bound(60, 1, 0.10)` returns **`0.132258 >= 0.10`** (`k = 1`, `threshold_counterexample.ok = True`).
- **Real Reduced-Tower `C6` Recursive Descent Campaign (`ADAPTIVE`, `c6_adaptive.json` replayed LIVE in `audit_bn_cm.py`):**
  To test `H5(i)-(iii)` directly on genuine recursive 12D descent trees of `C6-m184_w216-rec` (reduced pair `f_4(x) = x^4 + w x^2 + w^2`, `g_2(x) = 4 w x^2 - 1` over `O_h`, `B_0 = 32,768`, exact `O_h`-ideal primitivity `rank_mod_p_6x12 == 6`, exact degree-1 root matching, exact norm dedup, **ZERO censoring**), `c6_adaptive.json` records **`1,360` roots across 4 settings (`f` and `g` sides)** yielding **`959` after-failure attempts** across **`519` root clusters**:

| Setting | Roots OK | `mean(M)` | `P_emp(M<1) <= T_1` | After Fail | `OE_fail` | `ratio_1` | `OE_work` | Shared / Colls |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `f19_d1` (`f`, `19b`) | `339/400` | `3.73` | `0.0375 <= 0.3489` | `113/217` | `0.729` | `0.843` | `0.998` | `0` / `0` |
| `f20_d1` (`f`, `20b`) | `286/320` | `4.94` | `0.0031 <= 0.2882` | `169/270` | `0.847` | `1.093` | `1.004` | `0` / `0` |
| `f21_d1` (`f`, `21b`) | `301/320` | `6.53` | `0.0000 <= 0.2345` | `238/389` | `0.823` | `0.922` | `1.067` | `1` / `3` |
| `g19_d1` (`g`, `19b`) | `223/320` | `1.97` | `0.1469 <= 0.5035` | `40/83` | `0.832` | `1.008` | `0.954` | `0` / `0` |
| **POOLED (`959` after)** | **`1149/1360`** | — | **`12/12` `M<1,2,3` OK** | — | **`0.8095`** `[0.7685,0.8434]` | **`0.9607`** | **`1.0307`** `[0.9970,1.0632]` | **`1` / `3`** |

  `audit_bn_cm.py` replays all `959` raw after-failure rows live: `12/12` `M`-tails stay strictly below `T_k`, pooled `OE_fail = 0.8095` (95% CI upper `0.8434 < 1.00`) and `ratio_1 = 0.9607 <= 1.00` show zero post-failure degradation, `OE_work = 1.0307` (95% CI upper `1.0632 < CHILD_X = 1.25`), and shared child ideals across retries occur in only `1/959` (`0.104%`) with `3` total collisions across `1,360` roots. (Reduced-tower diagnostics support the plausibility of declared heuristic `H5` and are not substituted for production probabilities.)

---

## 5. Itemized 16-Row Ledger, Strict Additive Union Bound (`FIX-U`), & Exact Rational Run Gate

| Ledger Row | Computed (`log2` u) | Declared Upper Bound (`log2` u) |
|---|---:|---:|
| `collection` | `87.579176` | `87.579177` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `84.872875` | `84.872876` |
| `gen` | `78.666334` | `78.666335` |
| `lists` | `81.911150` | `81.911151` |
| `two_logs` | `83.472241` | `83.472242` |
| `canon` | `72.377324` | `72.377325` |
| `sm` | `73.833931` | `73.833932` |
| `misc` | `69.550634` | `69.550635` |
| `fb` | `65.038073` | `65.038074` |
| `bases` | `58.786036` | `58.786037` |
| `split` | `81.841117` | `81.841118` |
| `gtree` | `61.251028` | `61.251029` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| `output` | `11.761902` | `11.761903` |
| **Expected opened run (`exp_run`)** | **`87.903225260`** | **`87.903225261`** |
| **Gated claim (`E[W]`)** | **`86.548991569`** | **`86.548991570`** |

**Strict additive union bound on opened-run success (`FIX-U`, under declared heuristics `H1-H5`, `H-R`, `H-X`):**
- Main-phase conditional Schwartz-Zippel (`4 N_sw_cap / 2^64`): `4.875e-06`
- Filtering 4-parameter hard-cap abort under fixed-instance `H3`: `0.0`
- Split `2^112` budget abort: `8.342e-10`
- Log-phase `2^112` meter Markov bound: `2.584e-09`
- Descent node-budget abort (`H5`): `1.945e-03`
- Descent history-conditional repeated-ideal collision (`H5(iii)`, `C_peak = 4.0`): `9.714e-04`
- Group-order & numerical (`2/r + 2^-48`): `3.553e-15`
- **Opened-run success:** `1 - sum(failures) = 0.99707922 >= 0.997079`.

---

## 6. Declared Exploratory Heuristics

- **H1 / H1':** Standard NFS smoothness under full de Bruijn `rho_dB` with Murphy `alpha`, exact 2-adic split, live norm-histogram-weighted Euler-product score-only retention/survivor laws on `N_f`, `N_g`, and the 12D Gaussian volume heuristic for Band-3 `sigma`-orbit-folded LLL bases.
- **H2:** Uniform sampling over the certified 6-facet Voronoi fundamental domain (`200,000`-sample bundled pool quadrature, `SE = 0.0073`b).
- **H3 (Fixed-Instance Deterministic 4-Cap Conformance):** On the fixed deterministic public-seed production relation matrix (`U' = 2^53.64`, `exc = 0.82`), after filtering halts under the deterministic pre-merge guard `W_avg <= 161.0`, the filtered `Z[omega]` matrix satisfies all four deterministic hard caps (`N_sw <= 0.00160 * U'`, `N_l <= 0.00098 * U'`, `k5 <= 12.00`, `tab_r <= 0.85`). Supported qualitatively by the strictly monotonic `ALL4` 3-scale C ladder (`2^20 -> 2^21 -> 2^22`), live `2^13 -> 2^14` + `w_mid -> 0` monotone check, and `e_crit = 0.676 << 0.82`. LA, generator, and 1-owner RAM are priced deterministically at the 4-cap with `lazy = 0`.
- **H-R / H-X (Separated Rank Bridge, `RANK-BRIDGE.md`):** `(a)` `A = [M_prod | S_prod]` has exact 1D kernel `(w_true, c_true)` with `c_true != 0`, and `(b)` finite controllability `rank(C_d) = N_l` and observability `rank(O_e A) = rank(A)` at charged `d, e` on sparse `X`, yielding conditional Schwartz-Zippel bound `4 N_sw_cap / 2^64` (verified live in `hx_z6.py`).
- **H5 / H-EA / H-FFT / H-SEP (12D Descent Norm, Tail, Conditional Hazard, & Collision Premise):** At `LB = 59`, `(0)` pure `QBINS`/`BOOT_BINS` lower-bound paired smoothness and `QBINS_LO`/`BOOT_LO_OFF` upper-bound survivors; `(i)` usable candidate counts `M` obey the overdispersed Negative Binomial lower-tail bound (`mean = mu/2`, `r_sh in {0.5, 1.0}`); `(ii)` explicit conditional candidate-failure hazard `P(E_j | E_1, ..., E_{j-1}, M >= j) <= psi` giving `phi <= min_k (P(M < k) + psi^k)`; `(iii)` history-conditional orbit-collision hazard `<= 4.0 / (li(2^59)/6)`, together with `H-FFT` scaled remainder trees and `H-EA`/`H-SEP` stored-prefix multi-factor ECM isolation.
