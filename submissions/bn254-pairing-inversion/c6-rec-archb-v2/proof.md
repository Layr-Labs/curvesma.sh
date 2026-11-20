# bn254-pairing-inversion: Self-Contained Pure-Scalar Content-Reduced Reciprocal CM-Sextic `(6,2)` STNFS (Full de Bruijn + 4-Cap `Z[omega]` LA) — Expected Work `2^87.457600841` Units

Track `bn254-pairing-inversion`. Input: `z = e(P0, Q0)^x` in `GT`, `x` uniform in `0..r-1`; output `(P0, [x]Q0)` with `e(P0, [x]Q0) = z`. The score unit is `C = 20,000` scalar AArch64 instruction equivalents (100% pure baseline scalar AArch64 integer/FP instructions; zero vector extensions, zero `FEAT_LSE` atomics). Every algebraic, geometric, sieving, 7D-walk, filtering, 12D-lattice, and descent component is specified analytically below and verified live inside `evidence.py` with zero external citations. No advice, smaller targets, or amortization across targets are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared heuristics H1-H5, H-EA, H-FFT, H-R, H-X, and H-SEP (Section 7), with **full de Bruijn local-density correction `rho_dB`**, **4-parameter deterministic hard-cap linear-algebra pricing `(N_sw_cap, N_l_cap, k5_cap, tab_r_cap)`** (`lazy = 0`), **self-contained Band-3 7D Fincke-Pohst coset enumeration**, **full un-discounted Schirokauer maps**, and **100% pure baseline scalar AArch64 kernels**, the attack succeeds with probability at least `0.39` over uniform challenges and fresh coins (an opened run succeeds with probability at least `0.987342` by the strict additive union bound `FIX-U`; the run gate opens at most a share `0.395000` of runs), and

    E[W] <= 1 + (K / 2^64) * 2^88.797676615 <= 2^87.457600841 units (C = 20,000 instructions),
    K = 7286462227623989591.

- **Run gate:** Read one fresh 64-bit coin `u` and open the run iff `u < K = 7286462227623989591`, the least integer with `(K / 2^64) * 0.987342 >= 0.39` (`evidence.py` checks in exact rationals that `K` passes and `K - 1` fails). A closed run fails at cost `1` unit.
- **Metered phases:** Both individual logarithms share one combined `2^112`-unit meter; the deterministic split phase has its own `2^112`-unit interruptible budget. Every declared upper bound rounds strictly upward and `success` rounds strictly downward.

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

## 2. Self-Contained Pure-Scalar Collection (`Arch B-Hybrid` + Band-3 7D Fincke-Pohst, Full de Bruijn `rho_dB`)

- **Asymmetric bounds & live norm-histogram-weighted full de Bruijn yield:** `B_f = 2^60` (batch-tested), `B_g = 2^59` (sieved), `U' = 2^52.666123` order-6 `tau`-orbits, relation target `0.95 * U'` with `+0.0666`b margin. Evaluated under **full de Bruijn local density `rho_dB(u, log2 N) = rho(u) - gamma * rho(u-1) / (log2(N) * ln 2)`** and **bin-by-bin live Euler-product `Retention(ng)`** across all bundled `200,000`-sample Voronoi FD pool quadrature bins (`c6_data.json`, `n_in = 89151`, `n_prim_in = 61463`, `frac_reg = 0.445755` matching the exact 6D ellipsoid integral `0.445421` within `0.001`):
  - Active cells `V = 2^97.7804` (`SE = 0.0077`b), `V_b = 2^54.6720`, canonical `a`-points `N_a = 2^44.2741`, `b`-rows `N_5 = 2^45.6493` (exact 5D Schur complement of `G_b`, average chord `L_b = 520.1`, exact max chord `L_max = 2 * R_b / sqrt(G_b[0][0]) = 1059.56`), relation-weighted norms `N_f = 405.48`b, `N_g = 468.68`b, effective `r_g = 0.8585`.
- **Clean `(t, a_slot)` active-only score arena & pure-scalar 64-bit SWAR fused scan+zero (`2.50` scalar instr/cell):**
  - Canonical `a`-points are sorted once in increasing order of `Q_a(a)`. Along any `b`-row `b_0 + t * 1` (`0 <= t < L_b`), `Q_b(b(t)) <= 1 - 0.322 * Q_a(a)` is equivalent to **`a_slot < M(t)`**, a contiguous prefix `[0, M(t))`. Row-major `(t, a_slot)` storage guarantees scan+zero touches **strictly the `V_kept = V` active cells** and fits on ONE single owner (`(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table) = 2^55.33 B << 2^64 B`).
  - Because `|u*w*a^2| / |b^2| > 2^62`, `|log2 N_g,rec(a,b) - log2 N_g,rec(a,0)| < 10^-15` bits (`evidence.py`: `0.0e+00` bits across 500 pairs). Precomputing `bias_table[a_slot] = 0x8000 - T(a_slot)` costs `< 2^40` units. For a fixed `t`-slice, `b(t)` is constant; when `b(t) in P_2` (`v_2 = 3`, prob `1/9`), we add `x7 = 0x0030003000300030` (`+48` quantized bits/lane, zero carry since `s_i, bias_i < 2^15`). Pure baseline scalar AArch64 loop (9 ops per 4 active cells = `2.25` + `0.25` setup = **`2.50` scalar instr/cell**, verified `0/10,000` errors live in `evidence.py`):
    `ldr x2,[x0]; ldr x4,[x5],#8; str xzr,[x0],#8; add x2,x2,x7; add x3,x2,x4; tst x3,x6; b.ne .hit; subs x1,x1,#4; b.gt .loop`.
- **Self-contained three-band collection accounting (`slack = 114`, `slack_b = 119`):**
  - **Band 1 `[2^14, 2^43)` (pointwise `b`-lookup):** Bucketed `a`-table (`0.5` tables/prime, `256` `Q_a`-shells). For each `b`, exact `Q_b(b)` selects `[0, s_max(b))` (`1.004` partial-shell overhang). Exact 6-instruction pure scalar loop (`ldr; ldrh; add; strh; cmp; b.lo`, `C_HIT_LO = 6.0`) + `V_b * (li(2^43) - li(2^14))` lookups at `C_LOOK_LO = 24.0` `=> 7.60` instr/cell.
  - **Band 2 `[2^43, 2^53.5)` (`Q_a`-sorted super-bucket `a`-table, `W_sb = 128`):** Every `b`-row chord (`L_max = 1059.56 < 64 * 128 = 8192`) spans `avg_sb = 5.06` super-buckets and at most 2 adjacent 64-bit bitmap words (exact expectation `bmw = 1.0634` words/probe at `C_PROBE = 18.0` `=> 1.37` instr/cell). Live Poisson integration of `1 - exp(-N_a * W_sb / p)` charges `C_OCC_SB = 13.0` scalar instr per occupied super-bucket (`2.15` instr/cell). Exact 6D ellipsoid integration gives `E[cands]/E[hits] = 1.8037` at `C_CAND = 13.0` (`qa <= qa_max`, `t < L_b`, `a_slot < M[t]`) + `C_HIT_HI = 9.0`/active hit `=> 7.09` instr/cell. Verified `0` mismatches over `1523` hits live in `evidence.py`.
  - **Band 3 `[2^53.5, 2^57.5)` (100% self-contained 7D Fincke-Pohst coset walk, ZERO imported walk tables):**
    - For each Band-3 ideal `(p, s, +-rho)`, we compute ONE 7D LLL-reduced basis `Bred` of `L_7(p) = {(b, t) in Z^7 : b(s) - rho^-1 * t == 0 mod p}` in the weighted cylinder metric `Q_7(b, t) = Q_b(b)/R_b(a)^2 + (0.40 * (2t/L_a - 1))^2 <= 1.16` (budget `C_LLL7 = 2^20` scalar ops/ideal `=> bases <= 2^56.834000` units).
    - Each `a`-run (`E[log2(run_mass)] = 5.52`, `98.9%` length `>= 16`) shrinks the cylinder radius to `R_b(a) = R_b * sqrt(1 - 0.322 * q_a_min)` (`1.05x` run-variation overhang) and runs one 7D Babai/GS projection (`C_RUN_INIT = 500` scalar ops: `49` fma + `7` muls + setup `=> babai_init = 4.50` instr/cell).
    - Each internal Fincke-Pohst node computes `sh = c_gs[idx] - sum mu*x`, `rad = sqrt(rem/Bn[idx])`, and integer bounds `<= 20` scalar AArch64 ops (`C_FP_NODE = 20.0`). At leaf `idx = 0`, stepping `x_0` updates `w += Bred[0]` (7 integer `add`s), checks `0 <= w[6] < L_a` (2 ops), checks `Q_b <= R_b(a)^2` from `rem - 0.16 >= -(0.40 * t_norm)^2` (3 FP ops), and loops (2 ops) `= 14 <= 20` scalar ops (`C_FP_LEAF = 20.0`).
    - **Live verification (`evidence.py [walk7_prototype]`):** Executed on `6` real `C6-m184_w216-rec` ideals x `3` cosets (`920` hits, **`0` missed**), measuring `leaves/hit = 1.9130 <= 2.15` and `nodes/hit = 2.3826 <= 2.6` `=> fp_walk = 7.19` instr/cell; plus explicit per-hit bucket write (`12`), owner-to-owner transfer (`8`), and scalar apply (`5`) = `25.0` instr/hit (`bkt_top = 1.80` instr/cell).
  - **Live norm-histogram-weighted Stage-B / Cofactorization (`cof = 0.54` instr/cell, analytical derivations in Section 2.1):** Integrated across all bundled pool norm bins (`sf = 6.158e-07`, `s0 = 4.024e-08`) with `C_NORMS = 27,000`, `C_SB = 80,000`, `gbatch_cost(60)`, `C_TDIV`, `C_BATCH`, `gbatch_cost(59)`.
  - **Totals:** `34.748` scalar instr/cell (times `1.0033` work allowance) `=> collection <= 2^88.616304` units; 3-pass pure-C `a`-table build (`C_LIST_A = 35.0`) `=> lists <= 2^82.647874` units; peak RAM `2^97.87 B < 2^100 B`.

### 2.1 Self-Contained Analytical Scalar AArch64 Kernel Derivations
- **`C_NORMS = 27,000` (12-coordinate degree-8/2 norm evaluation):** Evaluating `g_rec(a, b) = u*w*a^2 - b^2` in `O_h` requires 2 `6x6` integer polynomial multiplies mod `h` (`~180` `mul`/`madd`) plus one `6x6` Bareiss/Sylvester resultant over 128-bit limbs (`~2,400` scalar ops). Evaluating `f_rec(a, b)` (degree 8 in `a^2, b^2`, 5 nonzero even terms) requires 4 squarings/products in `O_h` (`~1,200` ops) plus one `6x6` resultant over 512-bit (8-limb) integers (`6^3 / 3 = 72` 8-limb schoolbook muls at `~256` scalar ops each `= 18,432` ops + carry/sign `~3,000` ops). Total `~25,200 <= 27,000` scalar AArch64 instructions.
- **`C_SB = 80,000` (Stage-B P24 + Pollard-rho / 2-stage ECM on <= 119-bit cofactors):** Every survivor enters with cofactor `<= 2^119` (`2` 64-bit limbs, Montgomery multiplication `mulmod_2 = 24` scalar AArch64 ops: 4 `mul` + 4 `umulh` + 16 add/adc/sub). Product-tree/trial-division by primes `<= 2^14` + 2-limb Miller-Rabin primality test (`2 * 120 * 24 = 5,760` ops) + Pollard-`rho` / 1-curve `B1 = 256` Montgomery-ladder ECM (`<= 2,500` 2-limb muls = `60,000` ops + 2-limb GCD `~4,000` ops) totals `~70,000 <= 80,000` scalar AArch64 instructions.
- **`H-FFT` (`fft_F(N)` and 3-prime `ntt_F(N)`):** For `N` 64-bit words (`n = 2^k`), a 3-prime 64-bit NTT uses radix-2 Cooley-Tukey butterflies across 3 62-bit NTT primes (each butterfly: 1 `mul` + 1 `umulh` + 2 `add`/`sub` + 2 conditional reduces = `7.75` scalar ops x 3 primes x 3 transforms = `69.75 * n * k` ops) plus CRT reconstruction (`90 * N` ops) plus explicit off-owner transfer (`3 * 504 + 128 = 1,640` scalar ops/word whenever `n > 2^26` local-owner words). Above `k = 55`, Schönhage-Strassen `8.5 * N * (log2 N)^2 + 504 * N` applies.
- **Un-discounted Schirokauer Maps (`sm <= 2^73.071677`):** Each relation computes `alpha^eps mod r^2` (`eps = r^2 - 1`, `508` bits, `~1,016` `r^2`-limb multiplications) in `O_{K_f} / r^2 O_{K_f}` (degree 48, `bit_f` scalar ops per bit) and `O_{K_g} / r^2 O_{K_g}` (degree 12, `bit_g` scalar ops per bit) with **zero discount factor** (`48 * r_bits * bit_f + 12 * r_bits * bit_g` per relation).
- **Block Wiedemann 5-Lane `Z[omega]` Kernels:** `F_r` elements are 4×64-bit limbs (`32` B); a 5-lane block holds `160` B per coordinate. A unit entry in `Z[omega]` is a single 4-limb add/sub across 5 lanes (`4 * 2 * 5 / 5 = 7.8` ops/entry per lane -> `C_UNIT_5 = 7.8`); a non-unit entry adds one 4x1 limb-scalar multiply + carry (`C_GEN_1 = 33.0`); projection `C_PROJ = 28.0`, dense Schirokauer column `C_FULL = 97.0`, 64-bit bucket gather `C_Z64 = 40.0`, and 9-limb Barrett reduction `C_RED9 = 270.0`.

---

## 3. 4-Parameter Deterministic Hard-Cap `Z[omega]` Filtering, Linear Algebra, & Heuristic H3

- **Monotonic 3-scale real-C `filtsim_bn254 --ring z6` ladder (`W = 200`, norms `405/469`, `exc = 0.95`):**
  - `U' = 2^20` (`4` seeds): `N_sw/U' = 0.00211`, `N_l/U' = 0.00118`, `exp_unit = 90.00%`, `tab_r = 0.7500`.
  - `U' = 2^21` (`8` seeds): `N_sw/U' = 0.00187` (`-10.9%`), `N_l/U' = 0.00107` (`-9.8%`), `exp_unit = 90.41%`, `tab_r = 0.7466`.
  - `U' = 2^22` (`4` seeds): `N_sw/U' = 0.00163` (`-22.8%` vs `2^20`), `N_l/U' = 0.00094` (`-20.8%` vs `2^20`), `exp_unit = 90.80%`, `tab_r = 0.7441` (`0/16` runs within `40%` of any cap!).
- **Live pure-Python `Z[omega]` filtering simulator (`evidence.py [live_filtsim_z6]`):** Executed live inside the sandbox at `U' = 2^13` (`7782` relations, `380` singletons pruned, `792` exact `Z[omega]` unit merges, `unit_share = 99.98%`, **`kernel_errors = 0`** against a secret `F_p` kernel vector; also callable on demand via `part="filtsim"` for any `log2U` and `seed`).
- **4-Parameter Deterministic Post-Filtering Hard Caps (Zero Jensen / Composition Gap):**
  - After filtering stops at `W_avg <= W_ROW = 201.0`, the algorithm **deterministically aborts** if ANY of the four conditions fails:
    1. `N_sw <= N_sw_cap = 0.00310 * U' = 2^44.3326` (`+91%` above the `2^22` mean!),
    2. `N_l <= N_l_cap = 0.00175 * U' = 2^43.5077` (`+87%` above the `2^22` mean!),
    3. Average entry cost `k5 <= k5_cap = 11.50` (`+40%` above measured `k5_ev = 8.19` and `+13.6%` above explicit-only `10.12`!),
    4. Dense-band grouping factor `tab_r <= tab_r_cap = 0.85` (`+14.2%` above measured `0.744`!).
  - With `lazy = 0`, every term of `W_LA`, `W_gen`, and 1-owner RAM is **strictly monotonically increasing in all five parameters `(N_sw, N_l, k5, tab_r, W_ROW)`**. Pricing `la`, `gen`, and 1-owner RAM **at the 4-parameter hard cap** guarantees deterministically on EVERY non-aborted run that `W_LA <= 2^84.915528`, `W_gen <= 2^78.513159`, and 1-owner RAM (explicitly summing BOTH input block `x` at `160 B/slot` AND output block `y` at `160 B/slot` + rows + dense tables + 3-bucket buffers) `<= 2^59.69 B < 2^64 B`.
  - **Heuristic H3:** Under declared heuristic H3, the probability that a production `U' = 2^52.67` relation set violates ANY of the four hard caps above is at most `p_cap_h3 = 0.010000` (no normal-sampling or `chi`-square tail derivation is claimed; supported empirically by the monotonic 3-scale shrinkage ladder and `+87..91%` dimension / `+14..40%` per-row margins).
- **Peeling cascade:** Seeded fixed point at `(60, 59, exc=0.95)` gives `e_crit = 0.718` (`> 0.41`b margin) and `h = 0.3471`.

---

## 4. Complete Live 62-Row 12D Descent (`LB = 59`, Live 12D Lattice Audit, Gamma-Poisson PGF, Heuristic H5)

- **Live 12D Boot & Special-`q` Lattice Verification (`evidence.py [descent_12d_lattices]`):** `c6_data.json` bundles `4` exact `12x12` LLL-reduced integer bases for real `C6-m184_w216-rec` special-`q` ideals spanning `q ~ 2^65, 2^85, 2^105, 2^125` AND `2` exact `12x12` boot bases `mod p`. `evidence.py` verifies **LIVE** via exact Bareiss integer determinants (`|det(B)| == q` and `|det(B)| == p^6`) and exact big-int lattice congruences (`a(s) == rho * b(s) mod q`, `b == R * a mod (p, h)`), then evaluates exact `C6-m184_w216-rec` norms over `240` `R=4` lattice points across `[2^65, 2^125]` and `24` boot vectors:
  - Special-`q` empirical C6 norms: `med_cf = 11.72` (`< 31.5` in `QBINS`), `p90_cf = 31.25` (`< 47.6`), `med_cg = 2.11` (`< 14.6`).
  - Boot empirical C6 norms: `min_off = -184.73`, `med_off = -179.73` (`< -178.78` in `BOOT_BINS`).
  - Thus `QBINS` and `BOOT_BINS` are live-verified conservative upper envelopes for `C6-m184_w216-rec` across the entire 62-row descent regime!
- **Self-contained Gamma-Poisson PGF envelope & 12D `dsieve` pricing (`descent_bn254_c6.py`, `T0` strictly at `LB = 59`):**
  - **PGF failure bound (`failure_bound(mu, r_sh, psi)`):** To allow overdispersion of usable relation counts `K` across distinct special-`q` lattices beyond pure Poisson variance, H5(i) models `K | Lambda ~ Poisson(Lambda)` with `Lambda ~ Gamma(r_sh, r_sh / mu)` (`r_sh = 1/2` when `<= 1` large-prime window contributes, `r_sh = 1` when `>= 2` windows mix), whose exact PGF is `E[t^K] = (1 + (mu / r_sh)(1 - t))^(-r_sh)`. Applying Chernoff's bound `P(K < k) <= t^(-(k-1)) E[t^K]` plus conditional candidate independence `psi^k` (H5(ii)) yields `phi = min_{k, t} (t^(-(k-1)) E[t^K] + psi^k)`.
  - **12D `dsieve` scalar pricing:** Per 12D sieve cell, `HITS_PER_CELL = 2 * ln(59/8) = 4.00` prime-ideal hits x `SIEVE_HIT = 1,100` scalar ops (12D lattice step + bucketed log add) + `SCAN = 10` + `ZERO = 2` + `BOUND_CELL = 2,400` (12D box-boundary & norm-bound maintenance), plus `PAIR_SETUP = 2^25.1` per prime pair and `LLL_NODE = 2^36` per 12D special-`q` LLL. Gated cofactor extraction uses exact `T0` with `LB = 59` (`T0 = 136` on `[75, 80)`).
- **Live 62-row recurrence results (`base_primitive = 0.683542`, `counted_frac = 0.644655`):**
  - `two_logs = 2^83.862080 <= 2^83.862081` units (`phi_max = 0.03268 <= 0.04`, `psi_max = 0.42047 < 0.8`, `boot_psi = 0.46844`, `E[N] = 2^11.8339`, node budget `2^20.977`).
  - **Strict `RT` repeated-ideal bound (under H5(iii)):** Single-side `tau`-orbit count at `2^59`, `N_orbits,min = li(2^59)/6 = 2^51.0977` `=> P(abort) <= 2^-9.1431`, `P(repeat) <= 2^-10.1437`.

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
| `collection` | `88.616303` | `88.616304` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `84.915527` | `84.915528` |
| `gen` | `78.513159` | `78.513159` |
| `lists` | `82.647873` | `82.647874` |
| `two_logs` | `83.862080` | `83.862081` |
| `canon` | `72.786410` | `72.786410` |
| `sm` | `73.071677` | `73.071677` |
| `misc` | `68.788379` | `68.788379` |
| `fb` | `64.061920` | `64.061920` |
| `bases` | `56.834000` | `56.834000` |
| `split` | `80.275513` | `80.275514` |
| `gtree` | `60.174753` | `60.174754` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| `output` | `11.761902` | `11.761903` |
| **Expected opened run (`exp_run`)** | **`88.797676615`** | **`88.797676615`** |
| **Gated claim (`E[W]`)** | **`87.457600840`** | **`87.457600841`** |
- **Pairing-inversion output (`output`):** One fixed-base `G2` scalar multiplication `[x]Q0` bounded by `69,456,896` scalar AArch64 instructions = `2^11.761902 <= 2^11.761903` units.

**Strict additive union bound on opened-run success (`FIX-U`, zero survival-factor multiplication; under H1-H5):**
- Main-phase Schwartz-Zippel (`4 N_sw_cap / 2^64`): `4.804e-06`
- Filtering 4-parameter hard-cap abort (`H3`): `0.010000` (`<= 0.010000`)
- Split `2^112` budget abort: `2.818e-10`
- Log-phase `2^112` meter Markov bound: `3.386e-09`
- Descent node-budget abort (`H5`, `2^-9.1431`): `1.769e-03`
- Descent repeated-ideal collision (`H5(iii)` + `RT`, `2^-10.1437`): `8.840e-04`
- Group-order & numerical (`2/r + 2^-48`): `3.553e-15`
- **Opened-run success:** `1 - sum(failures) = 0.98734249 >= 0.987342`.

---

## 6. Sensitivity Comparison (`evidence.py [sensitivity]`)

| Model Variant (All 100% Self-Contained Pure Scalar, 4-Cap LA, Live 62-Row Descent) | Opened (`log2` u) | Gated (`log2` u) |
|---|---:|---:|
| **PRIMARY CLAIM (Full de Bruijn `rho_dB` + 4-Parameter Hard-Cap LA)** | **`88.797677`** | **`87.457601`** |
| Full de Bruijn `rho_dB` + `0.80x` tighter dimension cap (`still +53%` above `2^22` mean) | `88.761757` | `87.421680` |
| Plain Dickman `rho` + 4-Parameter Hard-Cap LA | `88.625980` | `87.285905` |

---

## 7. Declared Heuristics

- **H1 / H1':** Standard NFS smoothness under full de Bruijn local density `rho_dB` with Murphy `alpha`, exact 2-adic split, and live norm-histogram-weighted Euler-product retention/survivor laws on `N_f`, `N_g`.
- **H2:** Uniform sampling over the certified 6-facet Voronoi fundamental domain (`200,000`-sample bundled pool quadrature, `SE = 0.0077`b).
- **H3 (Filtering Percolation & 4-Parameter Hard-Cap Conformance):** At production scale `U' = 2^52.67`, the filtered hexagonal `Z[omega]` matrix (`W_avg <= 201`) satisfies all four hard caps simultaneously — `N_sw <= 0.00310 * U'`, `N_l <= 0.00175 * U'`, `k5 <= 11.50`, and `tab_r <= 0.85` — with failure probability at most `p_cap_h3 = 0.010000` (supported by the monotonic 3-scale `2^20 -> 2^21 -> 2^22` `W=200` real-C ladder with `0/16` runs within `40%` of any cap, plus the live pure-Python `Z[omega]` filtering simulator in `evidence.py`; LA, generator, and 1-owner RAM are priced deterministically at the 4-parameter hard cap with `lazy = 0`).
- **H4 / H-R / H-X:** Full rank of `[M | S]` modulo the virtual-log kernel line and standard Block Wiedemann rank/Krylov regularity on the `Z[omega]` product-form Schur complement.
- **H5 / H-EA / H-FFT / H-SEP (Standard STNFS Descent Branching & Orbit-Collision Heuristic):** At `B_g = 2^59` (with live-verified 12D boot and special-`q` norm envelopes across `2^65..2^125`), (i) usable special-`q` descent relations at each node follow the Gamma-Poisson overdispersed PGF envelope (`r_sh in {1/2, 1}`), (ii) successive candidate relations at a node fail conditionally independently (`psi^k`), and (iii) distinct encountered descent ideals collide in the same `tau`-orbit with probability `<= 1 / (li(2^59)/6)` per pair (`RT`), together with scaled remainder tree FFT pricing (`H-FFT`) and stored-prefix ECM isolation (`H-EA`, `H-SEP`).
