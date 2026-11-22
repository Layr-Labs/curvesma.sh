# `recover-P` (`bls12-381-fixed-g2-invert-g1`) through two `G_T` logarithms with `(4, 3)` `Q(zeta_15)` CM-STNFS behind a run gate: expected work `2^108.599317240806` units

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity `Q in G_2` and `T in G_T`; output `P in G_1` with `e(P, Q) = T` for the exact supplied `Q` (`P = O` for `T = 1`). The score unit is one dense `F_{p^12}` multiplication, `C = 20,000` pure baseline scalar ARMv8.0-A integer instruction equivalents (`simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64). No advice, smaller targets, or amortization are used. **No prior verdict or score is a premise.**

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is a run gate: draw a 64-bit word `u`; unless `u < K = 7201071206392798289`, stop as a failure (`< 1` u). `K` is the least integer with `K * 0.99905 >= 0.39 * 2^64`. Under `H1`–`H5`, `H-X`, and `H-FFT` (section 9), the algorithm returns the exact `P` with probability `>= 0.99905` once past the gate (section 7, even at `+10%` meter sensitivity); `u` is independent of the challenge and all other coins, so cold-run success is `>= 0.39` and

    E[W] <= 1 + (K / 2^64) * (2^109.9564 + 2^66) < 2^108.599317240806  dense-Fp12 units (C = 20,000).

`2^109.9564` is the sum of section 8's rows (every-run bounds for collection and linear algebra plus expected work for the two individual logarithms `2^103.0698` and ECM relation splitting). The score `108.599317240806` is an **expected** work bound, not a worst-case bound on every run: every run is bounded pointwise behind the hard aborting meter `T_run = 2^123` u plus wrapper:

    W <= 2^123 + 2^66 < 2^123.000000000000  units.

Analytic, exploratory-v1 claim. Every algebraic identity, kernel opcode count, filtering simulation, collection pool evaluation, linear-algebra memory/work bound, 384-mesh descent recurrence, and ledger/gate bound is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib).

## 1. Reduction (Fixed `Q`: Calibrated `gamma_0`-Blinded Two-Logarithm Inversion) and Field Map

Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md`) recovers a nonzero homogeneous virtual-logarithm functional `V(u) = lambda * log_g0(u) mod r` (with unknown nonzero scalar `lambda in (Z/rZ)^*`). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {1}` (`1` final exponentiation `<= 2^11` u, ZERO third descent!) with KNOWN nonzero virtual log `V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0` (`gcd(r, (p^12-1)/r) == 1` asserted LIVE in `audit_cm.py`). Given the challenge `(Q, T)` with `Q = [a] Q_0 in G_2 \ {O}` (`a in (Z/rZ)^*`) and `T = e(P, Q) = g_0^{a b} in G_T`:
  0. Run gate (`Claim`): draw 64-bit `u`; if `u >= K_gate`, stop as failure (`< 1` u).
  1. Check `Q` (on twist, `[r]Q = O`, `Q != O`) and `T` (`T^r = 1`). If `T == 1`, return `P = O`.
  2. Compute `h := e(P_0, Q) = g_0^a in G_T` (`1` pairing, `<= 2^11` u) and map `h, T` into the `Q(zeta_15)` tower (`field-map.md`, wrapper `<= 2^66` u, failure `< 2^-48`).
  3. For each target `z in {h, T}`: each fresh boot ball draws uniform `e in [0, r-1]`, blinds `z' := z * gamma_0^e` (exact uniform in `G_T`), runs individual-log descent on `z'` to obtain `V(z') = V(z) + e * V_0 mod r`, and unblinds `V(z) := (V(z') - e * V_0) mod r = lambda * log_g0(z) mod r`.
  4. Since `Q != O` and `lambda != 0`, `V(h) = lambda * a != 0 mod r`. Compute `b := V(T) * V(h)^(-1) mod r` — cancelling `lambda` identically without ever computing `V(g_0)` or a third descent! Output `P := [b] P_0 in G_1` if `e(P, Q) == T`, otherwise fail.

## 2. Tower, Polynomials, and `Aut = C_12` Symmetry

Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1` asserted live). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`, monogenic maximal order `O_H = Z[zeta_15]`, `h_K = 1`), `g(x) = x^3 - u omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = <sigma, omega> = Z/12Z`. Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U') = 63.317373`, target excess `exc = 1.09`. `audit_cm.py` verifies `algebra_bls12`, `murphy_alpha_exact_live` (`alpha_f = 2.74, alpha_g_total = 0.77`), `real_fb_joint_coll_live`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms`, `voronoi_exact`, and `srt_toy`.

## 3. Relation Collection (`P4_prod_2ad_prim_w10b_box102`, `Arch2-Prod-ZeroCopy-U32`)

- **2D Norm Pool & Finite-Region Margin:** `A2_prod = 313,000, r_ab = 1.02`, `V_kept = 2^120.118602`, `f_prim = 0.95518`, 100% deterministic `c_gbatch_full` (`P_{71.69}` product tree `<= 75.5` u, 256-owner RAM `2^61.22 B << 2^64 B`) replacing collection ECM -> explicit `det_credited = 0.99`. Live `audit_cm.py` values: `margin_net_rethi = +0.2256`b, **`margin_strict_hi = +0.1558`b (`>= 34.4` MC SE)**, `log2_rels_2se = 63.6583 >= 63.658 = LOG2_RELS` (descent input), `nf_mean = 555.75 >= 549.0`, `ng_mean = 680.35 >= 678.0`, `c_cof_2d_upper_live = 1.1558` Ir/cell (with live `Q(zeta_15)` split-prime census up to `2^14` verifying `kh_g, kh_f <= 0.10`).
- **Pure Baseline ARMv8.0-A Kernels (`KERNELS.md`, `STORAGE.md`):** `p_m4_s = 55.0`, `p_m4 = 56.5`, `p_lo = p_f = 59.0`, `p_g = 64.0`, `W_tile = 448 = 7 * 64`, `W_SB = 32`, `rev_b1a_apply_u16` (`3.625 <= 3.65`), 1-rec short buckets (`6.0`/`8.0`), 11-op `rev_b2_prod_emit_1` (`c_hit_hi = 14.85`), 64-cell `rev_scan_zero_u64` (`0.7031 -> 0.82`), `H2_ALLOW = 1.03` -> **`relation_collection <= 2^109.434964` u**, `replay_a_tables <= 2^106.01666` u, 1-owner build RAM `2^63.956 B < 2^64 B`, total tables `2^119.163 B < 2^120 B`.

## 4. Filtering and Linear Algebra (`ArchB-1Owner-3Vec-4Col` + 16 Read-Only Matrix Shards)

- **Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0` at every scale):** `nl_25 = 0.001682` vs `nl_cap = 0.00247` (`1.468x` headroom); `sw_25 = 0.003417` vs `sw_cap = 0.00500` (`1.463x` headroom); pure baseline ARMv8.0-A 3-vec 4-col `unit_run_3vec4_z12` (`92/12 = 7.6667 <= 7.70` Ir/entry-vec, `26 <= 29` GPRs, `0` spills) -> `k5_live = 9.07 <= k5_charged = 10.60` (`+16.9%` headroom, `hex6`, ZERO TABM: `tab_ratio_charged = 1.0`).
- **1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (`m = 16,384, n = 120 = 40 * 3, n_sm = 4`, ALL `< 2^64 B`):** Dense-vector owner `2 * N_sw * 96 B = 2^63.2585 B << 2^64 B`, `16` read-only sparse-matrix shard owners `2^62.325 B << 2^64 B` each, explicit priced cross-owner matrix streaming **`c_stream_mat = 268.0` Ir/row** (`4.0` Ir/8B entry / `3` vectors = `1.3333` Ir/vector-nonzero) charged inside `c_prod`, **`linear_algebra <= 2^107.822657` u**, `matrix_generator <= 2^95.6908` u.

## 5. Main-Phase Success & Missing-Log Reconstruction (`[6b]`)

Collection and filtering are coin- and challenge-independent (`rels_2se = 2^63.6583 >= 1.09 U'`). Under `H3` and `H-X` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx` with `34/34` unscreened two-sided `Q(zeta_15)` BSGS hits + `gamma0_blind_unblind_lambda_cancel_ok = True`, and `unplanted_two_sided_bw` with `pm_basis_verify`), Block Wiedemann fails with probability `<= 2^-66`, ECM relation split Chebyshev failure `< 2^-70`, misc `< 2^-80`, and Schirokauer injectivity `<= 2/r`.

**Reconstruction (`[6b]` `missing_recon_6b_live`):** Purged columns' logarithms are recovered by reverse back-substitution along the singleton-peel stack; columns lost even with raw count `>= 3` during singleton peeling are explicitly counted in `missing_recon_6b_live` (`2^14` `Q(zeta_15)` production profile at `exc = 1.09`): across every occurrence band `[1.5, 3), [3, 5), [5, 8), [8, 14), [14, 30)`, the descent's `P(Poisson(lambda) <= 2)` missing-log model exceeds the measured total unrecovered fraction (including all `raw_cnt >= 3` peeled losses!) by **`3.19x` to `6.73x`** (`meas_miss <= mod_miss` asserted LIVE in `audit_cm.py`, mirroring `c36p` `[6b]`).

## 6. Individual Logarithms

`DESCENT.md` (`verify_q15_descent.py`: `[7]`, `[7l]`, `[7e]`, `[11]`, `[11b]`, `[12]`): `76` uncapped sieve nodes (`bins (71, 76)..(256, 261)`) + whole-ball renewal boot (`BOOT16 = (261, 336, 14, 2)` -> `n1 = 261, nA = 237.25, kA = 14, R_b = 2, BOOT_RELS = 18`, exact `F_{p^6}`-subfield elimination `phi_12(a - b X) = (-phi_4(b) W_1^(-1)) * s * z'` verified live on `4` 381-bit `G_T` targets in `audit_cm.py`), `min_ret = 0.9056`, `max_loss = 9.844e-04`, `min_mu = 68.37`, `max_phi = 0.022509`, `boot_psi = 0.34432`, `boot_phi = 4.627e-09`. All descent LLL/setup constants (`PAIR_SETUP = 2^35.2`, `LLL_NODE = 2^40.0`, `LLL_BOOT = 2^46.5`, `BOUND_CELL = 3,600`, `SIEVE_HIT = 2,200`) are proved input-uniform in `KERNELS.md` and `verify_kernels_symbolic`.

Expected two-log work `E[D] <= 2^103.0698` u (`live = 2^103.065798` u, one log `2^116.3535` Ir). The unified hard meter caps the entire opened run at `T_run = 2^123` u (reached = failure). Metered work is `min(W_run, T_run) <= W_run`: its expectation is `<= 2^109.9564` u, and by Markov's inequality **`P(meter) <= E[W_run] / T_run = 0.00011844 <= 0.00012`** (`0.00013028` at `+10%` expected work).

## 7. Success

Core opened-run success:

    P(success | open) >= (1 - 2^-66 - 2^-70 - 2^-80 - 2/r) * (1 - 2 * 4.627e-09 - 0.00011844 - 2^-12.280 - 2^-10.660 - 4/r)
                      = 0.99906243

and with BOTH `+10%` meter sensitivity (`1.1 * p_meter = 0.00013028`) AND wrapper failure (`2^-48`), **`P(success | open) >= 0.99905059 >= 0.99905`**:

- **Two boots failing all 18 relations:** `2 * boot_phi = 2 * 4.627e-09 = 9.254e-09`.
- **Hard meter `T_run = 2^123` u:** `P(meter) <= 0.00011844` (`0.00013028` at `+10%`).
- **Node counter (`N_MAX = 2^27.8`):** The run aborts as failure before special-`q` node `N_MAX + 1`. `q15_counts16.py` evaluates `E[N]` from the exact same recurrence with unit node cost: **`E[N] = 2^15.5200`** for both logs -> `P(abort) <= 2^-12.280`.
- **Repeats:** `N_MAX` nodes in the smallest bin `(71, 76)` (`2^71.26` ideals), `2^6` for `1/p` weighting: `P(repeat) <= 2^-10.660`, charged.

Subtrees of distinct ideals are independent (`H5`). Past the gate `P(success | open) >= 0.99905`; every abort and a closed gate fail: `P(success) >= (K / 2^64) * 0.99905 >= 0.39`. `W <= 1 + 1[u < K] W_run`, `u` independent of `W_run`: `E[W] <= 1 + (K / 2^64) * (2^109.9564 + 2^66) < 2^108.599317240806`.

## 8. Ledger (`C = 20,000` Ir/u)

| Component | `log2(u)` |
|---|---:|
| `relation_collection` | `109.434964` |
| `replay_a_tables` | `106.01666` |
| `linear_algebra` | `107.822657` |
| `matrix_generator` | `95.6908` |
| `both_individual_logs_exp` | `103.0698` |
| `band_split` | `100.674` |
| `ecm_split_expected` | `96.4533` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Opened-Run Expected Sum (`E[W_run]`)** | **`109.9564`** |
| **Hard Pointwise Meter (`T_run + wrapper`)** | **`123.000000`** |
| **Gated Cold-Run Expected Claim (`E[W]`, `p_cold >= 0.39`)** | **`108.599317240806`** |

## 9. Heuristics (Explicit Exploratory-v1 Premises)

- **`H1` (Smoothness):** Dickman-de Bruijn / Murphy-alpha smoothness (`alpha_f = 2.74, alpha_g_total = 0.77`, exact 2-adic `v_2(N_g)` class law, sampled 200k `box102` norms; boot pair ratio `BALL_PAIR = 2^24` vs measured `1.00078`).
- **`H2` (Collection & Sieve Cost):** Uniform 200k `box102` pool with `H2_ALLOW = 1.03` (`>= 1 + 2 SE`), `det_credited = 0.99`, `kh_g = kh_f = 0.10`, and uncapped `dsieve16` mean per-hit cost `SIEVE_HIT = 2,200` Ir/hit.
- **`H3` (Matrix & Filtering):** `r` unramified and prime to class numbers, Schirokauer injectivity (`n_sm = 4`), `1`-dimensional kernel, and production zero-column-loss filtering caps `nl_cap = 0.00247`, `sw_cap = 0.00500` (`1.46x` above `2^25` ladder), `k5_charged = 10.60`, `c_stream_mat = 268.0`.
- **`H4` (Kernel Cost):** Pure baseline scalar ARMv8.0-A integer instruction counts (`KERNELS.md`, `STORAGE.md`, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64).
- **`H-X` (Block Wiedemann):** Full Krylov/Hankel ranks for block-diagonal `X` and `4` Schirokauer columns at charged `(d, e)` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx`, `unplanted_two_sided_bw`).
- **`H-FFT` (Local SSA):** `fft_F(N)` bounds local large-integer multiplication past `2^26` limbs.
- **`H5` (Descent):** `384`-mesh lower sums; per-ideal `BOUND_CELL` yield certificate (`DISC = 0.001`, backed by `0 / 3,116` top-edge ideals); lazy cofactorisation order independent of outcomes; NegBin(`r = 0.5 / 1.0`) recurrence with independent subtrees of distinct ideals (measured in `[7e]` (A)+(B)+(C)); Poisson missing-log occurrence at `2`-SE count `2^63.658`; exact `O_H` local omitted-mass laws (`q15_retain16.py`); random-integer ECM orders; `UDUP = 0.65` (measured `yr = 0.9739` at `R_boot = 2`).

## 10. Memory & Reproducibility

Peak storage `< 2^120` B (`2^119.163` B replay tables during collection; LA dense owner `2^63.2585 B << 2^64 B` + `16` matrix shards `2^62.325 B << 2^64 B` each; descent `< 2^105` B; wrapper `< 2^30` B). Reproduction: `python3 -I -B audit_cm.py <<< '{}'`.
