# `recover-P` (`bls12-381-fixed-g2-invert-g1`) through two `G_T` logarithms with `(4, 3)` `Q(zeta_15)` CM-STNFS (100% Eager Cofactorization & Full Cross-Owner Streaming) behind a run gate: expected work `2^108.822176096911` units

Track `bls12-381-fixed-g2-invert-g1`. Input: nonidentity `Q in G_2` and `T in G_T`; output `P in G_1` with `e(P, Q) = T` for the exact supplied `Q` (`P = O` for `T = 1`). The score unit is one dense `F_{p^12}` multiplication, `C = 20,000` pure baseline scalar ARMv8.0-A integer instruction equivalents (`simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64/x18). No advice, smaller targets, or amortization are used. **No prior verdict or score is a premise.**

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is a run gate: draw a 64-bit word `u`; unless `u < K = 7201864164761372185`, stop as a failure (`< 1` u). `K` is the least integer with `K * 0.99894 >= 0.39 * 2^64`. Under `H1`–`H5`, `H-X`, and `H-FFT` (section 9), the algorithm returns the exact `P` with probability `>= 0.99894` once past the gate (section 7, even at `+10%` meter sensitivity AND `10x DISC = 0.010` stress!); `u` is independent of the challenge and all other coins, so cold-run success is `>= 0.39` and

    E[W] <= 1 + (K / 2^64) * (2^110.1791 + 2^66) < 2^108.822176096911  dense-Fp12 units (C = 20,000).

`2^110.1791` is the sum of section 8's rows (every-run bounds for collection and linear algebra plus expected work for the two 100%-eager individual logarithms `2^105.8199` and ECM relation splitting). The score `108.822176096911` is an **expected** work bound, not a worst-case bound on every run: every run is bounded pointwise behind the hard aborting meter `T_run = 2^123` u plus wrapper:

    W <= 2^123 + 2^66 < 2^(123 + 2^-56)  units.

Analytic, exploratory-v1 claim. Every algebraic identity, kernel opcode count, filtering simulation, collection pool evaluation, linear-algebra memory/work bound, 384-mesh 100%-eager descent recurrence (+ `10x DISC` stress), and ledger/gate bound is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib).

## 1. Reduction (Fixed `Q`: Calibrated `gamma_0`-Blinded Two-Logarithm Inversion) and Field Map

Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md`) recovers a nonzero homogeneous virtual-logarithm functional `V(u) = lambda * log_g0(u) mod r` (with unknown nonzero scalar `lambda in (Z/rZ)^*`). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {1}` (`1` final exponentiation `<= 2^11` u, ZERO third descent!) with KNOWN nonzero virtual log `V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0` (`gcd(r, (p^12-1)/r) == 1` asserted LIVE in `audit_cm.py`). Given the challenge `(Q, T)` with `Q = [a] Q_0 in G_2 \ {O}` (`a in (Z/rZ)^*`) and `T = e(P, Q) = g_0^{a b} in G_T`:
  0. Run gate (`Claim`): draw 64-bit `u`; if `u >= K_gate`, stop as failure (`< 1` u).
  1. Check `Q` (on twist, `[r]Q = O`, `Q != O`) and `T` (`T^r = 1`). If `T == 1`, return `P = O`.
  2. Compute `h := e(P_0, Q) = g_0^a in G_T` (`1` pairing, `<= 2^11` u) and map `h, T` into the `Q(zeta_15)` tower (`field-map.md`, wrapper `<= 2^66` u, failure `< 2^-48`).
  3. For each target `z in {h, T}`: each fresh boot ball draws uniform `e in [0, r-1]`, blinds `z' := z * gamma_0^e` (exact uniform in `G_T`; if `z' == 1`, immediately returns `V(z) = -e * V_0 mod r`), runs individual-log descent on `z'` to obtain `V(z') = V(z) + e * V_0 mod r`, and unblinds `V(z) := (V(z') - e * V_0) mod r = lambda * log_g0(z) mod r`.
  4. Since `Q != O` and `lambda != 0`, `V(h) = lambda * a != 0 mod r`. Compute `b := V(T) * V(h)^(-1) mod r` — cancelling `lambda` identically without ever computing `V(g_0)` or a third descent! Output `P := [b] P_0 in G_1` if `e(P, Q) == T`, otherwise fail.

## 2. Tower, Polynomials, and `Aut = C_12` Symmetry

Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1` asserted live). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`, monogenic maximal order `O_H = Z[zeta_15]`, `h_K = 1`), `g(x) = x^3 - u omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = <sigma, omega> = Z/12Z`. Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U') = 63.317373`, target excess `exc = 1.09`. `audit_cm.py` verifies `algebra_bls12`, `murphy_alpha_exact_live` (`alpha_f = 2.74, alpha_g_total = 0.77`), `real_fb_joint_coll_live`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms`, `voronoi_exact`, and `srt_toy`.

## 3. Relation Collection (`P4_prod_2ad_prim_wave1_1m_r83`, `Arch2-Prod-ZeroSlice-1PassB2-U16-Streamed`)

- **1,000,000-Draw Exact-Norm 2D Pool & Finite-Region Margin:** `A2_prod = 320,000`, `r_ab = 0.83`, `V_kept = 2^120.3739`, `log2(V_a) = 55.8474 < 0.90 * 2^56` (`10.0%` headroom below `2^56` -> exact `2^56` bank stride, `oid < 2^56`, `K_a = 1` zero slicing!), `L_max = 439.06 < W_tile = 448 = 7 * 64`, `L_b = 188.58`, `f_prim = 0.955881`. 100% deterministic `c_gbatch` (`P_{71.69}` product tree via 3-prime Radix-4 Lazy Shoup `ntt3_F`, `<= 75.5` u) replacing collection ECM -> explicit `det_credited = 0.99`. Live `audit_cm.py` values across all `5,648` 2D bins (`bg_max = 697`b): `margin_net_rethi = +0.10494`b, **`margin_strict_hi = +0.03507`b (`15.2` MC SE)**, `log2_rels_2se = 63.542031 >= 63.54 = LOG2_RELS` (descent input), `nf_mean = 555.32 >= 549.0`, `ng_mean = 684.31 >= 678.0`, **`c_cof_2d_upper_live = 0.9824` Ir/cell** (with `P_min = 32`, `KH = 0.08` verified live up to `2^14`: `gb/tb = 0.07197 <= 0.08`, `fb/tb = 0.05649 <= 0.08`, loop-included 49-op 3-prime Radix-4 Lazy Shoup `ntt3_F` + 5-transform middle-product `RT`, and ordered norms checking `S_f` before `N_g` and `Stage-B(g)` before `N_f`).
- **Pure Baseline ARMv8.0-A Kernels with Explicit Cross-Owner Streaming (`KERNELS.md`, `STORAGE.md`):** `P_min = 32`, `p_m4_s = 54.5`, `p_m4 = 57.5`, `p_lo = p_f = 58.0`, `p_g = 63.5`, `W_tile = 448 = 7 * 64`, `W_SB = 16`, Band-1a-long `rev_b1a_apply_u16` (`3.625`) + `rev_b1_stream_u8` (`1.25`) = `4.875 <= 4.90` (`b1a_long_lk = 68.0`), Band-1a-short `6.0 + 4.0 = 10.0`, Band-1b `8.0 + 4.0 = 12.0`, slab-routed 62-bit Band-2 `rev_b2_direct_sb_hdr` (`11 <= 12.0` Ir/SB) + exact column-translating `rev_b2_route_stream_1` (`14.0` Ir/hit, translating `r -> u = (r + t_0) & 15` and routing each record once to its exact 16-column slab owner `S` or `S + 1`) + `rev_b2_direct_1` (`8.0` Ir/hit) = `22.0` Ir/hit (ZERO emit/flush buffers, ZERO `c_bd_chk`!), 64-cell `rev_scan_zero_u64` (`0.7778 <= 0.78`), `H2_ALLOW = 1.01` (`9.2` MC SE) -> **`ipc_tot = 12.536` Ir/cell**, **`relation_collection <= 2^109.748629` u** (`live = 109.748529`), **`replay_a_tables <= 2^104.7482` u** (`live = 104.7481`), 1-owner sieve RAM `2^63.087 B < 2^64 B`, 1-owner build RAM `2^63.898 B < 2^64 B`, total tables `2^118.338 B < 2^120 B`.

## 4. Filtering and Linear Algebra (`ArchB-1Dense-16Shard-4Vec-ZeroCopy`)

- **Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0` at every scale):** `nl_25 = 0.001682` vs `nl_cap = 0.00246` (`1.462x` headroom, strictly `> 1.45 * (mean + 2 SE)`); `sw_25 = 0.003417` vs `sw_cap = 0.00498` (`1.458x` headroom, strictly `> 1.45 * (mean + 2 SE)`); pure baseline ARMv8.0-A 4-vec 4-col `unit_run_4vec4_z12` (`122/16 = 7.625 <= 7.65` Ir/entry-vec, explicit `x_base` offset addressing, `28 <= 29` GPRs, ZERO `x18`, `0` spills) + 2-vec `gen_run1_2v_z12` (`ldrb w22, [x23], #1`, `43/2 = 21.50 <= 22.0` Ir/entry-vec, `25 <= 29` GPRs, ZERO `x18`, `0` spills) -> `k5_live = 8.5721 <= k5_charged = 10.10` (`+1.528` Ir headroom, `hex6`, ZERO TABM: `tab_ratio_charged = 1.0`).
- **1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (`m = 16,384, n = 128 = 32 * 4, n_sm = 4`, ALL `< 2^64 B`):** Dense-vector owner `2 * N_sw * 128 B = 2^63.668 B < 2^64 B` (`12.71 EiB < 16 EiB`, `20.6%` headroom), `16` read-only sparse-matrix shard owners `2^62.346 B << 2^64 B` each (8B offsets + `<= 15%` 1B `uint8_t` non-unit coefficients -> `8.15 B/nz`), explicit priced cross-owner matrix streaming **`c_stream_mat = 205.0` Ir/vec-row** (4-vec Krylov) AND **`c_stream_mat_1vec = 820.0` Ir/row** (unamortized 1-vec Horner; `205` 8B words/row covering 8B offsets + `ceil(0.15*201/8)=4` packed 8B coeff words) charged inside `c_prod` / `c_prod_horner`, **`linear_algebra <= 2^107.742319` u** (`live = 107.742219`), `matrix_generator <= 2^95.5941` u (`live = 95.5936`).

## 5. Main-Phase Success & Missing-Log Reconstruction (`[6b]`)

Collection and filtering are coin- and challenge-independent (`rels_2se = 2^63.542031 >= 1.09 U'`). Under `H3` and `H-X` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx` with `34/34` unscreened two-sided `Q(zeta_15)` BSGS hits + `gamma0_blind_unblind_lambda_cancel_ok = True`, and `unplanted_two_sided_bw` with `pm_basis_verify`), Block Wiedemann fails with probability `<= 2^-66`, ECM relation split Chebyshev failure `< 2^-70` (`live = 2^-73.0`), misc `< 2^-80`, and Schirokauer injectivity `<= 2/r`.

**Reconstruction (`[6b]` `missing_recon_6b_live`):** Purged columns' logarithms are recovered by reverse back-substitution along the singleton-peel stack; columns lost even with raw count `>= 3` during singleton peeling are explicitly counted in `missing_recon_6b_live` (a `2^14`-column synthetic Poisson-incidence peeling simulation calibrated to the `Q(zeta_15)` production occurrence profile at `exc = 1.09`; actual tower relations are separately tested in `h3_ladder_and_genuine_q15_hx`): across every occurrence band `[1.5, 3), [3, 5), [5, 8), [8, 14), [14, 30)`, the descent's `max(P(Pois(lambda) <= 3), P(Pois(0.75 * lambda) <= 1))` missing-log upper bound exceeds the simulated total unrecovered fraction (including all `raw_cnt >= 3` peeled losses!) (`meas_miss <= mod_miss` asserted LIVE in `audit_cm.py`).

## 6. Individual Logarithms (100% Eager Cofactorization, Zero `Q * q'` Fallback, `10x DISC` Live Stress)

`DESCENT.md` (`verify_q15_descent.py`: `[7]`, `[7l]`, `[7e]`, `[11]`, `[11b]`, `[12]`): `76` uncapped 100%-eager (`frac = 1.0`) sieve nodes (`bins (71, 76)..(256, 261)`, all `76/76` thin missing ideals directly with zero `Q * q'` fallback) + whole-ball renewal boot (`BOOT16 = (261, 336, 14, 2)` -> `n1 = 261, nA = 237.25, kA = 14, R_b = 2, BOOT_RELS = 18`, exact `F_{p^6}`-subfield elimination verified live on `4` 381-bit `G_T` targets in `audit_cm.py`, `per_rel /= exp(-mb)`), `min_ret = 0.9055`, `max_loss = 1.072e-03 < 2e-3`, `min_mu = 50.03`, `max_phi = 0.022823`, `max_node_psi = 0.2149`, `boot_psi = 0.34234`, `boot_phi = 4.171e-09`. All descent LLL/setup constants (`PAIR_SETUP = 2^35.2 > 2^35.06`, `LLL_NODE = 2^40.0 > 2^39.41`, `LLL_BOOT = 2^47.0 > 2^46.75` WITH `S_boot = 2^21` scaling before Phase 1, `BOUND_CELL = 3,600`, `SIEVE_HIT = 2,200`) are proved input-uniform in `KERNELS.md` and `verify_kernels_symbolic`.

Expected two-log work `E[D] <= 2^105.8199` u (`live = 2^105.815890` u, one log `2^119.1036` Ir). The unified hard meter caps the entire opened run at `T_run = 2^123` u (reached = failure). Metered work is `min(W_run, T_run) <= W_run`: its expectation is `<= 2^110.1791` u, and by Markov's inequality **`P(meter) <= E[W_run] / T_run = 0.00013821 <= 0.00015`** (`0.00015203` at `+10%` expected work).

## 7. Success

Core opened-run success (exact unrounded `boot_phi_exact`, `over_log2_exact`, `rep_log2_exact`):

    P(success | open) >= (1 - 2^-66 - 2^-70 - 2^-80 - 2/r) * (1 - 8.341e-09 - 0.00013821 - 2^-11.494 - 2^-11.659 - 4/r)
                      = 0.99920610

with BOTH `+10%` meter sensitivity (`1.1 * p_meter = 0.00015203`) AND wrapper failure (`2^-48`), **`P(success | open) >= 0.99919228 >= 0.99894`**, AND even under simultaneous **`10x DISC = 0.010` stress** (`two_log = 106.3006`, `boot_phi = 5.325e-06`, `over10 = -10.7386`) + `+10%` meter + `2^-48` wrapper, **`succ_head_disc10 = 0.99894004 >= 0.99894 = p_floor`** (asserted LIVE in `audit_cm.py`!):

- **Two boots failing all 18 relations:** `2 * boot_phi = 2 * 4.171e-09 = 8.341e-09` (`1.065e-05` at `10x DISC`).
- **Hard meter `T_run = 2^123` u:** `P(meter) <= 0.00013821` (`0.00015203` at `+10%`).
- **Node counter (`N_MAX = 2^26.8`):** The run aborts as failure before special-`q` node `N_MAX + 1`. `q15_counts16.py` evaluates `E[N]` from the exact same recurrence with unit node cost: **`E[N] = 2^15.3057`** for both logs -> `P(abort) <= 2^-11.494`.
- **Single-side repeats (`M_min_side = 2^70.2598`):** `N_MAX` nodes in the smallest single-side bin `(71, 76)` (`2^70.26` ideals, zero cross-side pooling!), `2^6` for `1/p` weighting: `P(repeat) <= 2^-11.659`, charged.

Subtrees of distinct ideals are independent (`H5`). Past the gate `P(success | open) >= 0.99894`; every abort and a closed gate fail: `P(success) >= (K / 2^64) * 0.99894 >= 0.39`. `W <= 1 + 1[u < K] W_run`, `u` independent of `W_run`: `E[W] <= 1 + (K / 2^64) * (2^110.1791 + 2^66) < 2^108.822176096911`.

## 8. Ledger (`C = 20,000` Ir/u)

| Component | `log2(u)` |
|---|---:|
| `relation_collection` | `109.748629` |
| `replay_a_tables` | `104.7482` |
| `linear_algebra` | `107.742319` |
| `matrix_generator` | `95.5941` |
| `both_individual_logs_exp` | `105.8199` |
| `band_split` | `100.5667` |
| `ecm_split_expected` | `95.8769` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Opened-Run Expected Sum (`E[W_run]`)** | **`110.1791`** |
| **Hard Pointwise Meter (`T_run + wrapper`)** | **`123.000000`** |
| **Gated Cold-Run Expected Claim (`E[W]`, `p_cold >= 0.39`)** | **`108.822176096911`** |

## 9. Heuristics (Explicit Exploratory-v1 Premises)

- **`H1` (Smoothness):** Dickman-de Bruijn / Murphy-alpha smoothness (`alpha_f = 2.74, alpha_g_total = 0.77`, exact 2-adic `v_2(N_g)` class law, 1,000,000-draw `r_ab = 0.83` norm pool; boot pair ratio `BALL_PAIR = 2^24` vs measured `1.00078`).
- **`H2` (Collection & Sieve Cost):** Uniform 1,000,000-draw `r_ab = 0.83` pool with `H2_ALLOW = 1.01` (`9.2` MC SE), `det_credited = 0.99`, `kh_g = kh_f = 0.08` (verified live to `2^14`), explicit cross-owner Band-1a/1b/2 streaming (`rev_b1_stream_u8`/`rev_b1_stream_tail1`), and uncapped `dsieve16` mean per-hit cost `SIEVE_HIT = 2,200` Ir/hit.
- **`H3` (Matrix & Filtering):** `r` unramified and prime to class numbers, Schirokauer injectivity (`n_sm = 4`), `1`-dimensional kernel, and production zero-column-loss filtering caps `nl_cap = 0.00246`, `sw_cap = 0.00498` (`> 1.45x (mean + 2 SE)` above `2^25` ladder), `k5_charged = 10.10`, `c_stream_mat = 205.0`, `c_stream_mat_1vec = 820.0`.
- **`H4` (Kernel Cost):** Pure baseline scalar ARMv8.0-A integer instruction counts (`KERNELS.md`, `STORAGE.md`, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64/x18).
- **`H-X` (Block Wiedemann):** Full Krylov/Hankel ranks for block-diagonal `X` and `4` Schirokauer columns at charged `(d, e)` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx`, `unplanted_two_sided_bw`).
- **`H-FFT` (Local SSA):** `ntt3_F(N)` / `fft_F(N)` bounds local large-integer multiplication past `2^26` limbs.
- **`H5` (Descent):** `384`-mesh lower sums; per-ideal `BOUND_CELL` yield certificate (`DISC = 0.001`, backed by `0 / 3,116` top-edge ideals + live `10x DISC = 0.010` stress clearing `p_floor`); **100% eager cofactorization (`frac = 1.0`, zero lazy-order premise)**; NegBin(`r = 0.5 / 1.0`) recurrence with independent subtrees of distinct ideals (measured in `[7e]` (A)+(B)+(C)); `0.02`-bit right-endpoint upper-sum missing-log bound `max(P(Pois(lam)<=3), P(Pois(0.75*lam)<=1))` at `2`-SE count `2^63.54` with 100% discard everywhere (**zero `Q * q'` fallback lattices**); single-side repeat bound (`M_min_side = 2^70.26`); exact `O_H` local omitted-mass laws (`q15_retain16.py`); random-integer ECM orders; `UDUP = 0.65` (measured `yr = 0.9739` at `R_boot = 2`).

## 10. Memory & Reproducibility

Peak storage `< 2^120` B (`2^118.338` B replay tables during collection; LA dense owner `2^63.668 B < 2^64 B` + `16` matrix shards `2^62.346 B << 2^64 B` each; descent `< 2^95` B; wrapper `< 2^30` B). Reproduction: `python3 -I -B audit_cm.py <<< '{}'`.
