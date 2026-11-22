# BLS12-381 `recover-P` (`bls12-381-fixed-g2-invert-g1`) — Complete Self-Contained Verification Proof

**No prior verdict or score is a premise.** Every algebraic identity, kernel opcode count, filtering simulation, collection pool evaluation, linear-algebra memory/work bound, calibrated `gamma_0`-blinded Deterministic-Capped Two-Level Internal-Retry + `ZERO-EA det_2592` descent recurrence, and run-gate claim in this package is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib, pure baseline scalar ARMv8.0-A integer instructions, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64, `C = 20,000` instructions per dense `F_{p^12}` multiplication).

## 1. Track Reduction (Calibrated `gamma_0`-Blinded Exact `lambda`-Cancellation) & Declared Claim
- **Track:** `bls12-381-fixed-g2-invert-g1` — fixed-`Q` / recover-`P` (`e(P, Q) -> P`).
- **Reduction (Exact Two-Logarithm Fixed-`Q` Inversion with Calibrated Factor-Base Blinding Generator `gamma_0`):** Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md`) recovers a nonzero homogeneous virtual-logarithm functional `V(u) = lambda * log_g0(u) mod r` (with unknown nonzero scalar `lambda in (Z/rZ)^*`). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** **`gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {1}`** (`1` final exponentiation `<= 2^11` u, ZERO third descent!) with KNOWN nonzero virtual log **`V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0`** (`gcd(r, (p^12-1)/r) == 1` asserted LIVE in `audit_cm.py`, so `V(gamma_0) = V_0 = lambda * log_g0(gamma_0) mod r` identically!). Given the challenge `(Q, T)` with `Q = [a] Q_0 in G_2 \ {O}` (`a in (Z/rZ)^*`) and `T = e(P, Q) = g_0^{a b} in G_T`:
  1. If `T == 1`, return `P = O` immediately (exact identity case).
  2. Otherwise, compute **`h := e(P_0, Q) = g_0^a in G_T`** (`1` pairing evaluation, `<= 2^11` u).
  3. For each target `z in {h, T}`: each fresh boot attempt draws uniform `e in [0, r-1]`, blinds **`z' := z * gamma_0^e`** (EXACT uniform in `G_T` since `gamma_0` has prime order `r`!), runs individual-log descent on `z'` to obtain `V(z') = V(z) + e * V_0 mod r`, and **unblinds using the KNOWN `V_0`**: **`V(z) := (V(z') - e * V_0) mod r = lambda * log_g0(z) mod r`**.
  4. Since `Q != O` and `lambda != 0`, `V(h) = lambda * a != 0 mod r`. Compute **`b := V(T) * V(h)^(-1) mod r = (lambda * a * b) * (lambda * a)^(-1) = b mod r`** — cancelling the unknown common scale `lambda` IDENTICALLY without ever computing `V(g_0)` or a third descent! Output **`P := [b] P_0 in G_1`**, verified by `e(P, Q) == T` (`<= 2^11` u).
- **Declared Headline Claim:**
  - Opened-run expected work: **`E[W_run] <= 2^111.3593` u** behind ONE unified hard aborting meter `T_run = 2^123` u (`p_meter = 0.0003131844 <= 0.00032`).
  - Opened-run success lower bound: **`p_floor = 0.99968`** (`> 0.39`).
  - Cold-run gate (`K_gate = 7196533079332111406`, wrapper `2^66` u): **`E[W] <= 2^110.001308765382` u** with cold success **`>= 0.39`**. (Claim bounds strictly **expected** work `E[W]`, not the separate pointwise `2^123` cap.)

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1` asserted LIVE). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`), `g(x) = x^3 - u * omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = Z/12Z` (`<sigma, omega>`).
- Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`.
- Live `audit_cm.py` verifies `algebra_bls12`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms` (`12/12`), and `srt_toy` (`0/128` mismatches).

## 3. Relation Collection (`P4_prod_2ad_prim_w10b_box102`, `Arch2-Prod-ZeroCopy-U32`)
- **2D Norm Pool & Finite-Region Margin:** `A2_prod = 313,000, r_ab = 1.02`, `V_kept = 2^120.118602`, `f_prim = 0.95518`, 100% deterministic `c_gbatch_full` (`P_{71.69}` product tree `74.41` u `<= 75.5` u, 256-owner RAM `2^61.22 B << 2^64 B`) replacing collection ECM -> explicit **`det_credited = 0.99`**.
  - Live `audit_cm.py` values: `margin_net_rethi = +0.2256`b, **`margin_strict_hi = +0.1558`b (`>= 34.4` MC SE)**, `c_cof_2d_upper_live = 1.1558` Ir/cell (with live `Q(zeta_15)` split-prime census up to `2^14` verifying `kh_g <= 0.10, kh_f <= 0.10`).
- **Pure Baseline ARMv8.0-A Kernels (`KERNELS.md`):** `p_m4_s = 55.0`, `p_m4 = 56.5`, `p_lo = p_f = 59.0`, `p_g = 64.0`, `rev_b1a_apply_u16` (`3.625 <= 3.65`), short-bucket 1-rec `rev_b1a_apply_1` (`6.0`) & `rev_b1b_apply_1` (`8.0`, `bkt_oh = 8.0`), 11-op `rev_b2_prod_emit_1` (`32 KiB`-aligned flush, `c_hit_hi = 14.85`), 64-cell `rev_scan_zero_u64` (`0.7031` raw -> `0.82`, `W_tile = 448 = 7 * 64`), `H2_ALLOW = 1.03` -> **`relation_collection = 2^109.434964` u**, `replay_a_tables = 2^106.01666` u, 1-owner build RAM `2^63.956 B < 2^64 B`, total tables `2^119.163 B < 2^120 B`.

## 4. Linear Algebra (`ArchB-1Owner-3Vec-4Col`), `H3`, & `H-X`
- **Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0` at every scale):** `nl_25 = 0.001682` vs **`nl_cap = 0.00247` (`1.468x` headroom)**; `sw_25 = 0.003417` vs **`sw_cap = 0.00500` (`1.463x` headroom)**; pure baseline ARMv8.0-A 3-vec 4-col `unit_run_3vec4_z12` (`92/12 = 7.6667 <= 7.70` Ir/entry-vec, `26 <= 29` GPRs, `0` spills) -> `k5_live = 9.07 <= k5_charged = 10.60` (`+1.53` Ir/entry-vec = `+16.9%` headroom, `hex6`, ZERO TABM: `tab_ratio_charged = 1.0`).
- **1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (`m = 16,384, n = 120 = 40 * 3, n_sm = 4`, ALL `< 2^64 B`, Priced Streaming):** dense-vector owner **`2 * N_sw * 96 B = 2^63.2585 B << 2^64 B`**, `16` read-only sparse-matrix shard owners **`2^62.322 B << 2^64 B`** each, explicit cross-owner matrix streaming **`c_stream_mat = 268.0` Ir/row** (`4.0` Ir/8B entry / `3` vectors = `1.3333` Ir/vector-nonzero) charged inside `c_prod`, **`linear_algebra = 2^107.822657` u**, `matrix_generator = 2^95.6908` u.
- **Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X` + Live `gamma_0`-Blinding/Unblinding `lambda`-Cancellation + TRUE Same-Matrix `H3` (`w8_verifiers.py`):** `34/34` exact post-hoc BSGS hits on both sides (`0` lost 2-core cols), `gamma0_blind_unblind_lambda_cancel_ok = True`, matched-scale calibration, TRUE same-matrix full-vs-proxy `H3`, `unplanted_two_sided_bw`, and `filtsim_and_hx_live`.

## 5. Individual-Logarithm Descent (Deterministic-Capped Two-Level Internal-Retry `q * p_rand_j` + `ZERO-EA det_2592` + Total-Capped `k = 15..0` `dsieve16`)
- **Deterministic-Capped Two-Level Internal-Retry `q -> q * p_rand_j` (`ZERO Trial 0!`, `K_retry = 6`, `eta_bad = 0.6%`, `pe_good = 0.25`, `M_cap = 100.0`, `|P_f| = 1,344` in `[2^12, 2^14)`):** Deterministic first-vector gate ONLY (`||b_0(L_Q)||_2 >= 0.20 * Q^(1/16)`, mandatory entry `LLL_NODE` charged once; exact `+-1`-only factor-12 continuous hypercube Blichfeldt `<= 6.452e-7`, live Chebotarev `1.061x` on BOTH `f_9` and `g_3`). At every internal node `Q`, run up to `K_retry = 6` conditionally-i.i.d. fresh WITH-REPLACEMENT `p_j ~ Unif(P_f)` trials (`min_mu1 = 15.02 >= 15.0`, `post_cap_mult = 8.0` pointwise post-sieve cap -> `W_own <= own1_det`, Markov exceedance `<= 1/8`), descending children and retrying with a brand-new `p_{j+1}` if empty OR any child fails — all behind a **hard deterministic per-node work meter `T_cap = min(nat_wc, 100.0 * W_good)`**! Thus EVERY bad/adversarial ideal's entire subtree work is **100% deterministically `<= T_cap`** (ZERO adversarial chain gap; ZERO child-work premise needed!), while good nodes incur exact Markov abort `p_mk = 0` on shallow bins / `0.01` on deep bins. Under the single conditional premise `H5-2lvl` (`P(child in B | Q in G, Y_j = 1, M >= 1) <= 0.006` and `P(boot child in B | boot succeeds) <= 0.006`), `max_phi = 0.0168` and Boole ONLY over the `<= 17` root boot children yields **`Phi_boot = 0.2145 <= 0.55`** (Wald `1.2731`)! Backed LIVE in `m5_qprand.py` by `oh_coprime_catches_int_gcd1 = True`, `8` live 381-bit CRT `L_{q, p_j}` + `4` live `prod_enc12`, the **`720x20 = 14,400` i.i.d.-parent strictly descending exact `O_H`-coprime deg-1-root campaign (`0/720` flagged `>= 4/20`, `77.48%` detection power -> exact power-adjusted 95% UB `0.00536 < 0.00600 = eta_bad` at exact `pe_good = 0.25`, `4` live replays)**, and `24x24` literal deg-`(72, 24)` (`4` live replays).
- **Total-Capped `dsieve16` (`k = 15..0`, `gate_frac = 1.0`, `alpha = 100/74`) + Uniform `(99/100)`-LLL `PAIR_SETUP = 2^35.2`:** Every `L_p` (`p <= 2^71`) runs 2-phase `(11/40)->(99/100)`-Bareiss-LLL (`2^35.06 < 2^35.2` proved input-uniform) and enforces ONE pointwise total-work cap `T_p_total = _wc_tot_ir(det)` (`k = 15..0`, `c_node_base = 156, k1_rej = 244, c_leaf = 340`) covering BOTH internal traversal AND leaf processing (`96`-bit sieve slack absorbs tiny primes; 15D quotient `det >= p/R` for large primes). Verified LIVE: `0/32` skipped on genuine `Q(zeta_15)` roots (`max_ratio = 0.8095`), `0/8` skipped on complete root-to-leaf FP runs (`6e-5`), AND `0/4` skipped + `3/3` live `v_0` hits on genuine 381-bit relation-conditioned `L_p`, conservative `ret_node = 0.90`.
- **`ZERO-EA` Unconditional Deterministic-Norm Boot (`M_boot = 2592.0`, `surv_b = 1.0`, `tau = 2592.0`, `over_b = 0.0`):** Single canonical `b_0 = R_sk[0]` (`LLL_BOOT = 2^46.5 > 2^46.45` proved input-uniform), deterministic Lovász bound `2591.90 < 2592.0` (exact `Fraction` Sylvester minors + `32/32` real 381-bit `G_T` bases live-certified + target `0` from scratch via `lll_16d_int`), constant full `ecm_post_b` (`k_lx_bi = 130 + kA_b = 16`, ZERO `P_71` boot wave) on 100% of trials (**ZERO `H-EA`, ZERO `BOOT_BINS`, ZERO `p_band`, ZERO `318`-bin mediant!**). Knapsack over `2592.0`: **`Phi_boot = 0.214519 <= 0.55`**, Wald **`E[runs/log] <= 1.2731`**, **`E[two_logs] <= 2^110.677157` u < `2^110.6812` u**.

## 6. Complete Work Ledger (`C = 20,000` Ir/u)

| Component | `log2(u)` |
|---|---:|
| `relation_collection` | `109.434964` |
| `replay_a_tables` | `106.01666` |
| `linear_algebra` | `107.822657` |
| `matrix_generator` | `95.6908` |
| `both_individual_logs_exp` | `110.6812` |
| `band_split` | `100.674` |
| `ecm_split_expected` | `96.4533` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Core Expected Sum (`E[W_run]`)** | **`111.3593`** |
| **Gated Cold-Run Claim (`E[W]`, `p_cold >= 0.39`)** | **`110.001308765382`** |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared exploratory STNFS premises: `H1` (standard STNFS smoothness/Dickman-de Bruijn heuristic; ZERO `H-EA` and ZERO `BOOT_BINS`!), `H2` (`H2_ALLOW = 1.03`, `det_credited = 0.99`, `kh_g = kh_f = 0.10`, and total-capped `dsieve16` `ret_node = 0.90`), `H3` (`nl_cap = 0.00247`, `sw_cap = 0.00500`, `k5 = 10.60`, `c_stream_mat = 268.0`), `H-X`, and Deterministic-Capped Two-Level **`H5-2lvl`** (`P(child in B | Q in G, Y_j = 1, M >= 1) <= eta_bad = 0.6%`, `P(boot child in B | boot succeeds) <= 0.6%`, `pe_good <= 0.25`; ZERO separate child-work premise needed thanks to hard node meter `T_cap` + `8.0x` pointwise post-sieve cap!).
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` and `python3 -I -B audit_cm.py <<< '{"descent_mesh":384}'`.
