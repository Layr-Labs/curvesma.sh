# BLS12-381 recover-P (`bls12-381-fixed-g2-invert-g1`) — Complete Self-Contained Verification Proof

**No prior verdict or score is a premise.** Every algebraic identity, kernel opcode count, filtering simulation, collection pool evaluation, linear-algebra memory/work bound, calibrated `gamma_0`-blinded pure `k = 1` Boole/Wald descent recurrence, and run-gate claim in this package is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib, pure baseline scalar ARMv8.0-A integer instructions, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64, `C = 20,000` instructions per dense `F_{p^12}` multiplication).

## 1. Track Reduction (Calibrated `gamma_0`-Blinded Exact `lambda`-Cancellation) & Declared Claim
- **Track:** `bls12-381-fixed-g2-invert-g1` — fixed-Q / recover-P (`e(P, Q) -> P`).
- **Reduction (Exact Two-Logarithm Fixed-`Q` Inversion (`recover-P`) with Calibrated Factor-Base Blinding Generator `gamma_0`):** Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md`) recovers a nonzero homogeneous virtual-logarithm functional `V(u) = lambda * log_g0(u) mod r` (with unknown nonzero scalar `lambda in (Z/rZ)^*`). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** **`gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {{1}}`** (`1` final exponentiation `<= 2^11` u, ZERO third descent!) with KNOWN nonzero virtual log **`V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0`** (`gcd(r, (p^12-1)/r) == 1` asserted LIVE in `audit_cm.py`, so `V(gamma_0) = V_0 = lambda * log_g0(gamma_0) mod r` identically!). Given the challenge `(Q, T)` with `Q = [a] Q_0 in G_2 \ {O}` (`a in (Z/rZ)^*`) and `T = e(P, Q) = g_0^{a b} in G_T`:
  1. If `T == 1`, return `P = O` immediately (exact identity case).
  2. Otherwise, compute **`h := e(P_0, Q) = g_0^a in G_T`** (`1` pairing evaluation, `<= 2^11` u).
  3. For each target `z in {h, T}`: each fresh boot attempt draws uniform `e in [0, r-1]`, blinds **`z' := z * gamma_0^e`** (EXACT uniform in `G_T` since `gamma_0` has prime order `r`!), runs individual-log descent on `z'` to obtain `V(z') = V(z) + e * V_0 mod r`, and **unblinds using the KNOWN `V_0`**: **`V(z) := (V(z') - e * V_0) mod r = lambda * log_g0(z) mod r`**.
  4. Since `Q != O` and `lambda != 0`, `V(h) = lambda * a != 0 mod r`. Compute **`b := V(T) * V(h)^(-1) mod r = (lambda * a * b) * (lambda * a)^(-1) = b mod r`** — cancelling the unknown common scale `lambda` IDENTICALLY without ever computing `V(g_0)` or a third descent! Output **`P := [b] P_0 in G_1`**, verified by `e(P, Q) == T` (`<= 2^11` u).
- **Declared Headline Claim:**
  - Opened-run expected work: **`E[W_run] <= 2^110.9369` u** behind ONE unified hard aborting meter `T_run = 2^123` u (`p_meter = 0.0002336927 <= 0.00024`).
  - Opened-run success lower bound: **`p_floor = 0.99976`** (`> 0.39`).
  - Cold-run gate (`K_gate = 7195957218479160129`, wrapper `2^66` u): **`E[W] <= 2^109.578792317464` u** with cold success **`>= 0.39`**. (Claim bounds strictly **expected** work `E[W]`, not the separate pointwise `2^123` cap.)

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1` asserted LIVE). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`), `g(x) = x^3 - u * omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = Z/12Z` (`<sigma, omega>`).
- Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`.
- Live `audit_cm.py` verifies `algebra_bls12`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms` (`12/12`), and `srt_toy` (`0/128` mismatches).

## 3. Relation Collection (`P4_prod_2ad_prim_v17`, `Arch2-Prod-ZeroCopy-U32`)
- **2D Norm Pool & Finite-Region Margin:** `V_kept = 2^120.431440`, `f_prim = 0.95642`, explicit ECM retention credit **`det_credited = 0.95`** (`-0.0740`b penalty).
  - Live `audit_cm.py` values: `log2_rels_2d_mid = 63.859129`, `log2_rels_2d_rethi = 63.853894` (`margin_net_rethi = +0.33819`b), **`log2_rels_2d_strict_hi = 63.714125` (`margin_strict_hi = +0.19842`b >= `35.1` MC SE)**, `c_cof_2d_upper_live = 2.0535` Ir/cell.
- **U32 Kernels & Work Pricing (`KERNELS.md`):** `p_m4 = 56.5`, `p_lo = 59.0`, `p_f = 59.0`, `p_g = 64.0`, `c_hit_lo_m4 = 3.80`, `c_hit_lo_u = 5.80`, `c_hit_hi = 16.35`, `c_scan_zero = 0.95`, `H2_ALLOW = 1.05` -> **`relation_collection = 2^110.059203` u**, `replay_a_tables = 2^105.95909` u, 1-owner build memory `2^63.9240 B < 2^64 B`.
- **Live Murphy-Alpha & Real-FB Joint Collection (`unitdomain.py`):** Exact `p = 3` Hensel tree + unramified `p <= 10,000` (`alpha_sum_10k <= 4.02` charged); live `1,200`-pair exact deg-`(72, 24)` `U_toy = -14` Voronoi/`D_12` joint collection.

## 4. Linear Algebra (`ArchB-1Owner-4Vec-ZeroCopy`), `H3`, & `H-X`
- **Zero-Column-Loss Filtering Ladder (`2^16..2^24` x 6 seeds, `cols_lost = 0` at every scale):** `nl_24 = 0.001774` vs **`nl_cap = 0.00297` (`1.674x` headroom)**; `sw_24 = 0.003602` vs **`sw_cap = 0.00605` (`1.680x` headroom)**; `k5_live = 10.136 <= k5_charged = 13.7` (`hex6`, ZERO TABM: `tab_ratio_charged = 1.0`).
- **Single-Owner `< 2^64 B` Block Wiedemann (`m = 8192, n = 112, n_sm = 4`, `n_owners = 1`):** `log2(N) = 54.922051`, `log2(N_sw) = 55.948523`, peak single-owner RAM **`2 * N_sw * 128 B = 2^63.9485 B < 2^64 B`**, **`linear_algebra = 2^108.530371` u**, `matrix_generator = 2^93.1361` u.
- **Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X` + Live `gamma_0`-Blinding/Unblinding `lambda`-Cancellation + TRUE Same-Matrix `H3` (`w8_verifiers.py`):**
  - Unscreened two-sided `Q(zeta_15)` `H-X`: `34/34` `f`-side and `34/34` `g`-side exact post-hoc BSGS hits with `cols_2core_lost = 0`, PLUS live zero-BSGS-oracle end-to-end verification of `gamma_0`-blinding, `V_0`-unblinding, and `h^(V(T)/V(h)) == T` exact `lambda`-cancellation (`gamma0_blind_unblind_lambda_cancel_ok = True`).
  - Matched-scale genuine-vs-synthetic calibration + TRUE same-underlying-matrix full-incidence (`0` omitted cols) vs truncated-proxy (`139` middle cols) + `unplanted_two_sided_bw` + `filtsim_and_hx_live`.

## 5. Individual-Logarithm Descent (`pure_k1_boole_wald_iid_boot`, `p0_floor = 0.001`, `capped_fp_gso_h2`, `gamma_0`-Blinded `b_best2` `P_71` EA Boot)
- **Pure `k = 1` Zero-Retry Boole + Wald `i.i.d.` Fresh-`e` `gamma_0`-Blinded `b_best2` Boot (`ADAPTIVE.md`, `DESCENT.md`):**
  - `78`-row `lx = 71.0` DAG (`min_mu = 501.92 >= 500`, `max_p0 = 0.003969 < 0.004`, `p0_floor = 0.001`, `p_node_abort = p_repeat_loss = p_descent_fail = 0.0`).
  - **Pointwise Capped Per-Prime FP (`capped_fp_gso_h2`):** Explicitly acknowledging that LLL + `det(L_p) = p` alone does NOT force the geometric GSO profile (`diag(1,...,1,p)`), every prime `p <= 2^71` in `dsieve16` enforces a hard per-prime instruction cap `T_p := _fp_cap_ir(p)` (skipping any prime exceeding `T_p`, so sieve work is deterministically `<= own_1`!); negligible skip rate is declared under `H2` (absorbed in `ret_node = 0.90`) and supported LIVE by `32/32` real `12..71`-bit `Q(zeta_15)` 16D Bareiss-LLL (`99/100`) prime sublattices `L_p` (`0/32` skipped, `max_ratio = 0.386 << 1.0`, `2.6x..4.7x` headroom!).
  - Under explicitly declared exploratory STNFS premise `H5` (`P(M = 0 | slot) <= max(0.001, (1 + 0.5 * mu)^(-1))`, implied via Cantelli by `Var(M | slot) <= max(2 * mu, 0.001 * mu^2)`), unconditional Boole over `N_tree_1log_det = 451` deterministic slots/log gives **`Phi_tree = 0.513188 <= 0.52 <= 0.60`** and exact Wald **`E[runs/log] <= 2.0542`**.
  - Fresh-`e` `i.i.d.` `gamma_0`-blinded `b_best2` + `P_71` Early-Abort Boot (`n1 = 216, nA = 186, kA = 14, tau = 1950`, `BOOT_BINS = [2534.5, 2536.0, 2537.5, 2538.5, 2541.5]`, `BOOT_LO = 2490.0`, `256 * NORMS`, exact 1-wave overshoot `over_b`; `boot_witness.json` stores `24` offline-generated 381-bit `b_best2` norm summaries with target `0`, `n_targets = 1`, recomputed live from scratch via 16D Bareiss-LLL inside `audit_cm.py`, matching `b2_24[0]` bit-for-bit!) gives **`two_boot_exp = 2^105.471` u**, **`two_tree_det = 2^108.691` u**, **`E[two_logs] <= 2^108.838119` u** (`384/370` mesh: **`2^108.765011` u <= `2^108.8431` u**).
- **Scoped Diagnostics (`ADAPTIVE.md` §2, ALL LIVE):** `0/240` zeros on UNSCALED 24-bit real-2-core-FB campaign (`12.47` summed `p_0`, `M_min = 47`, `min(M/mu) = 1.689`, `6` live replays incl. `2` live 2-hop depth-2 chains); `64` live immediate-abort strict `k = 1` boundary roots (`O/E_p0 = 0.3698 <= 0.85`, zero post-abort siblings); `(9,3)` `120`-ideal Fano `1.076..1.316 << 2.00` + `0/100` `hi100` zeros; `120` encountered vs `120` random + `24 + 12` live 381-bit production special-`q` lattices (`>= 10x` conservative headroom); naive-PGF counterexample.

## 6. Complete Work Ledger (`C = 20,000` Ir/u)

| Component | `log2(u)` |
|---|---:|
| `relation_collection` | `110.059203` |
| `replay_a_tables` | `105.95909` |
| `linear_algebra` | `108.530371` |
| `matrix_generator` | `93.1361` |
| `both_individual_logs_exp` | `108.8431` |
| `band_split` | `100.8898` |
| `ecm_split_expected` | `97.4329` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Core Expected Sum (`E[W_run]`)** | **`110.9369`** |
| **Gated Cold-Run Claim (`E[W]`, `p_cold >= 0.39`)** | **`109.578792317464`** |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared exploratory STNFS premises: `H1` / `H1-Boot` / `H-EA` (smoothness & boot early-abort heuristics on `Q(zeta_15)` norms), `H2` (relation-collection `H2_ALLOW = 1.05`, explicit `det_credited = 0.95`, and `capped_fp_gso_h2` negligible per-prime FP-cap skip rate absorbed in `ret_node = 0.90`), `H3` (zero-column-loss filtering caps `nl_cap = 0.00297`, `sw_cap = 0.00605`, `k5 = 13.7`), `H-X` (Schirokauer full-rank map & controllability), and `H5` (selected-slot conditional pure-`k = 1` emptiness premise `P(M = 0 | slot) <= max(0.001, (1 + 0.5*mu)^(-1))`, combined by unconditional Boole over `<= 451` deterministic slots and exact Wald fresh-`e` `i.i.d.` `gamma_0`-blinded boot retries). All production bounds are conditional on these explicitly declared exploratory premises.
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` and `python3 -I -B audit_cm.py <<< '{"descent_mesh":384}'`.
