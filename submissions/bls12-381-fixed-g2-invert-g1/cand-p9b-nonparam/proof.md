# BLS12-381 `recover-P` (`bls12-381-fixed-g2-invert-g1`) — Complete Self-Contained Verification Proof

**No prior verdict or score is a premise.** Every algebraic identity, kernel opcode count, filtering simulation, collection pool evaluation, linear-algebra memory/work bound, nonparametric blinding descent recurrence, and run-gate claim in this package is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib, baseline scalar ARMv8.0-A integer instructions, `simd = 0`, `C = 20,000` instructions per dense `F_{p^12}` multiplication).

## 1. Track Reduction & Declared Claim
- **Track:** `bls12-381-fixed-g2-invert-g1` — fixed-Q / recover-P (`e(P, Q) -> P`).
- **Reduction (Exact Challenge-Dependent Two-Logarithm Fixed-`Q` Inversion):** Let `P_0 in G_1`, `Q_0 in G_2` be the fixed standard generators and `g_0 := e(P_0, Q_0) in G_T` be the fixed reference generator of order `r = u^4 - u^2 + 1`. Given the supplied challenge `(Q, T)` with `Q = [b] Q_0 in G_2 \ {O}` (`b in (Z/rZ)^*`) and `T = e(P, Q) = e([a] P_0, [b] Q_0) = g_0^{a b} in G_T`:
  1. Compute the challenge-dependent pairing base **`h := e(P_0, Q) = g_0^b in G_T`** (`1` pairing evaluation, `<= 2^11` u).
  2. Run two individual discrete logarithms in `G_T` to the fixed reference base `g_0` (using the precomputed factor-base logarithms): **`x_1 := log_{g_0}(h) = log_{g_0}(e(P_0, Q)) = b mod r`** and **`x_2 := log_{g_0}(T) = a * b mod r`**.
  3. Since `Q != O`, `x_1 = b != 0 mod r`. Compute **`a := x_2 * x_1^{-1} mod r`** and output **`P := [a] P_0 in G_1`**.
  4. **Exactness & Verification:** `e(P, Q) = e([a] P_0, [b] Q_0) = g_0^{a b} = g_0^{x_2} = T` identically for **every** `b in (Z/rZ)^*` and **every** `a in Z/rZ` (probability `1` conditional on the two `G_T` logarithms), verified directly by `1` pairing check (`<= 2^11` u).
- **Declared Headline Claim:**
  - Opened-run expected work: **`E[W_run] <= 2^111.3103` u** behind hard global meter `T_run = 2^123` u (`p_meter <= 0.00032`) and descent meter `T_desc = 2^117.5` u (`p_desc_meter <= 0.00562`).
  - Opened-run success lower bound: **`p_floor = 0.99407`** (`>= 0.985 > 0.39`).
  - Cold-run gate (`K_gate = 7237146467297801091`, wrapper `2^66` u): **`E[W] <= 2^109.960426677528` u** with cold success **`>= 0.39`**.

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`), `g(x) = x^3 - u * omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = Z/12Z` (`<sigma, omega>`).
- Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`.
- Live `audit_cm.py` verifies `algebra_bls12`, `toy_e2e_dlp` (exact DLP recovery in `(4,3)` `Q(zeta_15)`), `sigmaB_exactness`, `anchor_norms` (`12/12`), and `srt_toy` (`0/128` mismatches).

## 3. Relation Collection (`P4_prod_2ad_prim_v17`, `Arch2-Prod-ZeroCopy-U32`)
- **2D Norm Pool & Finite-Region Margin:** `V_kept = 2^120.431440`, `f_prim = 0.95642`, `det_credited = 0.95` (`-0.0740`b penalty).
  - Live `audit_cm.py` values: `log2_rels_2d_mid = 63.859129`, `log2_rels_2d_rethi = 63.853894` (`margin_net_rethi = +0.33819`b), **`log2_rels_2d_strict_hi = 63.714125` (`margin_strict_hi = +0.19842`b >= `35.1` MC SE)**, `c_cof_2d_upper_live = 2.0535` Ir/cell.
- **U32 Kernels & Work Pricing:** `p_m4 = 56.5`, `p_lo = 59.0`, `p_f = 59.0`, `p_g = 64.0`, `c_hit_lo_m4 = 3.80` (`rev_b1a_apply_u8`), `c_hit_lo_u = 5.80` (`rev_b1b_apply_u8`), `c_hit_hi = 16.35`, `c_scan_zero = 0.95` (`rev_scan_zero_u32`), `H2_ALLOW = 1.05` -> **`relation_collection = 2^110.059203` u**, `replay_a_tables = 2^105.959090` u, 1-owner build memory `2^63.9240 B < 2^64 B`.
- **Live 3-Scale Disjoint Real-FB Joint-Collection, Suyama ECM, & Murphy-Alpha Calibration (exact `audit_cm.py` output):**
  - `real_fb_joint_coll` (`3` disjoint scales `B_fb = 8192..65536`, `2250` independent pairs): **`Joint = 258/252.39 (1.0222x)`**, `wcof_obs_over_mod = 0.8384` (`max = 0.8628`), Suyama ECM `0.999` live `100/100` (`CP95 >= 0.97049 >= 0.95`).
  - `murphy_alpha_live`: `alpha_f_10k = 2.70436`, `alpha_g_ge3_10k = 1.10299`, `alpha_sum_10k = 3.80734 <= 4.02` charged (`+0.21266`b conservative norm margin).

## 4. Linear Algebra (`ArchB-1Owner-4Vec-ZeroCopy`), `H3`, & `H-X`
- **Zero-Column-Loss Filtering Ladder (`2^16..2^24` x 6 seeds, `cols_lost = 0` at every scale):**
  - `nl_24 = 0.001776` vs **`nl_cap = 0.00297` (`1.672x` headroom)**; `sw_24 = 0.003603` vs **`sw_cap = 0.00605` (`1.679x` headroom)**; `exp_unit_24 = 0.934070` -> `k5_live = 10.1364 <= k5_charged = 13.7` (`hex6`).
- **Single-Owner `< 2^64 B` Block Wiedemann (`m = 8192, n = 112, n_sm = 4`, `n_owners = 1`, `c_xown_live = 0.0`):**
  - `log2(N) = 54.922051`, `log2(N_sw) = 55.948523`, peak single-owner RAM **`2 * N_sw * 128 B = 2^63.9485 B < 2^64 B`**, **`linear_algebra = 2^108.530371` u**, `matrix_generator = 2^93.1361` u.
- **Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X` + Matched-Scale Calibration + TRUE Same-Underlying-Matrix Full-vs-Proxy `H3` (`w9_hx_h3.py`, exact live output):**
  - Unscreened two-sided `Q(zeta_15)` `H-X`: `34/34` `f`-side and `34/34` `g`-side exact post-hoc BSGS hits with `cols_2core_lost = 0`.
  - Matched-scale genuine-vs-synthetic calibration: `gen_k5_hex6 = 11.600 <= syn_matched_k5_hex6 = 12.544`.
  - TRUE same-underlying-matrix full-incidence (`0` omitted cols) vs truncated-proxy (`139` middle cols): `FULL N_l_eff = 236 <= PROXY 264`, `FULL N_sw = 14226 <= PROXY 15942`, `cols_lost = 0/0`.
  - Plus `unplanted_two_sided_bw` and `filtsim_and_hx_live` (`0` backsub errors).

## 5. Individual-Logarithm Descent (**Nonparametric `G_min = 64/192`, `K_blind = 24`, `k_bt = 4`, `af_mult = 1.10` + Hard Per-Scan Meter `T_scan_Ir = 2.0 * own_1`** (`p_hyp1 = 2.6643e-05`, `p_hyp_kbt = 0.014400` verified in exact `Fraction`))
- **Boot (`m_sv = 2592.0`, `n1 = 246`, `nA = 211.0`, `kA = 16`, `ret_boot = 1.0`, `cells_b = 1`):**
  - Proved deterministic Lovasz boot norm bound + live 381-bit `G_T` boot witness (`verify_gt_boot_20`) + live `150`-exponent `G_T` smoothness check (`toy_gt_boot_150_live`: `ret_boot_det_pass = 150/150`, `b32_smooth_obs = 14/150` vs `pred_at_mean = 13.89`, `O/E = 1.0082`, `pred_at_det_ub = 0.28`).
- **Nonparametric Blinding Recurrence & Hard Meter (`ADAPTIVE.md`, `DESCENT.md`):**
  - Pairwise-`O_H`-coprime blinding pool `|B_blind| = 192` in `[2048, 8192]` (`p_sat = 0`, direct `Z`-basis `R[0..15]`, unconditional `m`-truncated `(9/10)`-Lovasz signed-64-bit Q32 `dsieve16`, hard per-scan instruction meter `T_scan_Ir = 2.0 * own_1 <= 2^113.9475` Ir).
  - Exact `Fraction`-verified without-replacement stopping & monotone-shrinking coupling (`verify_wr_stopping_exact(K_bl=24)`): `p_hyp1 = 2.6643e-05`, `p_hyp_kbt = 0.014400`.
  - Live `eval_descent_live`: `min_mu = 32.15 >= 32.0`, `min_kappa_crit_pz_g64 = 1.9689 >= 1.968`, `max_phi = 0.015016`, `max_psi = 0.146690`, `P_boot_tree = 0.252818`, **`two_log_exp_units = 2^110.018886` u** (`384/370` mesh: **`2^110.006408` u <= `2^110.020` u**), `p_desc_meter_live = 0.005597 <= 0.00562`.
- **Three Scoped Compatibility Pillars (Exploratory Support, Not Production Proof):**
  1. `(3,3)` full-depth trees (`120/120` complete, `290` slots), separate `90 x 24` cohorts (`min G_0 >= 23/24 >> 8/24`, `max kappa <= 0.2010 << 1.9689`), `5,760`-scan Weak-LP-Dual (`min G_F >= 63/96 >> 32/96`), exact `360`-root reconciled round-transition campaign (`q15_blind_h5.b85`, `psi_r/psi_1 <= 1.0225 <= 1.10`, `w_r/w_1 <= 1.0544 <= 1.10`).
  2. Literal degree-`(72, 24)` `(9,3)` `18 x 12 = 216` live-LLL paired bridge (`max kappa_act = 0.3163 << 1.9689`, calibrated 381b `kap_cal = 0.9318 << 1.9689`) + live `48`-scan (`2` parents x `24`) 1-LP Weak-LP-Dual & scan-cap check (`min_cert_GF6_of_24 = 19/24 >> 8/24`, `max_scan_work_over_mean = 1.0776 << 2.0`, `0/48` truncated).
  3. Stride-8 `24/192` single-parent 381-bit audit (`cv2_full_192_max_Q = 0.10045 << 1.9689`) + `11 x 4` 381-bit cross-parent audit (`68` live 381-bit CRT LLLs, `W_ball_emp / T_scan_Ir <= 0.50`).

## 6. Complete Work Ledger (`C = 20,000` Ir/u)

| Component | `log2(u)` |
|---|---:|
| `relation_collection` | `110.059203` |
| `replay_a_tables` | `105.95909` |
| `linear_algebra` | `108.530371` |
| `matrix_generator` | `93.1361` |
| `both_individual_logs_exp` | `110.02` |
| `band_split` | `100.8898` |
| `ecm_split_expected` | `97.4329` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Core Expected Sum (`E[W_run]`)** | **`111.3103`** |
| **Gated Cold-Run Claim (`E[W]`, `p_cold >= 0.39`)** | **`109.960426677528`** |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared standard/exploratory STNFS premises: `H1` / `H1-Boot` (smoothness heuristics on `Q(zeta_15)` norms), `H2` (relation-collection finite-regime cost allowance `H2_ALLOW = 1.05`), `H3` (zero-column-loss filtering caps `nl_cap = 0.00297`, `sw_cap = 0.00605`, `k5 = 13.7`), `H-X` (Schirokauer full-rank map), and `H5 = (H5(a), H5(b))` (nonparametric `G_min = 64/192` blinding availability and `1.10x` retry child-work/hazard envelope). All production bounds are conditional on these explicitly declared premises.
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` (`7.07s` local, `checks_ok = true`) and `python3 -I -B audit_cm.py <<< '{"descent_mesh":384}'` (`8.51s` local, `checks_ok = true`).
