# `DESCENT.md` — Pure `k = 1` Zero-Retry Boole + Wald `i.i.d.` Fresh-`e` `b_best2` `P_71` Early-Abort Descent (`lx = 71.0`, `p0_floor = 0.001`)

- **Architecture:** `pure_k1_boole_wald_iid_boot` behind ONE unified hard aborting run meter `T_run = 2^123` u (`p_meter = 0.0002336927 <= 0.00024`, ZERO separate `T_desc`, ZERO blinding pool, ZERO `af_mult`, ZERO internal retries: `p_node_abort = p_repeat_loss = p_descent_fail = 0.0`).
- **DAG (`78` rows, `39` 5-bit bins `[71,76)..[261,266)` x `2` sides, `B_16(R)` exact theta-series cells up to `R = 140`):**
  - Every node uses direct 16D Bareiss-LLL (`99/100`, proved input-uniform `< 2^40` Ir) + deterministic `dsieve16` worst-case pricing (`wc_lll99_det`, `udup_node = 0.4`, `ret_node = 0.9`, `p_ecm = 0.999`, `peel_h = 0.256845`).
  - **`p0_floor = 0.001`**, `p_0(b, s) = max(0.001, (1 + 0.5 * mu(b, s))^-1)`, **`min_mu = 501.92 >= 500`**, **`max_p0 = 0.003969 < 0.004`**.
  - Exact deterministic knapsack over post-`P_71` boot cofactor `<= 2^1950` (`n1 = 216, nA = 186, kA = 14`): **`N_tree_1log_det = 451`** slots/log (`902` across 2 logs), **`Phi_tree = 0.513188 <= 0.52 <= 0.60`**, exact Wald **`E[runs/log] <= 2.0542`**.
- **Fresh-`e` `i.i.d.` `b_best2` + `P_71` Early-Abort Boot:**
  - `LLL_BOOT = 2^48` Ir + `256 * NORMS` (`b_best2` over `72` short vectors `{R_sk[0..15]} U {R_sk[i] +- R_sk[j]}`) + `TDIV` + `P_71` strip wave (`2^28.04` Ir/leaf, `wave_n = 2^61.0`) + `surv_b = 2^-20.69` (computed by `2.0 * psurv_p71(2490.0, 1950.0)`) * `(14 * ECM(186) + MR + IDEAL)`.
  - Exact 1-wave overshoot `over_b = wave_n * trial_b` is explicitly included in `boot_own_1 = trial_b / p_boot_trial + over_b`.
  - 5-quintile `BOOT_BINS = [2534.5, 2536.0, 2537.5, 2538.5, 2541.5]` and `BOOT_LO = 2490.0` are certified across `24/24` live 381-bit `G_T` targets (target `0` recomputed live from scratch in `audit_cm.py`).
  - Live `eval_descent_live`: **`two_boot_exp = 2^105.4709` u**, **`two_tree_det = 2^108.6911` u**, **`E[two_logs] <= 2^108.838119` u** (`384/370` mesh: **`2^108.765011` u <= `2^108.8431` u**).
