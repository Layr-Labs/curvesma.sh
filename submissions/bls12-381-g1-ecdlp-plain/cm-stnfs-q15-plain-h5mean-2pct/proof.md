# BLS12-381 Plain G1 ECDLP (`bls12-381-g1-ecdlp-plain`): `(4,3)` `Q(zeta_15)` CM-STNFS, Zero-Gap `H5-Mean` Telescoping Internal-Retry `M1` (`p_bar = 2.0%`, `K = 6`, `mu_1 >= 20.15`) & `ZERO-EA` Deterministic-`2592` Boot (`E[W_G1] < 2^113.052349273957` G1 units)

## 1. Exact Problem Statement, Target Pairing Reduction, & Claimed Bound
- **Target (`tracks/bls12-381-g1-ecdlp-plain/rules.md`):** Given `T = [x]P_0 in G_1` on BLS12-381 (`P_0` standard generator, prime order `r = u^4 - u^2 + 1` of 255 bits, `u = -0xd201000000010000`, `p` 381 bits, uniform `x in {0, ..., r-1}`, no auxiliary powers), output 64-hex `x` satisfying `[x]P_0 = T`.
- **Trivial Zero-Point Check (`x = 0`):** If `T = O` (point at infinity), output `0` immediately (`<= 1` G1 unit).
- **Public Pairing Transfer `G_1 -> G_T` (`<= 2^12` core u `<< 2^66` wrapper):** Fix the standard public `G_2` generator `Q_0` (`ord(Q_0) = r`). Compute `h := e(P_0, Q_0) in G_T \ {1}` and `z := e(T, Q_0) = h^x in G_T` via the optimal Ate pairing + final exponentiation (two pairings = `<= 2^11` core u). Solve `z = h^x` in `G_T` using the self-contained `(4,3)` `Q(zeta_15)` CM-STNFS algorithm below, then verify `[x]P_0 == T` in `G_1` (`<= 2^9` G1 u).
- **Exact Reference-Unit Conversion (`COST_MODEL.md`):** This track prices 1 G1 unit at **`C = 2,000` baseline scalar ARMv8.0-A integer instructions** (no SIMD, no LSE atomics). All internal core rows in `cm_data.json` are priced in **`C_core = 20,000 = 10` G1 units**; final cold-run expected work is scaled by an exact factor of `10` (`+log2(10) = +3.32192809...` bits):
  - **Core opened-run expectation:** `E[W_run_core] <= 2^111.0885` core u (`exact = 2^111.088294`).
  - **Unified aborting meter:** `T_run = 2^123.0` core u (`10 * 2^123 < 2^126.3219282` G1 u), giving `p_meter = 2^{111.0885 - 123.0} <= 0.000261` and `p_succ_open >= 0.99974 = p_floor`.
  - **Deterministic cold-run gate (`K_gate = 7,196,101,175,052,238,713`):** `K_gate * p_floor / 2^64 >= 0.39` exact (`80`-digit `Decimal` & exact `Fraction` asserted LIVE in `audit_cm.py`).
  - **Final G1 claim:** **`E[W_G1] <= 1 + (K_gate / 2^64) * 10 * (2^111.0885 + 2^66) < 2^113.052349273957` G1 units** (`exact = 113.052349273457`, core `109.730421179070`, **`-1.6675` full bits below `#1 114.719881964658`**!).

## 2. Algebraic Tower, Automorphisms, & Declared Exploratory-v1 Premises
- **Tower (`field-map.md`):** `H = Q(zeta_15)` (`deg 8`, `Phi_15(y)`), `omega = y^5`, `g_3(X) = X^3 - u * omega^2` (`deg 24`), `f_9(X) = X^9 + (omega + 1) X^6 + X^3 - omega` (`deg 72`), quotient `F_{p^12} = F_p[y]/(h_pi(y))[X]/(X^3 - u * omega^2)` (`tau = Frob^1`, exact `<sigma, omega>` order-12 automorphism quotient verified LIVE in `check_algebra_bls12`).
- **Explicitly Declared Exploratory-v1 STNFS Premises (`(H1, H2, H3, H-R, H-SM, H-X, H5-Mean)`):**
  - `H1`: Standard Dickman / Canfield-Erdos-Pomerance smoothness with Murphy-`alpha` for collection (`2D` norm pool), `M1` descent trials, and unconditional worst-case `ZERO-EA` boot (`E[p_boot_actual] >= b(2592.0)`).
  - `H2`: Collection sieve overhead `<= H2_ALLOW = 1.05` and TOTAL-capped `dsieve16` retention `>= ret_node = 0.90` (supported by `0/32` live skips on genuine `Q(zeta_15)` roots + `0/12` live skips on real relation-conditioned primes `p | N_g(v)` + `<= 130` smooth-relation divisor bound).
  - `H3`: Structured-Gaussian zero-column-loss filtering caps (`nl_cap = 0.00297`, `sw_cap = 0.00605`, `>= 1.67x` headroom over `2^16..2^24` ladder).
  - `H-R` & `H-SM` (`RANK-BRIDGE.md` §§1–2): Full-rank 2-core orbit-reduced relation matrix (`rank = N_fb - 1` mod `r`, `H-R`) and non-degenerate 4-column Schirokauer-map compression into the `chi`-eigenspace (`H-SM`, verified live on our unscreened zero-oracle two-sided `Q(zeta_15)` tower).
  - `H-X`: Full-rank right kernel recovery of `(N + 4)`-column Schirokauer-augmented Block Wiedemann matrix (`m = 8192, n = 112`).
  - `H5-Mean`: Zero-Gap Telescoping Internal-Retry premise at modeled single-trial mean `mu_1 >= 20.15` over `K = 6` fresh i.i.d. WITH-REPLACEMENT draws from exact `|P_f| = 1,344` in `[2^12, 2^14)` behind a **`4.0x` pointwise survivor cap + deterministic single-factor ECM (`W_trial <= own_1` pointwise)**: for every reached ideal distribution `D` (any boot child under `e ~ Unif(Z/rZ)`, or any child `c_{j,r}` conditioned on `{Y_j = 1}` of a reached parent), **`E_{Q ~ D}[p_eff(Q)] <= p_bar = 0.02`** (`2.0%`), where `p_eff(Q) = 1` if `||b_0(L_Q)||_2 < 0.20 * q^(1/16)` else `p(Q) = Pr_{p ~ Unif(P_f)}(trial (Q, p) empty | Q)`. **NO good/bad ideal split, ZERO bad-parent child-work gap, ZERO `H-EA` premise, and ZERO boot-norm distribution premise!**

## 3. Relation Collection (`Arch2-Prod-ZeroCopy-U32`, `B_f = 2^71.29`, `B_g = 2^71.69`, `U' = 2^63.3174`)
- **2D Norm Pool & Finite-Region Margin:** `V_kept = 2^120.431440`, `f_prim = 0.95642`, explicit ECM retention credit **`det_credited = 0.95`** (`-0.0740`b penalty).
  - Live `audit_cm.py` values: `log2_rels_2d_mid = 63.859129`, `log2_rels_2d_rethi = 63.853894` (`margin_net_rethi = +0.33819`b), **`log2_rels_2d_strict_hi = 63.714125` (`margin_strict_hi = +0.19842`b >= `35.1` MC SE)**, `c_cof_2d_upper_live = 2.0535` Ir/cell.
- **U32 Kernels & Work Pricing (`KERNELS.md`):** `p_m4 = 56.5`, `p_lo = 59.0`, `p_f = 59.0`, `p_g = 64.0`, `c_hit_lo_m4 = 3.80`, `c_hit_lo_u = 5.80`, `c_hit_hi = 16.35`, `c_scan_zero = 0.95`, `H2_ALLOW = 1.05` -> **`relation_collection = 2^110.059203` core u**, `replay_a_tables = 2^105.95909` core u, 1-owner build memory `2^63.9240 B < 2^64 B`.
- **Live Murphy-Alpha & Real-FB Joint Collection (`unitdomain.py`):** Exact `p = 3` Hensel tree + unramified `p <= 10,000` (`alpha_sum_10k <= 4.02` charged); live `1,200`-pair exact deg-`(72, 24)` `U_toy = -14` Voronoi/`D_12` joint collection.

## 4. Linear Algebra (`ArchB-1Dense16Mat-4Vec-Streamed`), `H3`, `H-R`, `H-SM`, & `H-X`
- **Zero-Column-Loss Filtering Ladder (`2^16..2^24` x 6 seeds, `cols_lost = 0` at every scale):** `nl_24 = 0.001776` vs **`nl_cap = 0.00297` (`1.672x` headroom)**; `sw_24 = 0.003603` vs **`sw_cap = 0.00605` (`1.679x` headroom)**; `k5_live = 10.136 <= k5_charged = 13.7` (`hex6`, ZERO TABM: `tab_ratio_charged = 1.0`).
- **Explicit 1-Dense + 16-Matrix-Shard Streaming Block Wiedemann (`m = 8192, n = 112, n_sm = 4`, `STORAGE.md` §2):** `log2(N) = 54.922051`, `log2(N_sw) = 55.948523`; dense double-buffered owner `2 * N_sw * 128 B = 2^63.9485 B < 2^64 B`; read-only sparse matrix partitioned across **`16` matrix-shard owners (`2^62.597 B < 2^64 B` each, asserted LIVE in `audit_cm.py`)**; explicit cross-owner streaming transfer **`c_stream_mat = 201.0` Ir/row** priced inside `c_prod` -> **`linear_algebra = 2^108.604147` core u**, `matrix_generator = 2^93.1361` core u.
- **Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X` + Live `gamma_0`-Blinding/Unblinding `lambda`-Cancellation + TRUE Same-Matrix `H3` (`w8_verifiers.py`):** `34/34` `f`-side and `34/34` `g`-side exact post-hoc BSGS hits with `cols_2core_lost = 0`, live zero-BSGS `gamma0_blind_unblind_lambda_cancel_ok = True`, matched-scale genuine-vs-synthetic calibration, TRUE same-underlying-matrix full-incidence vs truncated-proxy, `unplanted_two_sided_bw`, and `filtsim_and_hx_live`.

## 5. Individual-Logarithm Descent (Zero-Gap `H5-Mean` Telescoping Internal-Retry `M1` + `ZERO-EA` Deterministic-`2592` Boot)
- **Zero-Gap `H5-Mean` Telescoping Internal-Retry `M1` Pure `q -> q * p_rand_j` (`ADAPTIVE.md`, `DESCENT.md`):**
  - `78`-row `lx = 71.0` DAG (`max_R = 86`, **`min_mu = 20.15 >= 20.0`**, `p_node_abort = p_repeat_loss = p_descent_fail = 0.0`).
  - Up to **`K = 6`** fresh i.i.d. WITH-REPLACEMENT known-log `f`-side prime ideals `p_j ~ Unif(P_f)` (`|P_f| = 1,344` in `[2^12, 2^14)`, certified LIVE in `m1_qprand.py`), retrying on fresh `p_{j+1}` whenever trial `j` is empty OR any recursive child fails!
  - Behind the **`4.0x` pointwise survivor cap + deterministic single-factor ECM (`W_trial <= own_1` pointwise)** and declared **`H5-Mean`** (`E_{Q ~ D}[p_eff(Q)] <= p_bar = 0.02`, NO `G/B` split!), the exact disjoint indicator identity `X_{j+1} = X_j * E_j + F_j = (X_j - Y_j) + F_j` yields `a_{j+1} <= cPhi * a_j + (1 - cPhi) * p_bar` and `(1 - cPhi) * sum_j b_j <= a_1 - a_{K+1} <= 1`, giving:
    - **`phi(s) = cPhi^6 + (1 - cPhi^6) * 0.02 <= 0.020063`** across ALL 78 rows,
    - **Exact zero-gap node work `W_node = LLL_NODE + (((1 - cPhi^6)/(1 - cPhi)) * 0.98 + 6 * 0.02) * own_1 + (1 / (1 - cPhi)) * cW`** (ZERO bad-parent child-work gap!),
    - Boot root (`<= 17` children): **`Phi_boot = 0.341066 <= 0.35`**, exact Wald **`E[runs/log] = 1.517603`**, **`two_tree_exp = 2^109.0138`** core u.
  - **Pointwise TOTAL-Capped Per-Prime FP (`k = 15..0`, `alpha = 100/74`):** Hard per-prime cap `T_p_total` bounds BOTH internal FP traversal (`k = 15..1`) AND `340 * hits` leaf processing (`k = 0`). Supported LIVE inside real 381-bit composite `L_{(Q, p_j)}` bases at exact priced `(9/10)`-LLL on `32` genuine `Q(zeta_15)` roots (`0/32` skipped, `max_ratio = 0.8143 < 1.0`) AND `12` real relation-conditioned primes `p | N_g(v)` (`0/12` skipped, `cond_div_rel_max_ratio = 0.8266 < 1.0`) + `<= 130` prime-factor (`> 2^{16}`) union bound (`ret_node = 0.90` under `H2`).
- **`ZERO-EA` Deterministic-`2592` Boot (`n1 = 246, nA = 211, kA = 16, tau = 2592.0, surv_b = 1.0`, `LLL_BOOT = 2^46.5 > 2^46.45` proved):**
  - Single first `(99/100)`-LLL vector `b_0 = R_sk[0]` with unconditional algebraic Lovasz bound `M <= 2591.90 < 2592.0` (`ret_boot = 1.0`; **ZERO `H-EA`, ZERO `318`-bin mediant, ZERO `BOOT_BINS`, ZERO `BOOT_LO`, ZERO `p_band`, ZERO `b_best2`, ZERO toy-EA**).
  - Full sequential Phase-2 deterministic single-factor ECM+MR+`IDEAL` factorization paid on **100% of boot trials (`surv_b = 1.0`)** (`over_b <= 2^93.2` core u). Constant `c_trial = c_pre + ecm_post_b` + live 4-bit monotonicity `b(M) >= b(2592.0)` on `[72, 2592]` directly give **`two_boot_exp = 2^107.1064`** core u and **`E[two_logs] <= 2^109.354735`** core u (`384/370` mesh = **`2^109.332680`** `<= 2^109.3597` core u ledger ceiling!). Live `verify_gt_boot_20` recomputes target `0` from scratch (`n_scratch = 1`) and certifies all `32/32` real 381-bit `G_T` bases (`max = 2547.66 < 2592.0`).
- **Scoped Live `M1` & `H5-Mean` Diagnostics (`ADAPTIVE.md` §2):** Live `m1_qprand_live` (`|P_f| = 1,344`, strict `+-1` `theta16`/Blichfeldt `1.631e-8` over full `[71, 265]`, `4` real 381-bit `prod_sq24` + `4` real parent-selected `prod_enc12` `||b_0||` & `L_{q,p_j}` geometry/headroom, planted adversarial ideal `0.1365 < 0.20`, **`300`-parent (`150 g + 150 f`) + `127`-reached-child + `60`-forced-retry-`j=2` (`0/28` retry-selected child empties!) exact `O_H`-coprime FULL END-TO-END RECURSIVE DESCENT campaign (`B_fb = 32768 < LP_max = 65536 < q in [2^16, 2^17)`, `mean_M = 7.89` [`2.55x` below production `min_mu = 20.15`!]): Parent Trial-1 empty = `1 / 300 = 0.00333` [exact CP 95% UB `0.01571 < 0.02`, 98% UB `0.01929 < 0.02 = p_bar` directly on `E[p_eff(Q)]`!], Reached-Child Trial-1 empty = `0 / 127 = 0.0000`, Combined `427` ideals Trial-1 empty = `1 / 427 = 0.00234` [`CP95 = 0.01106`], `0 / 300` full end-to-end recursive descent failures, and `3` live bit-exact parent+child end-to-end replays including `idx = 202` Trial-1-empty -> Trial-2-retry recovery!**, AND **literal degree-`(72, 24)` (`f_9, g_3`) exact `O_H`-coprime `q * p_rand_j` campaign [`28` ideals x `24` draws = `672` lattices, `emp_p1 = 0.2946`, `med_fano = 1.20` at `mean_M = 1.53`, + live bit-exact `2,496`-pt replay!]**); plus strict zero-omitted-children **real recovered 2-core-FB** 240-ideal 24-bit campaign (`0/240` zeros), `64` live boundary roots (`O/E_p0 = 0.3698`), exact deg-`(72, 24)` `120`-ideal dispersion + `hi100` (`0/100` zeros), `120` enc vs `120` rand + `24 + 12` live 381-bit lattices (`>= 10x` headroom), and naive-PGF negative control.

## 6. Complete Work Ledger (`C_core = 20,000 = 10` G1 units; Final Claim Scaled `10x` to `C = 2,000` G1 units)

| Component | `log2` core u (`C_core = 20,000`) |
|---|---:|
| `relation_collection` | `110.059203` |
| `replay_a_tables` | `105.95909` |
| `linear_algebra` (incl. 16-shard streaming) | `108.604147` |
| `matrix_generator` | `93.1361` |
| `both_individual_logs_exp` (`H5-Mean M1 p_bar=2%, K=6, ZERO-EA`) | `109.3597` |
| `band_split` | `100.8898` |
| `ecm_split_expected` | `97.4329` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_h3_u` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `subproduct_tree_sieve` | `75.0` |
| **Opened-run expectation `exp_run_claim`** | **`111.0885`** (`exact = 111.088294`) |
| **Unified aborting meter `T_run`** | **`123.0`** (`p_meter <= 0.000261`) |
| **Core gated claim (`C_core = 20,000`)** | **`109.730421179070`** (`K_gate = 7,196,101,175,052,238,713`, `p_floor = 0.99974`) |
| **Plain G1 ECDLP Claim (`C = 2,000`, exact `10x`)** | **`113.052349273957`** (`-1.6675` bits below `#1 114.719881964658`!) |
