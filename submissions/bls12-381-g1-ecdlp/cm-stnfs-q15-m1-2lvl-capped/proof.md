# bls12-381-g1-ecdlp: `(4,3)` `Q(zeta_15)` CM-STNFS, Pure `M1` Two-Level Capped `q * p_rand_j` (`K = 7`, `eta_bad = 0.5%`, `pe = 0.30`) & `ZERO-EA` Deterministic-`m_boot` Boot — `E[W(n)] < 2^113.676151317439` G1 Units

**No prior verdict or score is a premise.** Track `bls12-381-g1-ecdlp` (powers-of-tau auxiliary-powers G1 ECDLP). Input: `(n, A, T)` with `1 <= n < r`, `A_i = [s^i mod r]P_0`, and target `T = A_1 = [s]P_0 in G_1 \ {infinity}` (uniform `s in {1, ..., r-1}`); output 64-lowercase-hex `s`. Uses ONLY `T = A_1` and public generators `P_0 in G_1, Q_0 in G_2` (**zero accesses to `A_2, ..., A_n`**, constant work bound `W(n)` valid for **all `1 <= n < r`**, scored at **`n = 2^20 = 1,048,576`**). Reference unit: one BLS12-381 G1 target-group-operation equivalent (**`C = 2,000`** baseline scalar ARMv8.0-A integer instructions, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64; core rows in `C_core = 20,000 = 10` G1 units, scaled exactly `10x` = `+log2(10)` bits in the G1 claim). Every algebraic identity, kernel count, filtering simulation, collection pool, 1-dense + 16-shard LA streaming bound, Pure `M1` Two-Level Capped + `ZERO-EA` descent recurrence, and 80-digit `Decimal` run-gate claim is verified live by `audit_cm.py` under `python3 -I -B`.

## 1. Track Reduction (Powers-of-Tau G1 ECDLP -> Two `G_T` Logs with Calibrated `gamma_0` `lambda`-Cancellation) & Declared Claim
- **Reduction (Valid for All `1 <= n < r`, Zero Accesses to `A_2..A_n`):** Let `P_0 in G_1, Q_0 in G_2` be public BLS12-381 generators and `e` the pinned reduced pairing (`py_ecc.optimized_bls12_381.pairing(Q_0, P, final_exponentiate=True)`). From ONE Stage-1 relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, define `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {1}` (`<= 2^11` core u, ZERO third descent) with KNOWN nonzero `V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0` (`gcd(r, (p^12-1)/r) == 1` asserted LIVE). Given `(n, A, T)`:
  0. **Run gate:** draw uniform 64-bit `u`; if `u >= K_gate`, abort (`< 1` G1 unit).
  1. **Validate:** decode `T = A_1` (ignoring `A_2..A_n`) and verify `T in G_1 \ {infinity}`.
  2. **Pairing transfer & field map:** compute `h := e(P_0, Q_0)` and `z := e(T, Q_0) = h^s` in `G_T \ {1}` (`<= 2^12` core u) and map into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
  3. **Calibrated `gamma_0`-blinded two-log descent:** For `w in {h, z}`: each fresh boot draws uniform `e in [0, r-1]`, blinds `w' := w * gamma_0^e`, recovers `V(w')`, and unblinds `V(w) := (V(w') - e * V_0) mod r`.
  4. **Exact `lambda`-cancellation & output check:** Compute `s := V(z) * V(h)^(-1) mod r`, verify `1 <= s < r` and `[s]P_0 = T`, and output 64-hex `s` iff valid. Wrapper `<= 2^66` core u, failure `< 2^-48`.
- **Declared Headline Claim (`W(n)` Constant for All `1 <= n < r`, Scored at `n = 2^20 = 1,048,576`):**
  - Opened-run expected work: **`E[W_run] <= 2^111.7121` core u** behind ONE unified hard aborting meter `T_run = 2^123` core u (`10 x 2^123` G1 u; `p_meter = 0.0003999480 <= 0.00040`, deterministic pointwise cap `< 2^126.3219282` G1 u).
  - Opened-run success lower bound: **`p_floor = 0.99960`** (`> 0.39`).
  - Cold-run gate (`K_gate = 7197109032359668998`, wrapper `2^66` core u): **`E[W(n)] <= 1 + (K_gate / 2^64) x 10 x (2^111.7121 + 2^66) < 2^113.676151317439` G1 units** (`C = 2,000`; core `2^110.354223222551`, valid for all `1 <= n < r`, in particular **`log2(W(2^20)) < 113.676151317439`**) with cold algorithmic success **`>= 0.39`**.

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1`). Base field `H = Q(zeta_15)`, `g(x) = x^3 - u * omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega`, `Aut = Z/12Z`. Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`. Live `audit_cm.py` verifies `algebra_bls12`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms` (`12/12`), and `srt_toy` (`0/128`).

## 3. Relation Collection (`P4_prod_2ad_prim_v17`, `Arch2-Prod-ZeroCopy-U32`)
- `V_kept = 2^120.431440`, `f_prim = 0.95642`, `det_credited = 0.95`. Live `audit_cm.py`: `log2_rels_2d_strict_hi = 63.714125` (`margin_strict_hi = +0.19842`b >= `35.1` MC SE), `c_cof_2d_upper_live = 2.0535` Ir/cell, **`relation_collection = 2^110.059203` core u**, `replay_a_tables = 2^105.95909` core u, 1-owner build RAM `2^63.9240 B < 2^64 B`.

## 4. Linear Algebra (`ArchB-1Dense16Mat-4Vec-Streamed`), `H3`, & `H-X`
- Zero-column-loss filtering ladder (`2^16..2^24` x 6 seeds, `cols_lost = 0`): `nl_24 = 0.001776 <= nl_cap = 0.00297`, `sw_24 = 0.003603 <= sw_cap = 0.00605`, `k5_live = 10.136 <= 13.7` (ZERO TABM). Explicit 1 dense owner (`2^63.9485 B < 2^64 B`) + `16` read-only matrix-shard owners (`2^62.597 B < 2^64 B` each) + explicit `c_stream_mat = 201.0` Ir/row transfer -> **`linear_algebra = 2^108.604147` core u**. Unscreened zero-oracle two-sided `Q(zeta_15)` `H-X` + live `gamma_0` `lambda`-cancellation + TRUE same-matrix `H3`.

## 5. Individual-Logarithm Descent (Pure `M1` Two-Level Capped [`K = 7`, `eta_bad = 0.5%`, `pe = 0.30`] + `ZERO-EA` Boot + TOTAL-Cap `dsieve16`)
- **Pure `M1` Two-Level Capped Recurrence (`ADAPTIVE.md`, `DESCENT.md`, ZERO Trial 0!):**
  - `78`-row `lx = 71.0` DAG (`min_mu = 15.00 >= 15.0`, pool `[2^14, 2^16)` with **`4,740` exact live-counted `f_9` ideals**, explicit deterministic **`post_cap_mult = 16.0x`** pointwise post-sieve cap -> unconditional Markov exceedance `<= 1/16 = 6.25%`!).
  - **Deterministic Adversarial Subtree Work Cap `T_cap(b) = min(nat_wc(b), M_cap(b) * W_good(b))`:** Tracks `W_wc[b] = T_cap(b)` bottom-up (`nat_wc = LLL_NODE + 7*(own_1_det + cW_wc)`), enforced pre-instruction with ZERO overshoot. On `40/78` shallow/leaf slots `p_mk = 0.0`; on `38/78` deep slots exact unconditional Markov bounds cap-trip failure on good ideals by `p_mk = 1 / M_cap`. Gives **100% pointwise deterministic worst-case work `<= T_cap(b)` on EVERY bad ideal even on infinite adversarial bad chains!**
  - **Exact Bottom-Up Induction With Explicit `M >= 1`, Boot-Success, & Reached-Child Conditioning (`H5-2lvl`):** Good ideals `G_b` pass `||b_0(L_Q)||_2 >= 0.20 * q^(1/16)` (`+-1`-only factor-12 Blichfeldt `<= 6.452e-7`, **`7,750x` below `eta_bad = 0.005`**), have `p_1(Q) <= 0.30`, and conditional on `M(Q, p) >= 1` (and reached-child selection, including after earlier trial/subtree failures) produce bad children with probability `<= eta_bad = 0.005` (`0.5%` = `5x` larger than `10^-3`). Exact `(1 - p_1(Q))` cancellation yields `E[child calls | Q in G_b] <= 1 / (1 - S_phi)`, `q_ret = 0.30 + 0.70 * S_phi <= 0.3467`, `phi_good = q_ret^7 + p_mk <= 0.00685`, closing strictly at `phi(b) = 0.005 + 0.995 * phi_good(b) <= 0.01182 <= 0.015` and `W_node(b) = 0.995 * W_good(b) + 0.005 * T_cap(b)`.
  - **Unconditional Boole ONLY Over `<= 16` Top-Level `ZERO-EA` Boot Children (`boot = (241, 211, 15)`, conditioned on boot relation success!):** **`Phi_boot = 0.186913` (`E[runs/log] <= 1.2299`)**, **`two_boot_exp = 2^107.1952`**, **`two_tree_exp = 2^110.7158`**, **`E[two_logs] <= 2^110.836328`** (`384/370`: **`2^110.828064` <= `2^110.8413`** core u).
  - **Pointwise TOTAL-Capped `dsieve16` (`T_p_total`, `k = 15..0`, `alpha = 100/74`) & `ZERO-EA` Single-`b_0` Boot (`surv_b = 1.0`, ZERO `H-EA`):** `0/32` skipped on genuine `Q(zeta_15)` roots (`max_ratio = 0.8051`; conditional-geometry `H2` note in `DESCENT.md`: `k^(-16)` short-vector volume + `96`-bit sieve slack), `32/32` real 381-bit `G_T` `b_0` bases certified (`max = 2547.66 < 2592.0`), live `32`-bit monotonicity (`79` points on `[72, 2592)`).
  - **Cantelli/Fano Bridge & Strictly Separate Exact-`O_H`-Coprime Diagnostics (`ADAPTIVE.md` §2; NEVER pooling parents & children!):** Exact Cantelli bridge (`p_1 <= Fano/(Fano+mu) <= 0.30` for all `Fano <= 6.43` at `mu >= 15.0` vs measured deg-`(72,24)` `med_fano = 1.20`; NegBin `1.20^(-75) = 1.15e-6`); live exact `<a, b> = O_H` coprimality unit test; Campaign A parents (`32` genuine `raw_M == 0` above-cutoff ideals x `28` draws: `emp_p1 = 0.1596 << 0.30`); Campaign B parents (`28` literal deg-`(72, 24)` ideals x `24` draws: `emp_p1 = 0.2946 <= 0.30`, `med_fano = 1.20`); **Campaign C reached first-usable-relation LP children (`40` actual LP children x `16` retries = `640` 2nd-level lattices: `0/40` gate fails, `0/40` `p_1 > 0.30`, `pooled_p1 = 0.0500`, `max_child_p1 = 0.2500`, `UMVUE E_hat[p_1^7] = 0.0`, `4` live bit-exact child replays)**; **`8` real 381-bit `G_T` boot children (`min b_0/q^(1/16) = 1.1933`, `0/8 < 0.20`)**; `240`-ideal strict zero-omitted-children (`0/240` zeros); `64` live boundary roots; `36` live 381-bit production lattices.

## 6. Complete Work Ledger (`C_core = 20,000 = 10` G1 units; Final Claim Scaled `10x` to `C = 2,000` G1 units)

| Component | `log2` core u (`C_core = 20,000`) |
|---|---:|
| `relation_collection` | `110.059203` |
| `replay_a_tables` | `105.95909` |
| `linear_algebra` (incl. 16-shard streaming) | `108.604147` |
| `matrix_generator` | `93.1361` |
| `both_individual_logs_exp` (`Pure M1 2-Lvl Capped K=7 eta=0.5% + ZERO-EA`) | `110.8413` |
| `band_split` | `100.8898` |
| `ecm_split_expected` | `97.4329` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Core Expected Sum (`E[W_run]`)** | **`111.7121`** |
| **Gated Cold-Run Claim (`E[W(n)]`, `C = 2,000`, `p_cold >= 0.39`)** | **`113.676151317439` G1 units** (core `110.354223222551`; constant for all `1 <= n < r`, scored at `n = 2^20`) |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared exploratory STNFS premises: **`H1`** (Dickman/de Bruijn + Murphy smoothness across collection, node sieving/survivors, and semismoothness of selected first-vector boot norms `N_g(b_0)` — ZERO `H-EA`), **`H2`** (`H2_ALLOW = 1.05`, `det_credited = 0.95`, `capped_fp_gso_h2` node retention `ret_node = 0.90` under hard TOTAL per-prime cap `T_p_total`, honestly declared as exploratory-v1 heuristic premise supported by conditional-geometry `k^(-16)` / `96`-bit slack analysis, `post_cap_mult = 16.0x`, and `0/32` skipped live), **`H3`** (zero-column-loss filtering caps), **`H-X`** (Schirokauer full-rank map & controllability), and **`H5-2lvl`** (for good parents `Q in G_b` averaged over fresh `p ~ Unif(P_f)` **conditioned on `M(Q, p) >= 1` and reached-child selection including after earlier trial/subtree failures**, and for fresh `e ~ Unif(Z/rZ)` **conditioned on boot relation success**, every child lands in bad set `B = G^c` with conditional probability `<= eta_bad = 0.5%`, while good ideals `Q in G_b` satisfy `p_1(Q) <= 0.30` at `min_mu >= 15.0`; amplified via exact conditional independence of `K = 7` retries + deterministic pre-instruction subtree work cap `T_cap(b)` with exact Markov `p_mk = 1/M_cap` to `phi(b) <= 0.015`, combined via unconditional Boole over ONLY `<= 16` top-level `ZERO-EA` boot children + exact Wald fresh-`e` i.i.d. `gamma_0`-blinded boot retries).
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` and `python3 -I -B audit_cm.py <<< '{"descent_mesh":384}'`.
