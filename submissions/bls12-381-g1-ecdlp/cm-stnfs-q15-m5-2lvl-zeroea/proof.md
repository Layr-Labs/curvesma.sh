# bls12-381-g1-ecdlp: `(4,3)` `Q(zeta_15)` CM-STNFS, Two-Level Internal-Retry `M5 Hybrid` (`eta_bad = 1.0%`, `K = 8`) & `ZERO-EA` Deterministic-`m_boot` Boot — `E[W(n)] < 2^113.364850291911` G1 Units

**No prior verdict or score is a premise.** Track `bls12-381-g1-ecdlp` (powers-of-tau auxiliary-powers G1 ECDLP). Input: `(n, A, T)` with `1 <= n < r`, `A_i = [s^i mod r]P_0`, and target `T = A_1 = [s]P_0 in G_1 \ {infinity}` (uniform `s in {1, ..., r-1}`); output 64-lowercase-hex `s`. Uses ONLY `T = A_1` and public generators `P_0 in G_1, Q_0 in G_2` (**zero accesses to `A_2, ..., A_n`**, constant work bound `W(n)` valid for **all `1 <= n < r`**, scored at **`n = 2^20 = 1,048,576`**). Reference unit: one BLS12-381 G1 target-group-operation equivalent (**`C = 2,000`** baseline scalar ARMv8.0-A integer instructions, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64; core rows in `C_core = 20,000 = 10` G1 units, scaled exactly `10x` = `+log2(10)` bits in the G1 claim). Every algebraic identity, kernel count, filtering simulation, collection pool, 1-dense + 16-shard LA streaming bound, Two-Level Internal-Retry + `ZERO-EA` descent recurrence, and 80-digit `Decimal` run-gate claim is verified live by `audit_cm.py` under `python3 -I -B`.

## 1. Track Reduction (Powers-of-Tau G1 ECDLP -> Two `G_T` Logs with Calibrated `gamma_0` `lambda`-Cancellation) & Declared Claim
- **Reduction (Valid for All `1 <= n < r`, Zero Accesses to `A_2..A_n`):** Let `P_0 in G_1, Q_0 in G_2` be public BLS12-381 generators and `e` the pinned reduced pairing (`py_ecc.optimized_bls12_381.pairing(Q_0, P, final_exponentiate=True)`). From ONE Stage-1 relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, define `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {1}` (`<= 2^11` core u, ZERO third descent) with KNOWN nonzero `V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0` (`gcd(r, (p^12-1)/r) == 1` asserted LIVE). Given `(n, A, T)`:
  0. **Run gate:** draw uniform 64-bit `u`; if `u >= K_gate`, abort (`< 1` G1 unit).
  1. **Validate:** decode `T = A_1` (ignoring `A_2..A_n`) and verify `T in G_1 \ {infinity}`.
  2. **Pairing transfer & field map:** compute `h := e(P_0, Q_0)` and `z := e(T, Q_0) = h^s` in `G_T \ {1}` (`<= 2^12` core u) and map into `F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`).
  3. **Calibrated `gamma_0`-blinded two-log descent:** For `w in {h, z}`: each fresh boot draws uniform `e in [0, r-1]`, blinds `w' := w * gamma_0^e`, recovers `V(w')`, and unblinds `V(w) := (V(w') - e * V_0) mod r`.
  4. **Exact `lambda`-cancellation & output check:** Compute `s := V(z) * V(h)^(-1) mod r`, verify `1 <= s < r` and `[s]P_0 = T`, and output 64-hex `s` iff valid. Wrapper `<= 2^66` core u, failure `< 2^-48`.
- **Declared Headline Claim (`W(n)` Constant for All `1 <= n < r`, Scored at `n = 2^20 = 1,048,576`):**
  - Opened-run expected work: **`E[W_run] <= 2^111.4009` core u** behind ONE unified hard aborting meter `T_run = 2^123` core u (`10 x 2^123` G1 u; `p_meter = 0.0003223465 <= 0.00033`, deterministic pointwise cap `< 2^126.3219282` G1 u).
  - Opened-run success lower bound: **`p_floor = 0.99967`** (`> 0.39`).
  - Cold-run gate (`K_gate = 7196605068419303501`, wrapper `2^66` core u): **`E[W(n)] <= 1 + (K_gate / 2^64) x 10 x (2^111.4009 + 2^66) < 2^113.364850291911` G1 units** (`C = 2,000`; core `2^110.042922197024`, valid for all `1 <= n < r`, in particular **`log2(W(2^20)) < 113.364850291911`**) with cold algorithmic success **`>= 0.39`**.

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1`). Base field `H = Q(zeta_15)`, `g(x) = x^3 - u * omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega`, `Aut = Z/12Z`. Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`. Live `audit_cm.py` verifies `algebra_bls12`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms` (`12/12`), and `srt_toy` (`0/128`).

## 3. Relation Collection (`P4_prod_2ad_prim_v17`, `Arch2-Prod-ZeroCopy-U32`)
- `V_kept = 2^120.431440`, `f_prim = 0.95642`, `det_credited = 0.95`. Live `audit_cm.py`: `log2_rels_2d_strict_hi = 63.714125` (`margin_strict_hi = +0.19842`b >= `35.1` MC SE), `c_cof_2d_upper_live = 2.0535` Ir/cell, **`relation_collection = 2^110.059203` core u**, `replay_a_tables = 2^105.95909` core u, 1-owner build RAM `2^63.9240 B < 2^64 B`.

## 4. Linear Algebra (`ArchB-1Dense16Mat-4Vec-Streamed`), `H3`, & `H-X`
- Zero-column-loss filtering ladder (`2^16..2^24` x 6 seeds, `cols_lost = 0`): `nl_24 = 0.001776 <= nl_cap = 0.00297`, `sw_24 = 0.003603 <= sw_cap = 0.00605`, `k5_live = 10.136 <= 13.7` (ZERO TABM). Explicit 1 dense owner (`2^63.9485 B < 2^64 B`) + `16` read-only matrix-shard owners (`2^62.597 B < 2^64 B` each) + explicit `c_stream_mat = 201.0` Ir/row transfer -> **`linear_algebra = 2^108.604147` core u**. Unscreened zero-oracle two-sided `Q(zeta_15)` `H-X` + live `gamma_0` `lambda`-cancellation + TRUE same-matrix `H3`.

## 5. Individual-Logarithm Descent (Two-Level Internal-Retry `M5 Hybrid` [`eta_bad = 1.0%`, `K = 8`] + `ZERO-EA` Boot + TOTAL-Cap `dsieve16`)
- **Two-Level Conditionally-I.I.D. Internal-Retry `M5 Hybrid` (`ADAPTIVE.md`, `DESCENT.md`):**
  - `78`-row `lx = 71.0` DAG (`min_mu0 = 20.05 >= 20.0`, `min_mu1 = 20.06 >= 20.0`, explicit deterministic **`post_cap_mult = 4.0x`** pointwise post-sieve cap).
  - **Explicit `eta_bad = 1.0%` (`0.01`) Adversarial/Dead-Ideal Tolerance at EVERY Level:** Replaces flat `1,911`-slot Boole (`2e-5` tail) with internal per-node retries: if Trial 0 (`raw L_Q`) is empty OR any child fails, node `Q` draws up to `K = 8` conditionally-i.i.d. with-replacement `f`-side primes `p_j` from `[2^16, 2^18)` (`16,272` `O_H` ideals, 100% known logs) and retries fresh composite lattices `L_(Q, p_j)` and their children.
  - **Exact Pre-Draw Inductive Closure (`psi(Q) <= 0.00201` on good `G_b`, pre-draw mean `phi(b) <= 0.011991 <= 0.0120`):** For any fixed ideal `Q`, its subtree failure probability `psi(Q)` is deterministic over `Q`'s own fresh draws. Under declared exploratory premise **`H5-2lvl`** (pre-draw probability over a good parent's fresh `p ~ Unif(P_f)` — or the boot's fresh `e ~ Unif(Z/rZ)` — BEFORE revealing children that child `c_r` lands in bad set `B = G^c` is `<= eta_bad = 0.01`, **`10x` larger than `10^-3`** and **`15,500x` larger than our Blichfeldt bound `6.452e-7`**, plus `p_0(Q), p_1(Q) <= 0.35` on good `Q in G`), each fresh retry on fixed good `Q` fails with probability `p_(fail,1)(Q) <= 0.35 + S_phi <= 0.4602` (Boole over ONLY one relation's `<= m_ch` children!), giving `psi(Q) <= q_ret^8 <= 0.00201` and pre-draw mean `phi(b) = 0.01 + 0.99 * q_ret^8 <= 0.0120`!
  - **Exact Expected Work (Zero Correlation Gap, 100% Adversarial Bad-Ideal Pricing):** Pointwise-capped `own_0_det, own_1_det`, exact good-ideal child-call multiplier `(1 - 0.35)/(1 - 0.35 - S_phi)`, and worst-case adversarial `0.01 * (own_0_det + 8 * own_1_det + 9 * W_ch)` explicitly charged!
  - **Unconditional Boole ONLY Over `<= 16` Top-Level `ZERO-EA` Boot Children (`boot = (246, 211, 15)`):** **`Phi_boot = 0.188774` (`E[runs/log] <= 1.2327`)**, **`two_boot_exp = 2^106.8518`**, **`two_tree_exp = 2^110.0553`**, **`E[two_logs] <= 2^110.203972`** (`384/370`: **`2^110.193665` <= `2^110.209`** core u).
  - **Pointwise TOTAL-Capped `dsieve16` (`T_p_total`, `k = 15..0`, `alpha = 100/74`) & `ZERO-EA` Single-`b_0` Boot (`surv_b = 1.0`, ZERO `H-EA`):** `0/32` skipped on genuine `Q(zeta_15)` roots (`max_ratio = 0.8051`), `32/32` real 381-bit `G_T` `b_0` bases certified (`max = 2547.66 < 2592.0`), live `32`-bit monotonicity (`79` points on `[72, 2592)`).
  - **Scoped Exact-`O_H`-Coprime Diagnostics (`ADAPTIVE.md` §2):** Exact `<a, b> = O_H` coprimality unit-tested live (`(2-y)*(a_0,b_0)` caught despite `gcd = 1`); Campaign A (`32` genuine `raw_M == 0` above-cutoff ideals x `28` draws = `896` trials, `pool = 328`: **ALL `32/32` have `hat_p1 <= 0.3214 < 0.35`** even at `mean_hits = 6.21`, `3.2x` below `min_mu1 = 20.06`, with `4` live `raw_M == 0` checks + `3` live bit-exact fallback replays); Campaign B (`28` literal deg-`(72, 24)` ideals x `24` draws = `672` lattices, `emp_p1 = 0.2946 <= 0.35`, `med_fano = 1.20`, live bit-exact replay); plus `240`-ideal strict zero-omitted-children (`0/240` zeros!), `64` live boundary roots, and `36` live 381-bit production lattices.

## 6. Complete Work Ledger (`C_core = 20,000 = 10` G1 units; Final Claim Scaled `10x` to `C = 2,000` G1 units)

| Component | `log2` core u (`C_core = 20,000`) |
|---|---:|
| `relation_collection` | `110.059203` |
| `replay_a_tables` | `105.95909` |
| `linear_algebra` (incl. 16-shard streaming) | `108.604147` |
| `matrix_generator` | `93.1361` |
| `both_individual_logs_exp` (`Two-Level M5 K=8 eta=1% + ZERO-EA`) | `110.209` |
| `band_split` | `100.8898` |
| `ecm_split_expected` | `97.4329` |
| `canonical_a_and_sort` | `86.0` |
| `schirokauer_4col_crt` | `84.0` |
| `filtering_and_misc` | `80.0` |
| `pg_pf_product_trees` | `75.5` |
| `factor_base_roots_ideals` | `75.0` |
| `numerical_certificate` | `50.0` |
| `verification` | `11.0` |
| **Core Expected Sum (`E[W_run]`)** | **`111.4009`** |
| **Gated Cold-Run Claim (`E[W(n)]`, `C = 2,000`, `p_cold >= 0.39`)** | **`113.364850291911` G1 units** (core `110.042922197024`; constant for all `1 <= n < r`, scored at `n = 2^20`) |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared exploratory STNFS premises: **`H1`** (Dickman/de Bruijn + Murphy smoothness across collection, node sieving/survivors, and semismoothness of selected first-vector boot norms `N_g(b_0)` — ZERO `H-EA`), **`H2`** (`H2_ALLOW = 1.05`, `det_credited = 0.95`, `capped_fp_gso_h2` node retention `ret_node = 0.90` under hard TOTAL per-prime cap `T_p_total`), **`H3`** (zero-column-loss filtering caps), **`H-X`** (Schirokauer full-rank map & controllability), and **`H5-2lvl`** (pre-draw child bad-ideal probability `P_p(c_r(Q, p) in B | Q in G) <= eta_bad = 1.0%` over the parent's fresh uniform `p ~ Unif(P_f)` — and over the boot's fresh `e ~ Unif(Z/rZ)` — BEFORE revealing children, absorbing all gate failures `||b_0|| < 0.20*q^(1/16)` plus up to `1.0%` permanently dead ideals, and `p_0(Q), p_1(Q) <= 0.35` on good `Q in G` at `min_mu0, min_mu1 >= 20`, amplified via exact conditional independence of `K = 8` internal retries to `psi(Q) <= 0.00201`, pre-draw mean `phi(b) <= 0.0120`, and combined via unconditional Boole over ONLY `<= 16` top-level `ZERO-EA` boot children + exact Wald fresh-`e` i.i.d. `gamma_0`-blinded boot retries).
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` and `python3 -I -B audit_cm.py <<< '{"descent_mesh":384}'`.
