# Plain `G_1` ECDLP (`bls12-381-g1-ecdlp-plain`) through two `G_T` logarithms with `(4, 3)` `Q(zeta_15)` CM-STNFS behind a run gate: expected work `2^112.420017540845` G1 units

Track `bls12-381-g1-ecdlp-plain`. Input: `T = [x] P_0 in G_1` on BLS12-381 (`P_0` standard generator, uniform `x in {0, ..., r-1}`, zero auxiliary powers); output canonical 64-hex `x` with `[x] P_0 = T` (`x = 0` immediately for `T = O`). The score unit is one BLS12-381 `G_1` group-operation equivalent, **`C_G1 = 2,000`** pure baseline scalar ARMv8.0-A integer instruction equivalents (`simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64/x18; AAPCS64 `sub/add sp` `MONT5 = 279 <= 280`; zero undefined ASM labels and per-kernel disjoint GPR maps across all 30 hot-loop blocks; `mont2_cios = 47 <= 48 < 58.84` Ir for `L = 2`). Internal STNFS stages are accounted in dense `F_{p^12}` multiplication units (**`C_12 = 20,000 = 10` G1 units**) and converted by an exact factor of `10` (`+log2(10) = +3.321928094887...` bits, verified live in 80-digit `Decimal` inside `audit_cm.py`). No advice, smaller targets, or amortization are used. **No prior verdict or score is a premise.**

## Claim

`W_G1` is the total classical work (in `C_G1 = 2,000` Ir G1 units) of the one attack on the one target `T = [x] P_0 in G_1`. If `T = O`, return `x = 0` (`<= 1` G1 u). Otherwise Step 0 is a run gate: draw a fresh 64-bit word `u` from an independent entropy stream drawn before and isolated from the challenge `T` and all internal solver coins (blinding exponents `e`, Schirokauer/BW matrices, boot cells, `p_rand` draws, certificate samples); unless `u < K = 7201431620367092223`, stop as a failure (`<= 1` G1 u). `K` is the least integer with `K * 0.9990 >= 0.39 * 2^64`. Under exploratory premises `H1`–`H6`, `H-R`, `H-SM`, `H-X` and `H-FFT` (section 9), the algorithm returns the exact scalar `x` with probability `>= 0.9990` once past the gate (section 7); `u` is independent of the challenge and of all solver coins, so cold-run success is `>= 0.39` and

    E[W_G1] <= 1 + (K / 2^64) * 10 * (2^110.4551 + 2^66) < 2^112.420017540845  G1 units (C_G1 = 2,000),
    E[W_12] <= 1 + (K / 2^64) *      (2^110.4551 + 2^66) < 2^109.098089445958  dense-Fp12 core units (C_12 = 20,000).

`2^110.4551` core units (`= 2^113.776984` G1 units) is the sum of section 8's rows: every-run bounds for collection and linear algebra plus the expected work of the two individual logarithms (`2^107.5340` core u) and of ECM relation splitting; every major row carries `>= +0.0094`b slack above its live recomputed value. The score is an **expected** work bound; every run is bounded pointwise by the hard aborting meter `T_run = 2^123` core u plus wrapper:

    W_G1 <= 1 + 10 * (2^123 + 2^66) < 2^126.321928094888  G1 units.

Analytic, exploratory-v1 claim. Under `python3 -I -B` (pure Python 3.11+ stdlib, about 7 s), `audit_cm.py` recomputes from scratch the algebraic/CRT/NTT identities, Murphy-alpha, real-FB joint collisions and `(f, g)` survival coupling, a toy end-to-end DLP, `sigma/B` Band-2 exactness, a 400-draw live cross-check of the 1,000,000-draw `A^2 = 321,500` pool, Voronoi + 16D LLL + `4x` 381-bit `F_{p^6}`-subfield boot lifts, `srt_toy`, the 2D retention integration with the static pre-gates and collection-cut conditioning, the relation-split cap and collection ECM Chebyshev bound, the `2^12` filtering simulator + H-X, the `[6b]` missing-log reconstruction, the two-sided `Q(zeta_15)` H-X, the ARMv8.0-A kernel parser and 64-bit simulations, the `384`-mesh descent evaluation with the per-lattice yield certificate (the `912`-lattice production archive re-checked record by record, two lattices recomputed from their seeds, and the toy lattice-conditional archive re-checked), the ECM ladder, boot archive, `UDUP`, `BOUND_CELL` and retention checks, the `G_1`/`G_2` generator orders, and the 80-digit `Decimal` ledger, gate, exact `x10` G1 conversion and manifest.

## 1. Reduction (Plain `G_1` ECDLP -> Calibrated `gamma_0`-Blinded Two-Logarithm `G_T` Transfer) and Field Map

Let `P_0 in G_1`, `Q_0 in G_2` be the fixed public standard generators from `target-profile.json` (`P_0` on `y^2 = x^3 + 4`, `Q_0` on `y^2 = x^3 + 4(1 + u)`, `[r] P_0 = O`, `[r] Q_0 = O`, verified LIVE in `audit_cm.py`) and `g_0 := e(P_0, Q_0) in G_T \ {1}` of prime order `r = u^4 - u^2 + 1`. Linear algebra (`RANK-BRIDGE.md` §5) recovers a nonzero homogeneous virtual-logarithm functional `V_g(alpha) = lambda * log_g0(phi_iso^(-1)(alpha)^{e_cof}) mod r` (explicitly including `- c_g * Lambda_g(alpha)`, `e_cof := (p^12 - 1)/r`, `gcd(r, e_cof) == 1` asserted LIVE). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^{e_cof} in G_T \ {1}` (`1` final exponentiation `<= 2^11` core u, no third descent) with KNOWN nonzero virtual log `V_0 := e_cof * V_g(a_0 - b_0 X) mod r != 0`. Given `T = [x] P_0 in G_1`:
  0. Trivial identity check & run gate (`Claim`): if `T == O`, return `x = 0` (`<= 1` G1 u). Otherwise draw fresh independent 64-bit `u`; if `u >= K_gate`, stop as failure (`<= 1` G1 u).
  1. Check `T in G_1` (on `E(F_p)`, `[r] T = O`, `T != O`).
  2. Compute `h := e(P_0, Q_0) = g_0 in G_T \ {1}` and `z_T := e(T, Q_0) = g_0^x in G_T \ {1}` (`2` optimal Ate pairings + final exponentiations against public `Q_0`, `<= 2^11` core u) and map `h, z_T` into the `Q(zeta_15)` tower (`field-map.md`, wrapper `<= 2^66` core u, failure `< 2^-48`).
  3. For each target `z in {h, z_T}`: each fresh boot draw samples uniform `e in [0, r-1]`, blinds `z' := z * gamma_0^e` (if `z' == 1`, immediately returns `V(z) = -e * V_0 mod r`), lifts `z'` via Guillevic `F_{p^6}`-subfield elimination (`W_1 != 0` proved unconditionally in `DESCENT.md` §2), tests one uniform half-cell of `L_boot(z')` in `B_16(2)/{+-1}`, descends it (section 6) to obtain `V(z') = (lambda * e_cof) * log_g0(z) + e * V_0 mod r`, and unblinds `V(z) := (V(z') - e * V_0) mod r`.
  4. Since `lambda * e_cof != 0 mod r`, `V(h) != 0 mod r`. Compute `x := V(z_T) * V(h)^(-1) mod r` (the common factor `lambda * e_cof` cancels; no third descent). Verify `[x] P_0 == T` in `G_1` (`<= 2^9` G1 u) and output canonical 64-hex `x`, otherwise fail.

## 2. Tower, Polynomials, and `Aut = C_12` Symmetry

Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, lcm(p^4-1, p^6-1)) == 1` and `gcd(r, (p^12-1)/r) == 1` asserted live). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`, monogenic maximal order `O_H = Z[zeta_15]`, `h_K = 1`), `g(x) = x^3 - u omega^2`, `f(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = <sigma, omega> = Z/12Z`. Smoothness bounds: **`B_f = 2^70.79`, `B_g = 2^71.19`, `log2(U') = 62.827711`**, target excess `exc = 1.09`.

## 3. Relation Collection (`P5_prod_2ad_prim_v5_1m_r83`, `Arch2-Prod-ZeroSlice-Buf64K-B2-U16-Streamed`)

- **Region Volume (`A^2 = 321,500, r_ab = 0.83`):** `V_a = 2^55.8744 < 0.92 * 2^56` (`N_a^tbl = 4 V_a = 2^57.8744`), `V_b = 2^64.5535` (`L_max = 440.08 < W_tile = 448`, `L_b = 189.02`, `N_7 = 2^54.9910`, `V_kept = N_a^tbl * N_7 * L_b = 2^120.4279`, `f_prim = 0.955571`).
- **1,000,000-Draw 2D Pool & Static Pre-Gates:** 100% deterministic `c_gbatch` (`P_{71.19} + P_{70.79}` + retained `g`-batch candidate product trees across `2,048` owners at `3.091 EiB << 16.00 EiB` each, ALL THREE `L-11` subtrees `<= 1.053 EiB < 2^62 B`, `live = 71.39 <= 75.5` u) -> `det_credited = 0.99`. Static 2-adic class thresholds `T_g[c] = (299, 287, 275)` and constant `T_f_const = 118`b give `margin_net_rethi = +0.10531`b, `margin_strict_hi = +0.03474`b (`15.1` SE), `log2_rels_2se_live = 62.982177 >= 62.980 = LOG2_RELS`, conditioned odd-norm means `(547.280, 671.016) >= (546.0, 670.0)`, `c_cof_2d_upper_live = 1.0557` Ir/cell.
- **Per-Band `ipc_tot` Decomposition (`P_min = 32, p_m4_s = 52.0, p_m4 = 57.5, p_lo = p_f = 58.0, p_g = 63.5`, `H2_ALLOW = 1.01`):**

| Component | `g`-side (Ir/cell) | `f`-side (Ir/cell) | Notes (`KERNELS.md`, `STORAGE.md`) |
|---|---:|---:|---|
| Band-1a-long `[2^32, 2^52.0)` hits | `2.3848` | `2.3848` | `rev_b1_stream_u16` `1.125` + `rev_b1a_apply_u16` `3.625` <= `4.78` |
| Band-1a-short `[2^52.0, 2^57.5)` hits | `1.0332` | `1.0332` | `rev_b1_stream_tail1` `4.0` + `rev_b1a_apply_1` `6.0` = `10.0` |
| Band-1b `[2^57.5, 2^58.0)` hits | `0.1068` | `0.1068` | `rev_b1_stream_tail1` `4.0` + `rev_b1b_apply_1` `8.0` = `12.0` |
| Band-1 4-bank lookups + guards + row init | `1.0720` | `1.0720` | `4 * 96 = 384` (1a-long), `4 * 10 = 40` (1a-short), `12` (1b) |
| Band-2 `[2^58.0, 2^63.5)` row `k_0` + probe/SB hdr | `0.9767` | `0.0000` | `8.25 <= 10.0` `rev_b2_row_k0_step4` + `19 <= 20` probe + `cbz`-guarded `10` SB hdr |
| Sequential `64 KiB`-buffered Band-2 `(x4, w_lp)` | `2.2215` | `0.0000` | `15` emit `[x7],#16` + `6` recv `[x_rx],#16` + `0.004` flush = `21.004 <= 22.0` + `28` row-end drains |
| 64-cell `ORR` scan + blunt `1,000`-Ir cold | `0.8400` | *(combined)* | Live `sga_mean`-tied `0.8231 <= 0.84` (`47`-op `cbz`/`.Lscan_done: ret` + `1,000` Ir cold) |
| 2D Cofactorization (`c_cof_2d_upper_live`) | `1.0557` | *(combined)* | Static 3-class `T_g=(299,287,275)` + `T_f=118`b, `MONT5=280`, `ntt3_F` |
| **Total `ipc_tot` (exact sum, before `1.01` `H2_ALLOW`)** | **`14.2874`** | *(total)* | `V_kept + log2(ipc_tot * 1.01 / 20000)` = **`109.991165`** <= **`110.0006`** (`+0.0094`b slack) |

  `replay_a_tables <= 2^105.3423` u (`live = 105.3263`, `+0.0160`b slack), 1-owner sieve RAM `2^63.087 B`, 1-prime-per-owner Band-2 build peak `<= 14.981 EiB` & resident `<= 7.491 EiB`, total tables `2^118.349 B < 2^120 B`.

## 4. Filtering and Linear Algebra (`ArchB-1Dense-16Proj-16Shard-4Vec-ZeroCopy`, Stage-B-Capped Ladder, Conservative Caps `nl_cap = 0.00278`, `sw_cap = 0.00568`, `nonunit_ub = 0.25`)

- **Stage-B Row-Rejected Conditioned Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0`, conditioned norms `[547.0, 671.0]`, `sbf = 375`, `sbg = 315`, `W_max = 4,096`):** Raw smooth rows (`ret = 1.0`, zero prime relocation between rows) whose residual exceeds `row_bf > 375` or `row_bg > 315` are discarded outright (`stageb_row_keep_frac = 0.6483`, automatically inducing `ret_f_realized = 0.8246`, `ret_g_realized = 0.7770` at `U = 2^25`). At `U = 2^25`, `nl_25 = 0.0004961` vs conservative `nl_cap = 0.00278` (`5.60x` headroom), `sw_25 = 0.0008920` vs conservative `sw_cap = 0.00568` (`6.37x` headroom), `eu_25 = 0.918213`, exact expanded `avg_w_25 = 200.0182 < 201.0`, `max_row_w_25 = 1,596.05 <= 4,096`, `merge_stack_mean_w = 73.29 < 201.0`. Unit-fraction accounting: the Krylov kernel executes expanded 6-term hexagonal entries (`exp_unit_hex6 = 0.93443` at `U = 1024`, climbing to `0.9618` at `U = 4096`), and charging **`nonunit_ub = 0.25`** (`k5_charged = 12.00` for 4-vec Krylov vs `k5_live = 8.8284`; `k5_horner_charged = 14.00` for 1-vec Horner vs `k5_horner_live = 10.7825`) also covers the 4-tuple ring-entry non-unit fraction `1 - 0.7659 = 0.2341`.
- **Complete `< 16.00 EiB` Owner Map (`m = 32,768, n = 128 = 32 * 4, n_sm = 4`, ALL `< 2^64 B` asserted live):** Dense owner `2^63.368 B = 10.323 EiB` (`2.581 EiB` Horner reuse); `16` projection owners holding `Z + Y=[S,BZ] + X` at `2^63.359 B = 10.263 EiB`; `2,048` purge owners at `2^62.552 B = 5.866 EiB` (`1,040 B/rel`); `8,192` mutable filtering pools at `2^62.459 B = 5.498 EiB` (`nonunit_ub = 0.25`); `2,048` Krylov sequence owners at `2^63.354 B = 10.224 EiB`; `2^22` PM-Basis owners at `2^63.108 B = 8.623 EiB` (final `F(z)` `0.639 EiB`); `16` read-only matrix shard owners at `2^62.124 B = 4.358 EiB`; aggregate streaming `c_stream_mat = 73.0`, `c_stream_mat_1vec = 291.0`, `c_xproj = 70.0`.
- **Explicit LA Total Work (`audit_cm.py`):** With `N = 0.00278 * U' = 2^54.3370`, `N_sw = 0.00568 * U' = 2^55.3678`, `k5_charged = 12.00`, `k5_horner_charged = 14.00`, `1`-pass Horner justified by `PM-Basis` `[z^{s_out[j]}] q_j(z)[:s_sm] != 0`: **`linear_algebra = 107.232346 <= 107.2484` u (`+0.0160`b slack)**, **`matrix_generator = 98.2557 <= 98.2717` u (`+0.0160`b slack)**.

## 5. Main-Phase Success & Missing-Log Reconstruction (`[6b]`)

Collection ECM relation splitting deterministically truncates at `n_rels_upper = 2^63.530` (`Pr[N_rels > n_rels_upper] == 0`), using prime-exponent Newton perfect-power precheck and Cauchy-Schwarz correlated random-sum variance under post-`2^48` collection `p_r = 1e-3`, giving `ew_l2 = 94.5038 <= 95.3025` u and Chebyshev `< 2^-70.62 < 2^-70` (asserted live). Block Wiedemann fails with probability `<= 2^-66`, misc `< 2^-80`, Schirokauer injectivity `<= 2/r` (`+ 4/r` descent = `6/r`). `[6b]` `missing_recon_6b_live` (synthetic `2^14`-column Stage-B-capped peeling diagnostic) confirms `meas_miss <= mod_miss` across all occurrence bands.

## 6. Individual Logarithms (certified single trial per slot, Boole over slots, fresh-`e` Wald boot renewals)

`DESCENT.md` (`q15_descent16.py`, `q15_cert16.py`, `verify_q15_descent.py`, `q15_ecm16.py`, `q15_bound_cell.py`): `76` eager nodes (`71` `dsieve16` + `5` direct-ECM, `LB = LX = 68.5`, bins `[68.5, 75), [75, 80), .., [255, 260)` on both sides) and the boot `BOOT16 = (230, 312, 15, 2)`. Each slot draws `p_rand` uniformly from the `4,740` known-log `f_9` ideals in `[2^14, 2^16)`, reduces `L_{Q, p_rand}`, and keeps it only if a **per-lattice yield certificate** passes: the basis screen and `256` fresh uniform cells with exact norms, the search's admission rules (`v_Q >= 2`, `BOUND_CELL` cap) and exact unit-orbit sizes (rational Fincke-Pohst) give a layer-cake Clopper-Pearson lower bound `mu(L) >= 15` on the expected distinct usable relations, with draw `i` using error `1e-9 / (i (i + 1))` so that the sieved lattice violates `mu(L) >= 15` with probability `<= 1e-9` in total. After `64` uncertified draws the slot aborts (structural ideals such as the `(135, g)` example `t = 5/8 mod q` do). The slot sieves the certified ball, certifies the children of its first usable relation and descends them; a relation with an aborting child is replaced by the next usable relation of the same ball, and only a ball with `<= 3` usable relations or `4` rejected ones triggers a new attempt (fresh lattice, at most `6`). Under `H5` a certified ball is empty with probability `<= 1e-4` and has `<= 3` usable relations with probability `<= 0.01`; under `H6` a reached ideal is non-good with probability `<= 0.03` (`0.01` for boot children). Every reached slot then fails with probability `<= P_NODE / (1 - a) + a^6`, `a = 0.01 + (2 (kA + kB) * 0.03)^4`, and Boole over the expected `3,598.3` slot visits of a tree gives `Phi_tree <= 0.3634969` with no independence across slots. Failed trees are abandoned and the boot restarts with a fresh uniform `e` (i.i.d. renewals), so `p_descent_fail = 0`, `E[runs/log] = 1.57108` and, by Wald's identity, the expected two-log work is `E[D] <= 2^107.5340` core u (live `2^107.517970`, `+0.0160`b slack). The production campaign, run with the algorithm's exact trial law, certified all `28,000` first draws (`200` fresh ideals per plan entry, `600` per boot-child entry; minimum `74`, `>= 2.74x` the planned `mu`); a toy test with exact relation detection finds no empty lattice among `2,646` with per-lattice `mu >= 10` (`DESCENT.md` section 5).

The hard meter caps each opened run at `T_run = 2^123` core u: `E[min(W_run, T_run)] <= 2^110.4551` u and by Markov's inequality `P(meter) <= E[W_run] / T_run <= 0.00016735 <= 0.00017` (`<= 0.00018408` at `+10%`).

## 7. Success

Core opened-run success (80-digit `Decimal`; descent failures only trigger renewals):

    P(success | open) >= (1 - 2^-66 - 2^-70 - 2^-80 - 2/r) * (1 - 0.00016735 - 4/r) >= 0.99983265,

and with `+10%` meter sensitivity and the `2^-48` wrapper failure, `P(success | open) >= 0.99981592 >= 0.9990 = p_floor`. Past the gate `P(success | open) >= 0.9990`; cold-run success `>= (K / 2^64) * 0.9990 >= 0.39`, and `E[W_G1] <= 1 + (K / 2^64) * 10 * (2^110.4551 + 2^66) < 2^112.420017540845`.

## 8. Ledger (`C_12 = 20,000` Ir/core u; exact `x10` to `C_G1 = 2,000` Ir/G1 u)

| Component | Live `log2(u)` | Charged `log2(u)` | Slack |
|---|---:|---:|---:|
| `relation_collection` | `109.991165` | `110.0006` | `+0.0094`b |
| `replay_a_tables` | `105.3263` | `105.3423` | `+0.0160`b |
| `linear_algebra` | `107.232346` | `107.2484` | `+0.0160`b |
| `matrix_generator` | `98.2557` | `98.2717` | `+0.0160`b |
| `both_individual_logs_exp` | `107.517970` | `107.5340` | `+0.0160`b |
| `band_split` | `100.0057` | `100.0217` | `+0.0160`b |
| `ecm_split_expected` | `94.5038` | `95.3025` | `+0.7987`b |
| `canonical_a_and_sort` | `<= 86.0` | `86.0` | — |
| `schirokauer_4col_crt` | `<= 84.0` | `84.0` | — |
| `filtering_and_misc` | `<= 73.5` | `80.0` | `+6.5`b |
| `pg_pf_product_trees` | `71.3876` | `75.5` | `+4.1`b |
| `factor_base_roots_ideals` | `<= 75.0` | `75.0` | — |
| `numerical_certificate` | `<= 50.0` | `50.0` | — |
| `verification` | `<= 11.0` | `11.0` | — |
| **Opened-Run Expected Sum (`E[W_run]`, core u)** | **`110.455056`** | **`110.4551`** | — |
| **Same in G1 units (`x10`)** | **`113.776984`** | **`113.777029`** | — |
| **Hard Pointwise Meter (`T_run + wrapper`, core u)** | — | **`< 123 + 2^-56`** | — |
| **Gated Cold-Run Expected Claim (`E[W_12]`, core u)** | **`109.098089445957`** | **`109.098089445958`** | — |
| **Gated Cold-Run Expected Claim (`E[W_G1]`, G1 u, score)** | **`112.420017540844`** | **`112.420017540845`** | — |

## 9. Heuristics (Explicit Exploratory-v1 Premises)

- **`H1` (Smoothness & Split Density):** Dickman-de Bruijn / Murphy-alpha smoothness (`alpha_f = 2.74, alpha_g_total = 0.77`, exact 2-adic `v_2(N_g)` class law, 1,000,000-draw `A^2 = 321,500, r_ab = 0.83` norm pool cross-checked online), within-bin conditional `(f, g)` survival independence (accompanied at toy scale by `real_fb_joint_coll_live`), exploratory `0.99 / 1.01` `li(x)` ideal-density brackets, and split-prime fraction `<= KH = 0.08`; the one-cell-per-fresh-blinding boot has `BALL_PAIR = 0.0`.
- **`H2` (Collection & Sieve Cost):** Uniform 1,000,000-draw `r_ab = 0.83` pool with `H2_ALLOW = 1.01`, `det_credited = 0.99`, static 3-class `T_g[c] = (299, 287, 275)` + constant `T_f_const = 118`b pre-gates, explicit Band-1a/1b/2 streaming (`rev_b2_row_k0_step4` `8.25 <= 10.0`, `64 KiB`-buffered Band-2 `21.004 <= 22.0`), `47`-op scan + blunt `1,000`-Ir cold bound (`<= 0.84`), generous scalar build bounds (`32.0 / 10.0`), and `SIEVE_HIT = 2,600` Ir/hit.
- **`H3` (Matrix & Filtering):** `r` unramified and prime to class numbers; conservative zero-column-loss filtering caps `nl_cap = 0.00278`, `sw_cap = 0.00568` (`> 4.3x` above the Stage-B-capped `2^25` ladder; extrapolating from toy `Z[zeta_12]` scales `U = 2^16..2^25` to production `U' = 2^62.83`), deterministic `W_max = 4,096` guard, `nonunit_ub = 0.25`, `k5_charged = 12.00` (4-vec Krylov), `k5_horner_charged = 14.00` (1-vec Horner), `c_stream_mat = 73.0`, `c_stream_mat_1vec = 291.0`, `c_xproj = 70.0`.
- **`H-R` & `H-SM` (`RANK-BRIDGE.md` §1):** Filtered `A = [M | S]` has `rank(M) = N`, `rank(A) = D - 1`, and data-oblivious Schirokauer compression loses injectivity with probability `<= 4/r` (`6/r` total).
- **`H4` (Kernel Cost):** Pure baseline scalar ARMv8.0-A integer instruction counts (`KERNELS.md`, `STORAGE.md`, `30` hot-loop blocks incl. `rev_b2_row_k0_step4` and `47`-op `mont2_cios`; abstract 64-bit owner address model).
- **`H-X` (Block Wiedemann):** Full Krylov/Hankel ranks for block-diagonal `X` and `4` Schirokauer columns at charged `(d, e)` (supported live by `filtsim_and_hx_live`, `h3_ladder_and_genuine_q15_hx`, `unplanted_two_sided_bw`).
- **`H-FFT` (Local SSA):** Odd-`k`-exact `ntt3_F(N)` / `fft_F(N)` bounds local and 4-step distributed large-integer multiplication past `2^26` limbs.
- **`H5` (descent, certified ball):** a certified trial lattice `L_{Q, p_rand}` (`mu(L) >= 15` in the model, admission rules applied) has no usable relation with probability at most `P_EMPTY = 1e-4` and at most `3` with probability at most `P_LOW = 0.01`, over the algorithm's `p_rand` draw. Basis: the random-integer model of `H1` applied to the distinct orbits of one ball gives `3.1e-7` and `2.1e-4`; toy support in `DESCENT.md` section 5 (per-lattice dispersion `0.90`-`1.33`, empty lattices at or below the per-lattice prediction, none among `2,646` with `mu >= 10`, no `M <= 3` among `1,163` with `mu >= 15`).
- **`H6` (descent, reached ideals):** a reached ideal is *good* (each certificate draw passes with probability `>= 1/2`) except with probability `P_BAD = 0.03` (`0.01` for boot children), independently for distinct relations. Basis: goodness is a property of the ideal's own lattice; reached ideals are treated like uniform ideals of their side and norm (every degree-1 ideal of norm `q` contains a `1/q` fraction of the pair lattice); per entry, `200` (boot-child entries `600`) fresh production ideals at the hardest norms certify on their first draw, bounding the non-good fraction by `0.0299` (`0.00994`) at `95%`. Non-good ideals abort after `64` draws and their relation is replaced (`DESCENT.md` section 3). Nothing else is assumed about earlier failures or dependence across slots: the certificate (error `1e-9` over the node's coins), Boole over slots and i.i.d. fresh-`e` renewals give the rest. Cost premises: `UDUP = 0.40` and the plan's model `mu` set the radii (`912/912` production lattices certify), `BOOT_KEEP_FRAC = 0.90 < CP95 LB = 0.9395`, the `L_min = 2` `12`-bit ECM ladder (`sep = 0.5`, `p_r = 0.0043`, `2^26` Newton precheck) and random-integer ECM group-order smoothness (`q15_ecm16.py`).

## 10. Memory & Reproducibility

Peak storage `< 2^120` B (`2^118.349` B replay tables during collection; all per-owner spaces `< 2^64` B asserted live in `audit_cm.py`). Reproduction: `python3 -I -B audit_cm.py <<< '{}'`.
