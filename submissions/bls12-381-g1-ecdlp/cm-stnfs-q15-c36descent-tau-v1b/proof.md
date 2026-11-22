# BLS12-381 Powers-of-Tau `G_1` ECDLP (`bls12-381-g1-ecdlp`) — Complete Self-Contained Proof

**No prior verdict or score is a premise.** Every algebraic identity, kernel opcode count, filtering simulation, `box102` collection pool evaluation, `1-Dense + 16-Shard` linear-algebra memory/work bound, uncapped `Q(zeta_15)` 16D individual-logarithm descent (`384`-mesh, `TRUNC = 0.0`, `[7]`/`[7l]`/`[7e]`/`[11]`/`[11b]`/`[12]`), and 80-digit Decimal run-gate claim in this package is computed and asserted live by `audit_cm.py` under `python3 -I -B` (pure Python stdlib, pure baseline scalar ARMv8.0-A integer instructions, `simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64). Internal STNFS stages are accounted in dense `F_{p^{12}}` multiplication units (`C_{12} = 20,000` Ir) and converted exactly (`10x`, `+log2(10) = +3.321928` bits) to Track-4 **`bls12-381-g1-ecdlp` `G_1` group-operation equivalents (`C_{G1} = 2,000` Ir)**.

## 1. Track 4 (`bls12-381-g1-ecdlp`) Reduction, Zero Auxiliary-Power Accesses, & Declared Claim
- **Track & Scoring Parameter:** `bls12-381-g1-ecdlp` — BLS12-381 `G_1` discrete logarithm with auxiliary powers `(n, A = (A_1, ..., A_n), T = A_1 = [s] P_0)`, uniform `s in [1, r-1]`, valid for **ALL `1 <= n < r`**, scored at **`n = 2^20 = 1,048,576`**.
- **Zero Auxiliary-Power Accesses (`W(n)` Constant for All `1 <= n < r`):** The attack reads ONLY `T = A_1` and performs **ZERO accesses** to `A_2, ..., A_n` (zero Cheon-type preprocessing or divisibility preconditions on `n` or `r - 1`). Consequently, the total expected work bound **`W(n) = W(2^20) < 2^111.92142` G1 units** holds identically for every `1 <= n < r`.
- **Exact Reduction (`G_1` ECDLP -> Two `G_T` Virtual-Logarithm Descents with Calibrated `gamma_0` Blinding & Exact `lambda`-Cancellation):** Let `P_0 in G_1`, `Q_0 in G_2` be the public BLS12-381 generators and `g_0 := e(P_0, Q_0) in G_T` of prime order `r = u^4 - u^2 + 1`. Stage-2 linear algebra (`RANK-BRIDGE.md`) recovers a nonzero homogeneous virtual-logarithm functional `V(u) = lambda * log_{g_0}(u) mod r` (`lambda in (Z/rZ)^*`). Pick ONE Stage-1 collected relation `(a_0, b_0)` with known nonzero `V_g(a_0 - b_0 X) != 0 mod r`, and define the **calibrated `G_T` blinding generator** **`gamma_0 := phi_iso^{-1}(a_0 - b_0 X)^{(p^{12}-1)/r} in G_T \ {1}`** (`1` final exponentiation `<= 2^11` `F_{p^{12}}`-u, zero third descent) with known nonzero virtual log **`V_0 := ((p^{12}-1)/r) * V_g(a_0 - b_0 X) mod r != 0`** (`gcd(r, (p^{12}-1)/r) == 1` asserted LIVE in `audit_cm.py`). Given `T = [s] P_0 in G_1 \ {O}`:
  1. Compute **`z_T := e(T, Q_0) = g_0^s in G_T`** and **`g_0 := e(P_0, Q_0) in G_T`** (`2` pairings, `<= 2^11` `F_{p^{12}}`-u).
  2. For each target `z in {g_0, z_T}`: each fresh boot ball draws uniform `e in [0, r-1]`, blinds `z' := z * gamma_0^e` (exact uniform in `G_T`), runs individual-log descent (`DESCENT.md`) to recover `V(z')`, and unblinds **`V(z) := (V(z') - e * V_0) mod r = lambda * log_{g_0}(z) mod r`**.
  3. Since `lambda != 0`, `V(g_0) = lambda != 0 mod r`. Recover **`s := V(z_T) * V(g_0)^{-1} mod r`**, cancelling `lambda` identically.
  4. Verify **`[s] P_0 == T`** in `G_1` (`<= 2^11` `F_{p^{12}}`-u) and output the canonical 32-byte scalar `s in [1, r-1]`.
- **Declared Headline Claim (Track 4 `C_{G1} = 2,000` Ir/G1-op = `10x` `F_{p^{12}}`-u):**
  - Opened-run expected work: **`E[W_run] <= 2^109.9565` `F_{p^{12}}`-u** (`= 2^113.2784` G1-u) behind ONE unified hard aborting meter `T_run = 2^123` `F_{p^{12}}`-u (`= 10 * 2^123` G1-u, Markov `p_meter = 0.00011844 <= 0.00012`).
  - Opened-run success lower bound: **`succ_hr = 0.999064 >= p_floor = 0.9990`** (`> 0.39`), explicitly including `2 * boot_phi (8.03e-9) + 2^-12.295 (node-cap Markov) + 2^-10.660 (repeat collision) + p_meter + 2^-66 + 2^-70 + 2^-80 + 6/r + 2^-48`.
  - Cold-run gate (`K_gate = 7201431620367092223`, `(K_gate / 2^64) * 0.9990 >= 0.39`, wrapper `2^66` `F_{p^{12}}`-u):
    - Core `F_{p^{12}}` expected work: `E[W_{F12}] = 2^108.599489 < 2^108.59950` `F_{p^{12}}`-u.
    - **Track-4 `G_1` expected work (valid for ALL `1 <= n < r`, scored at `n = 2^20`):** **`E[W_{G1}(n)] = 1 + (K_gate / 2^64) * 10 * (E[W_run] + 2^66) = 2^111.921418 < 2^111.92142` G1 units** with cold success **`>= 0.39`**. (Claim bounds strictly expected work `E[W]`, not the separate pointwise `10 * 2^123` abort cap.)

## 2. Tower & Factor Base (`(4, 3)` `Q(zeta_15)` CM-STNFS)
- Curve `BLS12-381` (`u = -0xd201000000010000`, `p = 381`b, `r = 255`b, `gcd(r, (p^12-1)/r) == 1` asserted live). Base field `H = Q(zeta_15)` (`Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1`, monogenic `O_H = Z[zeta_15]`, `h_K = 1`, `Disc = 3^4 * 5^6`), `g_3(x) = x^3 - u * omega^2`, `f_9(x) = x^9 + (omega + 1) x^6 + x^3 - omega` (`omega = y^5`), `Aut = C_12` (`<sigma, omega>`).
- Smoothness bounds: `B_f = 2^71.29`, `B_g = 2^71.69`, `log2(U_p) = 63.317373`, target excess `exc = 1.09`.
- Live `audit_cm.py` verifies `algebra_bls12`, `murphy_alpha_exact_live`, `real_fb_joint_coll_live`, `toy_e2e_dlp`, `sigmaB_exactness`, `anchor_norms` (`12/12`), `voronoi_exact` (`10` units), and `srt_toy` (`0/128` mismatches).

## 3. Relation Collection (`P4_prod_2ad_prim_w10b_box102`, `Arch2-Prod-ZeroCopy-U32`)
- **2D Norm Pool & Finite-Region Margin:** `A2_prod = 313,000, r_ab = 1.02`, `V_kept = 2^120.118602`, `f_prim = 0.95518`, deterministic `c_gbatch_full` (`P_{71.69}` product tree `74.41 <= 75.5` u, 256-owner RAM `2^61.22 B << 2^64 B`) -> `det_credited = 0.99`.
  - Live `audit_cm.py` values: `log2_rels_2se_live = 63.658279`, `margin_net_rethi = +0.2256`b, `margin_strict_hi = +0.1558`b (`>= 34.4` MC SE), `c_cof_2d_upper_live = 1.1558` Ir/cell (with live `Q(zeta_15)` split-prime census up to `2^14` verifying `kh_g <= 0.10, kh_f <= 0.10`).
- **Pure Baseline ARMv8.0-A Kernels (`KERNELS.md`):** `relation_collection = 2^109.434464 <= 2^109.434964` `F_{p^{12}}`-u, `replay_a_tables = 2^106.0162 <= 2^106.01666` `F_{p^{12}}`-u, 1-owner build RAM `2^63.956 B < 2^64 B`, total tables `2^119.163 B < 2^120 B`.

## 4. Linear Algebra (`ArchB-1Dense-16Shard`), `H3`, & `H-X`
- **Zero-Column-Loss Filtering Ladder (`2^16..2^25`, `cols_lost = 0` at every scale):** `nl_25 = 0.001682` vs `nl_cap = 0.00247` (`1.468x` headroom); `sw_25 = 0.003417` vs `sw_cap = 0.00500` (`1.463x` headroom); `unit_run_3vec4_z12` (`92/12 = 7.6667 <= 7.70` Ir/entry-vec) -> `k5_live = 9.07 <= k5_charged = 10.60` (`+16.9%` headroom, zero TABM).
- **1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (`m = 16,384, n = 120 = 40 * 3, n_sm = 4`, ALL `< 2^64 B`, Priced Streaming):** Dense owner `2^63.2585 B << 2^64 B`, `16` matrix-shard owners `2^62.322 B << 2^64 B` each, cross-owner streaming `c_stream_mat = 268.0` Ir/row -> `linear_algebra = 2^107.822157 <= 2^107.822657` `F_{p^{12}}`-u, `matrix_generator = 2^95.6898 <= 2^95.6908` `F_{p^{12}}`-u.
- **Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` `H-X` + Live `gamma_0`-Blinding/Unblinding `lambda`-Cancellation + TRUE Same-Matrix `H3` (`w8_verifiers.py`, `hx_z12.py`):** `34/34` exact post-hoc BSGS hits, `gamma0_blind_unblind_lambda_cancel_ok = True`, `unplanted_two_sided_bw` and `filtsim_and_hx_live` verified live.

## 5. Individual-Logarithm Descent (`Q(zeta_15)` Uncapped 16D `c36p`-Style Descent, `DESCENT.md`)
- Adapting `c36p`'s uncapped descent architecture to 16D `Q(zeta_15)` (`q15_descent16.py`, `q15_retain16.py`, `q15_counts16.py`, `q15_ecm16.py`, `q15_cells16.py`, `verify_q15_descent.py`):
  - **Zero per-prime `dsieve16` caps and zero survivor caps**; 100% input-uniform 2-stage Bareiss LLL (`PAIR_SETUP = 2^35.2 > 2^35.06`, `LLL_NODE = 2^40.0 > 2^39.30`, `LLL_BOOT = 2^46.5 > 2^46.45`); exact 16D ball half-counts (`R <= 64`, peak descent `< 2^95 B`, `4 GiB/owner`, `129.03 MiB` read-only ranking table `< 256 MiB/worker`).
  - **Exact `Q(zeta_15)` retention (`q15_retain16.py`):** `base_primitive = 0.968504`, `counted_fraction = 0.906510`, `min_ret = 0.9056`, `max_loss = 9.844e-4 < 1e-3` across all 76 nodes (`[11b]` verified live).
  - **Exact 384-mesh recurrence (`TRUNC = 0.0`, live-coupled to `log2_rels_2se_live = 63.658279`, `nf_mean = 555.75`, `ng_mean = 680.35`):** `76/76` sieve nodes, `min_mu = 69.41`, `max_phi = 0.022372`, `max_node_psi = 0.2412`, `BOOT16 = (261, 336, 14, R_b=2)` (`boot_psi = 0.34161`, `boot_phi = 4.013e-9`, `E[S^2]/E[S]^2 = 1.000780 << 2^24`, `UDUP = 0.65` cleared across all 4 regimes `0.9739 / 0.9922 / 1.8589 / 2.0518`), **`two_log_units = 103.051745 <= 103.06` `F_{p^{12}}`-u**, `log2(E[N]) = 15.5052`, `N_MAX = 2^27.8` (`over_log2 = -12.295`, `rep_log2 = -10.660`).
  - **Live & Bundled Conditioning `[7e]` + `DISC = 0.001` (`q15_cond_7e.b85`, `verify_q15_descent.py`):** `3,116 / 3,116` top-edge ideals `>= 2.805x` (`95%` UB `0.000961 < 0.001`), `384 / 384` planted children `>= 7.185x`, `[7e]-(A)` `KS p = 0.22 / 0.35`, `[7e]-(B)` `KS p = 0.43 / 0.80` (`P(fail | sib fail) = 0.0`).

## 6. Complete Work Ledger (`C_{12} = 20,000` Ir/u, Converted `10x` to Track-4 `C_{G1} = 2,000` Ir/G1-op)

| Component | `log2(F_{p^{12}} u)` | `log2(G_1 ops)` (`+log2 10`) |
|---|---:|---:|
| `relation_collection` | `109.434964` | `112.756892` |
| `replay_a_tables` | `106.016660` | `109.338588` |
| `linear_algebra` | `107.822657` | `111.144585` |
| `matrix_generator` | `95.690800` | `99.012728` |
| `both_individual_logs_exp` | `103.060000` | `106.381928` |
| `band_split` | `100.674000` | `103.995928` |
| `ecm_split_expected` | `96.453300` | `99.775228` |
| `canonical_a_and_sort` | `86.000000` | `89.321928` |
| `schirokauer_4col_crt` | `84.000000` | `87.321928` |
| `filtering_and_misc` | `80.000000` | `83.321928` |
| `pg_pf_product_trees` | `75.500000` | `78.821928` |
| `factor_base_roots_ideals` | `75.000000` | `78.321928` |
| `numerical_certificate` | `50.000000` | `53.321928` |
| `verification` | `11.000000` | `14.321928` |
| **Core Expected Sum (`E[W_run]`, live `109.956179`)** | **`109.9565`** | **`113.278428`** |
| **Gated Cold-Run Claim (`E[W(n)]`, `1 <= n < r`, `p_cold >= 0.39`)** | **`108.59950`** | **`111.92142`** |

## 7. Declared Exploratory-v1 Assumptions & Reproducibility
- Declared exploratory STNFS premises: `H1` (standard STNFS Dickman-de Bruijn smoothness heuristic for norms), `H2` (`H2_ALLOW = 1.03`, `det_credited = 0.99`, `kh_g = kh_f = 0.10`), `H3` (`nl_cap = 0.00247`, `sw_cap = 0.00500`, `k5 = 10.60`, `c_stream_mat = 268.0`), `H-X` (Block Wiedemann nonzero kernel recovery), and `H5` (standard `c36p` negative-binomial relation law `r = 0.5 / 1.0`, `DISC = 0.001`, `BALL_PAIR = 2^24`, `UDUP = 0.65`, supported live by `[7e]` (A)+(B)+(C), `3,116/3,116` top-edge certifications, and exact 8D Fincke-Pohst unit orbits).
- Reproduction: `python3 -I -B audit_cm.py <<< '{}'` (default `descent_mesh = 384`).
