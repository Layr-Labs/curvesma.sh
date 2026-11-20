# `ADAPTIVE.md` — Pure `k = 1` Boole + Stage-1 Exact-Integer `F_Vor / mu_30` Runtime Norm Certificate `S_cert(node) >= 10.0` (`H1-cert`) & Full-Production 381-Bit Calibration

## 1. Replacing `H5-Empty` with Standard `H1-cert` + Runtime Norm Certificate `S_cert(node) >= T_cert = 10.0`

Production BLS12-381 `Q(zeta_15)` descent uses **pure `k = 1` zero-retry Boole descent** per root attempt, wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting work meter `t_run = 2^{123.0}` units (`p_node_abort = p_repeat_loss = 0` proved), **WITHOUT any `H5-Empty` / Cantelli conditional-variance / zero-inflation premise**:

1. **Stage-1 Exact-Integer `F_Vor / mu_30` Runtime Certificate (`KERNELS.md` §5, `BOUND_CELL_g = 1,920` Ir/cell, ZERO Stage-2 `f_9` Evals!):**
   - For every reached descent node `w_j`, Stage 1 (`BOUND_CELL_HI_g`) ALREADY computes 60-bit certified embeddings `(r_{a,k}, i_{a,k}, r_{b,k}, i_{b,k}, e_k)` (`k = 0..3`) and unconditional `L_hi_g(x) >= log2|N_g(x)|` on `100%` of cells in `S_16(R)`.
   - Reusing those exact registers (`bound_cell_scert_sector/emb4/fvor/acc`, `<= 216 <= 240` Ir/cell worst-case pointwise -> `1,879 <= BOUND_CELL_g = 1,920` Ir/cell), each node accumulates:
     $$S_{\text{cert}}(\text{node}) := \sum_{x \in S_{16}(R) \cap \mathcal{C}_{\text{int}}} 0.95 \cdot \mathrm{ss}_f\bigl(L_{\text{cert},f}(x)\bigr)\cdot \mathrm{ss}_g\bigl(L_{\text{hi},g}(x)\bigr)\cdot \mathrm{det}_n\cdot \mathrm{ret}_{\text{node}}\cdot e^{-\mathrm{miss}_n},$$
     where:
     - `L_cert_f(x) := 9 * sum_{k=0}^3 log2_ceil(X_k^{hi} + Y_k^{hi}) - lq_f >= log2|N_f(x)|` is an **unconditional polynomial upper bound** (exact because `|w_k + 1| = |-w_k^2| = 1` and `(X^3+Y^3)(X^6+Y^6) <= (X^2+Y^2)^{9/2}` for all `X, Y >= 0`, verified live in `verify_kernels_symbolic`);
     - `\mathcal{C}_{\text{int}} \subset F_{\text{Vor}} / \mu_{30}` is the **strict integer inner cone** (`r_{a,0} > e_0`, `i_{a,0} > e_0`, `10000*(i_{a,0}+e_0) < 2125*(r_{a,0}-e_0)` since `tan(12 deg) > 0.2125`, AND all 10 Voronoi linear inequalities `sum_k c_{u,k} X_k^{lo/hi} > 0` proved in `unitdomain.py` + Sinnott full cyclotomic unit group check). **Consequence (ZERO Unit-Duplicate Overcounting):** Every cell `x in \mathcal{C}_{\text{int}}` represents a **100% pairwise distinct `O_H^*`-orbit**!
     - `0.95 < 1 / zeta_H(2) = 0.956` accounts explicitly for `O_H`-coprimality (`min N(P) = 16 => sum 1/N(P)^2 < 0.045`).
   - **Runtime Certificate Rule:** If `S_cert(node) < T_cert = 10.0`, the node deterministically treats the ball as failed, aborts the current root attempt, and draws a fresh uniform boot exponent `e in {0..r-1}`.
2. **Standard `H1-cert` Independence (`NO` Cantelli, `NO` Variance / Zero-Inflation Premise!):**
   - **Declared Exploratory Premise `H1-cert` (standard NFS `H1` independence on distinct `O_H^*`-orbits given certified norms):** Conditional on the node's exact lattice norms and `S_cert(node) >= T_cert = 10.0`, the events `{x yields a usable O_H-coprime relation}` across the pairwise-distinct `O_H^*`-orbits `x in S_16(R) \cap \mathcal{C}_{\text{int}}` are independent with probability `>= 0.95 * ss_f(L_cert_f(x)) * ss_g(L_hi_g(x)) * det_n * ret_node * exp(-miss_n)`.
   - **1-Line Unconditional Consequence:** Whenever `S_cert(node) >= 10.0`,
     $$P(M_{w_j} = 0 \mid \text{norms},\; S_{\text{cert}} \ge 10) = \prod_{x} (1 - p_x) \le \exp\!\Bigl(-\sum_x p_x\Bigr) \le e^{-S_{\text{cert}}} \le e^{-10} = 4.540 \times 10^{-5}.$$
3. **Exploratory Norm-Certificate Premise `H5-Cert` (Directly Measurable at Full Production 381-Bit Scale!):**
   - For every reached DFS slot `(j, b)` (`|T| <= 538` deterministically via bit-sum knapsack), we declare exploratory norm-stability premise **`H5-Cert`**: `P(S_cert(w_j) < 10.0 | A_{j,b}, Q_boot) <= eps_cert = 0.0005`.
   - Combining `H1-cert` and `H5-Cert` by exact Boole gives **`p_0(b) = exp(-10.0) + 0.0005 = 0.00054540`** per slot, **`P(E_fail | Q_boot) <= P_boot_tree = 538 * 0.00054540 = 0.29343 << 0.60`**, and exact `i.i.d.` Wald fresh-`e` boot gives **`E[two_logs] <= 2^{110.776344} <= 2^{110.777344}`** u (default `32/185` == `96/185`!).
4. **Constructive Exact `O_H`-Submodule Saturation `|N_{H/Q}(det_{O_H}(v_1, v_2))| == q` (`obpair.py`, `obasis.py`, `unitdomain.py`) & Zero Fixed Divisors:**
   - `obpair.ob_pair(q, ry, rt, R)` selects the shortest `R[j]` (`j >= 1`) with `|N_{H/Q}(a_1 b_j - b_1 a_j)| == q` (succeeding on `36/36` live production lattices!), and if none matches, executes the constructive PID completion `obasis.prime_generator` + `obasis.basis` (proven Fincke–Pohst generator search in `F_Vor` cone `T_2 <= 2 q^{1/4}/m_0`, `m_0 = 0.028378`, exact content division, integer Hermite complement `x_1 y - y_1 x = 1`, and Babai size-reduction against `{y^i b_1}` — live-verified on real 72-bit production prime `enc12[6]` inside `audit_cm.py`!).
   - Because `M = O_H v_1 + O_H v_2 == L_q` identically (`[L_q : M] = 1`), CRT makes `M -> (O_H/P)^2` surjective for EVERY prime ideal `P != Q`; exhaustive live proof in `check_algebra_bls12` verifies `min N(P) = 16 > 9 = deg(f_9) > 3 = deg(g_3)` (with unit coefficients `lead(f_9) = 1`, `const(g_3) = -1`), unconditionally ruling out ANY fixed prime-ideal divisor `P != Q` on `M = L_q`!

## 2. Pillar 1 — Full-Production 381-Bit `S_cert` Calibration (`1,000/1,000` Fresh + `36/36` Live-Verified!)

Unlike `M`-count lower tails, `S_cert(node)` is a pure function of exact lattice norms and is **measured directly on 1,036 genuine full-production 381-bit `Q(zeta_15)` special-`q` lattices** (`72b..263b`, `u_f in [9.33, 24.04]`, `100%` exact saturated `|N(det)| == q`):
- **`1,000` fresh full-production lattices across all `39` bins (`prod_scert1000.b85`, verified live in `audit_cm.py`):** **`1,000/1,000` pass `S_cert >= 10.0` (`0` fails!)**; `min(S_cert) = 29.7` (`2.97x` above `10.0`, even with 256-pt MC noise!), `p0.1% = 133.0` (`13.3x`), `p1% = 751.2` (`75x`), **`median = 9,968.8` (`997x` above `10.0`!)**, `max = 1,069,476.9`. Across all `1,000` lattices, `log2(S_cert)` has `mean = 13.266`, `sd = 1.582` -> **`T_cert = 10.0` (`log2 = 3.322`) lies `6.29` standard deviations below the production mean!**
- **`12/12` adaptively encountered production children (`prod_enc12`, recomputed LIVE from scratch in `audit_cm.py`):** **`12/12` exact saturated (`obpair`), `min(S_cert) = 1,000.0` (`100x` above `T_cert = 10.0`!), `median(S_cert) = 8,087.0` (`809x`!)**, `36.39..52.15`b sextile headroom (`183.49x..595.8x`).
- **`24/24` fresh production lattices (`prod_sq24`, recomputed LIVE from scratch in `audit_cm.py`):** `24/24` exact saturated (`obpair`), `35.97..52.92`b sextile headroom (`192.89x..668.18x`), live `sp_rho_mean = -0.0024 / -0.0227`.

## 3. Pillar 2 — Strictly-Decreasing Zero-Recount Multi-Level `Q(zeta_15)` Adaptive Tree (`q15_adaptive_h5.b85`, `q15_adaptive_h5_live`)

Companion multi-level adaptive descent diagnostic on `(9,3)` `Q(zeta_15)` (`u_t = -5`, `lll_16d_int(99,100)`, `B_16(sqrt(3))/{+-1}`, `B_FB = 8,192`, `3LP <= 2^{24}`):
- **Strict norm descent & ZERO ancestor recounting:** Every child prime satisfies **`2,048 < p_child < q_parent`**, and every ancestor relation `(nf_anc, ng_anc)` is **explicitly excluded** when counting child hits `M` (making parent-relation recounting or backtracking algebraically impossible!).
- **Fresh-ideal calibration:** `437` fresh uniform split ideals calibrated across all `12` `(side, bin)` strata (`b = 11..16`, `Fano in [0.56, 1.19]` across `11/12` bins).
- **Multi-level pure `k = 1` adaptive tree:** `500` roots (`209` ok) -> **`286` strictly-smaller depth-1 child slots + `124` strictly-smaller depth-2 child slots = `410` total adaptively reached child slots** (`232` distinct ideals, `90` real observed failures), ALL recomputed from raw rows in `audit_cm.py`:

| Stratum | Reached slots `n` | Unique ideals | Obs. fails (`rate`) | **`O/E` vs `p_0(b)`** | Reached / fresh fail | Reached / fresh mean | Exact **`C_bday`** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `d1_all` | `286` | `193` | `73` (`0.255`) | **`0.4756`** (`< 1`) | `1.1055` | `0.9472` | **`0.905`** |
| `d1_g` | `106` | `84` | `44` (`0.415`) | **`0.6136`** (`< 1`) | `1.0331` | `0.8687` | **`0.885`** |
| `d1_f` | `180` | `109` | `29` (`0.161`) | **`0.3546`** (`< 1`) | `1.2372` | `0.9657` | **`0.908`** |
| `d2_all` | `124` | `81` | `17` (`0.137`) | **`0.3112`** (`< 1`) | `1.2231` | `0.8651` | **`0.859`** |
| **Pooled `d1+d2`** | **`410`** | **`232`** | **`90` (`0.220`)** | **`0.4324`** (`< 0.50`) | **`1.1260`** | **`0.9175`** | **`0.899`** |

- **Key diagnostic observations:** `4/4` strata clear `O/E_p0 <= 0.614 << 1.0` (pooled `0.4324`, `2.3x` headroom); reached-slot mean is `91.8%` of matched fresh controls (`94.7%` at depth-1); exact child-ideal birthday collision **`C_bday in [0.859, 0.908]`** (`~ 1.0` uniform, ZERO bad-ideal concentration); and **3 LIVE bit-exact strictly-decreasing depth-2 witnesses** are replayed from scratch via `lll_16d_int(99,100)` + exact `norm_q15_z` in `audit_cm.py`: `[g, 19381] > [f, 15061] > [f, 7591]` (`4` new hits), `[g, 23911] > [f, 6121] > [f, 5881]` (`2` new hits), `[g, 40591] > [g, 7951] > [f, 3571]` (`2` new hits).

## 4. Pillar 3 — Direct Real-FB `Q(zeta_15)` `(9,3)` Balls (`220` Low- & High-Yield Balls, `q15_93_b16_live`)

All 6 subsets (`120` low-yield `B_FB = 2^{13}, 3LP <= 2^{24}` + `100` higher-yield `B_FB = 2^{15}, 3LP <= 2^{30}`) are live-verified in `audit_cm.py` (plus 1 live regenerated 9-hit ball):

| Group | `n` | `mean(M)` | `min(M)` | Empty (`obs` vs `p_0`) | `Fano` (`Var/E`) | `k_hat` (`max(0,(V-E)/E^2)`) |
|---|---:|---:|---:|---:|---:|---:|
| `lo_g40` | `40` | `1.93` | `0` | `7/40` (`0.175 < 0.537`) | `1.316` | `0.1640` |
| `lo_f40` | `40` | `3.58` | `1` | `0/40` (`0.000 < 0.406`) | `1.117` | `0.0328` |
| `lo_enc40` | `40` | `2.02` | `0` | `8/40` (`0.200 < 0.562`) | `1.076` | `0.0375` |
| `hi_enc40` | `40` | `14.65` | `2` | **`0/40`** (`0.000`) | `2.743` | `0.1190` |
| `hi_g30` | `30` | `15.63` | `4` | **`0/30`** (`0.000`) | `3.911` | `0.1862` |
| `hi_f30` | `30` | `25.33` | `6` | **`0/30`** (`0.000`) | `3.475` | `0.0977` |

## 5. Pillar 4 — Full-Production 381-Bit `Q(zeta_15)` Norm & Yield-Headroom Calibration (`prod_sq24_live`, `prod_cond_q15_live`)

All `36/36` real production BLS12-381 `(4,3)` `Q(zeta_15)` special-`q` lattices (`24` fresh controls `q in [2^{73}, 2^{258}]` spanning `u_f in [9.33, 24.04]` + `12` adaptively encountered production children from real parent norms) run live `lll_16d_int(99, 100)` + **exact `|N_{H/Q}(det_{O_H}(v_1, v_2))| == q` saturation (`36/36`)** inside `audit_cm.py`:
- Every exact-saturated lattice beats all 6 `qbins` sextiles by `35.97..52.92` bits (`joint_pr_over_charged_min >= 192.89x`, median `550x..668x` on fresh controls; **`>= 183.49x`**, median `595.8x` on all `12` encountered children!), with live within-ball small-prime co-divisibility `sp_rho_mean = -0.0024` (`g`) and `-0.0227` (`f`).
- **Scope:** Production diagnostics (`1,000` fresh + `12` live encountered `S_cert` lattices + `24` fresh controls) calibrate exact `S_cert(node)` norm-certificate stability directly at full 381-bit production scale (`u_f in [9.33, 24.04]`, `6.29` sigma above `T_cert = 10.0`), while small-scale `Q(zeta_15)` real-FB campaigns (`u ~ 5..8`) test multi-level selection and empirical empty rates; neither unconditionally proves the full-production independence or norm-certificate tail, which are explicitly declared as exploratory premises **`H1-cert`** and **`H5-Cert`** (§1).
