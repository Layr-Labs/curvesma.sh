# `ADAPTIVE.md` — Pure `k = 1` Boole + Exact One-Sided Cantelli `H5-Empty` Premise & 3-Pillar Live `Q(zeta_15)` Support

## 1. Pure `k = 1` Zero-Retry Boole + Exact One-Sided Cantelli `H5-Empty` Declaration

Production BLS12-381 `Q(zeta_15)` descent uses **pure `k = 1` zero-retry Boole descent** per root attempt, wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting work meter `t_run = 2^{123.0}` units (`p_node_abort = p_repeat_loss = 0` proved):

1. **Zero-Retry `k = 1` DAG & Exact Coarse-Bin Linearity Bridge (`DESCENT.md` §3):** Each of the `78` nodes (`39` bins) scans `S_16(R) = B_16(R)/{+-1}` once (`ZERO` skipped primes in `dsieve16`, certified upper-bound `BOUND_CELL_HI`, `p_nc = 0`) and takes its first usable relation (`k = 1`); `M_{w_j} = 0` aborts the root attempt immediately and draws a fresh uniform `e in {0..r-1}` (ZERO internal retries, ZERO `THRESHOLD-REPAIR` `psi^k`, ZERO Gate E). Indexing realized DFS slots `j = 1..|T|` (`|T| <= 538` deterministically via bit-sum knapsack `sum lo(b_i) <= M`) and coarse bin events `A_{j,b} := {j <= |T| and w_j in bin b}`, exact Boole + linearity of expectation gives:
   $$P(E_{\text{fail}} \mid Q_{\text{boot}}) \le \sum_{j=1}^{538} \sum_{b=1}^{78} P(A_{j,b} \mid Q_{\text{boot}}) \cdot P(M_{w_j} = 0 \mid A_{j,b}, Q_{\text{boot}}).$$
2. **Explicitly Declared Exploratory Premise `H5-Empty` (Exact One-Sided Cantelli Form — NO Gamma–Poisson Overclaim!):**
   - For any non-negative integer count `M >= 0` with `E[M] >= mu` and variance `Var(M) <= v`, the **unconditional one-sided Cantelli inequality** (requiring NO Gamma–Poisson distributional assumption and immune to zero-inflation whenever `Var(M) <= v`) gives `P(M = 0) <= Var(M) / (Var(M) + E[M]^2) <= (1 + mu^2 / v)^{-1}`.
   - Accordingly, we **explicitly declare `H5-Empty` as an exploratory STNFS descent premise**: for every reached `(j, b)` with `P(A_{j,b} | Q_boot) > 0`, the zero-relation probability obeys
     $$P(M_{w_j} = 0 \mid A_{j,b}, Q_{\text{boot}}) \le p_0(b) := \max\bigl(0.001,\; (1 + 0.5\cdot \mu(b))^{-1}\bigr)$$
     (implied unconditionally via one-sided Cantelli by `Var(M_{w_j} | A_{j,b}, Q_boot) <= max(2*mu(b), 0.001*mu(b)^2)` — **explicitly stated as a declared exploratory STNFS premise rather than deduced from mean yield or empirical `CV^2` alone**).
   - Reverse linearity + deterministic pointwise bit-sum knapsack yields **`P(E_fail | Q_boot) <= P_boot_tree = 0.56454 < 0.60`** (`min_mu = 838.21`, `max_p0 = 0.00238`, `N_tree_max_1log = 538`), and exact `i.i.d.` Wald fresh-`e` boot gives `E[two_logs] <= 2^{111.357575} <= 2^{111.358575}` u (`96/185`: `2^{111.335467}`).
3. **Why Algebraic Zero-Intensity Obstructions Do Not Occur on `Q(zeta_15)` Split Prime Ideals:**
   - Every split degree-1 prime ideal `q` (`det(L_q) = q`) undergoes deterministic `lll_16d_int(99, 100)` reduction, obeying uniform Lovász basis-length upper bounds across all `|S_16(R)| >= 2^{80}` primitive lattice points (and exact norm caps `|N_g| <= 2^{ng_max}, |N_f| <= 2^{nf_max}` enforced inside `NORMS`).
   - Because `f_9(X), g_3(X)` are irreducible and primitive over PID `O_H = Z[zeta_15]` with fixed small discriminants, no split degree-1 prime ideal `q > 2^{70}` induces a local fixed-divisor obstruction on `N_f(L_q)/q` or `N_g(L_q)/q` beyond fixed small-prime content (already accounted for in `alpha_f, alpha_g`). Ruling out a hypothetical production-only zero-inflated subpopulation at `u_f ~ 15..23` remains part of the explicitly declared exploratory premise `H5-Empty`, supported diagnostically (not proved) by the three live pillars below.

## 2. Pillar 1 — LIVE Strictly-Decreasing, Zero-Ancestor-Recount Multi-Level `Q(zeta_15)` `(9,3)` Adaptive Tree (`q15_adaptive_h5.b85`, `q15_adaptive_h5_live`)

Directly executes pure `k = 1` multi-level adaptive descent on the exact `(9,3)` `Q(zeta_15)` tower (`u_t = -5`, `lll_16d_int(99,100)`, `B_16(sqrt(3))/{+-1}`, `O_H`-coprime & `O_H^*`-orbit deduplicated, `B_FB = 8,192`, `3LP <= 2^{24}`):
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

## 3. Pillar 2 — Direct Real-FB `Q(zeta_15)` `(9,3)` Balls (`220` Low- & High-Yield Balls, `q15_93_b16_live`)

All 6 subsets (`120` low-yield `B_FB = 2^{13}, 3LP <= 2^{24}` + `100` higher-yield `B_FB = 2^{15}, 3LP <= 2^{30}`) are live-verified in `audit_cm.py` (plus 1 live regenerated 9-hit ball):

| Group | `n` | `mean(M)` | `min(M)` | Empty (`obs` vs `p_0`) | `Fano` (`Var/E`) | `k_hat` (`max(0,(V-E)/E^2)`) |
|---|---:|---:|---:|---:|---:|---:|
| `lo_g40` | `40` | `1.93` | `0` | `7/40` (`0.175 < 0.537`) | `1.316` | `0.1640` |
| `lo_f40` | `40` | `3.58` | `1` | `0/40` (`0.000 < 0.406`) | `1.117` | `0.0328` |
| `lo_enc40` | `40` | `2.02` | `0` | `8/40` (`0.200 < 0.562`) | `1.076` | `0.0375` |
| `hi_enc40` | `40` | `14.65` | `2` | **`0/40`** (`0.000`) | `2.743` | `0.1190` |
| `hi_g30` | `30` | `15.63` | `4` | **`0/30`** (`0.000`) | `3.911` | `0.1862` |
| `hi_f30` | `30` | `25.33` | `6` | **`0/30`** (`0.000`) | `3.475` | `0.0977` |

## 4. Pillar 3 — Full-Production 381-Bit `Q(zeta_15)` Norm & Yield-Headroom Calibration (`prod_sq24_live`, `prod_cond_q15_live`)

All `36/36` real production BLS12-381 `(4,3)` `Q(zeta_15)` special-`q` lattices (`24` fresh controls `q in [2^{73}, 2^{258}]` spanning `u_f in [9.33, 24.04]` + `12` adaptively encountered production children from real parent norms) run live `lll_16d_int(99, 100)` inside `audit_cm.py`:
- Every lattice beats all 6 `qbins` sextiles by `35.97..51.59` bits (`joint_pr_over_charged_min >= 87.94x`, median `475x..668x` on fresh controls; `>= 40x` on all `12` encountered children), with live within-ball small-prime co-divisibility `sp_rho_mean = -0.0091` (`g`) and `-0.0227` (`f`).
- **Scope:** Production diagnostics calibrate norm bounds and modeled mean yields at `u_f in [9.33, 24.04]`, and small-scale `Q(zeta_15)` real-FB campaigns (`u ~ 5..8`) test multi-level selection and empirical empty rates; neither unconditionally proves the production lower tail, which is explicitly declared as exploratory premise **`H5-Empty`** (§1).
