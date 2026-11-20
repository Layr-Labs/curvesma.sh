# `ADAPTIVE.md` — Quantitative Distributional, Dispersion, & Multi-Level Reached-Slot Selection Support for `H5-Empty`

## 1. Variance-Decomposition & Overdispersed Negative-Binomial Model (Why `p_0(b) = max(0.001, (1 + 0.5*mu)^{-1})` Is Ultra-Conservative)

Large modeled means alone do not bound `P(M = 0)`; the empty-ball tail depends on **dispersion** (`Var(M)/mu^2`) and **reached-slot selection bias** under `A_{j,b}`. Decompose the relation count `M(q) = sum_{v in B_16(R)/O_H^*} 1_{E_v}` at a reached ideal `q` in bin `(j, b)` by the Law of Total Variance over the reduced lattice `L_q` (conditional mean `Lambda(L_q) := E[M(q) | L_q]`):

$$\frac{\operatorname{Var}(M)}{\mu^2} = \frac{1}{\mu} + \underbrace{\operatorname{CV}^2(\Lambda)}_{\text{ideal-to-ideal lattice/}\alpha\text{ dispersion}} + \underbrace{\rho_{\text{within}}}_{\text{within-ball pairwise correlation}}.$$

Under the Gamma–Poisson (Negative Binomial) overdispersed mixture `M | Lambda ~ Poisson(Lambda)`, `Lambda ~ Gamma(shape r = 1/k_tot, scale k_tot * mu)` (`Var(M)/mu^2 = 1/mu + k_tot`), the exact empty-ball probability is:

$$P(M = 0) = (1 + k_{\text{tot}} \cdot \mu)^{-1/k_{\text{tot}}}.$$

- **Extreme conservatism of charged `p_0(b) = max(0.001, (1 + 0.5*mu(b))^{-1})`:** Our charged bound sets **`r = 1`** (equivalent to **`k_tot = 1.0`**, i.e., **100% relative standard deviation `CV(Lambda) = 1.0`** with an exponential density all the way down to `0`!), applies an additional **`2x` mean haircut (`0.5 * mu`)**, AND imposes a hard **`0.001` floor** on every single slot (`min_mu = 838.21`, `max_p0 = 0.00238`, `P_boot_tree = 0.56454 < 0.60` over `<= 538` deterministic bit-sum knapsack slots).
- **Sensitivity check:** Even at `k_tot = 0.25` (`CV = 50%`, `r = 4`, far above our measured production `CV^2`), `(1 + 0.25 * 838.21)^{-4} = 5.1e-10` — **`4.6 x 10^6`x smaller** than our charged `p_0 = 0.00238`! Even at `k_tot = 0.50` (`r = 2`), `(1 + 0.5 * 838.21)^{-2} = 5.6e-6`, still **`420x` smaller** than `0.00238`.
- **Honest scope:** Transferring bounded dispersion (`k_tot << 1`) and bounded reached-slot selection from the three diagnostic pillars below to every production reached slot `(j, b)` at `u_f ~ 15..23` is the **explicitly declared exploratory STNFS heuristic `H5-Empty`**, NOT an unconditional theorem.

## 2. Pillar 1 — LIVE Strictly-Decreasing, Zero-Ancestor-Recount Multi-Level `Q(zeta_15)` `(9,3)` Adaptive Tree (`q15_adaptive_h5.b85`, `q15_adaptive_h5_live`)

Directly tests reached-slot selection under `A_{j,b}` (`P(M_{w_j}=0 | A_{j,b}, Q_boot) <= p_0(b)`) on the exact `(9,3)` `Q(zeta_15)` tower (`u_t = -5`, `lll_16d_int(99,100)`, `B_16(sqrt(3))/{+-1}`, `O_H`-coprime & `O_H^*`-orbit deduplicated, `B_FB = 8,192`, `3LP <= 2^{24}`):
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

- **Key conclusions from Pillar 1:**
  1. **Every reached-slot stratum clears `O/E_p0 <= 0.614 << 1.0` (pooled `0.4324`)** — our charged `p_0(b)` is **`2.3x` more conservative** than the actual observed failure rate across `410` adaptively reached depth-1 & depth-2 slots!
  2. **Negligible reached-slot selection bias:** Reached-slot mean is `91.8%` of fresh uniform controls (`94.7%` at depth-1, with the tiny `-8%` gap including the deliberate exclusion of the 1 ancestor relation), and reached/fresh fail is `1.126` — completely swamped by our `2x` mean haircut and `r = 1` tail.
  3. **Zero bad-ideal concentration:** Exact birthday collision factor **`C_bday in [0.859, 0.908]`** (measured against the exact count `N_bin_exact` of split `Q(zeta_15)` ideals in every bin) matches uniform independent draws (`~ 1.0`), ruling out concentration on bad ideals.
  4. **3 LIVE bit-exact strictly-decreasing depth-2 witnesses** replayed from scratch via `lll_16d_int(99,100)` + exact `norm_q15_z` in `audit_cm.py`:
     `[g, 19381] > [f, 15061] > [f, 7591]` (`4` new hits), `[g, 23911] > [f, 6121] > [f, 5881]` (`2` new hits), `[g, 40591] > [g, 7951] > [f, 3571]` (`2` new hits).

## 3. Pillar 2 — Direct Real-FB Empirical Dispersion (`Fano` & `k_hat`) Across 220 Low- & High-Yield `Q(zeta_15)` `(9,3)` Balls (`q15_93_b16_live`)

All 6 subsets (`120` low-yield `B_FB = 2^{13}, 3LP <= 2^{24}` + `100` higher-yield `B_FB = 2^{15}, 3LP <= 2^{30}`) are live-verified in `audit_cm.py` for empirical Fano factor `Var(M)/E[M]` and empirical NB overdispersion **`k_hat = max(0, (Var - E)/E^2)`** (plus 1 live regenerated 9-hit ball):

| Group | `n` | `mean(M)` | `min(M)` | Empty (`obs` vs `p_0`) | `Fano` (vs `NB(0.5)`) | **`k_hat` (vs charged `1.0`)** |
|---|---:|---:|---:|---:|---:|---:|
| `lo_g40` | `40` | `1.93` | `0` | `7/40` (`0.175 < 0.537`) | `1.316` (`1.96`) | **`0.1640`** |
| `lo_f40` | `40` | `3.58` | `1` | `0/40` (`0.000 < 0.406`) | `1.117` (`2.79`) | **`0.0328`** |
| `lo_enc40` | `40` | `2.02` | `0` | `8/40` (`0.200 < 0.562`) | `1.076` (`2.01`) | **`0.0375`** |
| `hi_enc40` | `40` | `14.65` | `2` | **`0/40`** (`0.000`) | `2.743` (`8.32`) | **`0.1190`** |
| `hi_g30` | `30` | `15.63` | `4` | **`0/30`** (`0.000`) | `3.911` (`8.82`) | **`0.1862`** |
| `hi_f30` | `30` | `25.33` | `6` | **`0/30`** (`0.000`) | `3.475` (`13.67`) | **`0.0977`** |

Every single subset — both fresh and adaptively encountered, at both low and high yields — exhibits **`k_hat in [0.0328, 0.1862] << 0.50 < 1.0`** (`5x..30x` tighter than our charged `r = 1` dispersion!).

## 4. Pillar 3 — Full-Production `u_f = 10..24` BLS12-381 `(4,3)` `Q(zeta_15)` Lattice Dispersion (`CV^2(Lambda_prod)`) & Within-Ball Co-Divisibility (`prod_sq24_live`, `prod_cond_q15_live`)

To bridge the `u`-regime from toy real-FB counting (`u ~ 5..8`) to full production (`u_f ~ 15..23`), we measure BOTH variance components of the decomposition in §1 directly on **real 381-bit production BLS12-381 `(4,3)` `Q(zeta_15)` special-`q` lattices**:
1. **Live inside `audit_cm.py` (`prod_sq24_live` & `prod_cond_q15_live`, `36/36` production lattices — `24` fresh controls + `12` adaptively encountered children from parent norms, `q in [2^{73}, 2^{258}]`):**
   - Evaluates directly in the production `u_f` range: `g_side` **`u_f in [10.56, 24.04]`**, `f_side` **`u_f in [9.33, 20.25]`**!
   - All `36/36` lattices beat every `qbins` sextile by `35.97..51.59` bits (`joint_pr_over_charged_min >= 87.94x`, median `475x..668x`), with zero degradation on the `12` adaptively encountered production children (`min >= 40x`!).
   - **Within-ball small-prime co-divisibility (`rho_within`, `p in [7..97]`):** Live pairwise correlation across primitive ball vectors in `prod_sq24_live` is **`sp_rho_mean = -0.0091` (`g_side`)** and **`-0.0227` (`f_side`)** (asserted `|sp_rho_mean| < 0.04`, confirming negligible positive pairwise smoothness correlation `rho_within ~ 0`!).
2. **Ideal-to-ideal lattice-yield dispersion `CV^2(Lambda_prod)` at production `u_f` (`120` real BLS12-381 `(4,3)` `Q(zeta_15)` lattices, `48` primitive vectors/lattice):**
   - At actual child large-prime bound `B_ch = max(70, 0.72 * log2 q)`: `g-100b` (`R=16`): `CV^2 = 0.1193` (`Gamma r = 8.38`); `f-100b` (`R=16`): `CV^2 = 0.1910` (`r = 5.24`); **`g-180b` (`R=48, u_f = 17.02`)**: **`CV^2 = 0.0435` (`Gamma r = 22.98`, `min/mean = 0.592`, `sp_rho = +0.0007`)**; **`f-180b` (`R=48, u_f = 14.45`)**: **`CV^2 = 0.0555` (`Gamma r = 18.02`, `min/mean = 0.548`, `sp_rho = -0.0015`)**.
   - Even at strict `lx = 70.0` (zero large primes): `CV^2(Lambda_prod) in [0.1307, 0.2414]` -> **`Gamma r in [4.14, 7.65] >> 1.0`**!
   - **Scope:** Production diagnostics evaluate modeled smoothness densities and small-prime co-divisibility on exact production norms (`u_f` up to `24`), not full production relation counts; together with Pillars 1–2 they provide consistent quantitative heuristic support for declared coarse-bin `H5-Empty` across both dispersion and reached-slot selection.
