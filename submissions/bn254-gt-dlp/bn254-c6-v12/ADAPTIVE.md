# `ADAPTIVE.md` — Structural Lower-Tail Non-Degeneracy & Multi-Regime Tail-Sensitive Support for `H5-Empty`

## 1. Honest Scope & Structural Lower-Tail Non-Degeneracy (Why Zero-Inflated Point Masses Have No Mechanism in Reduced Special-`q` Lattices)

For a reached special-`q` ideal $\mathfrak{q}$ in bin $(j, b)$, the charged empty-ball bound $p_0(b) = \max(0.001, (1 + 0.5\,\mu(b))^{-1})$ (`mu_min = 761.30`, `max_p0 = 0.002620`, `Phi_tree = 0.443208 <= 0.60` over `c_det_one = 427` deterministic knapsack slots, pure `k = 1` zero-retry Boole + exact Wald fresh-`e` boot) cannot be deduced from mean yield $\mu$ or second-moment dispersion $\operatorname{CV}^2(\Lambda)$ alone: an abstract zero-inflated mixing law putting mass $0.2$ at $\Lambda = 0$ and mass $0.8$ at $1.25\,\mu$ has $\operatorname{CV}^2 = 0.25$ yet empty-ball probability $\ge 0.2$. We therefore **never** invert $\operatorname{CV}^2$ into an unproved Gamma shape parameter.

Instead, three exact structural properties of reduced production BN254 `C6` special-`q` lattices rule out degenerate zero-inflated mechanisms while framing the remaining production tail transfer honestly as the **declared exploratory STNFS heuristic `H5-Empty`**:

1. **Deterministic Minkowski / Bareiss-LLL Pointwise Basis-Norm Ceiling (Zero Degenerate-Lattice Mass):** For *every* degree-1 split prime ideal $\mathfrak{q}$ with $q \in [2^{59}, 2^{220}]$, our $\delta = 9/10$ Bareiss 12D LLL (`bareiss_lll_12d`, input-uniform worst-case $\le 2^{32.89} < 2^{33}$ Ir) deterministically bounds $\prod_{i=0}^{11} \|b_i\|_G \le (4/(4\delta - 1))^{33}\,\det(L_\mathfrak{q})$. Every vector $(a, b) = B x \in B_{12}(R)$ therefore obeys a deterministic pointwise upper bound on $|N_f(a,b)|, |N_g(a,b)|$ — ruling out pathological lattices whose ball norms explode or collapse conditional smoothness density to zero.
2. **Exact Small-Prime CRT Independence ($\gcd(q, P_{\text{FB}}) = 1$):** Every production special-`q` satisfies $q > 2^{59} = B_{\text{FB}}$ and $\gcd(q, 2 \cdot 3 \cdot 19) = 1$. By the Chinese Remainder Theorem, reducing $L_\mathfrak{q}$ modulo every factor-base prime $p \le 2^{59}$ is surjective onto $(\mathcal{O}_H / p\mathcal{O}_H)^2$, so no factor-base prime can be systematically excluded or locally obstructed by $\mathfrak{q}$.
3. **Why Naive PGF Threshold Bounds Are NOT Used (`verify_threshold_counterexample`):** Plugging $t = \psi$ into $\mathbb{E}[t^M] \le (1 + 0.5\mu(1-t))^{-1}$ fails on bimodal counts ($0.10 t + 0.90 t^{100} \le G(t)$ yet $0.10 > G(0.10) = 0.03571$). Our production descent uses **pure `k = 1` zero-retry Boole** ($p_0(b)$ ONLY) + **exact Wald identity** over fresh independent-`e` boot retries (`two_logs = 83.798544`), NEVER PGF-at-$\psi$.

## 2. Pillar 1 — Tail-Sensitive Low-`p_0` Empirical Validation Across `324` Real-FB `C6` Ideals (`174` Adaptively Reached; `c6_pgf_campaign` & `deg8_dispersion_and_prod_yield`)

Directly tests the **lower tail** of real-FB relation counts ($M = 0$, $M \le 0.5\,\mu$, $M \le 0.8\,\mu$, and exponential lower-tail PGF $\mathbb{E}[0.8^M]$) in the low-$p_0$ regime ($p_0 \in [0.034, 0.25]$), all recomputed LIVE in `audit_bn_cm.py`:

1. **Real-2-core-FB `C6` usable-relation campaign (`u_t = 13, p_t = 79`, `224` ideals: `100` fresh ctrl + `100` adaptively reached depth-1 + `24` depth-2 + `16` live replays + `4` live 2-hop depth-2 chains, `mu in [27.64, 57.67]`, `p_0 in [0.0335, 0.0675]`):**

| Group | `n` | `M_min` | **`n_zero`** | `sum p_0` | **`O/E_p0`** | **`min(M/mu)`** | `p05(M/mu)` | **`M <= 0.5*mu`** | **`M <= 0.8*mu`** | `E_emp[0.8^M]` (vs `r=1`) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ctrl_100` | `100` | `26` | **`0`** | `4.594` | **`0.000`** | `0.7060` | `0.844` | **`0 / 100`** | `3 / 100` | `1.60e-4` (`<< 0.192`) |
| `slot_d1_100` | `100` | `28` | **`0`** | `4.962` | **`0.000`** | **`0.8527`** | `0.9792` | **`0 / 100`** | **`0 / 100`** | `1.40e-4` (`<< 0.205`) |
| `slot_d2_24` | `24` | `32` | **`0`** | `1.119` | **`0.000`** | **`0.8117`** | `0.8886` | **`0 / 24`** | **`0 / 24`** | `1.09e-4` (`<< 0.195`) |
| **`all_224`** | **`224`** | **`26`** | **`0`** | **`10.675`** | **`0.000`** | **`0.7060`** | **`0.8823`** | **`0 / 224`** | `3 / 224` | **`1.45e-4`** (`1,364x` below `0.198`) |

   - Across **ALL `124` adaptively reached depth-1 & depth-2 `C6` slots**, not a single slot fell even below **`0.80 * mu`** (`min(M/mu) = 0.8117`, `M_min = 28`, `0 / 124` zeros vs `6.081` expected), ruling out any heavy lower-tail or zero-inflated mass in real-FB `C6` counts!

2. **Exact deg-8 `f_rec` + deg-2 `g` mini-BN254 `C6` tower (`u_t = 5, p_t = 27631, r_t = 27481`, `B_FB = 1024`, `100` ideals: `50` fresh + `50` strict `k = 1` depth-1 reached + `6` live 2-hop depth-2 replays):**
   - `deg8_f_ctrl`: `n=25, mean=8.20, min=4, Fano=0.742, k_hat=0.0000, n_zero=0`; `deg8_f_d1_k1`: `n=25, mean=7.88, min=2, Fano=0.902, k_hat=0.0000, n_zero=0` (`d1/ctrl = 0.9610`).
   - `deg2_g_ctrl`: `n=25, mean=4.92, min=1, Fano=0.846, k_hat=0.0000, n_zero=0`; `deg2_g_d1_k1`: `n=25, mean=5.28, min=2, Fano=1.381, k_hat=0.0722, n_zero=0` (`d1/ctrl = 1.0732`).
   - **Pooled across all `324` low-`p_0` real-FB `C6` ideals (`174` adaptively reached):** **`0 / 324` zeros** vs **`34.72` expected** (`0 / 174` reached vs `18.10` expected -> **`O/E_p0 = 0.0000`**), empirical `k_hat <= 0.0722` across every stratum!

## 3. Pillar 2 — High-`p_0` Boundary Selection Stress-Test (`c6_adaptive.py`, `c6_adaptive_h5.b85`)

Complements the low-`p_0` tail check in §2 by testing reached-slot selection right at the high-failure boundary (`p_0 ~ 0.64`) on the exact `C6` `O_h` real-2-core-FB tower (`u_t = 2, p_t = 13, B_FB = 2^{14}`, `bareiss_lll_12d(9/10)`):
- **LIVE strict `k = 1` zero-retry traversal (`_desc_k1_slot`):** Runs `160` 17-bit `C6` roots LIVE inside `audit_bn_cm.py` -> **`76` nonempty roots == `49` completed + `27` empty-slot terminations** (exact partition asserted live), reaching **`48` strictly-smaller child slots** with **`27` observed slot zeros** vs **`30.621`** sum of fresh calibrated `p_cal(side, bit)` -> **direct per-slot `O/E_p0 = 0.8818 < 0.90 < 1.0`** (confirming zero adverse selection amplification even right at the empty-ball boundary!), plus `3` live bit-exact `k = 1` witnesses.
- **Companion `9,870`-root Boole-union & birthday-collision check (`c6_adaptive_h5.b85`):** `1,856` post-fail attempts (`1,946` total children) satisfy the attempt-level Boole union check (`boole_union_attempt_oe = 0.7935`, `pooled_oe_work = 0.9985`) with exact birthday collision factor **`max_C_bday = 2.62 < 3.0`** across all 10 concentration bins (`shared_prev_total = 3, shared_run_total = 4`).

## 4. Pillar 3 — LIVE Production Lower-Tail Quantiles & Within-Ball Correlation Across All `24/24` Real 254-Bit BN254 `C6` Lattices (`deg8_dispersion_and_prod_yield`)

To verify that neither lattice-yield lower-tail outliers nor positive within-ball small-prime clustering appear at full 254-bit BN254 scale (`q in [2^{59}, 2^{220}]`), `audit_bn_cm.py` evaluates **all `24/24` real production BN254 `C6` Bareiss-LLL (`delta = 9/10`) special-`q` lattices** LIVE:
1. **Production Lower-Tail Quantiles of Modeled Yield Ratio $\Lambda(L_\mathfrak{q}) / \mu_{\text{model}}$:**
   - Conservative (`alpha_f = 0, alpha_g = +3.4`): **`min = 2.46x`**, **`p10 = 3.22x`**, **`p25 = 4.38x`**, **`median = 6.07x`** (`prod_cv2_cons = 0.2357`) — **every single lattice exceeds `2.46x` the conservative modeled mean** (zero low-tail outliers!).
   - Actual (`alpha_f = -19.93, alpha_g = -0.04`): **`min = 7.9x`**, **`median = 25.0x`** (`prod_cv2_act = 0.2574`).
2. **LIVE Within-Ball Pairwise Small-Prime Co-Divisibility (`prod_sp_rho_mean`):**
   - Across `48` primitive ball vectors (`1,128` distinct pairs) and all `20` split primes $p \in [5, 400]$ ($p \ne 19$) on all `24/24` real production BN254 `C6` lattices: **`prod_sp_rho_mean = -0.0026`** (asserted `|sp_rho_mean| < 0.04`), confirming negligible positive within-ball smoothness clustering ($\rho_{\text{within}} \approx 0$)!
   - Minimum per-lattice mean small-prime smoothness score across all `24` lattices: **`prod_sp_score_min = 4.35 > 3.0`** (confirming §1.2 CRT local-density invariance with zero obstructed lattices!).
- **Honest Scope:** Production diagnostics evaluate modeled smoothness densities and exact small-prime co-divisibility on real 254-bit production norms rather than full production relation counts; together with the structural LLL+CRT non-degeneracy in §1 and the multi-regime real-FB lower-tail validation in §§2–3, they provide consistent exploratory support for declared heuristic `H5-Empty`.
