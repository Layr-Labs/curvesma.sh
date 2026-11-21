# `ADAPTIVE.md` — Pairwise-`O_H`-Coprime `FB`-Blinded CRT Descent (`p_sat = 0`), Pointwise-in-`E` Single-Cap Tree Closure (`T_tree = 4.5 * W_tree_1`, `R_tree = 44`), & Multi-Regime Empirical Diagnostics

## 1. Protocol & Exact Pointwise-in-`E` Single-Cap Tree Closure (NO Local Per-Scan Cutoff)

Production BLS12-381 `Q(zeta_15)` descent (`39` bins / `78` nodes, `lx = 71.0`, `boot = (265, 220.0, 16)`, `p_bcap = 0.50`, `K_blind = 24`, `N_tree_max_1log = 3,376` deterministic knapsack slots) uses without-replacement `FB`-blinded CRT first-hit (`k = 1`, zero sibling retries) descent wrapped in `R_tree = 44` conditionally `i.i.d.` fresh-`e'` restarts given `E in G_pre` and unified meter `t_run = 2^123.0` u:

1. **Preprocessing Environment `E in G_pre`:** Let `G_pre` denote the good preprocessing event (`H1/H3/H4` succeed; `Pr[E not in G_pre] <= p_pre < 2^-50`, included in `fail_core`). Because BLS12-381 `(p, r, f, g, B_blind)` is fixed, relation collection is deterministic, and `H3` enforces zero column loss (`post_peel_cols_lost == post_merge_cols_lost == 0`), on `G_pre` the solved factor base `F_known = F_full` is uniquely determined (`E = E*` deterministic singleton).
2. **Pairwise-`O_H`-Coprime Blinding Pool `B_blind` (`|B_blind| = 192` in `[2048, 8192]`) & Exact CRT `Z`-Basis (`p_sat = 0` Proved):**
   - `check_algebra_bls12` verifies via exact polynomial GCD that `[2048, 8192]` (`q' == 1 mod 15`) has `0` ramified primes and `192` distinct `O_H`-prime ideals admitting `f`-side degree-1 ideals.
   - For any reached ideal `Q` (`q >= 2^71`) and draw `Q'_t in B_blind`, `gcd(q, q'_t) = 1` unconditionally; exact Bareiss-LLL (`lll_16d_int(99, 100)`, `assert dsq[16] == (q*q')^2`) spans `L_{Q, Q'_t} = L_Q^(s) cap L_{Q'_t}^(f)` with index `1` (`p_sat = 0` proved).
3. **Full-Ball Scan (NO Local Per-Scan Cutoff) & Single Deterministic Tree Cap `T_tree = 4.5 * W_tree_1`:**
   - Every blinding draw `L_{Q, Q'_t}` scans its **entire** deterministic box `S_16(R)` with zero local per-scan truncation. Thus any draw with `M(Q, Q'_t, E) >= 1` is guaranteed to return a qualifying relation unless total cumulative tree work exceeds `T_tree = 4.5 * W_tree_1`.
4. **Declared Per-Environment Conditional MODEL Inequality `H5-Blind(a, b)` (Supported, Not Proved):**
   - **`H5-Blind(a)` (Pointwise Per-Environment Exceptional-Tree Fraction):** For **every** good preprocessing environment `E in G_pre` and any fixed non-identity target `h in G_T \ {1}`, over a fresh uniform rerandomizer `e' ~ Unif(Z_r^*)` (which makes `h' = h * g^(e') ~ Unif(G_T \ {1})` conditionally independent of prior restarts given `E`), the conditional probability that the descent tree on `(h', E)` encounters any exceptional parent (violating `H5-Blind(b)`) satisfies:
     `f_exc(E) := Pr_(e' ~ Unif(Z_r^*))[E_exc(h * g^(e'), E) | E] <= eps_exc = 0.12`  for all `E in G_pre`.
   - **`H5-Blind(b)` (Non-Exceptional Slot Availability):** On non-exceptional trees (`not E_exc`), every reached slot `(j, b)` satisfies `p_blind(Q; E) := P_{Q' ~ Unif(B_blind)}(M(Q, Q', E) = 0 | Q, E) <= 0.50` (by Paley-Zygmund, sufficient whenever `mu >= 32.15` and `kappa_blind <= 0.9689`).
5. **Conditional Proof of Expected Work & Restart Failure Pointwise for Every `E in G_pre` (ZERO `E[f(E)^R]` vs `E[f(E)]^R` Gap!):**
   - *Newton-Maclaurin & Conditional Scan-Work Theorem (proved unconditionally in every fixed `E in G_pre`):* Sampling `K_blind = 24` distinct primes without replacement from `B_blind` proves `P(all K_blind draws empty | Q, E) <= 0.50^24 <= 2^-24` (`P_NM_tree = 0.0002012` by union bound over `<= 3,376` slots) AND `E[draws | Q, E] <= 1 / (1 - 0.50) = 2.0000`. Since each scan traverses the deterministic box `S_16(R)` with cost `<= own_1(b)`, linearity of expectation over the knapsack fanout **proves** `E[W_tree | not E_exc, E] <= W_tree_1` given `H5-Blind(b)`.
   - *Single Tree-Budget Markov Bound (pointwise in `E`):* On non-exceptional trees, exact Markov's inequality gives `P(W_tree > 4.5 * W_tree_1 | not E_exc, E) <= 1/4.5`. On exceptional trees (`f_exc(E) <= eps_exc`), `T_tree` deterministically aborts at work `<= 4.5 * W_tree_1`.
   - *Pointwise-in-`E` Restart Closure:* For **every** `E in G_pre`, single-attempt conditional failure satisfies `f(E) <= f_exc(E) + P_NM_tree + 1/4.5 <= P_fail_1 = 0.34242 < 0.45`, and `E[W_attempt | E] <= boot_own_exp + (1 + 3.5 * eps_exc) * W_tree_1`. Because the `R_tree = 44` restarts draw fresh independent `e'_1, ..., e'_{R_tree} ~ Unif(Z_r^*)` **given fixed `E`**, the attempts are conditionally `i.i.d.` given `E`, so **`P(all R_tree fail | E) = f(E)^44 <= P_fail_1^44 <= 1.0000e-18 < 2^-57`** and **`E[W_desc | E] <= E[W_attempt | E] / (1 - f(E)) <= E[W_attempt | E] / (1 - P_fail_1)`** hold **pointwise for every `E in G_pre`** (eliminating any `E[f(E)^R]` vs `(E[f(E)])^R` or geometric-expectation conditioning gap!). Expected two-log cost is `2^108.970155 <= 2^108.971200` u (`384/370`: `2^108.961237` u).

## 2. Scoped Multi-Regime Empirical Support (Explicitly Distinguished from 381-bit Production Proof)

Production 381-bit `71`-bit-smooth retained counts `M` cannot be enumerated directly on classical verification hardware. Accordingly, we connect `H5-Blind(a, b)` via Paley-Zygmund (`p_blind(Q) <= (mu^-1 + kappa) / (1 + mu^-1 + kappa)`, requiring mean stability `mu >= 32` and overdispersion `kappa <= 0.9689`) using **four distinct, non-conflated empirical layers**, all verified live inside `pkg/audit_cm.py`:

### Layer 1 — Direct Per-Parent `p_hat(Q)` Audit of `H5-Blind(a)` Across Adaptive `(3,3)` Trees (`q15_adaptive_h5_live`)

Built on the exact BLS12-381-polynomial-structured `(3,3)` tower (`f_3 = a^3 - (2+w)b^3`, `g_3 = a^3 + 2 w^2 b^3`, `B_fb = 524,288`, `|B_blind_toy| = 96`, direct `16 x 16` CRT `Z`-basis `assert dsq[16] == (q*q')^2`, strict squarefree factorization, ancestor exclusion):
- **Live Depth-5 8-Root Run + 100% Per-Parent `p_hat(Q)` Audit (`live_8root_strict_fb_blind.parent_p_hat_audit`):** Runs `8` fresh roots starting above `2^24` live (`8/8` completed, `0` aborts, `max_depth = 5`, `10 d1 + 25 d2+` slots, `mean_scans = 1.023`), AND then **directly audits EVERY SINGLE unique adaptively encountered parent `Q` (`43/43` parents across depths `0..5`)** across `6` evenly spaced blinding primes from `B_blind_toy`: **`0/8` trees violate `max_Q p_hat(Q) > 0.50` or `> 0.35`** (`global_max_parent_p_hat <= 0.1667`). Note that raw tree completion alone does not test `p_blind <= 0.50` on every encountered parent; this per-parent audit directly measures `H5-Blind(a)`'s exact parent-level condition!
- **120-Root Campaign & `90 x 24 = 2,160` Observed CRT Balls:** `120/120` roots completed (`max_depth = 4`, `109 d1 + 61 d2+` slots); per-stratum dispersion:

| Stratum | Ideals x Draws | `mean(M)` | `p_hat` mean (`max_Q`) | Ideals Failing All 24 | Observed `kappa_hat` mean (`max_Q`) | PZ Upper at `mu = 32` |
|---|---:|---:|---|---:|---|---:|
| `d1` (reached depth 1) | `30 x 24` | `15.57` | `0.0000` (`0.0000`) | `0 / 30` | `0.0772` (`0.1602`) | `0.1607` |
| `d2+` (reached depth `2..4`) | `30 x 24` | `12.40` | `0.0014` (`0.0417`) | `0 / 30` | `0.1119` (`0.1957`) | `0.1849` |
| `ctrl` (fresh control) | `30 x 24` | `21.22` | `0.0000` (`0.0000`) | `0 / 30` | `0.0703` (`0.2010`) | `0.1884` |

### Layer 2 — LIVE Degree-`(72, 24)` `(9,3)` CRT-Blinded Exact-Counting WITH ACTUAL PRODUCTION `U_BLS = -0xd201000000010000` (`q15_93_ubls_live`)

- **Exact Production `U_BLS` Live Verification (`4` Parents x `8` Blinding Draws = `32` Live CRT Balls):** Directly evaluates the **unmodified 64-bit production BLS12-381 polynomials** `f(a,b)` (degree `72`) and `g(a,b) = a^3 - U_BLS * omega^2 * b^3` (degree `24`, `U_BLS = -0xd201000000010000`) on `32` live `16 x 16` CRT-blinded lattices (`B_fb = 262,144`, `LP_g = 2^536`, `LP_f = 2^80`, `1,024` points/ball, exact `O_H`-coprimality check): **`0 / 32` empty balls (`p_hat_max_Q = 0.0000`)**, **`kappa_hat_max_Q = 0.0000` (sub-Poisson on all 4 ideals!)**, `pz_worst_case_at_mu32 = 0.0303 << 0.50`!

### Layer 3 — High-Mean `(9,3)` Degree-`(72, 24)` CRT Balls at `u = -5` (`q15_93_b16_live`)

- **`18` Parents x `24` Blinding Draws = `432` Observed CRT Balls (`crt_blind_93_18x24`):** **`0 / 432` empty draws (`p_hat_max_Q = 0.0000`)**, `mean_M = 42.1` (`min_Q_mean_M = 16.38`), observed `kappa_hat_mean = 0.2098` (`max_Q = 0.3163 < 0.9689`), `pz_worst_case_at_mu32 = 0.2579 < 0.50`, plus `0 / 100` empty on `hi100`.

### Layer 4 — 381-Bit CRT Lattice-Geometry & Modeled-Yield Stability (`prod_sq24_live`, `48` Live 381-Bit Bareiss-LLL Lattices)

- **Scope Clarification:** `verify_prod_sq24` runs `48` live 381-bit Bareiss-LLL CRT reductions (`8` parents spanning `73..258` bits x `6` blinding primes in `[2048, 8192]`) with actual `U_BLS` and evaluates **modeled** smoothness yield from the reduced basis norms (NOT direct `71`-bit retained counts `M`, which cannot be enumerated).
- **What It Confirms:** `det_Z(R) == q * q'` (`p_sat = 0`), `>= 42.04x` (`g`) and `>= 153.62x` (`f`) modeled yield headroom over the charged plan, and stable modeled lattice geometry across blinding primes (`cv2_modeled_mean = 0.0634`, `cv2_modeled_max_Q = 0.1274`). Together with Layers 1–3, this supports — without unconditionally proving — the declared per-environment exploratory model inequality `H5-Blind(a, b)`.
