# `ADAPTIVE.md` — Pairwise-`O_H`-Coprime `FB`-Blinded CRT Descent (`p_sat = 0`), Single-Cap Tree Closure (`T_tree = 8.0 * W_tree_1`, `R_tree = 40`), & Multi-Regime Empirical Diagnostics

## 1. Protocol & Exact Single-Cap Tree Closure (NO Local Per-Scan Cutoff)

Production BLS12-381 `Q(zeta_15)` descent (`39` bins / `78` nodes, `lx = 71.0`, `boot = (246, 211.0, 16)`, `N_tree_max_1log = 2,729` deterministic knapsack slots) uses without-replacement `FB`-blinded CRT first-hit (`k = 1`, zero sibling retries) descent wrapped in `R_tree = 40` `i.i.d.` fresh-`e'` restarts and unified meter `t_run = 2^123.0` u:

1. **Pairwise-`O_H`-Coprime Blinding Pool `B_blind` (`|B_blind| = 192` in `[2048, 8192]`) & Exact CRT `Z`-Basis (`p_sat = 0` Proved):**
   - `check_algebra_bls12` verifies via exact polynomial GCD that `[2048, 8192]` (`q' == 1 mod 15`) has `0` ramified primes and `192` distinct `O_H`-prime ideals admitting `f`-side degree-1 ideals.
   - For any reached ideal `Q` (`q >= 2^71`) and draw `Q'_t in B_blind`, `gcd(q, q'_t) = 1` unconditionally; exact Bareiss-LLL (`lll_16d_int(99, 100)`, `assert dsq[16] == (q*q')^2`) spans `L_{Q, Q'_t} = L_Q^(s) cap L_{Q'_t}^(f)` with index `1` (`p_sat = 0` proved).
2. **Full-Ball Scan (NO Local Per-Scan Cutoff) & Single Deterministic Tree Cap `T_tree = 8.0 * W_tree_1`:**
   - Every blinding draw `L_{Q, Q'_t}` scans its **entire** ball `S_16(R)` with zero local per-scan truncation. Thus any draw with `M(Q, Q'_t, E) >= 1` is guaranteed to return a qualifying relation unless total cumulative tree work exceeds the single deterministic budget `T_tree = 8.0 * W_tree_1`.
3. **Declared Production Model Inequality `H5-Blind(a, b)` (Supported, Not Proved):**
   - **`H5-Blind(a)` (Exceptional-Tree Fraction):** Over uniform production challenges and fresh rerandomizers `e'`, the fraction of adaptively reached descent trees that encounter any exceptional parent (violating `H5-Blind(b)`) is at most `eps_exc = 0.15`.
   - **`H5-Blind(b)` (Non-Exceptional Slot Availability & Work):** On non-exceptional trees, every reached slot `(j, b)` satisfies `p_blind(Q, E) := P_{Q' ~ Unif(B_blind)}(M(Q, Q', E) = 0 | Q, E) <= 0.50` (consistent with Paley-Zygmund at `mu >= 32.15` whenever `kappa_blind <= 0.9689`) and expected full-ball scan work `<= own_1(b)` (`E[draws] <= 2.0`, total expected tree work `<= W_tree_1`).
4. **Conditional Proof of Success & Cost (Given `H5-Blind(a, b)`):**
   - *Newton-Maclaurin (proved unconditionally in every fixed `E`):* Sampling `K_blind = 24` distinct primes without replacement from `B_blind` proves `P(all 24 draws empty | Q, E) <= 0.50^24 = 2^-24` and `E[draws] <= 2.0` (verified in exact `Fraction` arithmetic in `pm_adaptive.py`). Union bound over `<= 2,729` slots gives `P_NM_tree = 0.0001627`.
   - *Single Tree-Budget Markov Bound:* On non-exceptional trees (`>= 1 - eps_exc`), exact Markov's inequality gives `P(W_tree > 8.0 * W_tree_1 | non-exc) <= 1/8.0`. On exceptional trees (`<= eps_exc`), `T_tree` deterministically aborts at work `8.0 * W_tree_1`.
   - *Combined Attempt Bound & `R_tree = 40` Restarts:* Single-attempt failure `P_fail_1 = eps_exc + P_NM_tree + 1/8.0 = 0.27516 < 0.45`; expected attempt work `<= boot_own_exp + (1 + 7.0 * eps_exc) * W_tree_1`; after `R_tree = 40` `i.i.d.` fresh-`e'` restarts, `p_desc_loss <= P_fail_1^40 <= 1.0000e-18 < 2^-57` (included in `fail_core`) and expected two-log cost is `2^109.311372 <= 2^109.312400` u (`384/370`: `2^109.290306` u).

## 2. Scoped Multi-Regime Empirical Support (Explicitly Distinguished from 381-bit Production Proof)

Production 381-bit `71`-bit-smooth retained counts `M` cannot be enumerated directly on classical verification hardware. Accordingly, we provide **three distinct, non-conflated empirical layers**, all verified live inside `pkg/audit_cm.py`:

### Layer 1 — Observed Retained-Relation Counts `M` in Full-Depth Adaptive `(3,3)` `Q(zeta_15)` Descent (`q15_blind_h5.b85` + Live Depth-5 8-Root Run)

Built on the exact BLS12-381-polynomial-structured `(3,3)` tower (`f_3 = a^3 - (2+w)b^3`, `g_3 = a^3 + 2 w^2 b^3`, `B_fb = 524,288`, `|B_blind_toy| = 96` solved `f`-side ideals, direct `16 x 16` CRT `Z`-basis `assert dsq[16] == (q*q')^2`, strict squarefree factorization, explicit ancestor exclusion):
- **120-Root Full-Depth Adaptive Campaign:** **`120 / 120` roots (`100.0%`, `0` aborts -> observed exceptional-tree rate `0.00% << eps_exc = 0.15`)** complete to `<= 524,288` FB leaves (`max_depth = 4`, `109` `d1` + `61` `d2+` slots, `mean_scans = 1.0`).
- **Observed Retained-Count `M` Dispersion (`90` Ideals x `24` Draws = `2,160` Observed CRT Balls):**

| Stratum | Ideals x Draws | `mean(M)` | `p_hat` mean (`max_Q`) | Ideals Failing All 24 | Observed `kappa_hat` mean (`max_Q`) | PZ Upper at `mu = 32` (`vs 0.50`) |
|---|---:|---:|---|---:|---|---:|
| `d1` (reached depth 1) | `30 x 24` | `15.57` | `0.0000` (`0.0000`) | `0 / 30` | `0.0772` (`0.1602`) | `0.1607` |
| `d2+` (reached depth `2..4`) | `30 x 24` | `12.40` | `0.0014` (`0.0417`) | `0 / 30` | `0.1119` (`0.1957`) | `0.1849` |
| `ctrl` (fresh control) | `30 x 24` | `21.22` | `0.0000` (`0.0000`) | `0 / 30` | `0.0703` (`0.2010`) | `0.1884` |

- **Live Depth-5 8-Root Adaptive Run Inside `audit_cm.py` (`q ~ 2^24`, `3` Live Depth-2 CRT Witnesses + Fresh `8`-Root Run):** Re-executes `lll_16d_int(99, 100)` on `3` depth-2 CRT witnesses AND runs `8` fresh roots starting above `2^24` live (`8 / 8` completed, `0` aborts, **`max_depth = 5`**, **`10` `d1` + `25` `d2+` slots**, `mean_scans = 1.023`).

### Layer 2 — Observed Retained-Relation Counts `M` on Exact BLS12-381 Degree-`(72, 24)` `(9,3)` CRT Lattices (`q15_93_b16_live`)

- **`18` Parents x `24` Blinding Draws = `432` Observed CRT Balls (`crt_blind_93_18x24`):** **`0 / 432` empty draws (`p_hat_max_Q = 0.0000`)**, `mean_M = 42.1` (`min_Q_mean_M = 16.38`), observed `kappa_hat_mean = 0.2098` (`max_Q = 0.3163 << 0.9688`), **`pz_worst_case_at_mu32 = 0.2579 << 0.50`**.
- **High-Regime `(9,3)` (`hi_enc40 + hi_g30 + hi_f30`):** `0 / 100` empty balls (`mean_M in [14.65, 25.33]`, `k_hat in [0.0977, 0.1862]`).

### Layer 3 — 381-Bit CRT Lattice-Geometry & Modeled-Yield Stability (`prod_sq24_live`, `48` Live 381-Bit Bareiss-LLL Lattices)

- **Scope Clarification:** `verify_prod_sq24` runs `48` live 381-bit Bareiss-LLL CRT reductions (`8` parents spanning `73..258` bits x `6` blinding primes in `[2048, 8192]`) and evaluates **modeled** smoothness yield from the actual reduced basis norms (NOT direct `71`-bit retained counts `M`, which cannot be enumerated).
- **What It Confirms:** `det_Z(R) == q * q'` (`p_sat = 0`), `>= 42.04x` (`g`) and `>= 153.62x` (`f`) modeled yield headroom over the charged plan, and stable modeled lattice geometry across blinding primes (`cv2_modeled_mean = 0.0634`, `cv2_modeled_max_Q = 0.1274`). Together with Layers 1–2 (observed retained counts `M`), this supports — without proving — the declared production exploratory model inequality `H5-Blind(a, b)`.
