# `ADAPTIVE.md` — Pure `k = 1` Boole + Wald `i.i.d.` Boot & Basis-Conditioned Paley–Zygmund Decomposition for `H5`

## 1. Pure `k = 1` Zero-Retry Boole, Exact Wald Boot, & Basis-Conditioned Paley–Zygmund (`pure_k1_boole_wald_iid_boot`)

Production BN254 `C6` `v12b` keeps the exact same `62`-row `LB = 59` descent geometry and parameters as accepted `c8699d82`: **pure `k = 1` zero-retry Boole** per root attempt wrapped in **unbounded `i.i.d.` fresh-`e` boot retries** (Wald's identity) and ONE unified hard aborting run meter `t_run = 2^{112}` units:

1. **Exact Basis-Conditioned Law of Total Probability + Unconditional Paley–Zygmund:** At any reached special-`q` slot `(j, b)`, after running `bareiss_lll_12d(9/10)`, let `B_q` denote the observed reduced 12D basis and `Lambda(B_q) := E[M(q) | B_q]` its conditional expected count of usable relations in `B_{12}(R)`. Conditioning on the event `{Lambda(B_q) >= mu(b)}` vs `{Lambda(B_q) < mu(b)}` and applying the **unconditional one-sided Cantelli / Paley–Zygmund inequality** on `{Lambda(B_q) >= mu(b)}` (requiring ZERO Gamma, ZERO Poisson, and ZERO higher-moment assumption!) gives the exact identity:

$$P(M = 0 \mid \text{slot}) \le \underbrace{P(\Lambda(B_q) < \mu(b) \mid \text{slot})}_{\text{H5-geom: }\le\,\varepsilon_{\text{basis}} = 2.5 \times 10^{-4}} + \underbrace{\frac{1/\mu(b) + \rho_{\text{within}}^+}{1 + 1/\mu(b) + \rho_{\text{within}}^+}}_{\text{Unconditional Paley–Zygmund; H5-corr: }\rho_{\text{within}}^+ \le 2.5 \times 10^{-4}}.$$

   - **Exact 62-Node Verification (`pz_62_nodes_ok = True` LIVE in `evaluate_descent`):** For **every** real `mu > 0` — and asserted live across all `62/62` production nodes (`mu_min = 761.30`, `p0_max = 0.002620`) — `0.00025 + (1/mu + 0.00025)/(1 + 1/mu + 0.00025) <= max(0.001, (1 + 0.5*mu)^{-1}) = p_0(b)` (minimum slack `+6e-8` at `mu = 1,998`; at `mu_min = 761.30`, `0.001811 << 0.002620`).
   - **Honest Exploratory Scope (Zero Overclaim):** We explicitly declare both halves — `H5-geom` (`P(Lambda(B_q) < mu(b) | slot) <= 2.5e-4`) and `H5-corr` (`rho_within^+ <= 2.5e-4` once `B_q` is fixed) — as **declared exploratory STNFS heuristics**, NOT deduced unconditionally from `CV^2`, mean yield, LLL norm ceilings, or CRT surjectivity alone.
   - **Unconditional Boole + Wald Propagation:** `p_node_abort = p_repeat_loss = p_descent_fail = 0`, **`Phi_tree = 0.443208 <= 0.60`** over `c_det_one = 427` deterministic slots, `boot = [180, 76, 12, 4]` (`x_EA = 65.3021`, `2*P_x = 2^{62.094} B << 2^{64} B`), `E[runs/log] <= 1.796003`, `E[two_logs] <= 2^{83.798544}`.
   - **Naive-PGF Counterexample (`verify_threshold_counterexample`):** `0.10*t + 0.90*t^{100} <= (1+30(1-t))^{-1}` on `[0, 0.80]`, yet `0.10 > G(0.10) = 0.03571`. Pure `k = 1` Boole avoids `G(psi)` and `psi^k`.

## 2. Scoped Multi-Regime Empirical Diagnostics (All Recomputed LIVE in `audit_bn_cm.py`)

Finite prototypes and modeled production checks are scoped exploratory diagnostics, not unconditional proofs:

1. **Low-`p_0` Real-FB Lower-Tail Campaigns (`324` `C6` ideals, `174` adaptively reached):**
   - **Real-2-core-FB `C6` usable relations (`u_t = 13, p_t = 79`, `224` ideals: `100` ctrl + `100` depth-1 + `24` depth-2 + `16` live replays + `4` live 2-hop depth-2 chains, `mu in [27.64, 57.67]`):**
     - `all_224`: `M_min = 26`, **`0 / 224` zeros** vs **`10.675` expected** (`O/E_p0 = 0.000`), **`0 / 224` below `0.5 * mu`**, `min(M/mu) = 0.7060`, `p05 = 0.8823`, `E_emp[0.8^M] = 1.45e-4 << 0.1983`.
     - **Adaptively reached `124` slots (`slot_d1_100 + slot_d2_24`):** **`0 / 124` zeros** vs `6.081` expected, **`0 / 124` even below `0.8 * mu`** (`min(M/mu) = 0.8117`, `p05 = 0.942`, `M_min = 28`, `d1/ctrl = 1.0757`, `d2/ctrl = 1.0074`).
   - **Exact deg-8 `f_rec` + deg-2 `g` mini-BN254 `C6` (`u_t = 5, p_t = 27631`, `100` ideals: `50` ctrl + `50` strict `k = 1` depth-1 reached + `6` live 2-hop replays):** **`0 / 100` zeros** vs **`24.04` expected** (`O/E_p0 = 0.000`), `Fano in [0.742, 1.381]`, `k_hat in [0.0000, 0.0722]`, `f_d1/ctrl = 0.9610`, `g_d1/ctrl = 1.0732`.
2. **High-`p_0` Boundary Strict `k = 1` Selection Check (`c6_adaptive.py`, `c6_adaptive_h5.b85`, `u_t = 2, p_t = 13, B_FB = 2^{14}`):**
   - `160` live strict `k = 1` 17-bit roots (`_desc_k1_slot`, zero retry): **`76` nonempty == `49` completed + `27` empty-slot terminations** (exact partition asserted live), reaching **`48` strictly-smaller children** with `27` observed slot zeros vs `30.621` sum of fresh calibrated `p_cal(side, bit)` -> **direct per-slot `O/E_p0 = 0.8818 < 0.90`**, plus `3` live bit-exact `k = 1` witnesses.
   - Companion `9,870`-root Boole-union & concentration dataset (`1,856` post-fail attempts, `1,946` total children): `6/6` strata `boole_union_attempt_oe <= 0.871` (pooled `0.7935`, `pooled_oe_work = 0.9985`), `10/10` child bins **`max_C_bday = 2.6198 < 3.0`** (`shared_prev = 3, shared_run = 4`).
3. **Full 254-Bit Production Lattice Diagnostics (`24/24` real BN254 `C6` lattices, `deg8_dispersion_and_prod_yield`):**
   - Supporting `H5-geom`: across all `24/24` real Bareiss-LLL (`delta = 9/10`) production lattices (`q in [2^{59}, 2^{220}]`), conservative modeled yield ratio `Lambda_cons(B_q) / mu_model` has **`min = 2.46x, p10 = 3.22x, p25 = 4.38x, med = 6.07x`** (`0/24` below `1.0x`), and actual Murphy-`alpha` ratio has **`min = 7.9x, med = 25.0x`**.
   - Supporting `H5-corr`: cross-basis small-prime co-divisibility across `48` primitive ball vectors (`1,128` pairs evaluated over the `24` bases, `20` split primes `5 <= p <= 400`) gives **`prod_sp_rho_mean = -0.0026`** (asserted `|sp_rho_mean| < 0.04`), with `prod_sp_score_min = 4.35 > 3.0` across all `24/24` lattices (a scoped small-prime diagnostic over 24 bases, not a proof of fixed-basis full-relation covariance).
