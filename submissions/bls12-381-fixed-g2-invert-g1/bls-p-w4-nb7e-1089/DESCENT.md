# `DESCENT.md` — Standard Negative-Binomial Sibling-Retry Special-`q` Descent (`lx = 71.0`, `Q(zeta_15)`)

## 1. Special-`q` Lattice, Norm Model, & Sieve/ECM Node Pricing

Production individual logarithms descend every large prime ideal `Q` above `lx = 71.0` (`39` five-bit bins `[71, 76), ..., [261, 266)` on both `f` and `g` sides = `78` nodes) down to the solved 2-core factor base `F_known = F_2core subset F_full` (`B_f = 2^71.29`, `B_g = 2^71.69`, `U' = 2^63.317373`).

- **Single Special-`q` Lattice `L_Q` (`det = q`):** For a degree-1 prime ideal `Q` of norm `q in [2^lo, 2^hi)` (`q == 1 mod 15`), the exact integer Bareiss-LLL reduction (`lll_16d_int(99, 100)`, `det(R) = q`) produces a 16-vector `Z`-basis of `L_Q`. The node scans a deterministic 16D rank-2 ball `S_16(R_n)` (`R_n in [38, 82]`, exact cell counts in `cells_rank2`, verified live in `verify_b16_dsieve_live`), applying `g`-first threshold gating (`sieve_gateG` or `ecm_gateG`, whichever is cheaper) and exact unit-duplicate thinning (`udup_node = 0.50`, verified by `unitdomain.py`).
- **Norm Sextiles & Missing-Log Thinning (`F_known = F_2core`):** For each node, norms are evaluated conservatively at the upper bin endpoint `c_hi = (hi + log2(cells))/16` over the 6 upper sextiles and 24 sub-sextiles (`qbins`, `qbins_sub`, live-verified on 24 + 12 production 381-bit Bareiss-LLL lattices in `prod_sq24_live` and `prod_cond_q15_live`, with `>= 10x` modeled smoothness headroom and all 6 sextiles strictly below `qbins`). Relations touching any prime ideal in `F_full \ F_2core` (post-singleton-purge missing density `1 - peel_h`, `peel_h = 0.256845`) are discarded; this thinning factor `exp(-miss_n)` is explicitly charged at every node and at boot via `miss_side(m, u, lb)`.
- **Cofactorization & `P_71` Product Tree:** Each surviving candidate strips primes `<= 2^71` via Bernstein's product-remainder tree `P_71` (`128` owners `< 2^62 B`/owner, build `<= 2^73.369` u `<= 2^75.5` u), followed by ECM (`ecm_success = 0.999`) and Miller-Rabin tests.

## 2. Negative-Binomial Sibling-Retry Recurrence (`H5`)

Following the standard STNFS negative-binomial sibling-retry model (`H5`, explicitly declared as an exploratory heuristic, not proved at 381-bit scale):

- **Per-Node Distinct Relation Count:** Each special-`q` node yields a distinct retained relation count `M ~ NB(mu, r)` with shape parameter:
  - `r = 0.5` for nodes with `<= 1` large-prime child (`kA + kB <= 1`, including empty-window base nodes `[71, 76)`), and
  - `r = 1.0` for nodes with `kA + kB >= 2`.
- **Sibling-Retry Failure & Expected-Work Recurrence:** A node tries distinct retained relations one by one until one relation's entire child subtree succeeds. Let `phi(b, s)` denote the upper bound on subtree failure probability for bin `b`, side `s`, and `psi = sum_children phi` (worst-case knapsack sum over the node's permitted child windows `(nA, n, kA, kB)` subject to `sum log2(q_child) <= m_max`). Then:
  - `phi(node) = max((1 + mu * (1 - psi) / r)^(-r) + DISC, 0.02)` (`DISC = 0.001` per-node certified yield discard charge) <= 0.05`,
  - `psi <= 0.80`, `mu >= 18.0`,
  - Expected node work `W(node) = own_1 + (1 / (1 - psi)) * W_children`,
  - Expected special-`q` node count `E[N(node)] = 1 + (1 / (1 - psi)) * N_children`.
- **Live Quadrature Results (`audit_cm.py`):**
  - Default mesh `32/185`: `min_mu = 34.08 >= 18.0`, `max_psi = 0.3519 <= 0.80`, `max_phi = 0.0364 <= 0.05`, `psi_boot = 0.4654 <= 0.70`, `E[N_two] = 2^14.696`, **`two_log_exp_units = 107.189790`** (`boot = 105.0999`, `tree = 106.8035`).
  - Refined 12x mesh `384/370` (lower sums): `min_mu = 35.70`, `max_psi = 0.3499`, `max_phi = 0.0354`, `psi_boot = 0.4600`, `E[N_two] = 2^14.677`, **`two_log_exp_units = 107.125952`**.
  - Ledger charge: `<= 107.190790` log2 units.

## 3. Boot (`G_T` Initial Splitting), Node Counter `N_MAX`, & Repeat Bound

- **Initial Splitting (Boot):** Parameters `(n1 = 265, nA = 245.0, kA = 12, kB = 1, cells_b = 1, m_sv = 2592.0)`. Under `H1/H5`, each independent randomized boot trial (`e' ~ Unif(Z_r^*)`) succeeds with modeled semismoothness + `F_2core` retention probability `pr_b = 2^-39.5120` (`p_boot_trial = 1.257e-12`). Each log tries up to `K_boot = 48` independent boot relations sequentially until one subtree completes:
  - Per-log boot failure probability: `p_boot_fail_1 = psi_boot^48 = 1.142e-16 < 2^-30` (both logs `<= 2.283e-16`).
  - Expected two-log descent work: `2 * (boot_own_1 + W_boot_children) / (1 - psi_boot) = 2^107.189790 <= 2^107.190790` units.
- **Explicit Node Counter (`N_MAX = 2^26.5`):** A single global counter tracks all special-`q` nodes visited across both logarithms and aborts (counted failure) if it would start node `N_MAX + 1`. By Markov's inequality, `P(counter abort) <= E[N_two] / N_MAX = 0.000280`.
- **Conservative Special-`q` Repeat Bound:** Because missing-log relations are discarded (`exp(-miss_n)` thinning), all special-`q` nodes lie in the `39` bins above `2^71`. Conservatively placing all `N_MAX = 2^26.5` nodes into the smallest single-side bin `[71, 76)` (`|bin_1side| = li(2^76) - li(2^71) = 2^70.26` prime ideals) and allowing `2^6` for `1/p` weighting inside a bin gives `P(repeat) <= C(N_MAX, 2) * 2^6 / |bin_1side| = 2^(2*26.5 - 1 + 6 - 70.26) = 0.000204`.
- **Total Descent Failure Charge:** `2 * p_boot_fail_1 + P(counter abort) + P(repeat) = 0.000484 <= 0.000508` (included in `fail_core` alongside the unified work meter `t_run = 2^123.0`).
