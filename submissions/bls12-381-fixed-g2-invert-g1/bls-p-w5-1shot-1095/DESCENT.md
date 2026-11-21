# `DESCENT.md` — 1-Shot Zero-Sibling-Retry Pure-Boole Special-`q` Descent (`lx = 71.0`, `Q(zeta_15)`)

## 1. Design Principle: Zero Sibling Retries, Zero Independence Assumptions, Pure Boole Union Bound

Instead of adaptive sibling retries (`k > 1`), per-node discard tuning, or multi-attempt boot powers (`K_boot > 1`), our individual-logarithm descent uses the **simplest possible 1-shot architecture (`k = 1`, `K_boot = 1`)**:

1. **Single-Pass Nodes (`k = 1`):** For every special-`q` prime ideal `Q` above `lx = 71.0` (`39` five-bit bins `[71, 76), ..., [261, 266)` on both `f` and `g` sides = `78` nodes), the algorithm performs **one** deterministic scan of its 16D Bareiss-LLL ball `S_16(R_n)`, takes the **first** surviving relation whose large-prime children lie in the node's permitted windows, and **never retries** a sibling relation if a child fails.
2. **Deterministic Worst-Case Tree Size (`N_tree_1log = 3,157`, `N_tree_2log = 6,314`):** Because there are zero sibling retries, every node's worst-case subtree size `N(b, s)` and worst-case tree work `W(b, s)` satisfy exact deterministic knapsack recurrences:
   - `W(b, s) = own_1(b, s) + max_children sum W(child)`,
   - `N(b, s) = 1 + max_children sum N(child)` -> at most **`3,157` special-`q` nodes per log (`6,314` across both logs)**, deterministically! No `N_MAX` node counter is needed.
3. **Pure Boole Union Bound (ZERO Subtree or Sibling Independence!):** Let `p_0(b, s) = (1 + mu(b, s) / r)^(-r)` be the marginal `H5` empty-node probability (`M ~ NB(mu, r)`, `r = 2.0` -> dispersion `1/r = 0.5000`). By **Boole's union bound** (which holds under **arbitrary dependence** across siblings, parent-child pairs, ancestors, and targets!), the probability that *any* node in either descent tree is empty is at most:
   - `P_sub(b, s) = p_0(b, s) + max_children sum P_sub(child)`,
   - **`P_2log_Boole = 2 * P_boot_children = 0.055100`** (`P_1log_Boole = 0.027550`).
4. **Zero Repeat-Collision Loss:** In a 1-shot descent without retries, if any prime ideal `Q` appears twice across the tree(s), its single already-found relation is simply reused (memoized) with zero additional failure probability and zero cycle risk (since every child norm is strictly smaller than its parent norm `q`).

## 2. Special-`q` Lattice, Norm Model, & Sieve/ECM Node Pricing

- **Single Special-`q` Lattice `L_Q` (`det = q`):** Exact integer Bareiss-LLL (`lll_16d_int(99, 100)`, `det(R) = q`) produces a 16-vector `Z`-basis of `L_Q`. Each node scans `S_16(R_n)` (`cells_rank2`, live-verified in `verify_b16_dsieve_live`), applying `g`-first threshold gating (`sieve_gateG` or `ecm_gateG`, whichever is cheaper) and exact unit-duplicate thinning (`udup_node = 0.50`, verified by `unitdomain.py`).
- **Norm Sextiles & `F_known = F_2core` Missing-Log Thinning:** Norms are evaluated conservatively at the upper bin endpoint `c_hi = (hi + log2(cells))/16` over `qbins` / `qbins_sub` (verified live on `24 + 12` production 381-bit Bareiss-LLL lattices in `prod_sq24_live` and `prod_cond_q15_live`, plus `120` production encountered children in `q15_enc_bias_live`: `156/156` strictly below `qbins`, `>= 10x` modeled yield headroom). Relations touching `F_full \ F_2core` (missing density `1 - peel_h`, `peel_h = 0.256845`) are discarded via explicit `exp(-miss_n)` thinning at every node and at boot.
- **Cofactorization & `P_71` Product Tree:** Surviving candidates strip primes `<= 2^71` via Bernstein's `P_71` tree (`128` owners `< 2^62 B`/owner, build `<= 2^73.369` u `<= 2^75.5` u), followed by ECM (`ecm_success = 0.999`) and Miller-Rabin.

## 3. Boot (`G_T` Initial Splitting, `K_boot = 1`) & Live Quadrature

- **Boot (`K_boot = 1`, `n1 = 265`, `nA = 245.0`, `kA = 12`, `kB = 1`, `cells_b = 1`, `m_sv = 2592.0`):** For each of the two `G_T` logs, the algorithm draws fresh i.i.d. uniform random exponents `s in Z_r^*` until the **first** semismooth boot relation is found (capped at `T_boot = ceil(20 / p_boot_trial)` trials; because the `s` draws are i.i.d. uniform before any descent tree is entered, failing to find a boot relation has probability `(1 - p_boot_trial)^T_boot <= e^(-20) < 2.1e-9` per log, `<= 4.2e-9` for both logs). Once that single boot relation is found, the algorithm runs **exactly one (`K_boot = 1`)** 1-shot descent tree!
- **Expected Two-Log Descent Work (`audit_cm.py`):**
  - Default mesh `32/185`: `min_mu = 48.19 >= 30.0`, `max_p0 = 0.001588`, `N_tree_1log = 3,157`, `P_2log_Boole = 0.055100`, **`two_log_exp_units = 106.638482`** (`boot = 104.1963`, deterministic `tree = 106.3451`).
  - Refined 12x mesh `384/370` (lower sums): `min_mu = 106.43`, `max_p0 = 0.000340`, `P_2log_Boole = 0.044367`, **`two_log_exp_units = 106.602643`**.
  - Total descent failure charge (pure Boole + `2*e^(-20)`): **`0.055100 <= 0.055300`**; ledger work charge `<= 106.639482` log2 units.
