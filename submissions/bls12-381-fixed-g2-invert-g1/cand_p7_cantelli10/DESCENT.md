# `DESCENT.md` — Pure `k = 1` Boole + Unbounded Fresh-`e` i.i.d. Boot Retries (`Wald` Identity, `lx = 71.0`)

## 1. Architecture: Zero Sibling Retries (`k = 1`), Pure Boole Union Bound, & Exact Wald Identity

Our individual-logarithm descent uses the **`pure_k1_boole_wald_iid_boot`** architecture:

1. **Strict Single-Pass Nodes (`k = 1`, Zero Sibling Retries):** For every special-`q` prime ideal `Q` above `lx = 71.0` (`39` five-bit bins `[71, 76), ..., [261, 266)` on `f` and `g` sides = `78` nodes), the algorithm performs **one** deterministic scan of its 16D Bareiss-LLL ball `S_16(R_n)`, takes the **first** surviving relation whose large-prime children lie in the node's permitted windows, and **never retries** a sibling relation if any child fails (instead aborting the current attempt immediately).
2. **Deterministic Worst-Case Tree Size (`N_tree_1log = 2,348`, `N_tree_2log = 4,696`):** Because there are zero sibling retries, every node's worst-case subtree size `N(b, s)` and worst-case tree work `W(b, s)` satisfy exact deterministic knapsack recurrences (`W(b, s) = own_1(b, s) + max_children sum W(child)`, `N(b, s) = 1 + max_children sum N(child)` -> at most **`2,348` slots per log**). No `N_MAX` node counter is needed, and repeat-encountered ideals simply reuse their single memoized relation with zero additional failure probability.
3. **Pure Boole Union Bound Over Selected First-Survivor Slots (ZERO Subtree or Sibling Independence!):** Under our **explicitly declared selected-slot conditional `H5` premise** (`P(M(I) = 0 | reached first-choice slot (b, s)) <= p_0(b, s)` using one-sided Cantelli `p_0(b, s) = max(0.0001, (1 + 1 * mu)^(-1))`; see `ADAPTIVE.md`), **Boole's union bound** (which holds under **arbitrary dependence** across siblings, parent-child pairs, and ancestors!) bounds the per-attempt tree failure probability by:
   - `P_sub(b, s) = p_0(b, s) + max_children sum P_sub(child)`,
   - **`Phi_tree = P_boot_children = 0.467995 <= 0.60`**.
4. **Unbounded Fresh-`e` i.i.d. Boot Retries (Exact Wald Identity -> `p_descent_fail = 0.0`!):** For each of the two `G_T` logs, the outer loop draws a **fresh independent uniform `e_i in Z_r^*`**, runs the boot splitter on `target^(e_i)` until a semismooth boot relation is found (expected cost `boot_own_1 = (1/p_boot_trial)(LLL_boot + trial_b)`), and runs ONE strict `k = 1` tree (worst-case cost `W_tree`). If any slot is empty (probability `<= Phi_tree`), it discards the attempt and draws a fresh independent `e_(i+1) in Z_r^*`. Because fresh-`e_i` attempts are **genuinely i.i.d.** and unbounded until success, **Wald's exact identity** gives:
   - Expected attempts per log `<= 1 / (1 - Phi_tree) = 1.8797`,
   - Expected two-log descent work `<= (2 / (1 - Phi_tree)) * (boot_own_1 + W_tree) / 20000`,
   - **`p_descent_fail = p_node_abort = p_repeat_loss = 0.0`**! The only termination cutoff is the unified execution meter `t_run = 2^123.0` u, whose Markov failure probability `p_meter <= E[W_open] / t_run <= 0.000330284183` is charged directly in §6!

## 2. Special-`q` Lattice, Norm Model, & Sieve/ECM Node Pricing

- **Single Special-`q` Lattice `L_Q` (`det = q`):** Exact integer Bareiss-LLL (`lll_16d_int(99, 100)`, `det(R) = q`) produces a 16-vector `Z`-basis of `L_Q`. Each node scans `S_16(R_n)` (`cells_rank2`, live-verified in `verify_b16_dsieve_live`), applying `g`-first threshold gating (`sieve_gateG` or `ecm_gateG`, whichever is cheaper) and exact unit-duplicate thinning (`udup_node = 0.50`, verified by `unitdomain.py`).
- **Norm Sextiles & `F_known = F_2core` Missing-Log Thinning:** Norms are evaluated conservatively at the upper bin endpoint `c_hi = (hi + log2(cells))/16` over `qbins` / `qbins_sub` (verified live on `24 + 12` production 381-bit Bareiss-LLL lattices in `prod_sq24_live` and `prod_cond_q15_live`, plus `120` production encountered children in `q15_enc_bias_live`: `156/156` strictly below `qbins`). Relations touching `F_full \ F_2core` (missing density `1 - peel_h`, `peel_h = 0.256845`) are discarded via explicit `exp(-miss_n)` thinning at every node and at boot.
- **Cofactorization & `P_71` Product Tree:** Surviving candidates strip primes `<= 2^71` via Bernstein's `P_71` tree (`128` owners `< 2^62 B`/owner, build `<= 2^73.369` u `<= 2^75.5` u), followed by ECM (`ecm_success = 0.999`) and Miller-Rabin.

## 3. Boot (`G_T` Initial Splitting) & Live Quadrature

- **Boot (`n1 = 265`, `nA = 225.0`, `kA = 12`, `kB = 1`, `cells_b = 1`, `m_sv = 2592.0`):** Deterministic 16D Bareiss-LLL Lovász bound `2591.9b < m_sv = 2592.0b` verified live in `voronoi_and_boot_live` (`pr_b = 2^-44.8138` under `H1`).
- **Expected Two-Log Descent Work (`audit_cm.py`):**
  - Default mesh `32/185`: `min_mu = 103.49`, `max_p0 = 0.009570`, `N_tree_1log = 2,348`, `Phi_tree = 0.467995 <= 0.60`, `Wald mult = 1.8797x`, **`two_log_exp_units = 110.296322`** (`boot = 107.5824`, `tree = 110.0577`).
  - Refined 12x mesh `384/370` (lower sums): `min_mu = 255.25`, `max_p0 = 0.003902`, `Phi_tree = 0.448683`, **`two_log_exp_units = 110.202740`**.
  - Descent failure charge: **`0.0`** (cut off ONLY by unified meter `t_run = 2^123.0` u); ledger work charge `<= 110.297322` log2 units.
