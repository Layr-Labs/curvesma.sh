# Count Tails and Sequential Availability (`H5` Repair)

This repairs recursive node-failure accounting without assuming count/descendant independence or using the invalid `PGF(psi)` substitution.

## 1. Half-Scale NB-PGF Envelope and Separate Sequential Hazard

At each reached special-`q` node, let `M` be the number of distinct qualifying relations in the node's lattice ball and `mu` its lower-sum modeled expectation (`DESCENT.md`). Under `H5`:
1. **Count-PGF envelope (`r = 2`, `alpha = 0.5`, `beta = 0.5 mu / r = 0.25 mu`):** For the pre-entry reached-ideal ensemble,
   `E[t^M] <= G(t) = (1 + beta (1 - t))^(-r)` for all `t in [0, 0.8]`.
   `DISC` and pooled per-type confidence arguments are removed (`DISC = 0`).
2. **Explicit sequential availability hazard:** Let `F_i` be failure of candidate relation `i`'s children and `A_i = {M >= i, F_1, ..., F_(i-1)}`. We explicitly assume `P(F_i | A_i, prehistory) <= psi` whenever `P(A_i) > 0`, where `psi` is the additive sum of the candidate's child-failure bounds. This conditional hazard and the separate conditional per-attempt child-work bound remain **disclosed transfer heuristics** — neither is inferred from `G(t)` alone.

## 2. Elementary Union Bound `T_k + psi^k`

Since `A_(i+1) subset A_i intersect F_i`, induction gives `P(A_(k+1)) = P(M >= k, F_1, ..., F_k) <= psi^k`. Any failing node either has `M < k` or fails its first `k` relations, so for every integer `k >= 1`:

    P(node fails) <= P(M < k) + psi^k.

Since `t^M >= t^(k-1)` on `{M < k}`, Markov's inequality gives `P(M < k) <= G(t) / t^(k-1)` for every `t in (0, 0.8]`. Setting `T_1 = G(0)` and, for `k > 1`, `t_k = min(0.8, (k - 1)(1 + 1/beta) / (k - 1 + r))`, `T_k = min(1, G(t_k) / t_k^(k-1))`, we evaluate live (`cm_models.phi_nb_pgf_union`):

    phi_node = max(0.02, min(1, min_(k=1..64) (T_k + psi^k + 1e-10))).

Expected relation attempts obey `E[attempts] = sum_(i>=1) P(A_i) <= 1 / (1 - psi)`, preserving `own_full_ball + child_work / (1 - psi)` with **100% full survivor work** (`1.0x`, zero lazy-ordering discount).

## 3. Refutation of `PGF(psi)` and Live Verification

The old substitution `P(fail) <= PGF(psi)` is invalid: take `M = 1` with probability `0.1` (and sole relation failing) and `M = 100` with probability `0.9` (immediate success). For `r = 1, beta = 30`, `E[t^M] <= G(t)` holds for all `t in [0, 0.8]`, yet `P(fail) = 0.1 > G(0.1) = 1/28 = 0.035714`. The repaired union bound `min_(k=1..64)(T_k + psi^k) = 0.132258 >= 0.1` separates the correlated low-count event and holds strictly. `audit_cm.py` (`threshold_and_adaptive_live.counterexample`) replays and verifies this counterexample live.
