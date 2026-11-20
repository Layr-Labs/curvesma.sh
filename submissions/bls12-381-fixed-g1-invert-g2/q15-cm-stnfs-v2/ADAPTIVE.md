# Conditional Next-Relation Evidence (`H5`), `THRESHOLD-REPAIR`, & Adaptive Child Anti-Concentration

`THRESHOLD-REPAIR` separates low-count events from sequential child failures without assuming count/solvability independence: let `F_i` be failure of relation `i`'s children and `A_i = {M >= i, F_1, ..., F_(i-1)}`. `H5` explicitly declares the sequential conditional hazard `P(F_i | A_i, pre-entry history) <= psi` (whenever `P(A_i) > 0`, where `psi` is the additive union bound over the relation's child failure envelopes) and a matching conditional per-attempt child-work envelope. Neither premise is inferred from `H-PGF` alone. By induction `P(A_(i+1)) <= psi P(A_i)`, so `P(M >= k, F_1..F_k) <= psi^k` and for every `k >= 1`:

    P(node fails) <= P(M < k) + psi^k <= min(1, G(t_k) / t_k^(k-1)) + psi^k,

where `G(t) = (1 + beta*(1 - t))^(-r)` (`beta = 0.5*mu/r`, `t in [0, 0.80]`), `T_1 = G(0)`, and for `k >= 2` the closed-form stationary point clipped to `[0, 0.80]` is `t_k = min(0.80, (k - 1)*(1 + 1/beta)/(k - 1 + r))`. Minimizing over `k = 1..64` and flooring at `0.01` gives `phi` (`cm_models.eval_descent_live`). Expected attempts satisfy `E[attempts] = sum_i P(A_i) <= 1/(1 - psi)`. `pm_adaptive.verify_threshold_counterexample` checks LIVE that naive `G(psi)` fails on a correlated count/solvability distribution (`M = 1` with prob `0.10` failing, `M = 100` with prob `0.90` succeeding: satisfies `H-PGF(0.5*60)` for all `t <= 0.80`, yet `P(fail) = 0.10 > G(0.10) = 0.03571`) while `THRESHOLD-REPAIR` gives `0.13226 >= 0.10`.

## 1. Deterministic-Order Real-2-Core-FB Adaptive Campaign over `Q(zeta_15)` (`1,178` After-Failure Attempts)

We run uncensored multi-level recursive descents over the EXACT production base ring `O_H = Z[zeta_15]` using the reduced `(f_3, g_3)` tower (`f_t = X^3 - (2 + w)`, `g_t = X^3 + 2 w^2`), exact `Z[omega]` `<tau>`-orbit column separation (`c0^2 - c0 c1 + c1^2`), the **real 2-core factor base `ret_f, ret_g` (`B_FB = 2^14`)**, pure-integer 2-stage Bareiss `lll_16d_int(99, 100)`, and the complete deterministic `P^1`-unique, zero-unit-duplicate ball `alpha = 1, beta in B_8(sqrt(3))` (`576` cells in fixed lexicographic order). Relations are tried in deterministic ball order, stopping at the first relation whose children ALL succeed; every recursive descent runs to completion (uncensored).

Root ideals exhaust all eligible degree-1 `<tau>`-orbits of bit-length `b_root in {20, 21}` (`27,819` roots across 6 settings) with root window `(14, b_root - 1, k = 2)`; a child of `b_c` bits uses strictly smaller window `(14, max(14, b_c - d), k = 1 if b_c - d > 14 else 0)` for shrink `d in {1, 2, 3}`. Child-row failure rates `p_row(setting, side, b)` and fresh mean subtree work `w_row(setting, side, b)` are calibrated **EXCLUSIVELY on children of first attempts (`18,023` first attempts with children)**; after-failure attempts (`1,178` attempts, `194` failures) NEVER enter calibration. For each after-failure attempt, `psi_hat = sum_c p_row(c)` and `W_hat = sum_c w_row(c)` (full charge across all children even if the first child fails).

| Setting | Roots | First `n` (`fail`) | After-fail `n` (`fail`, rate) | `O/E fail` (95% boot CI) | Ratio to first | `Work O/E` (95% CI) | Shared prev/run |
|---|---:|---:|---:|---:|---:|---:|---:|
| `f20_d1` | 3,199 | 1,967 (`0.148`) | 100 (`10`, `0.100`) | `0.605` (`0.294-0.992`) | `0.615` | `1.007` (`0.943-1.070`) | `0 / 0` |
| `g20_d1` | 3,239 | 2,133 (`0.144`) | 90 (`6`, `0.067`) | `0.447` (`0.146-0.773`) | `0.456` | `1.058` (`0.956-1.160`) | `0 / 0` |
| `f21_d1` | 6,054 | 4,256 (`0.158`) | 251 (`40`, `0.159`) | `0.862` (`0.630-1.112`) | `0.890` | `1.021` (`0.979-1.069`) | `0 / 0` |
| `g21_d1` | 6,074 | 4,392 (`0.165`) | 294 (`42`, `0.143`) | `0.734` (`0.518-0.945`) | `0.762` | `0.997` (`0.954-1.041`) | `0 / 0` |
| `f20_d2` | 3,199 | 1,781 (`0.177`) | 123 (`15`, `0.122`) | `0.630` (`0.334-0.915`) | `0.644` | `1.002` (`0.957-1.050`) | `0 / 0` |
| `f21_d3` | 6,054 | 3,494 (`0.229`) | 320 (`81`, `0.253`) | `0.959` (`0.789-1.138`) | `0.998` | `1.016` (`0.989-1.047`) | `0 / 0` |
| **Pooled** | **27,819** | **18,023 (`0.173`)** | **1,178 (`194`, `0.165`)** | **`0.802` (`0.703-0.905`)** | **`0.828`** | **`1.013` (`0.993-1.033`)** | **`0 / 0`** |

- **Conditional failure hazard:** Every single setting (6/6) has `O/E fail < 1.00` and `ratio to first < 1.00`; pooled `O/E fail = 0.8022` (95% root-cluster bootstrap CI `[0.703, 0.905]` strictly below `1.00`) and `ratio to first = 0.8279`.
- **Conditional child work:** Pooled `work O/E = 1.0132` (95% CI `[0.993, 1.033]`, consistent with fresh mean `1.0`; and well below production's conservative full-ball/upper-edge/`1/(1-psi)` node-work envelope).
- **Zero shared descendants:** Across all `1,178` after-failure attempts, **0** share any child `<tau>`-orbit with an earlier relation of the same node (`sh_prev = 0`) and **0** reuse any ideal visited earlier in the same descent run (`sh_run = 0`).

- **History-stratified breakdown (recomputed LIVE in `adaptive_h5_live`):** No history stratum exhibits runaway failure escalation: failure depth `d2` (1 prior fail, `n=1,104`): `O/E = 0.800`, `d3+` (`>=2` prior fails, `n=74`): `O/E = 0.828`; current relation child count `c1` (`n=843`): `0.862`, `c2` (`n=335`): `0.729`; hardest top-bin child (`n=265`): `0.746` vs subtop (`n=913`): `0.834`; prior failed relation had `1` child (`n=828`): `0.673` vs `2` children (`n=350`): `1.112` (disclosed `+11%` above raw first-attempt empirical mean, well beneath production `THRESHOLD-REPAIR` formula slack where `36/36` production `(f_9, g_3)` lattices achieve `>= 16.96x` joint-cell yield headroom above charged `pr`).

## 2. Adaptive Child-Ideal Anti-Concentration Diagnostic (`C_rep = 256`)

Across all `23,474` adaptively selected encountered children in all 12 `(side, bit)` bins (`b = 15..20`), we measure both the empirical max-hit concentration `C_emp = (max_hit / tot) * N_elig` and the exact unbiased pairwise birthday collision factor `C_bday = [sum_q c_q(c_q - 1) / (tot*(tot - 1))] * N_elig` against the exact count `N_elig` of eligible degree-1 `<tau>`-orbits in `[2^(b-1), 2^b)`:
- `f`-side (`b = 15..20`): `C_emp = 1.74, 2.88, 3.33, 4.10, 10.36, 19.11`; **`C_bday = 1.070, 1.116, 1.182, 1.237, 1.462, 2.258`**.
- `g`-side (`b = 15..20`): `C_emp = 2.18, 2.74, 3.25, 5.21, 7.63, 15.36`; **`C_bday = 1.081, 1.139, 1.141, 1.363, 1.723, 2.815`**.

Every bin exhibits near-uniform pairwise collision concentration (`C_bday <= 2.815`, `C_emp <= 19.11`), over an order of magnitude below production's charged `C_rep = 256` against the conservative `N_lx = 2^70/(12*70*ln 2) = 2^60.81` floor (which is itself `59x` below `[70, 76)`'s `2^66.69` orbits, giving `> 15,000x` effective anti-concentration allowance).

## 3. Live Verification & Exploratory Scope

`audit_cm.py` (`adaptive_h5_live` via `pm_adaptive.py`) runs LIVE on every default invocation: (i) verifies `threshold_counterexample`; (ii) replays all `1,178` stored after-failure rows in `adaptive_h5.b85`, recomputing all 6 settings' `O/E fail`, `rat1`, `work O/E`, `0/0` shared descendants, and all 12 concentration bins; (iii) **regenerates 3 complete multi-level after-failure witness descents from scratch** (`f20_d1: q=526531`, `g20_d1: q=531871`, `f21_d3: q=1063351` — exercising both succeeding and failing post-failure attempts!) using pure-integer `lll_16d_int(99, 100)` and the real 2-core FB, asserting bit-exact match. **Exploratory scope:** This campaign uses degree `(3, 3)` polynomials, `B_FB = 2^14`, and `15..21`-bit ideals; transfer of `H-PGF`, sequential `psi`, conditional child work, and `C_rep` to production `(f_9, g_3)` and `>= 70`-bit ideals remains an explicitly declared exploratory `H5` hypothesis (not a proof).
