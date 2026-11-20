# Real-2-Core-FB `Q(zeta_15)` Multi-Level Descent Diagnostic (Supplementary to Pure `k = 1` Boole Descent)

> **Scope note:** Our production descent (`DESCENT.md` §3) uses **pure `k = 1` zero-retry Boole-union-bound descent** with fresh `i.i.d.` boot exponent retries — it does **NOT** use `THRESHOLD-REPAIR`, does **NOT** assume any sequential post-failure hazard `P(F_i | M >= i, F_1..F_(i-1), history) <= psi`, and does **NOT** assume any `C_rep` adaptive repeat-collision bound (`p_node_abort = 0` and `p_repeat_loss = 0` are proved deterministically by strict-DAG memoization and Boole's inequality!). `pm_adaptive.verify_threshold_counterexample` demonstrates LIVE why naive `G(psi)` fails under correlated count/solvability (`M = 1` failing prob `0.10`, `M = 100` succeeding prob `0.90`: `P(fail) = 0.10 > G(0.10) = 0.03571`), motivating our strict `k = 1` zero-retry Boole-union-bound architecture.

## 1. Real-2-Core-FB Multi-Level Campaign over `Q(zeta_15)` (`27,819` Roots, `18,023` First Attempts)

We run uncensored multi-level recursive descents over the EXACT production base ring `O_H = Z[zeta_15]` using the reduced `(f_3, g_3)` tower (`f_t = X^3 - (2 + w)`, `g_t = X^3 + 2 w^2`), exact `Z[omega]` `<tau>`-orbit column separation (`c0^2 - c0 c1 + c1^2`), the real 2-core factor base `ret_f, ret_g` (`B_FB = 2^14`), pure-integer 2-stage Bareiss `lll_16d_int(99, 100)`, and the complete deterministic `P^1`-unique, zero-unit-duplicate ball `alpha = 1, beta in B_8(sqrt(3))` (`576` cells). Root ideals exhaust all `27,819` eligible degree-1 `<tau>`-orbits of bit-length `b_root in {20, 21}` across 6 settings (`d in {1, 2, 3}`), yielding `18,023` first attempts (`21,961` first-relation `k = 1` encountered children) plus `1,178` supplementary after-failure attempts (`194` failures):

| Setting | Roots | First `n` (`fail`) | After-fail `n` (`fail`) | `O/E fail` | `Work O/E` | Shared prev/run |
|---|---:|---:|---:|---:|---:|---:|
| `f20_d1` | 3,199 | 1,967 (`0.148`) | 100 (`10`) | `0.605` | `1.007` | `0 / 0` |
| `g20_d1` | 3,239 | 2,133 (`0.144`) | 90 (`6`) | `0.447` | `1.058` | `0 / 0` |
| `f21_d1` | 6,054 | 4,256 (`0.158`) | 251 (`40`) | `0.862` | `1.021` | `0 / 0` |
| `g21_d1` | 6,074 | 4,392 (`0.165`) | 294 (`42`) | `0.734` | `0.997` | `0 / 0` |
| `f20_d2` | 3,199 | 1,781 (`0.177`) | 123 (`15`) | `0.630` | `1.002` | `0 / 0` |
| `f21_d3` | 6,054 | 3,494 (`0.229`) | 320 (`81`) | `0.959` | `1.016` | `0 / 0` |
| **Pooled** | **27,819** | **18,023 (`0.173`)** | **1,178 (`194`)** | **`0.802`** | **`1.013`** | **`0 / 0`** |

## 2. Live Verification & Exploratory Scope

`audit_cm.py` (`adaptive_h5_live` via `pm_adaptive.py`) runs LIVE on every default invocation: (i) verifies `threshold_counterexample`; (ii) replays all `1,178` stored rows and 12 child-concentration bins (`C_bday <= 2.815`, `C_emp <= 19.11`, `0/1,178` shared descendants) in `adaptive_h5.b85`; (iii) **regenerates 3 complete multi-level witness descents from scratch** (`f20_d1: q=526531`, `g20_d1: q=531871`, `f21_d3: q=1063351`) using pure-integer `lll_16d_int(99, 100)` and the real 2-core FB, asserting bit-exact match. Together with `q15_samefield_descent_live` (`3,728`-cell `Q(zeta_15)` balls: `1/120 = 0.0083` empty vs `1/(1 + 0.5*3.53) = 0.362`, `43.4x` below `H5-Empty`), `prod_sq24_live` (`>= 16.96x / 55.01x` joint-cell yield headroom), and `prod_cond_q15_live` (`>= 41.60x` encountered joint-cell headroom), this provides scoped exploratory support (NOT proof) for `H5-Empty` in `DESCENT.md`.
