# Real-2-Core-FB `Q(zeta_15)` Multi-Level Descent Diagnostic (Supplementary to `k = 1` Boole Descent)

> **Scope:** Production descent (`DESCENT.md` §3) uses **pure `k = 1` zero-retry Boole descent** with fresh `i.i.d.` boot retries — NO `THRESHOLD-REPAIR`, NO `F_i | F_1..F_{i-1}`, NO `C_rep` (`p_node_abort = p_repeat_loss = 0` proved). `verify_threshold_counterexample` shows LIVE why naive `G(psi)` fails under correlated `(M, fail)` (`0.10 > 0.03571`).

## 1. Real-2-Core-FB `Q(zeta_15)` Campaign (`27,819` Roots, `18,023` First Attempts)

Uncensored multi-level descents over `O_H = Z[zeta_15]` (`f_t = X^3 - (2+w)`, `g_t = X^3 + 2 w^2`, `Z[omega]` `<tau>`-orbit separation, real 2-core FB `B_FB = 2^14`, `lll_16d_int(99, 100)`, `576`-cell `P^1`-unique ball) exhausting all `27,819` degree-1 `<tau>`-orbits of `20..21`b across 6 settings (`18,023` first attempts / `21,961` `k=1` children + `1,178` after-fail attempts):

| Setting | Roots | First `n` (`fail`) | After-fail `n` (`fail`) | `O/E fail` | `Work O/E` | Shared prev/run |
|---|---:|---:|---:|---:|---:|---:|
| `f20_d1` | 3,199 | 1,967 (`0.148`) | 100 (`10`) | `0.605` | `1.007` | `0 / 0` |
| `g20_d1` | 3,239 | 2,133 (`0.144`) | 90 (`6`) | `0.447` | `1.058` | `0 / 0` |
| `f21_d1` | 6,054 | 4,256 (`0.158`) | 251 (`40`) | `0.862` | `1.021` | `0 / 0` |
| `g21_d1` | 6,074 | 4,392 (`0.165`) | 294 (`42`) | `0.734` | `0.997` | `0 / 0` |
| `f20_d2` | 3,199 | 1,781 (`0.177`) | 123 (`15`) | `0.630` | `1.002` | `0 / 0` |
| `f21_d3` | 6,054 | 3,494 (`0.229`) | 320 (`81`) | `0.959` | `1.016` | `0 / 0` |
| **Pooled** | **27,819** | **18,023 (`0.173`)** | **1,178 (`194`)** | **`0.802`** | **`1.013`** | **`0 / 0`** |

## 2. Live Verification & Scope

`adaptive_h5_live` verifies `threshold_counterexample`, replays all `1,178` rows & 12 concentration bins (`C_bday <= 2.815`, `C_emp <= 19.11`, `0/1,178` shared descendants), and **regenerates 3 complete multi-level witnesses from scratch** (`q = 526531, 531871, 1063351`) bit-exactly via `lll_16d_int(99, 100)`. Together with `q15_samefield_descent_live` (`1/120 = 0.0083` empty vs `0.362`), `prod_sq24_live` (`>= 91.84x / 195.48x` in `B_16(8)`), `prod_cond_q15_live` (`>= 48.57x`), and live complete-orbit `yr >= 89.46x`, this gives exploratory support (not proof) for `H5-Empty`.
