# `ADAPTIVE.md` — Pairwise-`O_H`-Coprime `FB`-Blinded CRT Descent (`p_sat = 0`), Quantitative Paley–Zygmund + Newton–Maclaurin Bridge (`H5-Blind`), & Direct-`R` Multi-Regime `Q(zeta_15)` Diagnostics

## 1. Pairwise-`O_H`-Coprime `FB`-Blinded CRT First-Hit (`k = 1`) Descent + Quantitative Bridge (`H5-Blind`)

Production BLS12-381 `Q(zeta_15)` descent replaces unblinded `p_0 ~ 0.001` empty-ball assumptions (and avoids sibling-correlated `THRESHOLD-REPAIR` `psi^k`) with **without-replacement `FB`-blinded CRT first-hit (`k = 1`) descent** (`39` bins / `78` nodes, `R in [52, 92]`, `N_tree_max_1log = 2,729` deterministic slots), wrapped in `i.i.d.` fresh-`e` boot retries and unified meter `t_run = 2^{123.0}` u (`R_tree = 40 i.i.d. restarts, p_desc_loss <= 1.0000e-18 < 2^-57`):

1. **Pairwise-`O_H`-Coprime Blinding Pool `B_blind` (`|B_blind| = 192` in `[2048, 8192]`) & Exact Direct `Z`-Basis `R[0..15]` (`p_sat = 0` Proved):**
   - `check_algebra_bls12` verifies live that `[2048, 8192]` (`q' == 1 mod 15`) has `0` ramified roots, `684` relative-degree-1 `f`-side ideals, and **`192 >= 192` distinct `O_H`-prime ideals `q'_OH`** admitting an `f`-side degree-1 ideal. Post-LA setup fixes `B_blind` = `192` solved `f`-side ideals above `192` distinct `O_H`-primes (`gcd_{O_H}(q'_{OH,i}, q'_{OH,j}) = (1)` for `i != j`; expected 2-core multiplicity `~2^{50} >> 2`).
   - For any reached ideal `Q` (`q >= 2^{71}`) on side `s` and draw `Q'_t in B_blind` (`q'_t <= 8192`), `gcd(q, q'_t) = 1` unconditionally; exact CRT HNF `B` (`det_Z(B) = M = q * q'_t < 2^{278}`, exact initial GSO `d_1 = ... = d_{16} = M^2 < 2^{560}`) reduced by `lll_16d_int(99, 100)` (`assert dsq[16] == M * M`, `<= 2^{39.28} < 2^{40}` Ir) spans `L_{Q, Q'_t} = L_Q^{(s)} cap L_{Q'_t}^{(f)}` with index `1` (`p_sat_abort_proved = 0.0`).
   - Every `(a, b) = v * R` has `q'_t | N(F(a, b))` and `q | N(S(a, b))`, while `Q'_t` already has a solved virtual log from BW LA (`0` extra child slots!).
2. **Quantitative Three-Step Bridge from `H3` Mean Yield `mu >= 32.15` to Slot Failure `<= 2^{-24}`:**
   - **Step A — Conditional Second-Moment / Paley–Zygmund Bridge to `p_blind <= 0.50` (`H5-Blind`):** For fixed post-LA environment `E` and reached ideal `Q` in bin `b`, let `M = M(Q, Q', E)` be the retained qualifying relation count in `S_16(R)` on `L_{Q, Q'}` for `Q' ~ Unif(B_blind)`, with `mu = E[M | Q, E] >= 32.15` (`384/370`: `32.76`, already discounted by `ret_node = 0.90` and `udup_node = 0.40`) and cross-blinding overdispersion `kappa_blind(Q, E) := (Var(M | Q, E) - mu) / mu^2`. By Cauchy–Schwarz (`E[M]^2 <= P(M >= 1) E[M^2]`):
     $$p_{\text{blind}}(Q, \mathcal{E}) := \mathbb{P}_{Q' \sim \mathrm{Unif}(\mathcal{B}_{\text{blind}})}(M = 0 \mid Q, \mathcal{E}) \le \frac{\frac{1}{\mu} + \kappa_{\text{blind}}}{1 + \frac{1}{\mu} + \kappa_{\text{blind}}}.$$
     At our audited floor `mu >= 32.15`, **`p_blind(Q, E) <= 0.50` holds whenever `kappa_blind(Q, E) <= 1 - 1/32.15 = 0.9689`** (unlike unblinded `p_0 <= 0.001` at `mu = 889`, which required `kappa <= 0`!). We explicitly declare exploratory STNFS premise **`H5-Blind`** that every reached slot `(j, b)` satisfies `p_blind(Q_{w_j}, E) <= p_blind_cap = 0.50` (`50%` empty balls allowed!).
   - **Step B — Unconditional Newton–Maclaurin Without-Replacement Amplification (`0.50 -> 2^{-24}`):** Sampling `K_blind = 24` distinct `(Q'_1, ..., Q'_{24})` uniformly without replacement from `B_blind` (`|B_blind| = 192`) using fresh slot coins in ANY fixed environment `E`, let `f_i = 1_{M(Q, Q'_i, E) = 0} in {0, 1}` (`1 <= i <= 192`) with `e_1(f)/192 = p_blind(Q, E) <= 0.50`. By **Newton–Maclaurin's deterministic algebraic inequality** (verified in exact `Fraction` arithmetic in `pm_adaptive.py`, requiring **ZERO independence across draws in fixed `E`**):
     $$\mathbb{P}(\text{all } t \text{ draws empty} \mid Q, \mathcal{E}) = \frac{e_t(f)}{\binom{192}{t}} \le \Bigl(\frac{e_1(f)}{192}\Bigr)^t \le 0.50^t \quad (1 \le t \le 24),$$
     proving `P(slot fail | A_{j,b}, E) <= 2^{-24} = 5.9605e-8` (exact binary hypergeometric `C(96,24)/C(192,24) = 2.05e-8 < 2^{-25.5}`) and expected draws `sum_{t=0}^{23} P(fail first t) <= 1/(1 - 0.50) = 2.0` (`own = 2.0 * own_1` charged!).
   - **Step C — First-Hit `k = 1` Coarse-Bin Boole + Wald Boot:** Taking the first qualifying relation (`k = 1`, zero sibling retries if a descendant fails!) over `|T| <= 2,729` deterministic bit-sum knapsack slots gives **`P_NM_tree = 0.0001627, P_fail_1 = eps_exc + P_NM_tree + 0.125 = 0.27516 < 0.35`** and default `32/185` **`E[two_logs] <= 2^{109.311372} <= 2^{109.312400}`** u (`384/370`: **`2^{109.290306}`** u; `boot = (246, 211.0, 16)`).

## 2. Scoped Multi-Regime Quantitative Diagnostics (All Recomputed LIVE in `audit_cm.py`, 100% Direct `16`-Vector `Z`-Scan!)

### Pillar 1 — Direct-`R` Full-Depth Without-Replacement `FB`-Blinded CRT `Q(zeta_15)` `(3,3)` Descent + `kappa_blind` Paley–Zygmund Bridge (`q15_blind_h5.b85` + Live 8-Root Run)

Built on the exact BLS12-381-polynomial-structured `(3,3)` tower (`f_3 = a^3 - (2+w)b^3`, `g_3 = a^3 + 2 w^2 b^3`, `B_FB = 8,192`, `|B_blind_toy| = 64` solved `f`-side ideals, `K_blind_toy = 8` without-replacement draws vs production `K_blind = 24`, direct `16 x 16` CRT `Z`-basis `Rr[0..15]` with `assert dsq[16] == (q * q')^2`, strict squarefree `_fch_strict_sf`, unique `_ch_root_unique`, explicit `anc` exclusion):
- **120-Root Full-Depth Campaign:** **`119 / 120` roots (`99.17%`) complete to `100%` `<= 8,192` FB leaves** (`max_depth = 7`, `202` `d1` + `425` `d2+` slots, **`mean_scans = 1.5248 < 2.0`** even with `K_blind_toy = 8`!).
- **Cross-Blinding Retained-Relation `kappa_blind` & Paley–Zygmund Bridge (`40` `d1` + `40` `d2+` + `40` `ctrl` ideals x `8` draws = `960` CRT balls):**

| Stratum | Ideals x Draws | `mean(M)` | `p_hat` mean (med) | Measured **`kappa_blind` med** | **PZ Upper at `mu = 32.15`** (`vs 0.50` cap) |
|---|---:|---:|---|---:|---:|
| `d1` (reached depth 1) | `40 x 8` | `1.37` | `0.2625` (`0.250`) | **`0.0000`** (`<< 0.9689`) | **`0.0322`** (`15.5x` below `0.50`) |
| `d2p` (reached depth `2..7`) | `40 x 8` | `1.06` | `0.3875` (`0.375`) | **`0.0635`** (`<< 0.9689`) | **`0.0883`** (`5.7x` below `0.50`) |
| `ctrl` (fresh control) | `40 x 8` | `2.58` | `0.1187` (`0.125`) | **`0.0794`** (`<< 0.9689`) | **`0.1013`** (`4.9x` below `0.50`) |

- **3 Live Depth-2 CRT Witnesses + Live 8-Root Full-Depth `FB`-Blinded CRT Run Inside `audit_cm.py`:** Re-executes `lll_16d_int(99, 100)` (`dsq[16] == (q*q')^2`), exact divisibility `nf % q' == 0` & side `% q == 0`, and strict squarefree resolution into `<= 8,192` FB on 3 witnesses, PLUS runs a fresh 8-root full-depth without-replacement `FB`-blinded CRT `k = 1` descent live (`8/8` completed, `0` aborts, `max_depth = 5`, `mean_scans < 2.0`).

### Pillar 2 — Literal Degree-`(72, 24)` `Q(zeta_15)` `(9,3)` Direct `16`-Vector `Z`-Basis Scans (`q15_93_b16_live`, ZERO 2-Vector `O_H`-Span!)

- **100% Direct `16`-Vector `Z`-Basis Regeneration:** `verify_q15_93_b16_live` regenerates BOTH unblinded (`directR_f0_cnt = 7`, `dsq[16] == p^2`) AND CRT-blinded (`blind93[0][0] = 6`, `dsq[16] == (p * q')^2`) literal `(9,3)` retained-relation counts **exclusively via the production direct 16-vector scan `(a, b) = v * Rr[0..15]`** (`direct_16vec_Z_scan_only = true`, ZERO 2-vector `O_H`-span!).
- **`18` Literal `(9,3)` Ideals (`6` `f`, `6` `g`, `6` `enc`) x `8` CRT-Blinding Draws = `144` Retained-Relation CRT Balls:** Only `14 / 144` empty draws (**`p_hat_mean = 0.0972 << 0.50`**), **`kappa_blind_med = 0.0000`**, **`kappa_blind_p90 = 0.2857 < 0.9689 = kappa_crit_PZ`**, giving Paley–Zygmund upper bound **`0.0322 << 0.50`** at `mu = 32.15` (alongside the `220` unblinded `(9,3)` balls with `k_hat in [0.0328, 0.1862] << 0.9689` and `0/100` empty in `hi100`).

### Pillar 3 — `36/36` Live-LLL 381-Bit Production CRT Lattices `L_{Q, Q'}` (`prod_cond_q15_live` + `prod_sq24_live`, `0` Presupplied Blinding Ideals)

- **`prod_enc12` (`12/12` CRT `det == q * q'`):** Live parent divisibility (`12/12`), live-discovered `f`-side blinding ideals in `[2048, 8192]` (`presupplied_blinding_ideals = 0`), live Bareiss `99/100` LLL (`dsq[16] == (q*q')^2`), exact CRT divisibility `nf % q' == 0` & side `% q == 0`, sextile headroom **`38.25`b `f` / `30.32`b `g`**, and modeled joint yield **`121.77x` min / `366.77x` med** above charged `pr`.
- **`prod_sq24` (`24/24` CRT `det == q * q'`, `73b..258b` x `2048..8192`):** Live Bareiss `99/100` LLL (`dsq[16] == (q*q')^2`), exact CRT divisibility across all `384` sampled vectors, sextile headroom **`38.53..41.27`b `f` / `27.28..27.50`b `g`**, `sv_emp / sv_plan <= 0.121`, and modeled joint yield **`107.57x..169.94x` min / `414.45x..516.13x` med** above charged `pr`, plus pure-`Decimal(50)` `unitdomain.py`.
- **Honest Scope:** Toy retained-relation `kappa_blind` campaigns and 381-bit CRT modeled-yield checks are scoped compatibility support for the quantitative Paley–Zygmund bridge, NOT unconditional proofs of `H5-Blind`.


---
### Active H5-Blind(a, b) Operating Parameters & Deterministic Tail Closure (`lx = 71.0`, `R_tree = 40`)

- **Blinding Prime Pool `B_blind`:** Exact poly-GCD check in `audit_cm.py` confirms `192` pairwise-`O_H`-coprime unramified primes (`ram = 0`) in `[2048, 8192]` (`N_blind_pool = 192`).
- **Declared MODEL Inequality (Supported, Not Proved):** `H5-Blind(a)` bounds the exceptional-tree fraction by `eps_exc = 0.15`; `H5-Blind(b)` gives bin-wise smoothness on non-exceptional trees (`48` live 381-bit CRT LLLs confirm `CV^2_{Q'} < 0.15`).
- **Deterministic Caps & Restarts (No `k_bt=4` sibling retry, No separate environment-abort fudge):** Per-node scan cap `T_scan(b) = 2.5 * own_1(b)`, per-tree budget `T_tree = 8.0 * W_tree_1`, `R_tree = 40` `i.i.d.` fresh-`e'` restarts. Single-attempt failure `P_fail_1 = eps_exc + P_NM_tree + 0.125 = 0.27516 < 0.35`, residual failure `p_desc_loss <= P_fail_1^40 <= 1.0000e-18 < 2^-57` (explicitly included in `fail_core`), expected two-log work `2^109.311372 <= 2^109.312400` u (`384/370`: `2^109.290306` u).
