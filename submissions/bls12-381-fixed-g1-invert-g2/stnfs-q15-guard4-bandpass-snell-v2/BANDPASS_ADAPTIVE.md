# Per-Column 2D Adaptive Thresholding, Exact 2D Panjer DP, 4M-MC & 200k-Chord Certification (`BANDPASS_ADAPTIVE.md`)

This companion document details the Stage 1a dual band-pass geometry, exact floor-quantized 2D Panjer smooth-recall DP, orthogonal 2D column feature $(N_g(b), \text{nf}_0(b))$, deterministic zero-wraparound proof (`max_byte(c) <= 255` for every column), zero-interpolation $N = 4{,}000{,}000$ 16D Monte Carlo yield certification, `200,000`-chord 16B-granular adaptive slab timing, and explicit kill-gates for pruned alternatives (`Guard2`, $q \le 0.60$, $K = 3$).

---

## 1. Exact 2D Column Feature $(N_g(b), \text{nf}_0(b))$ & Per-Column Adaptive Thresholding

In Mode-A $(4, 3)$ CM-STNFS over $K = \mathbb{Q}(\zeta_{15})$, the algebraic norms of $(a, b) \in \mathcal{A}_k \times \mathcal{B}_k$ along a fixed chord column $b$ satisfy:
1. **Column-Constant $g$-Norm:** Since $|x| = 2^{63.7143}$ dominates $|a|^3 / |b|^3$ by $\sim 10^{12}$, $N_g(a, b) = |N_{K/\mathbb{Q}}(a^3 - x b^3)|$ is **deterministic to $< 10^{-7}$ bits across all $V_{a,\text{sh}}$ rows `a_id` in column $b$**! Over `4,000,000` 16D MC draws,
   $$\max_{(a, b)} \bigl|\log_2|N_g(a, b)| - \log_2|N_g(0, b)|\bigr| = 4.327 \times 10^{-8}\text{ bits},$$
   re-verified live on real Mode-A $(a, b)$ vectors in `experiment.py`.
2. **Orthogonal Column Proxy $\text{nf}_0(b, k)$ for $N_f(a, b)$:** On the non-convex $L_{0.65}(b)$ ball, spiky boundary vectors $b$ have **small** $N_g(b) \propto \prod_{m=0}^3 |\sigma_m(b)|^6$ yet **large** $N_f(a, b) = \prod_{m=0}^3 |Q_{\sigma_m}(a^3, b^3)|^2 \propto \max_m |\sigma_m(b)|^{18}$. Defining the deterministic column-only proxy (using the shell RMS radius $A_{\text{rms},k}$):
   $$\text{nf}_0(b, k) := \sum_{m=0}^3 \log_2\!\bigl(A_{\text{rms},k}^{18} + A_{\text{rms},k}^{12}|\sigma_m(b)|^6 + A_{\text{rms},k}^6|\sigma_m(b)|^{12} + |\sigma_m(b)|^{18}\bigr)$$
   achieves $\text{RMSE}(\log_2|N_f(a, b)| - \text{nf}_0(b, k)) = 9.20$ bits and exposes a **`4.05x`** smooth-yield-to-volume contrast orthogonal to $N_g(b)$:

| $\text{nf}_0(b)$ Bin (bits) | Volume Fraction | Smooth-Yield Fraction | Smooth / Volume Ratio |
|---|---:|---:|---:|
| `[538, 542)` | `0.03470` | `0.10316` | **`2.973`** |
| `[542, 546)` | `0.06899` | `0.10485` | `1.520` |
| `[546, 550)` | `0.10685` | `0.12336` | `1.155` |
| `[550, 554)` | `0.16578` | `0.17028` | `1.027` |
| `[554, 558)` | `0.20161` | `0.17613` | `0.874` |
| `[558, 562)` | `0.17145` | `0.13507` | `0.788` |
| `[562, 566)` | `0.12587` | `0.09399` | `0.747` |
| `[566, 570)` | `0.07517` | `0.05512` | **`0.733`** |
| `[570, 574)` | `0.03486` | `0.02601` | `0.746` |
| `[574, 586)` | `0.01467` | `0.01200` | `0.818` |

3. **Neyman–Pearson Lagrange Optimal Per-Column Policy:** For each chord column $c \in [0, L)$ with deterministic features $(k, N_g(b_c), \text{nf}_0(b_c, k))$, precompute the exact wrap-safe integer byte pair $(T_g(c), T_f(c))$ maximizing
   $$\text{ret}(c; T_g, T_f) - \lambda \cdot s_g(T_g)\,s_f(T_f)$$
   using the exact 0.5-bit 2D Panjer floor-quantized smooth-recall table (computed once in Stage 0; $< 2^{25}$ u). At fixed Candidate #1 bands (`g = [32, 61.5]`, `f = [36, 61.0]`, `q = 3.8`, matched `ret = 0.685`), exact ROC comparison proves:
   - Best Global Constant $(T_g, T_f) = (69, 43)$: `log2(s1) = -23.4421` (`ret = 0.68624`).
   - 1D Adaptive $(T_g, T_f)(\text{shell}, N_g(b))$: `log2(s1) = -23.6599` (**`-0.2178` b** on `s1`, `-0.0136` b net cold-run).
   - **2D Adaptive $(T_g, T_f)(\text{shell}, N_g(b), \text{nf}_0(b))$:** `log2(s1) = -23.7680` (**`-0.3260` b** on `s1`, `-0.0238` b net cold-run — **1.75x larger** than 1D $N_g$ alone!).

---

## 2. Deterministic Per-Column Zero-Wraparound Proof (`max_byte(c) <= 255`)

Every sieved prime of norm $p$ in the band $[P_{\min}, P_{\max}]$ adds the **floor-quantized** byte weight $w_p = \lfloor (\log_2 p) / q_{\text{step}} \rfloor$ to a column-specific pre-biased byte `bias_u8[c] = 128 - T_u8(c)`. Because $\lfloor \cdot \rfloor$ is subadditive ($\sum_i \lfloor x_i \rfloor \le \lfloor \sum_i x_i \rfloor$), every cell $(a, b_c)$ in column $c$ unconditionally satisfies:
$$\text{byte}(a, b_c) \le \text{bias\_u8}[c] + \left\lfloor \frac{N_{\max}(c)}{q_{\text{step}}} \right\rfloor = 128 - T_{\text{u8}}(c) + \left\lfloor \frac{N_{\max}(c)}{q_{\text{step}}} \right\rfloor.$$

Enforcing the hard per-column deterministic floor on every lookup table entry:
- **$g$-side (`q_step_g = 3.4`, `g_band = [30.0, 61.3]`):** Since $N_g(a, b_c) \le N_g(0, b_c) + 10^{-6} \le N_{g,\max,\text{col}}(c)$ (with global domain AM–GM ceiling $N_g(a, b) \le \mathbf{686.839} < 705.0$ at $R_0 = 1405.74$), enforcing
  $$T_g(c) \ge \max\!\left(1,\; \left\lfloor \frac{N_{g,\max,\text{col}}(c) + 1.0}{3.4} \right\rfloor - 127\right)$$
  guarantees `0 <= bias_g[c] <= 127` and **`max_byte_g(c) <= 255` deterministically for every column $c$**! Even at the global worst-case column ($N_{g,\max} = 705.0$), $T_g \ge \lfloor 705.0 / 3.4 \rfloor - 127 = 207 - 127 = 80$, while our mean optimal threshold is $\overline{T}_g = 266.13 / 3.4 = 78.27$ ticks (columns with $N_g(b_c) \approx 650$–$670$ bits have per-column wrap floor $\lfloor 671/3.4 \rfloor - 127 = 70 \ll T_g^*(c)$, so $q_g = 3.4$ stays completely unconstrained at the Neyman–Pearson optimum!).
- **$f$-side (`q_step_f = 3.8`, `f_band = [34.0, 61.0]`):** Using the unconditional global corner ceiling $N_f(a, b) \le \mathbf{617.300} < 630.0 = N_{f,\max,\text{charged}}$ (re-verified over the full $(k, f_a, f_b)$ simplex in `experiment.py`), enforcing
  $$T_f(c) \ge \left\lfloor \frac{630.0}{3.8} \right\rfloor - 127 = 165 - 127 = \mathbf{38}$$
  guarantees `0 <= bias_f[c] <= 90` and **`max_byte_f(c) <= 255` deterministically for every cell in the entire sieve domain**! (Our mean optimal threshold is $\overline{T}_f = 177.23 / 3.8 = 46.64$ ticks $\gg 38$.)

---

## 3. Zero-Interpolation `N = 4,000,000` MC & `200,000`-Chord Certification (`CAND2_FINAL_ETA85_7000_NL205`)

At the locked operating point `CAND2_FINAL_ETA85_7000_NL205` (`K = 2, q = 0.65, rb = [1.0, 0.82], V_a = 0.9351 * 2^{56}`, `B_f = 70.00, B_g = 70.36, eta = 0.85`, `g_band = [30.0, 61.3] (q_g = 3.4)`, `f_band = [34.0, 61.0] (q_f = 3.8)`, $\lambda = 6.7108864 \times 10^5$):

| Metric | Direct `4M`-MC / `200k`-Chord Value | Charged / Padded Value | Notes |
|---|---:|---:|---|
| 2D Panjer Policy Retention `ret` | `0.678727` | `0.678727` | Direct 0.5-bit table inside 4M MC (zero interpolation!) |
| Survival $\log_2 s_g$ / $\log_2 s_f$ / $\log_2 s_1$ | `-14.5751` / `-9.0419` / `-23.6169` | `-14.5751` / `-9.0419` / `-23.6169` | Exact 2D-adaptive mean over all columns |
| Mean Adaptive Thresholds $(\overline{T}_g, \overline{T}_f)$ | `(266.13, 177.23)` bits | `(78.27, 46.64)` ticks | Deterministic `max_byte <= 255` everywhere |
| Direct `4M`-MC Volume $\log_2 V_{\text{kept,MC}}$ | `119.653105` (`SE = 0.001185` b) | `119.724637` | Includes `3*SE (+0.003554)` + `gap (+0.002538)` + `anchor (+0.065440)` |
| Pure First-Principles $\log_2 V_{\text{kept,pure}}$ | — | `119.659197` | Excludes only the `+0.065440` b comparator anchor pad |
| Retained vs Required Orbits | $R_{\text{ret}} = 2^{61.767486}$ | $R_{\text{req}} = 2^{61.767101}$ | `eta * U = 0.85 * 2^{62.0016}` |
| Outer $L_{0.65}$ Radius $R_0$ | `1405.7365` | `1405.7365` | Shell 1 radius $0.82\,R_0 = 1152.7040$ |
| Conditional Mean Norms $(\bar{N}_f, \bar{N}_g)$ | `(543.46, 671.02)` bits | `(543.46, 671.02)` bits | Observed maxima `(596.91, 686.82)` bits |
| Deterministic Domain Bounds $(N_{f,\max}, N_{g,\max})$ | `(617.30, 686.84)` bits | `(630.0, 705.0)` bits | Unconditional AM–GM & simplex corner proofs |
| `200k`-Chord Effective Window $L_{\text{win,eff}}$ | `111.5343` active cols/win | `111.5343` | `Guard4 (W_B = 4, L <= 248)` |
| Non-Convex $L_{0.65}$ Interior Gap Fraction | `0.001657` measured | `0.001757` (`gap_ch`) | `gap_pad = log2(1/(1 - gap_ch)) = +0.002538` b |
| 16B-Granular Adaptive Prebias (2 passes) | `0.192400` Ir/active cell | Part of `0.4950` | `_slab_prebias_adaptive_sl` (`0` spills, `k16 = 1..16`) |
| 16B-Granular Masked $S_g$ Scan (fast path) | `0.277542` Ir/active cell | Part of `0.4950` | `_slab_scan_g_active_16` (`248/248` bit-exact w/ `0xFF` poison) |
| Slow-Path Extraction + Direct $S_f$ Probe | `0.007743` (`189 * 2^{-14.5751}`) | Part of `0.4950` | Total live `0.477685 <= 0.4950 Ir/active cell` (`+0.0173` pad!) |

---

## 4. Explicit Kill-Gates for Pruned Stage 1a Alternatives

Every alternative below was benchmarked under identical Critic-proof rules (`2.00x` uncached Bernstein, exact 2D Panjer, `W_pointwise = 2^{124}` u) and **explicitly pruned**:

1. **`Guard2` (`W_B = 2, L <= 252`) — KILLED (`+0.0291` b loss):** Although `W_B = 2` reduces wasted guard hits (`-0.063` Ir/cell), doubling bucket count `n_bkt` inflates `replay_a_tables` (`104.79 -> 105.16` u) and peak RAM (`2^{117.26} -> 2^{117.97}` B), causing a net loss of **`+0.0291` bits**. `Guard4 (W_B = 4, L <= 248)` strictly dominates.
2. **$q \in \{0.55, 0.60\}$ Ball Exponent — KILLED (`+0.008` to `+0.022` b loss):** Tightening $q$ from `0.65` to `0.60` gains only `-0.0029` b on $V_{\text{kept}}$ while doubling interior chord gaps (`0.00167 -> 0.00336`) and shrinking $L_{\text{win,eff}}$ (`111.8 -> 106.2`), losing **`+0.008` bits** net (`+0.022` b at $q = 0.55$). $q = 0.65$ is optimal.
3. **$K = 3$ Equal-$V_a$ Shells — KILLED (`+0.0237` b loss):** Multiplying replay table construction by `1.5x` (`+0.585` b on `replay_a_tables`) forces narrower bands that sacrifice 2D Panjer recall and shrink $L_{\text{win,eff}}$, leaving $K = 3$'s certified optimum `+0.0237` bits worse than $K = 2$.
4. **$q_f = 3.4$ Quantization on $f$-Side — KILLED:** Deterministic wrap-safety at $N_{f,\max} = 630$ forces $T_f \ge \lfloor 630/3.4 \rfloor - 127 = 58$ (`197.2` bits), severely clipping $f$-side smooth recall. Keeping $q_f = 3.8$ ($T_f \ge 38$) and $q_g = 3.4$ (per-column $N_g(b)$ floor) is strictly optimal.
