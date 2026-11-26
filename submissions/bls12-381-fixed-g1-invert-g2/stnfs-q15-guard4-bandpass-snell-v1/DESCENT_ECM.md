# Stage 4 12D Subfield $M_g$ Hilbert-90 Descent and Torsion-16 Edwards ECM (`DESCENT_ECM.md`)

This companion document specifies the complete Stage 4 individual-logarithm recovery pipeline for both the challenge target $z \in G_T$ and the fixed-argument base $g_0 = e(P, Q_0) \in G_T$. Every entry in the 12D subfield norm proof, the $q = 2^{70}..2^{275}$ ECM calibration table, the 5-layer descent table, and the explicit `ecm_split_expected` derivation is reproduced to $< 10^{-12}$ relative agreement by `experiment.py` and grounded in the compiled scalar ARMv8.0-compatible kernels `_m256_mont_mul` (`212 Ir`) and `_ecm_mont_xdbladd_11m` (`2,704 Ir/bit`) in `KERNELS_LA_ECM.md`.

---

## 1. Exact 12D Subfield $M_g = \mathbb{Q}(\zeta_5)[Y]/(Y^3 - x)$ Hilbert-90 Coset Lift

### 1.1 Exact Field Isomorphism $\mathcal{O}_{M_g} / (p) \cong \mathbb{F}_{p^{12}}$
Consider the degree-12 totally complex subfield of $L_g = \mathbb{Q}(\zeta_{15})[Y]/(Y^3 - x)$:
$$M_g := \mathbb{Q}(\zeta_5)[Y] / (Y^3 - x), \qquad [M_g : \mathbb{Q}] = \varphi(5)\cdot 3 = 12, \qquad [L_g : M_g] = 2.$$
**Theorem (Exact Residue-Field Isomorphism, verified live in `experiment.py`):**
1. Since $p \equiv 2 \pmod 5$ has multiplicative order $\operatorname{ord}_5(p) = 4$, $\Phi_5(t) = t^4 + t^3 + t^2 + t + 1$ is irreducible over $\mathbb{F}_p$, so $(p)$ is **inert** in $\mathbb{Z}[\zeta_5]$: $\mathbb{Z}[\zeta_5] / (p) \cong \mathbb{F}_{p^4}$.
2. Since $x^{(p-1)/3} \not\equiv 1 \pmod p$ and $(p^4 - 1)/(p - 1) \equiv 4 \equiv 1 \pmod 3$, $g(Y) = Y^3 - x$ is **irreducible** over $\mathbb{F}_{p^4}$, yielding
$$\mathcal{O}_{M_g} / (p) \cong \mathbb{F}_{p^4}[Y] / (Y^3 - x) \cong \mathbb{F}_{p^{12}}.$$
3. Because $M_g \hookrightarrow L_g$ is a quadratic extension ($L_g = M_g(\omega)$), every degree-1 prime ideal $\mathfrak{q}_{M}$ of $M_g$ above a rational prime $q \equiv 1 \pmod 5$ behaves according to $q \bmod 3$:
   - **Split primes ($q \equiv 1 \pmod{15}$, ~50%):** $\mathfrak{q}_M \mathcal{O}_{L_g} = \mathfrak{Q}_1 \mathfrak{Q}_2$ splits into two degree-1 prime ideals of $L_g$ with **exact norm $q$**. Those with $q \le B_g$ are immediately in the factor base; those with $q \in (B_g, B_{\text{boot}}]$ enter Layers 1–5 directly.
   - **Inert primes ($q \equiv 11 \pmod{15}$, i.e., $q \equiv 2 \pmod 3$, ~50%):** $\mathfrak{Q} = \mathfrak{q}_M \mathcal{O}_{L_g}$ remains inert of residue degree $2$ over $\mathbb{Q}$ with **$L_g$-norm $\operatorname{Norm}_{L_g/\mathbb{Q}}(\mathfrak{Q}) = q^2$**! We handle EVERY inert prime explicitly:
     - $q \in (2^{137.5}, 2^{275}]$ has $L_g$-norm $q^2 \in (2^{275}, 2^{550}]$: handled by an explicit **Layer-0 special-$\mathfrak{Q}$ descent** ($2^{550} \to 2^{275}$, $\le 5$ nodes, each costing `2^{61.825}` u $\implies$ Layer-0 total `2^{64.147}` u, decomposing $\mathfrak{Q}$ into split degree-1 ideals $\le 2^{275}$).
     - $q \in (2^{34.98}, 2^{137.5}]$ has $L_g$-norm $q^2 \in (B_g = 2^{69.96}, 2^{275}]$ (notably including $q \in (2^{34.98}, 2^{69.96}]$ whose squared norm $q^2 \in (2^{69.96}, 2^{139.92}]$ lies above $B_g$!): handled by adding **`+20` extra descent nodes** conservatively charged at the **most expensive Layer-5 cost** (`180 + 20 = 200` Layer-5 nodes, strictly upper-bounding any split across Layers 3–5!).
     - $q \le 2^{34.98}$ has $L_g$-norm $q^2 \le 2^{69.96} = B_g$, already inside the factor-base norm ceiling (or trivially descended in $< 2^{35}$ u if degree-2).

### 1.2 Unconditional Deterministic 12D LLL + AM–GM Bootstrapping Norm Proof ($N_{\text{boot}} \le 2{,}572$ bits)
Given a randomized target $h_k = h \cdot g_{\text{pub}}^{e_k} \in G_T \subset \mathbb{F}_{p^{12}}^\times$, applying Hilbert's Theorem 90 across $\mathbb{F}_{p^{12}} / \mathbb{F}_{p^6}$ expresses $h_k$ up to $\mathbb{F}_{p^6}^\times$ (which vanishes upon exponentiation by the cyclotomic cofactor $p^6 - 1$, since $\gcd(p^6 - 1, r) = 1$) via **6 linear equations over $\mathbb{F}_p$** on the 12 integer coordinates of $A(Y) = \sum_{j=0}^2 \sum_{i=0}^3 c_{i,j}\,\zeta_5^i Y^j \in \mathcal{O}_{M_g}$. These 6 linear congruences define an index-$p^6$ sublattice $\Lambda_{12} \subset \mathcal{O}_{M_g}$ of rank $12$:

1. **Exact Weighted Basis Determinant:** Scaling $\{1, Y, Y^2\}$ by $\{1, |x|^{1/3}, |x|^{2/3}\}$ multiplies $p^6$ by $(1 \cdot |x|^{1/3} \cdot |x|^{2/3})^4 = |x|^4$:
$$\log_2\!\bigl(p^6 |x|^4\bigr) = 6\log_2 p + 4\log_2 |x| = \mathbf{2{,}539.061119\text{ bits}}.$$
2. **Exact Canonical $\mathbb{R}^{12}$ Discriminant & Covolume:** Under the Minkowski canonical embedding $\sigma: M_g \hookrightarrow \mathbb{R}^{12}$ with coordinates $(\sqrt{2}\operatorname{Re}\sigma_k, \sqrt{2}\operatorname{Im}\sigma_k)_{k=1}^6$ (so $\|A\|_2^2 = 2\sum_{k=1}^6 |\sigma_k(A)|^2$ and $\operatorname{covol}_{\mathbb{R}^{12}}(\mathcal{O}_{M_g}) = \sqrt{|\Delta_{M_g}|}$ identically),
$$|\Delta_{M_g}| = |\Delta_{\mathbb{Q}(\zeta_5)}|^3 \cdot \bigl|N_{\mathbb{Q}(\zeta_5)/\mathbb{Q}}(\operatorname{disc}(Y^3 - x))\bigr| = 5^9 \cdot 3^{12} \cdot |x|^8,$$
$$\operatorname{covol}_{\mathbb{R}^{12}}(\Lambda_{12}) = \sqrt{|\Delta_{M_g}|}\,p^6 = 5^{9/2}\cdot 3^6 \cdot |x|^4 \cdot p^6 = \mathbf{2^{2559.019571}}.$$
3. **Unconditional Deterministic LLL + AM–GM Norm Upper Bound:** Even standard deterministic LLL ($\delta = 3/4$) in dimension $d = 12$ (without requiring BKZ-20!) guarantees a nonzero vector $A \in \Lambda_{12}$ with $\|A\|_2 \le 2^{(12-1)/4}\,\operatorname{covol}_{\mathbb{R}^{12}}(\Lambda_{12})^{1/12} = 2^{11/4}\,\operatorname{covol}_{\mathbb{R}^{12}}(\Lambda_{12})^{1/12}$. By AM–GM on the 6 conjugate pairs $|\sigma_k(A)|^2$:
$$|N_{M_g/\mathbb{Q}}(A)| = \prod_{k=1}^6 |\sigma_k(A)|^2 \le \left(\frac{\|A\|_2^2}{12}\right)^6 \le \frac{2^{33}}{12^6}\sqrt{|\Delta_{M_g}|}\,p^6 = \mathbf{2^{2570.509796} < 2^{2572.0} = 2^{N_{\text{boot}}}}!$$
*(Using BKZ-20 in 12D yields even shorter norms near $2^{2543}$ bits; our ledger conservatively charges the unconditional worst-case ceiling $N_{\text{boot}} = 2{,}572.0$ bits!)*

---

## 2. Strictly-Stage-2 Torsion-16 Edwards ECM Specification & Exact Calibration Table ($N_{\text{boot}} = 2{,}572$ bits)

We use twisted Edwards curves with rational torsion subgroup $\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/8\mathbb{Z}$ ($16 \mid \#E(\mathbb{F}_q)$, parametrized in birational Montgomery form with $x$-only differential doubling/addition `_ecm_mont_xdbladd_11m`):
- **Curve generation:** `40` modular multiplications per curve (`ir_gen = 40 * m_ir`), with exact compiled 4-limb CIOS `m_ir = 212 Ir` ($L = 4$), `13.25 L^2 + 80` ($4 < L \le 16$), and $\min(13.25 L^2 + 80,\; 2.5 \cdot \text{ntt\_mul\_ir}(L))$ ($L > 16$, using compiled Shoup61 NTT from `KERNELS_LA_ECM.md`).
- **Stage 1 cost:** $(B_1 / \ln 2)\cdot(11\,m_{\text{ir}} + 93 L)$ Ir/curve (`2,704 Ir/bit` at $L = 4$, exact match to `_ecm_mont_xdbladd_11m`).
- **Stage 2 cost:** $2.5 \cdot (\pi(B_2) - \pi(B_1)) \cdot (m_{\text{ir}} + 40 L)$ Ir/curve, `1.15x` loop/control/GCD overhead.
- **Torsion-16 Hasse-interval smoothness credit:** Conservatively `3.5` bits ($q_{\text{eff}} = q_{\text{bits}} - 3.5$), using the exact Dickman–Buchstab Stage-1 + Stage-2 probability:
$$p_{\text{curve}}(B_1, B_2; q) = \rho\!\left(\frac{q_{\text{eff}}}{\log_2 B_1}\right) + \int_{\log_2 B_1}^{\log_2 B_2} \rho\!\left(\frac{q_{\text{eff}} - t}{\log_2 B_1}\right)\frac{dt}{t}.$$

### 2.1 Self-Consistent Unconstrained ECM Calibration Table ($N_{\text{boot}} = 2{,}572$ bits, $q = 2^{70}..2^{275}$)

Optimizing over the unconstrained grid $\log_2 B_1 \in [8.0, 62.0]$ (step `0.25`) and $B_2 / B_1 \in \{10, 25, 60, 120, 250\}$ gives the exact table below. Every single row satisfies $\mathbb{E}[\text{cost}] = (\text{1-curve cost}) \cdot (1/p_{\text{curve}})$ to $< 10^{-12}$ relative error (verified live in `experiment.py`):

| Prime $q$ | $\log_2 B_1$ | $B_1$ | $B_2$ | $1/p_{\text{curve}}$ | 1-Curve Cost ($\log_2 u$) | **$\mathbb{E}[\text{cost}]$ ($\log_2 u$)** |
|---|---:|---:|---:|---:|---:|---:|
| $2^{70}$ | `14.25` | `1.948e+04` | `1.169e+06` | `1.409e+02` | `19.387` | **`26.525`** |
| $2^{85}$ | `16.25` | `7.794e+04` | `4.676e+06` | `3.779e+02` | `21.332` | **`29.894`** |
| $2^{100}$ | `18.25` | `3.117e+05` | `1.870e+07` | `8.535e+02` | `23.286` | **`33.023`** |
| $2^{120}$ | `20.50` | `1.483e+06` | `8.897e+07` | `2.734e+03` | `25.490` | **`36.907`** |
| $2^{140}$ | `22.75` | `7.054e+06` | `4.232e+08` | `7.274e+03` | `27.700` | **`40.529`** |
| $2^{160}$ | `24.75` | `2.822e+07` | `1.693e+09` | `1.978e+04` | `29.669` | **`43.941`** |
| $2^{185}$ | `27.25` | `1.596e+08` | `9.577e+09` | `5.827e+04` | `32.134` | **`47.965`** |
| $2^{215}$ | `29.75` | `9.029e+08` | `5.417e+10` | `2.446e+05` | `34.604` | **`52.504`** |
| $2^{245}$ | `32.50` | `6.074e+09` | `3.644e+11` | `7.209e+05` | `37.324` | **`56.784`** |
| **$2^{275}$** | **`34.75`** | **`2.889e+10`** | **`1.734e+12`** | **`2.572e+06`** | **`39.553`** | **`60.848`** |

---

## 3. Complete 6-Layer (Layer-0..5) Special-$\mathfrak{q}$ Descent Tree & Explicit `ecm_split_expected` Derivation

After stripping small primes $\le 2^{26}$ via $T_{26}$ ($< 2^{45}$ u), we test $N_{\text{boot}} = 2{,}572$ bits for $B_{\text{boot}} = 2^{275}$ smoothness with **zero early-abort** ($\rho(2572/275) = 2^{-31.68}$, survival factor `0.918`, charging `2.5x` full $2^{275}$-ECM discovery work per candidate), then descend every prime ideal above $B_g = 2^{69.96}$ — explicitly including **Layer-0** for inert $q \equiv 11 \pmod{15}$ with $q^2 \le 2^{550}$ AND **`+20` extra inert-prime nodes** (`180 + 20 = 200` conservatively charged at Layer-5 cost) for inert $q \in (2^{34.98}, 2^{69.96}]$ with $q^2 > B_g$:

| Stage 4 Component ($B_f = 2^{69.60}, B_g = 2^{69.96}$, $N_{\text{boot}} = 2572$ b, $B_{\text{boot}} = 2^{275}$) | Max Nodes (1 log) | Optimal `(deg, s)` | Enum $T$ ($\log_2$) | Per-Node ($\log_2 u$) | Layer Cost ($\log_2 u$) |
|---|---:|---|---:|---:|---:|
| **Bootstrapping (1 log, 12D $M_g$, ZERO early-abort)** | `1` | 12D Hilbert-90 LLL | — | `93.978` | **`93.978`** |
| **Layer 0** (inert $q \equiv 11 \bmod 15$: $2^{550.0} \to 2^{275.0}$) | `<= 5` | `(2, 0.0b)` | `54.50` | `61.825` | `64.147` |
| Layer 1 ($2^{275.0} \to 2^{155.0}$) | `<= 5` | `(2, 0.0b)` | `40.25` | `43.36` | `45.68` |
| Layer 2 ($2^{155.0} \to 2^{130.0}$) | `<= 12` | `(1, 0.0b)` | `35.84` | `38.33` | `41.91` |
| Layer 3 ($2^{130.0} \to 2^{106.0}$) | `<= 28` | `(1, 0.0b)` | `46.23` | `42.91` | `47.72` |
| Layer 4 ($2^{106.0} \to 2^{86.0}$) | `<= 65` | `(1, 0.0b)` | `61.47` | `58.15` | `64.17` |
| Layer 5 ($2^{86.0} \to 2^{69.96}$, incl. `+20` extra inert nodes!) | `<= 200` | `(1, 0.0b)` | `86.03` | `82.71` | **`90.35`** |
| **1-Log Descent Tree Subtotal ($W_{\text{tree},1}$, incl. Layer 0 & `+20` inert)** | — | — | — | — | **`90.350`** |
| **BOTH Logs ($z$ and $g_0$), Live Zero-Early-Abort (no Snell)** | — | — | — | — | **`95.090`** |
| **BOTH Logs with Snell Missing-Log Penalty $(1.04094 / 0.82)^2 = 1.61147\times$** | — | — | — | — | **`95.778`** |
| **CHARGED `both_individual_logs_exp` (`+0.72` b extra headroom)** | — | — | — | — | **`96.5001`** |

### 3.1 Explicit Derivation of `ecm_split_expected` (`2^{94.0707}` u)
To cover all residual multi-prime cofactor splitting inside the same Bernstein leaf bucket (Stage 1b) and secondary composite splits inside Stage 4, we derive and recompute live in `experiment.py`:
$$W_{\text{ecm\_split,live}} = 2 \cdot \bigl(W_{\text{tree},1} + 0.25\,W_{\text{boot},1}\bigr) \cdot \left(\frac{1.04094}{0.82}\right)^2 = 2 \cdot \bigl(2^{90.350} + 0.25 \cdot 2^{93.978}\bigr) \cdot 1.61147 = \mathbf{2^{94.070586}\text{ u}},$$
charged in `ledger_data.json` as **`ecm_split_expected = 2^{94.0707}` u** (exceeding Stage 1b same-bucket splitting work by $> 10\times$!).

---

## 4. H5 Torsion-16 Edwards ECM Live Witnesses

To verify the birational Torsion-16 ($\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/8\mathbb{Z}$) Edwards/Montgomery parametrization ($x = (m^2 - 2m + 2)/(2 - m^2)$, $d = (2x^2 - 1)/x^4 = ((m^2 - 4m + 2)(2 - m^2)/(m^2 - 2m + 2)^2)^2$, $A_{24} = 1/(1 - d)$), the $x$-only differential doubling/addition ladder (`_ecm_mont_xdbladd_11m`), and strictly-Stage-2 BSGS factor recovery directly in pure Python, `experiment.py` executes three live witnesses on every audit run:
1. **Torsion-16 Group-Order Witness (`21/21` primes in `[101, 200]`):** Verifies $16 \mid \#E(\mathbb{F}_p)$ for 100% (`21/21`) of primes $p \in [101, 200]$ ($p \le 200$) on the rational family $m = 3$.
2. **Witness 1 (`50023`, Strictly-Stage-2 Recovery on $N_1 = 50{,}023 \times 50{,}033 = 2{,}502{,}800{,}759$):** Running the Torsion-16 curve $m = 3$, base point $X_0 = 3$, $B_1 = 30$, $B_2 = 350$, BSGS step $D = 30$, Stage 1 exits with $\gcd(Z_{B_1}, N_1) = 1$, and Stage 2 giant step $270$ vs baby step $1$ recovers the exact prime factor **`50023`** strictly in Stage 2.
3. **Witness 2 (`35 = 5 \times 7`, Stage-2 & Zero-Leaf Isolation Witness):** Running the $x$-only Montgomery ladder at $A_{24} = 2$ ($A = 6$, good reduction mod $5$ and mod $7$), $X_0 = 4$, $B_1 = 2$, Stage 1 exits with $\gcd(Z_2, 35) = 1$, and Stage 2 BSGS step $D = 4$ isolates the nontrivial split **`5`** (and cofactor **`7`**) from $N_2 = 35$ cleanly.
