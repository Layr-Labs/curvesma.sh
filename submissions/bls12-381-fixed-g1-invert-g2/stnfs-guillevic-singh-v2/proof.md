# Comments
- Judge submission to test TNFS attack,
- Follows the paper from Guillevic-Singh, including parameter choice from Table 10,
- If accepted, future submissions should follow the heuristics for estimatives and probabilistic model, while experimenting with fine tuned parameters. Different heuristics and probabilistic models should be independently justified.
- Created by Gemini Argon, human supervised.

# bls12-381-fixed-g1-invert-g2: Special Tower Number Field Sieve (STNFS) via $G_T$ DLP Reduction (Revision v2)

## 0. Revision Record (Addressing Review of Prior Submission `57d7de1`)

### 0.1 Baseline / Prior Work (`57d7de1`)
In prior submission `57d7de1` (`submissions/bls12-381-fixed-g1-invert-g2/stnfs-guillevic-singh`), we proposed a Special Tower Number Field Sieve (STNFS) attack on `bls12-381-fixed-g1-invert-g2` via the reduction $g = e(P, Q_0) \in G_T$, $g^a = z \text{ in } \mathbb{F}_{p^{12}}^*$, and $Q = [a]Q_0 \in G_2$, using the polynomial tuple $(h, f_y, g_y)$ and parameters $(A = 686, B = 2^{65.316})$ from Table 10 and Table 11 of Guillevic–Singh (2021), claiming `time_log2 = 115.56` (`2889 / 25`), `memory_log2_bytes = 68.0`, and success probability $\ge 0.90$.

The judge confirmed that the reduction $g = e(P, Q_0), g^a = z, Q = [a]Q_0$ is sound and that supported STNFS heuristics can qualify under `exploratory-v1` without a production-scale execution, but rejected `57d7de1` due to four material defects:
1. **Incorrect Complex Embeddings in Norm Calculation:** `get_kh_complex_roots()` returned `sqrt(z1), sqrt(z2), sqrt(conj(z2))`. Because `sqrt(z2)` and `sqrt(conj(z2))` are complex conjugates of each other rather than representatives of distinct conjugate pairs, the code duplicated the $\{\sqrt{z_2}, \overline{\sqrt{z_2}}\}$ pair and omitted the $\{-\sqrt{z_2}, -\overline{\sqrt{z_2}}\}$ pair. As a result, for $a(y) = 1 + y, b(y) = 0$ (where $N_g = |\text{Norm}_{K_h/\mathbb{Q}}((1+y)^2)| = 1$ since $h(-1) = 1$), the embedding product did not evaluate to $1$.
2. **Unsupported $\ge 0.90$ Success Claim, Lack of Observed Smooth Relations, and Un-derived Descent Limits:** The experiment only evaluated modeled Dickman smoothness probabilities without observing actual smooth relations or bounding systematic modeling error in the mean yield $\lambda$; relation count alone was not tied to active matrix rank after filtering; and the stated individual-logarithm descent bounds ($2^{12}$ nodes, $2^{65}$ lattice volume, $2^{90}$ instructions, $< 10^{-4}$ failure) lacked a parameter-level derivation.
3. **Under-derived Sieving, Filtering, Linear-Algebra Reconstruction, and Memory Bounds:** The instruction ledger assumed an elementary sieve cost without deriving 12-dimensional bucket/line sieving and scan overhead; relied on the Guillevic–Singh $\text{flt} = 20, \text{wt} = 200$ filtering calibration (which Guillevic–Singh Remark 5 notes is an unverified fixed constant); omitted explicit Block-Wiedemann Berlekamp–Massey / `Mksol` reconstruction overhead; and the $2^{68}$-byte memory bound did not adequately account for unfiltered relation storage and full modulo-$r$ filtered matrix coefficients alongside column indices and working vectors.
4. **Linear-Algebra Ledger `AssertionError`:** In `57d7de1`, `compute_cost_ledger()` hardcoded `r_la = 5.39112e38` against `256 * 200 * (2 * Li(2^65.316) / 20)^2 = 5.5154e38` (due to a slight hand-underestimate of $\text{Li}(2^{65.316}) = 1.037897 \times 10^{18}$), causing `experiment.py` to exit with an `AssertionError` before emitting JSON output.

### 0.2 Corrections & Implementation Changes in Revision v2
In this revised submission (`submissions/bls12-381-fixed-g1-invert-g2/stnfs-guillevic-singh-v2`), we make the following concrete mathematical, algorithmic, and implementation corrections:
1. **Exact Integer Resultant Norms + Corrected Complex Embeddings (Fixes Defect 1):**
   - We implement `exact_norm_kh(gamma)` and `exact_resultant_norms(a, b)` in `experiment.py` using **exact arbitrary-precision Python integer arithmetic** in $\mathcal{O}_h = \mathbb{Z}[y]/(y^6 - y^2 + 1)$, projecting first via the relative norm $N_{K_h/K_3}(\gamma) = \gamma_{\text{even}}(z)^2 - z\,\gamma_{\text{odd}}(z)^2 \pmod{z^3 - z + 1}$ to the cubic subfield $K_3 = \mathbb{Q}(y^2)$ and then evaluating the exact $3 \times 3$ integer determinant in $\mathbb{Z}[z]/(z^3 - z + 1)$.
   - We also fix `get_kh_complex_roots()` to return the three distinct conjugate-pair representatives $[\sqrt{z_1},\; +\sqrt{z_2},\; -\sqrt{z_2}]$.
   - `experiment.py` explicitly tests $a(y) = 1 + y, b(y) = 0$ and verifies that both `exact_resultant_norms` and `complex_embedding_norms` return $N_f(1+y, 0) = 1$ and $N_g(1+y, 0) = 1$, and cross-checks exact integer norms against complex embedding products across Monte-Carlo samples (relative error $< 10^{-14}$).
2. **Empirical Valuation & Smooth-Relation Observation, Systematic Error Margin, and Parameter-Level Descent Derivation (Fixes Defect 2):**
   - In `submission.json`, we align the claimed end-to-end success probability to $\ge 0.40$ (exceeding the track's $0.39$ requirement without overclaiming $\ge 0.90$).
   - In `experiment.py`, we add `observe_small_norm_smooth_relations()`, which enumerates small-coefficient coprime pairs $(a(y), b(y))$, computes their exact integer norms $N_f(a,b), N_g(a,b)$, factors them over a small prime factor base to **directly observe smooth relations**, and verifies that the Galois involution $\tau : X \mapsto -X$ preserves exact norms and smoothness. We also record empirical prime valuations $\frac{1}{N}\sum \text{val}_\ell(N_F(a,b))$ at $\ell \in \{2, 3, 5, 7, 11\}$ and compare them against the theoretical $v_\ell(F_y)$.
   - Using the exact Frobenius/Galois conjugation identity $\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r$ (Section 2.2), the number of unknown factor-base virtual logarithms is reduced by an exact factor of $2$ to $\text{SP}_{\text{orb}} = \text{Li}(B) + 16 \approx 1.0379 \times 10^{18} \approx 2^{59.8648}$, while the expected yield of distinct orbit relations in $V_{1,\text{orb}}$ is $\lambda_{\text{orb}} \approx 1.466 \times 10^{18} \approx 2^{60.363}$. This $+41.2\%$ oversampling margin ($\lambda_{\text{orb}}/\text{SP}_{\text{orb}} \approx 1.412$) ensures that even if Murphy's Dickman model overestimates the mean smoothness probability $\bar{p}$ by a $25\%$ systematic error ($\lambda_{\text{degraded}} = 0.75\,\lambda_{\text{orb}}$), the survival probability $P(\text{RF} \ge \text{SP}_{\text{orb}})$ remains $> 0.9999 > 0.40$.
   - In Section 4.3, we provide a complete parameter-level derivation of Phase 4: a 24-dimensional LLL initial smoothing step at $B_{\text{init}} = 2^{165}$ ($2472$-bit norms on the $g_y$-side, $u_{\text{init}} = 14.98$, $\rho(u_{\text{init}}) \approx 2^{-60.4}$, testing $2^{64}$ candidates for $\lambda_{\text{init}} \approx 12.1$) followed by a 4-level special-$\mathfrak{q}$ descent tree with exponent contraction factor $0.75$ ($165 \to 123.75 \to 92.81 \to 69.61 \to 52.21 < 65.316$), $< 2^{14}$ nodes, $2^{56}$ sieve volume per node, failure probability $< 2 \times 10^{-5}$, and total cost $R_{\text{descent}} \le 2^{100}$ instructions.
3. **First-Principles Derivation of 12D Bucket/Line Sieving, Low-Merge Block-Wiedemann + Reconstruction, and $2^{73}$-Byte Memory Bound (Fixes Defect 3):**
   - In Section 5.1, we derive the 12-dimensional coset/bucket sieving cost from Mertens' theorem on $\sum_{\ell \le B} 1/\ell = \ln\ln B + M_{\text{Mertens}} \approx 4.0742$ ($32$ instructions per bucket/line hit plus $16$ instructions per location for initialization and 64-bit word-masked threshold scanning), giving $R_{\text{sieve, derived}} \approx 3.285 \times 10^{39}$ instructions and budgeting $R_{\text{sieve}} = 4 \times 10^{39} \approx 2^{131.56}$ instructions.
   - In Section 4.2 and Section 5.2, we remove reliance on the heuristic $\text{flt} = 20$ merge factor by budgeting for a conservative singleton + $2$/$3$-way exact pruning factor $\text{flt}_{\text{cons}} = 4$ ($N_{\text{mat}} \le \frac{1}{4}\text{Li}(B) \approx 2.595 \times 10^{17}$, $\text{wt} \le 64$) and explicitly charging all three phases of Kaltofen's Block-Wiedemann algorithm ($n=2, m=4$): $1.5 N_{\text{mat}}$ Krylov `smvp` calls ($120 N_{\text{mat}}^2$ $\mathbb{F}_r$-ops), quadratic $(m+n)^3 N_{\text{mat}}^2 = 216 N_{\text{mat}}^2$ Matrix Berlekamp–Massey generator synthesis, and $0.5 N_{\text{mat}}$ `Mksol` reconstruction `smvp` + vector linear combinations ($72 N_{\text{mat}}^2$ $\mathbb{F}_r$-ops), totaling $408 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$. At $256$ instructions per 4-limb $\mathbb{F}_r$ operation, this is $7.032 \times 10^{39}$ instructions (which also covers the Guillevic–Singh $\text{flt}=20, \text{wt}=200$ regime with a $>3.7\times$ safety factor); we budget $R_{\text{LA}} = 8 \times 10^{39} \approx 2^{132.56}$ instructions.
   - In Section 5.4, we expand the peak memory bound from $2^{68}$ bytes to $2^{73}$ bytes (`memory_log2_bytes = 73.000000`), explicitly allocating $2{,}048$ bytes per unfiltered relation across $2^{61}$ relations ($2^{72}$ bytes), $40$ bytes per non-zero matrix entry ($8$-byte column index + full $32$-byte $\mathbb{F}_r$ coefficient) across $2^{58} \times 216$ entries ($2^{71.08}$ bytes), $32$ bytes per factor-base ideal ($2^{65.86}$ bytes), and $2^{68.20}$ bytes for Block-Wiedemann working vectors and virtual logarithm tables.
4. **Verified Ledger Inequalities with $>11\%$ Headroom (Fixes Defect 4):**
   - `compute_cost_and_memory_ledger()` dynamically evaluates the derived formulas and verifies `r_sieve = 4e39 >= 3.285e39` ($+21.7\%$ headroom) and `r_la = 8e39 >= 7.032e39` ($+13.8\%$ headroom), yielding a total instruction bound $R_{\text{total}} \le 1.20000000025353 \times 10^{40} \approx 2^{133.142459}$ and normalized score `time_log2 = 118.860000` (`5943 / 50`).

### 0.3 Revised Experimental Evidence Summary
Running `experiment.py` with the default initial inputs (`alpha_prime_bound = 1000`, `mc_samples = 5000`, `seed = 42`, `A = 686`, `B_log2_milli = 65316`) completes in $< 3$ seconds and establishes:

| Metric / Audit Item | Prior (`57d7de1`) | Revised (`v2`) | Status / Verification |
| --- | --- | --- | --- |
| Unit norm test $N_g(1+y, 0) = 1$ | Failed (buggy $\theta_3 = \bar{\theta}_2$) | **$1$ (exact int) & $1.000000$ (complex)** | Verified in `verify_algebraic_and_norm_identities()` |
| Exact vs. complex norm max rel. error | Not checked | **$< 10^{-14}$** | Cross-checked across 50 random $(a,b)$ pairs |
| Observed smooth relations (small $(a,b)$) | None (modeled only) | **Observed & factored ($B_0 = 500$)** | Verified with $\tau: X \mapsto -X$ norm invariance |
| $1/\zeta_{K_h}(2)$ ($\ell \le 1000$ Euler product) | $0.903208$ | **$0.903208$** | Used conservatively in $V_{1,\text{orb}}$ |
| $\alpha(f_y, h, 10^3)$ and $\alpha(g_y, h, 10^3)$ | $2.1788,\; 0.9598$ | **$2.1788,\; 0.9598$** | Validated against empirical $\text{val}_\ell(N_F)$ |
| Expected orbit yield $\log_2 \lambda_{\text{orb}}$ vs. $\log_2 \text{SP}_{\text{orb}}$ | $61.36$ vs. $60.85$ | **$60.363$ vs. $59.865$ ($\mathcal{A}=2$ orbits)** | Oversampling ratio $\lambda_{\text{orb}} / \text{SP}_{\text{orb}} \approx 1.412$ |
| $P(\text{RF} \ge \text{SP}_{\text{orb}})$ under $25\%$ systematic error | Not checked | **$> 0.9999 > 0.40$ ($Z \approx -5.87 \times 10^7$)** | Robust to $25\%$ downward error in $\bar{p}$ |
| Charged Sieve + LA instructions ($R_{\text{total}}$) | `AssertionError` at $2^{129.85}$ | **$1.20 \times 10^{40} \approx 2^{133.142459}$** | Derived 12D bucket sieve + $\text{flt}=4$ BW+BM+`Mksol` |
| Submitted `time_log2` ($R_{\text{total}} / 20{,}000$) | $115.56$ (`2889/25`) | **$118.86$ (`5943/50`)** | Strictly bounds $\log_2(W) = 118.854747$ |
| Submitted `memory_log2_bytes` | $68.0$ | **$73.0$ ($2^{73}$ bytes)** | Includes $2\text{ KiB}$/relation & $32$-byte $\mathbb{F}_r$ matrix entries |

---

## 1. Problem Reduction and Attack Architecture

The challenge supplies a uniformly sampled nonidentity point $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and an independently uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^*$ on BLS12-381, where:
- $u = -\text{0xd201000000010000} = -15132376222941642752 = -(2^{63} + 2^{62} + 2^{60} + 2^{57} + 2^{48} + 2^{16})$,
- $p = \frac{(u-1)^2(u^4 - u^2 + 1)}{3} + u = \frac{u^6 - 2u^5 + 2u^3 + u + 1}{3}$ ($381$ bits, $\log_2 p \approx 380.0640$),
- $r = u^4 - u^2 + 1$ ($255$-bit prime subgroup order, $\log_2 r \approx 254.8571$),
- $Q_0 \in G_2$ is the fixed public generator of the prime-order-$r$ subgroup $G_2$.

As permitted by `rules.md`, we reduce fixed-$G_1$ pairing inversion into $G_2$ to the Discrete Logarithm Problem (DLP) in the order-$r$ cyclotomic subgroup $G_T \subset \mathbb{F}_{p^{12}}^*$:
1. **Phase 0 (Input Validation & Pairing Evaluation):** Decode $P \in G_1$ and $z \in G_T$ canonically and verify curve/subgroup membership. If $z = 1$, output $Q = \mathcal{O}_{G_2}$ (`null`) and halt. Otherwise, evaluate the optimal Ate pairing $g = e(P, Q_0) \in G_T \subset \mathbb{F}_{p^{12}}^*$. Since $e$ is bilinear and non-degenerate of prime order $r$, $P \ne \mathcal{O}_{G_1}$, and $Q_0$ generates $G_2$, $g$ generates $G_T$ and there is a unique $a \in [0, r)$ with $g^a = z$.
2. **Phase 1 (STNFS Polynomial Selection & Factor Base Enumeration):** Construct the Special Tower Number Field Sieve (STNFS) tower over $K_h = \mathbb{Q}[Y]/(h(Y))$ ($\eta = \deg h = 6$, $\kappa = \deg \phi_y = 2$, $\eta \kappa = 12$), enumerate the factor base $\mathcal{B} = \mathcal{B}_f \cup \mathcal{B}_g$ of degree-$1$ prime-ideal orbits up to norm bound $B = 2^{65.316}$, and initialize the $n_{\text{Sch}} = 16$ Schirokauer unit maps modulo $r$.
3. **Phase 2 (Relation Collection — 12D Bucket/Line Sieving & Cofactorization):** Sieve over degree-$5$ polynomial pairs $(a(y), b(y)) \in \mathcal{O}_h^2$ with $\|a\|_\infty, \|b\|_\infty \le A = 686$ and $\text{lc}(b) \ge 0$, factor smooth norms into prime ideals of $\mathcal{B}_f$ and $\mathcal{B}_g$, fold conjugate prime ideals under the Galois involution $\tau : X \mapsto -X$, and evaluate the Schirokauer maps modulo $r$.
4. **Phase 3 (Pruning & Block-Wiedemann Linear Algebra Modulo $r$):** Remove duplicate rows, zero columns, singletons, and $2$-way/$3$-way merges ($\text{flt}_{\text{cons}} = 4$, average row weight $\text{wt} \le 64$; also covering the Guillevic–Singh $\text{flt} = 20, \text{wt} = 200$ regime with $>3.7\times$ instruction headroom), and solve for the right kernel vector modulo $r$ via Kaltofen's Block-Wiedemann algorithm (charging Krylov `smvp`, dense Schirokauer columns, quadratic Matrix Berlekamp–Massey generator synthesis, and `Mksol` solution reconstruction in full).
5. **Phase 4 (Individual Discrete Logarithm — LLL Smoothing & 4-Level Special-$\mathfrak{q}$ Descent):** Apply 24-dimensional LLL smoothing ($B_{\text{init}} = 2^{165}$) and a 4-level special-$\mathfrak{q}$ descent tree in $\mathbb{F}_{p^{12}}^*$ for $g$ and $z$ down to the factor base bound $B = 2^{65.316}$, recovering $a \equiv \log(z)\cdot\log(g)^{-1} \pmod r$.
6. **Phase 5 (Output Construction & Verification):** Compute $Q = [a]Q_0 \in G_2$, verify $e(P, Q) = z$ by one full pairing evaluation, and output $Q$.

All challenge-independent phases (Phases 1–3) are charged 100% in full inside total work $W$ with zero amortization.

---

## 2. STNFS Tower Construction and Exact Resultant Norm Formula

### 2.1 Base Field $K_h$, Cubic Subfield $K_3$, and Maximal Order $\mathcal{O}_h$
We select the monic irreducible polynomial of degree $\eta = 6$:
$$h(Y) = Y^6 - Y^2 + 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 1.$$
Let $y = Y \bmod h(Y)$ and $K_h = \mathbb{Q}(y)$. Because $h(Y) = g_3(Y^2)$ with $g_3(Z) = Z^3 - Z + 1$ ($\text{Disc}(g_3) = -23$), $K_h$ contains the cubic subfield $K_3 = \mathbb{Q}(z)$ where $z = y^2$ satisfies $z^3 - z + 1 = 0$.
- **Monogenic Maximal Order and Class Number $1$:** $\text{Disc}(h) = 2^6 \cdot 23^2 = 33856$. By Dedekind's criterion at $\ell = 2$ and $\ell = 23$, $[\mathcal{O}_h : \mathbb{Z}[y]] = 1$, so $\mathcal{O}_h = \mathbb{Z}[y]$. Since $K_h$ is totally complex ($r_{\mathbb{R}} = 0, r_{\mathbb{C}} = 3$), its Minkowski bound is $M_{K_h} = \frac{6!}{6^6}(4/\pi)^3 \sqrt{33856} \approx 5.8706 < 6$. Because $h(Y) \pmod \ell$ has no degree-$1$ factors for $\ell \in \{2, 3, 5\}$ (smallest prime ideal norm is $N(\mathfrak{l}_7) = 7 > 5$), $\mathcal{O}_h$ is a Principal Ideal Domain ($h_{K_h} = 1$). Its torsion unit subgroup is $U = \{\pm 1\}$ ($w = |U|/2 = 1$).
- **Exact Integer Norm via $K_3$:** Any element $\gamma(y) = \sum_{i=0}^5 \gamma_i y^i \in \mathcal{O}_h$ decomposes into even and odd parts $\gamma(y) = \gamma_{\text{even}}(y^2) + y\,\gamma_{\text{odd}}(y^2)$, where $\gamma_{\text{even}}(z) = \gamma_0 + \gamma_2 z + \gamma_4 z^2$ and $\gamma_{\text{odd}}(z) = \gamma_1 + \gamma_3 z + \gamma_5 z^2$. Its relative norm to $K_3 = \mathbb{Q}[z]/(z^3 - z + 1)$ is
  $$N_{K_h/K_3}(\gamma(y)) = \gamma(y)\gamma(-y) = \gamma_{\text{even}}(z)^2 - z\,\gamma_{\text{odd}}(z)^2 \equiv c_0 + c_1 z + c_2 z^2 \pmod{z^3 - z + 1},$$
  and its absolute norm to $\mathbb{Q}$ is the exact $3 \times 3$ integer determinant of multiplication by $c_0 + c_1 z + c_2 z^2$ in $\mathbb{Z}[z]/(z^3 - z + 1)$:
  $$N_{K_h/\mathbb{Q}}(\gamma(y)) = \det \begin{pmatrix} c_0 & -c_2 & -c_1 \\ c_1 & c_0 + c_2 & c_1 - c_2 \\ c_2 & c_1 & c_0 + c_2 \end{pmatrix} \in \mathbb{Z}.$$
- **Three Distinct Complex Conjugate Pairs:** The cubic $g_3(Z) = Z^3 - Z + 1$ has one negative real root $z_1 \approx -1.324717957244746$ and one complex conjugate pair $z_2, \bar{z}_2 \approx 0.662358978622373 \pm 0.5622795120623013\,i$. Taking square roots gives the six roots of $h(Y) = 0$ partitioned into three distinct complex conjugate pairs $\{\theta_j, \bar{\theta}_j\}_{j=1}^3$:
  1. $\theta_1 = +\sqrt{z_1} \approx 1.1509639252968597\,i$ (with $\bar{\theta}_1 = -\sqrt{z_1} = -\theta_1$),
  2. $\theta_2 = +\sqrt{z_2} \approx +0.8748557533622111 + 0.3213589724840035\,i$ (with $\bar{\theta}_2 = +\sqrt{\bar{z}_2}$),
  3. $\theta_3 = -\sqrt{z_2} \approx -0.8748557533622111 - 0.3213589724840035\,i$ (with $\bar{\theta}_3 = -\sqrt{\bar{z}_2}$).
  For any $\gamma(y) \in \mathcal{O}_h$, $\prod_{j=1}^3 |\gamma(\theta_j)|^2 = |N_{K_h/\mathbb{Q}}(\gamma(y))|$. In particular, for $a(y) = 1 + y$ and $b(y) = 0$, both the exact determinant formula and $\prod_{j=1}^3 |(1 + \theta_j)^2|^2$ evaluate to $1$, as verified in `experiment.py`.

### 2.2 Top Polynomials $f_y(X)$ and $g_y(X)$ and Exact Automorphism Reduction ($\mathcal{A} = 2$)
Since $3p = P_6(u)$ where $P_6(U) = U^6 - 2U^5 + 2U^3 + U + 1$ and $y^6 = y^2 - 1$ in $\mathcal{O}_h$, setting
$$g_y(X) = X^2 - u y = X^2 + 15132376222941642752\, y \in \mathcal{O}_h[X]$$
and clearing denominators in $y^6 P_6(X^2/y)$ yields:
$$f_y(X) = y^6 P_6\!\left(\frac{X^2}{y}\right) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + y^2 - 1 \in \mathcal{O}_h[X].$$
Modulo $\mathfrak{p} = (p, h(y))$, $f_y(X) \equiv (y^2 - 1)\cdot 3p \equiv 0 \pmod{\mathfrak{p}, g_y(X)}$, so $\phi_y(X) = X^2 - u y$ is a common irreducible factor of degree $\kappa = 2$ over $\mathbb{F}_{p^6} = \mathbb{F}_p[y]/(h(y))$.

**Exact Automorphism Action on Virtual Logarithms ($\mathcal{A} = 2$):**
Both $f_y(X)$ and $g_y(X)$ are polynomials in $X^2$ and are invariant under the involution $\tau : X \mapsto -X$. In $\mathbb{F}_{p^{12}} = \mathbb{F}_{p^6}[X]/(X^2 - u y)$, since $u y \in \mathbb{F}_{p^6}$ is a quadratic non-residue, the $p^6$-Frobenius automorphism of $\mathbb{F}_{p^{12}} / \mathbb{F}_{p^6}$ acts on $X$ as $X^{p^6} = (X^2)^{(p^6-1)/2} X = (u y)^{(p^6-1)/2} X = -X = \tau(X)$. Because $r \mid (p^4 - p^2 + 1) \mid (p^6 + 1)$, we have $p^6 \equiv -1 \pmod r$. Therefore, for every degree-$1$ prime ideal $\mathfrak{q} \in \mathcal{B}_f \cup \mathcal{B}_g$, its conjugate prime ideal $\tau(\mathfrak{q})$ has virtual logarithm:
$$\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r.$$
Each conjugate pair $\{\mathfrak{q}, \tau(\mathfrak{q})\}$ thus introduces only **one** unknown virtual logarithm modulo $r$. During relation collection, whenever $a(y) - b(y)X$ is $B$-smooth, we fold each prime-ideal factor $\tau(\mathfrak{q})$ into $-\log(\mathfrak{q}) \pmod r$, reducing the number of unknown factor-base virtual logarithms from $2\,\text{Li}(B)$ to $\text{Li}(B)$:
$$\text{SP}_{\text{orb}} = \text{Li}(B) + n_{\text{Sch}} \approx 1.037897 \times 10^{18} \approx 2^{59.864815}.$$
Furthermore, because $N_f(a, -b) = N_f(a, b)$ and $N_g(a, -b) = N_g(a, b)$, $(a(y), b(y))$ and $(a(y), -b(y))$ yield the same orbit relation up to the involution $\tau$, so we restrict sieving to pairs where both $\text{lc}(b) > 0$ and $\text{lc}(a) \ge 0$ (volume $V_{0,\text{orb}} = (2A+1)^{12}/4 = V_0 / 2$).

---

## 3. Exact Computation of $1/\zeta_{K_h}(2)$, $\alpha(f_y, h, 10^3)$, $\alpha(g_y, h, 10^3)$, and Empirical Verification

### 3.1 Euler Product $1/\zeta_{K_h}(2)$ and TNFS Murphy $\alpha$-Values
As computed from first principles in `experiment.py` over all primes $\ell \le B_0 = 1000$:
- **Coprime Ideal Density:**
  $$\frac{1}{\zeta_{K_h}(2)} = \prod_{\ell \le 1000} \prod_{\mathfrak{l} \mid \ell} \left(1 - \frac{1}{N(\mathfrak{l})^2}\right) = 0.903208,$$
  which matches empirical coprime-ideal sampling on $[-A, A]^6 \times ([-A, A]^5 \times [0, A])$ within sampling variance (Guillevic–Singh Table 10 reports $0.9192$; we use the more conservative exact Euler product $0.903208$ for our primary success check).
- **Exact Prime-Ideal Root Valuations:**
  For $F_y \in \{f_y, g_y\}$ (both monic in $X$, so $n_{\mathfrak{l}}^{\text{pro}} = 0$), at every good prime ideal $\mathfrak{l} = (\ell, h_i(y))$ of norm $q = \ell^{\deg h_i}$, $v_{\mathfrak{l}}(F_y) = \frac{n_{\mathfrak{l}}^{\text{aff}} q}{q^2 - 1}$. At the bad prime ideals ($\ell \in \{2, 3, 7, 13\}$ for $f_y$, where $\text{Disc}_U(P_6) = 3^3 \cdot 7 \cdot 13 \cdot 50933$, and $\ell = 2$ for $g_y$), recursive root lifting (`average_val_affine_TNFS_Pr`) yields:
  $$\alpha(f_y, h, 10^3) = 2.1788, \qquad \alpha(g_y, h, 10^3) = 0.9598 \quad (\text{in base } e).$$
- **Empirical Verification in `experiment.py`:**
  1. *Empirical Local Valuations:* Across $N = 5{,}000$ random coprime pairs $(a(y), b(y))$, `experiment.py` computes the exact integer resultants $N_f(a,b), N_g(a,b) \in \mathbb{Z}$ and verifies that the empirical mean valuations $\frac{1}{N}\sum \text{val}_\ell(N_F(a,b))$ at small primes $\ell \in \{2, 3, 5, 7, 11\}$ agree closely with the theoretical $v_\ell(F_y)$.
  2. *Observed Smooth Relations on Small-Norm Pairs:* `experiment.py` enumerates coprime pairs $(a(y), b(y))$ with coefficients in $\{-1, 0, 1\}$, computes their exact integer norms $N_f(a,b)$ and $N_g(a,b)$, factors them by trial division to **directly observe simultaneous smooth relations**, and verifies that $\tau : X \mapsto -X$ preserves smoothness and norm factorizations.

---

## 4. Relation Yield, Systematic Modeling Margin, Matrix Rank, and Descent Analysis

### 4.1 Expected Yield (Algorithm 6.1) and Robustness to Systematic Modeling Error
At $A = 686$ and $B = 2^{65.316}$:
- Unreduced sieving volume: $V_0 = \frac{(2A+1)^{12}}{2} = \frac{1373^{12}}{2} \approx 2.243967 \times 10^{37} \approx 2^{124.077346}$.
- $\mathcal{A}=2$ orbit sieving volume ($\text{lc}(b) > 0, \text{lc}(a) \ge 0$): $V_{0,\text{orb}} = \frac{V_0}{2} \approx 1.121983 \times 10^{37} \approx 2^{123.077346}$.
- Valid coprime orbit space (using the conservative $1/\zeta_{K_h}(2) = 0.903208$):
  $$V_{1,\text{orb}} = V_{0,\text{orb}} \cdot \frac{1}{\zeta_{K_h}(2)} \approx 1.013384 \times 10^{37} \approx 2^{122.930479}.$$
- Using exact integer resultant norms $N_f(a,b), N_g(a,b) \in \mathbb{Z}$ across $N = 5{,}000$ coprime samples in `experiment.py`, the mean norm bit-lengths are $\mathbb{E}[\log_2 N_f] \approx 724.9\text{ bits}$ and $\mathbb{E}[\log_2 N_g] \approx 497.0\text{ bits}$, and the mean Murphy smoothness probability is $\bar{p} \approx 2^{-62.567} \approx 1.447 \times 10^{-19}$.
- The expected number of distinct orbit relations collected in $V_{1,\text{orb}}$ is:
  $$\lambda_{\text{orb}} = V_{1,\text{orb}} \cdot \bar{p} \approx 1.466 \times 10^{18} \approx 2^{60.363},$$
  compared against the required factor-base orbit threshold:
  $$\text{SP}_{\text{orb}} = \text{Li}(2^{65.316}) + 16 \approx 1.037897 \times 10^{18} \approx 2^{59.8648}.$$
- **Systematic Modeling Error Margin:**
  Notice that the oversampling ratio is
  $$\frac{\lambda_{\text{orb}}}{\text{SP}_{\text{orb}}} \approx \frac{1.466 \times 10^{18}}{1.0379 \times 10^{18}} \approx 1.412 \quad (+0.498\text{ bits of excess relations}).$$
  Even if the true mean smoothness probability $\bar{p}_{\text{true}}$ is systematically lower than Murphy's Dickman approximation by $25\%$ ($\bar{p}_{\text{true}} = 0.75\,\bar{p}$), the degraded expected yield is $\lambda_{\text{degraded}} = 0.75\,\lambda_{\text{orb}} \approx 1.0995 \times 10^{18} > \text{SP}_{\text{orb}} \approx 1.0379 \times 10^{18}$, for which the Poisson–Normal z-score is:
  $$Z_{\text{degraded}} = \frac{\text{SP}_{\text{orb}} - \lambda_{\text{degraded}}}{\sqrt{\lambda_{\text{degraded}}}} \approx \frac{-6.16 \times 10^{16}}{1.049 \times 10^9} \approx -5.87 \times 10^7 \ll -10,$$
  so $P(\text{RF} \ge \text{SP}_{\text{orb}}) = 1 - \Phi(Z_{\text{degraded}}) > 0.9999 > 0.40$. (Moreover, because our sieving budget $R_{\text{sieve}}$ below charges the full $V_0 = 2 V_{0,\text{orb}}$ volume, if needed the sieve can also scan $\text{lc}(a) < 0$ for asymmetric sub-lattices.)

### 4.2 Conservative Pruning ($\text{flt}_{\text{cons}} = 4$) and Sparse Matrix Rank
Rather than relying solely on the heuristic $\text{flt} = 20$ reduction criticized in Guillevic–Singh Remark 5, we analyze and budget for a **conservative singleton + $2$/$3$-way exact merge pruning** ($\text{flt}_{\text{cons}} = 4$):
1. **Active Ideals and Singleton/2-Way Pruning:**
   In an unfiltered relation row, since the smallest prime ideal of $\mathcal{O}_h$ has norm $7$, $N_f(a,b) \approx 2^{725}$ and $N_g(a,b) \approx 2^{497}$ have on average $\approx \ln\ln N_f + \ln\ln N_g \approx 12.1$ prime-ideal factors (and at most $\lfloor 725/\log_2 7 \rfloor + \lfloor 497/\log_2 7 \rfloor \le 435$ factors). By Mertens' theorem and Poisson hit distribution across prime ideals of norm $\ell \in [7, B]$, ideals near $B$ have mean hit count $\lambda(\ell) = \lambda_{\text{orb}} / \ell$; ideals with $0$ hits do not appear in the matrix (reducing the active column count), ideals with $1$ hit (singletons) are removed along with their unique relation, and ideals with $2$ or $3$ hits are eliminated by exact $2$-way/$3$-way row additions. This standard exact pruning reduces the active unknown count to $N_{\text{mat}} \le \frac{1}{4}\text{Li}(B) \approx 2.5947 \times 10^{17} \approx 2^{57.8648}$ while keeping average row weight $\text{wt} \le 64$.
2. **Useful Matrix Rank and Block-Wiedemann Kernel Recovery:**
   With $\lambda_{\text{orb}} \approx 1.466 \times 10^{18}$ collected relations, the excess of relations over active ideals ensures that after singleton/2-way pruning, we retain $N_{\text{mat}} + C_{\text{excess}}$ rows ($C_{\text{excess}} \ge 10^6$). Taking $N_{\text{mat}}$ random linear combinations over $\mathbb{F}_r$ of the excess rows yields a square $N_{\text{mat}} \times N_{\text{mat}}$ sparse matrix $M$ over $\mathbb{F}_r$ whose right kernel contains the true virtual-logarithm vector of the active factor-base ideals. By Kaltofen [4, Theorem 7], Block-Wiedemann with blocking parameters $n = 2, m = 4$ over $\mathbb{F}_r$ ($r > 2^{254}$) recovers a non-zero kernel vector with failure probability at most $2 N_{\text{mat}} / r < 2^{-195}$.

### 4.3 Parameter-Level Derivation of Phase 4 (Individual Logarithm Descent)
We derive the exact parameters and cost of Phase 4 for $\tau \in \{g, z\} \subset G_T \subset \mathbb{F}_{p^{12}}^*$:
1. **Step 4a — 24-Dimensional LLL Smoothing in $\mathbb{F}_{p^{12}}^*$ ($B_{\text{init}} = 2^{165}$):**
   For a random exponent $e \in [1, r-1]$, set $\tilde{\tau} = \tau^e \in \mathbb{F}_{p^{12}}^*$. Consider the $\mathbb{Z}$-module $\Lambda_{\tilde{\tau}} = \{(a(y) + b(y)X,\; c(y) + d(y)X) \in \mathcal{O}_h[X]^2 : (a + bX) \equiv \tilde{\tau}(c + dX) \pmod{\mathfrak{p}, X^2 - u y}\} \subset \mathbb{Z}^{24}$, which has determinant $p^{12}$ in dimension $24$. Scaling the $b, d$ coordinates by $|u|^{1/2} \approx 2^{31.855}$ to balance $a^2$ and $u y b^2$ in $N_g$ and applying LLL reduction produces a reduced basis $v_1, \dots, v_{24} \in \mathbb{Z}^{24}$ in which short vectors have $\|a\|_\infty, \|c\|_\infty \approx p^{1/2}|u|^{1/4} \approx 2^{206}$ and $\|b\|_\infty, \|d\|_\infty \approx p^{1/2}|u|^{-1/4} \approx 2^{174}$.
   Consequently, $a(y)^2 - u y b(y)^2$ and $c(y)^2 - u y d(y)^2$ have coefficients of size $\approx 2^{412}$, so their norms $N_g(a,b)$ and $N_g(c,d)$ in the degree-$6$ field $K_h$ have bit-length:
   $$\log_2 N_g \approx 6 \times 412 = 2472\text{ bits}.$$
   We test linear combinations of the LLL basis until both $N_g(a,b)$ and $N_g(c,d)$ are $B_{\text{init}}$-smooth with:
   $$B_{\text{init}} = 2^{165}.$$
   The Dickman smoothness parameter for each $2472$-bit norm is $u_{\text{init}} = \frac{2472}{165} = 14.9818$, giving $\rho(14.9818) \approx 2^{-60.4}$. Requiring $N_g(a,b)$ and $N_g(c,d)$ to be $B_{\text{init}}$-smooth (by testing $a+bX$ and $c+dX$ via a bipartite birthday/sieve search or testing $2^{64}$ one-sided smooth numerators/denominators against $g^e$) has success probability $\approx 2^{-60.4}$ per candidate. Testing $N_{\text{init}} = 2^{64}$ candidates yields an expected $\lambda_{\text{init}} = 2^{64 - 60.4} \approx 12.1$ valid $B_{\text{init}}$-smooth representations per target (failure probability $\exp(-12.1) < 6 \times 10^{-6}$). At $\le 2^{34}$ scalar AArch64 instructions per candidate (for LLL/norm evaluation and early-abort ECM smoothness testing), Step 4a for both $g$ and $z$ costs at most $2 \times 2^{64} \times 2^{34} = 2^{99}$ instructions.
2. **Step 4b — 4-Level Special-$\mathfrak{q}$ Descent Tree from $B_{\text{init}} = 2^{165}$ to $B = 2^{65.316}$:**
   Each prime ideal $\mathfrak{q}$ dividing $N_g(a,b)$ or $N_g(c,d)$ with norm $Q = N(\mathfrak{q}) \in (2^{65.316}, 2^{165}]$ is eliminated by special-$\mathfrak{q}$ lattice sieving:
   - The divisibility condition $\mathfrak{q} \mid (a'(y) - b'(y)X)$ defines an index-$Q$ sublattice $L_{\mathfrak{q}} \subset \mathbb{Z}^{12}$ whose LLL-reduced basis vectors have infinity norm $\approx Q^{1/12}$.
   - Sieving linear combinations with coefficients bounded by $A_{\text{desc}} = 2^{4.5}$ (sieve volume $V_{\text{desc}} = 2^{56}$ per node) produces pairs $(a', b') \in L_{\mathfrak{q}}$ with $\|a'\|_\infty, \|b'\|_\infty \approx 2^{4.5} Q^{1/12}$.
   - Dividing out $Q$, the residual norms on the $f_y$ and $g_y$ sides have bit-lengths:
     $$\log_2 N_f' \approx 6(12 \log_2(2^{4.5} Q^{1/12}) + 2) = 6 \log_2 Q + 336\text{ bits},$$
     $$\log_2 N_g' \approx 6(2 \log_2(2^{4.5} Q^{1/12}) + 63.71) - \log_2 Q = 436\text{ bits}.$$
   - Enforcing that all prime factors of $N_f'$ and $N_g'$ are $\le Q_{\text{next}} = Q^{0.75}$ (and discarding any cofactor split with more than $k_{\text{branch}} = 6$ prime ideals above $B$):
     - At Level 0 ($Q_0 \le 2^{165}$): $\log_2 N_f' \le 1326$, $\log_2 N_g' \le 436$, $Q_1 \le 2^{123.75}$, $u_f + u_g \le \frac{1326 + 436}{123.75} = 14.24 \implies \rho(u_f)\rho(u_g) \ge 2^{-47}$.
     - At Level 1 ($Q_1 \le 2^{123.75}$): $Q_2 \le 2^{92.81}$, $u_f + u_g \le \frac{1078.5 + 436}{92.81} = 16.32 \implies \rho(u_f)\rho(u_g) \ge 2^{-49}$.
     - At Level 2 ($Q_2 \le 2^{92.81}$): $Q_3 \le 2^{69.61}$, $u_f + u_g \le \frac{892.9 + 436}{69.61} = 19.09 \implies \rho(u_f)\rho(u_g) \ge 2^{-51}$.
     - At Level 3 ($Q_3 \le 2^{69.61}$): $Q_4 \le 2^{52.21} < B = 2^{65.316}$ (all resulting prime ideals lie in the factor base $\mathcal{B}$, and any inactive ideal outside the pruned matrix is rejected during node selection).
   - With sieve volume $V_{\text{desc}} = 2^{56}$ per node and joint smoothness probability $\ge 2^{-51}$, each node yields on average $\ge 2^5 = 32$ smooth relations (per-node failure probability $\le e^{-32} < 10^{-13}$).
   - Starting from at most $\lceil 2472 / 65.316 \rceil \le 38$ prime ideals at Level 0 and branching by at most $k_{\text{branch}} = 6$ over $4$ levels, the total number of descent tree nodes per target is at most $38 \times (1 + 6 + 6^2 + 6^3) = 9{,}842 < 2^{14}$ nodes.
   - At $\le 2^{64}$ instructions per node (covering $V_{\text{desc}} = 2^{56}$ sieve steps and ECM cofactorization), Step 4b for both targets costs at most $2 \times 2^{14} \times 2^{64} = 2^{79}$ instructions.
   - Combining Step 4a ($2^{99}$ instructions) and Step 4b ($2^{79}$ instructions), **total Phase 4 work is bounded by $R_{\text{descent}} = 2^{100}$ scalar AArch64 instructions**, with total descent failure probability $< 2 \times 6 \times 10^{-6} < 2 \times 10^{-5}$.

---

## 5. Complete Scalar AArch64 Instruction Ledger and Memory Accounting

Per `COST_MODEL.md` and `classical-group-work-v2.json`, total computation is $W = H + R / C$ where $C = 20,000$ scalar AArch64 instructions per reference dense $\mathbb{F}_{p^{12}}$ multiplication ($H = 0$, all operations charged in $R$).

### 5.1 First-Principles Derivation of Phase 2 (12D Bucket/Line Sieving)
For each degree-$1$ prime ideal $\mathfrak{q} = (\ell, y - r_y, X - r_X) \in \mathcal{B}_f \cup \mathcal{B}_g$ of norm $\ell \le B$, divisibility $\mathfrak{q} \mid (a(y) - b(y)X)$ is the linear congruence
$$a_0 \equiv r_X \sum_{i=0}^5 b_i r_y^i - \sum_{i=1}^5 a_i r_y^i \pmod \ell.$$
Because the coefficient of $a_0$ is $1 \not\equiv 0 \pmod \ell$, $\mathfrak{q}$ hits exactly a $1/\ell$ fraction of all points in $V_0 = \frac{1}{2}(2A+1)^{12}$. By Mertens' Second Theorem and the Prime Ideal Theorem, the total number of prime-ideal hits across both sides ($\mathcal{B}_f$ and $\mathcal{B}_g$) over $V_{0,\text{orb}} = V_0 / 2$ (or per side over $V_0$) is:
$$N_{\text{hits}} = 2 V_{0,\text{orb}} \sum_{\ell \le B} \frac{1 + o(1)}{\ell} = V_0 (\ln\ln B + M_{\text{Mertens}}) \approx 4.0742 \cdot V_0.$$
In cache-blocked bucket sieving (for primes $\ell > 2A+1$) and line sieving (for $\ell \le 2A+1$):
- Each hit requires $12$ instructions for modular offset stepping and 32-bit `(offset, logp)` bucket insertion, $12$ instructions for bucket-array log-norm accumulation (`ldrb`/`add`/`strb`), and $8$ instructions amortized for Gray-code outer hyperplane stepping ($32$ scalar AArch64 instructions per hit).
- Initializing, threshold-scanning (via 64-bit word mask `0x8080808080808080`), and clearing the $f$-side and $g$-side sieve byte arrays across $V_0$ locations takes at most $16$ scalar AArch64 instructions per location.
- Total sieving instructions:
  $$R_{\text{sieve, exact}} = \big(32(\ln\ln B + 0.2615) + 16\big) V_0 \approx 146.375 \cdot V_0 \approx 3.2846 \times 10^{39}\text{ instructions}.$$
  We allocate the explicit upper bound:
  $$R_{\text{sieve}} = 4 \times 10^{39} \approx 2^{131.557503}\text{ scalar AArch64 instructions}.$$

### 5.2 First-Principles Derivation of Phase 3 (Low-Merge Block-Wiedemann + Reconstruction)
For the pruned sparse matrix of dimension $N_{\text{mat}} \le \frac{1}{4}\text{Li}(B) \approx 2.59474 \times 10^{17}$ and average row weight $\text{wt} \le 64$ (plus $n_{\text{Sch}} = 16$ dense Schirokauer unit columns) over $\mathbb{F}_r$ ($4$ limbs of $64$ bits):
- **Krylov Stage ($n = 2, m = 4$):** $n(N_{\text{mat}}/m + N_{\text{mat}}/n) = 1.5 N_{\text{mat}}$ sparse matrix-vector products (`smvp`), each performing $(\text{wt} + n_{\text{Sch}}) N_{\text{mat}} = 80 N_{\text{mat}}$ multiply-accumulates in $\mathbb{F}_r$, totaling $120 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$.
- **Matrix Berlekamp–Massey / Generator Synthesis:** Computing the $(m+n) \times (m+n) = 6 \times 6$ matrix polynomial generator of degree $\lceil N_{\text{mat}}/n \rceil$ via the deterministic quadratic algorithm takes at most $(m+n)^3 N_{\text{mat}}^2 = 216 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$.
- **`Mksol` Solution Reconstruction & Verification:** $0.5 N_{\text{mat}}$ `smvp` evaluations ($40 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$) plus $N_{\text{mat}}$ dense vector linear combinations and pruning back-substitution ($32 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$), totaling $72 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$.
- **Total $\mathbb{F}_r$ Operations:** $(120 + 216 + 72) N_{\text{mat}}^2 = 408 N_{\text{mat}}^2 \approx 2.7469 \times 10^{37}$ operations in $\mathbb{F}_r$.
- Charging $256$ scalar AArch64 instructions per $4$-limb $\mathbb{F}_r$ multiply-accumulate + sparse index/vector load gives $256 \times 408 N_{\text{mat}}^2 \approx 7.0321 \times 10^{39}$ scalar AArch64 instructions (which also exceeds the Guillevic–Singh $\text{flt}=20, \text{wt}=200$ cost of $256 \times 680 \times (2\,\text{Li}(B)/20)^2 \approx 1.875 \times 10^{39}$ instructions by $>3.7\times$). We allocate:
  $$R_{\text{LA}} = 8 \times 10^{39} \approx 2^{132.557503}\text{ scalar AArch64 instructions}.$$

### 5.3 Summary Instruction Ledger and Score Normalization ($C = 20,000$)

| Phase / Component | Scalar AArch64 Instructions ($R_i$) | $\log_2(R_i)$ | Target-Group-Op Equivalents ($R_i / 20{,}000$) |
| --- | ---: | ---: | ---: |
| Phase 0: Input check & $g = e(P, Q_0)$ | $16,777,216$ | $24.000000$ | $838.8608$ |
| Phase 1: Poly selection & factor base | $75,557,863,725,914,323,419,136$ | $76.000000$ | $3.777893 \times 10^{18}$ |
| Phase 2a: 12D Bucket/Line Sieving | $4,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $131.557503$ | $2.000000 \times 10^{35}$ |
| Phase 2b: ECM cofactorization & Schirokauer | $1,267,650,600,228,229,401,496,703,205,376$ | $100.000000$ | $6.338253 \times 10^{25}$ |
| Phase 3: Pruning, BW Krylov, BM & `Mksol` | $8,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $132.557503$ | $4.000000 \times 10^{35}$ |
| Phase 4: LLL smoothing & 4-level descent | $1,267,650,600,228,229,401,496,703,205,376$ | $100.000000$ | $6.338253 \times 10^{25}$ |
| Phase 5: $Q = [a]Q_0$ & $e(P, Q) = z$ check | $16,777,216$ | $24.000000$ | $838.8608$ |
| **Total Attack Work** | **$12,000,000,002,535,301,200,532,016,666,718,622,875,648$** | **$133.142459$** | **$6.000000 \times 10^{35}$** |

Dividing $R_{\text{total}}$ by $C = 20,000$ gives:
$$W = \frac{R_{\text{total}}}{20,000} \approx 6.00000000126765 \times 10^{35} \approx 2^{118.854747} < 2^{118.860000}.$$
We therefore submit the exact rational upper bound:
$$\text{time\_log2} = 118.860000 = \frac{5943}{50}, \qquad \text{preprocessing\_log2} = 118.860000.$$

### 5.4 Complete Peak Memory Ledger (`memory_log2_bytes = 73.000000`)
We account for all storage simultaneously resident at peak:
1. **Factor Base Table ($\mathcal{B}_f \cup \mathcal{B}_g$):** $2.0758 \times 10^{18}$ prime-ideal records at $32$ bytes each (`norm_l`, `r_y`, `r_X`, flags) $\le 6.643 \times 10^{19}\text{ bytes} \approx 2^{65.86}\text{ bytes}$.
2. **Unfiltered Relation Table ($\text{RF} \le 2^{61}$ relations):** Each unfiltered relation stores its $12$ `int16` polynomial coefficients $(a_0, \dots, a_5, b_0, \dots, b_5)$ ($24$ bytes), up to $64$ factor-base ideal orbit indices ($64 \times 8 = 512$ bytes) and signed multiplicities ($64 \times 2 = 128$ bytes), and $16$ precomputed $32$-byte Schirokauer map evaluations modulo $r$ ($512$ bytes), padded to $2{,}048 = 2^{11}$ bytes per relation:
   $$M_{\text{unfiltered}} = 2^{61} \times 2^{11}\text{ bytes} = 2^{72.00}\text{ bytes}.$$
3. **Pruned Matrix Storage (Column Indices + Full $32$-Byte $\mathbb{F}_r$ Coefficients):** For $N_{\text{mat}} \le 2^{58}$ rows and up to $216$ non-zero entries per row (including $16$ Schirokauer columns), storing an $8$-byte column index **and** a full $32$-byte $\mathbb{F}_r$ coefficient ($40$ bytes per entry) takes:
   $$M_{\text{matrix}} = 2^{58} \times 216 \times 40\text{ bytes} \approx 2.490 \times 10^{21}\text{ bytes} \approx 2^{71.08}\text{ bytes}.$$
4. **Block-Wiedemann Working Vectors, Virtual Log Table & Sieve Buckets:** $32$ vectors of length $N_{\text{mat}} \le 2^{58}$ over $\mathbb{F}_r$ ($32 \times 2^{58} \times 32 = 2^{68}$ bytes) plus the virtual logarithm table ($\text{SP}_{\text{orb}} \times 32 \approx 2^{64.86}$ bytes) and cache/bucket buffers ($< 2^{40}$ bytes) $\le 2^{68.20}\text{ bytes}$.
5. **Total Peak Memory:**
   $$M_{\text{total}} = M_{\text{FB}} + M_{\text{unfiltered}} + M_{\text{matrix}} + M_{\text{working}} < 2^{72.60}\text{ bytes} < 2^{73.000000}\text{ bytes}.$$
   Nonuniform advice is $0$ bytes ($\text{nonuniform\_advice\_log2\_bytes} = 0$).

---

## 6. Summary of Heuristics and References
1. **Murphy $\alpha$/$E$ Smoothness Model & Systematic Error Margin:** Validated locally in `experiment.py` against exact integer resultant valuations and small-norm smooth relation counts; the $+41.2\%$ oversampling margin ($\lambda_{\text{orb}} / \text{SP}_{\text{orb}} \approx 1.412$) ensures $P(\text{RF} \ge \text{SP}_{\text{orb}}) > 0.40$ even under a $25\%$ systematic downward error in $\bar{p}$.
2. **Low-Merge Pruning ($\text{flt}_{\text{cons}} = 4$) and Sparse Rank:** Requires only singleton removal and $2$/$3$-way merges rather than aggressive $20\times$ merge filtering, with full Block-Wiedemann reconstruction charged.
3. **Special-$\mathfrak{q}$ Descent Contraction:** Uses a 4-level descent tree with contraction factor $0.75$ from $B_{\text{init}} = 2^{165}$ to $B = 2^{65.316}$ following Guillevic [3].

### References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. E. Kaltofen, "Analysis of Coppersmith's Block Wiedemann Algorithm for the Parallel Solution of Sparse Linear Systems," *Mathematics of Computation*, 64(210):777–806, 1995.