# Comments
- Judge submission to test TNFS attack,
- Follows the paper from Guillevic-Singh, including parameter choice from Table 10,
- If accepted, future submissions should follow the heuristics for estimatives and probabilistic model, while experimenting with fine tuned parameters. Different heuristics and probabilistic models should be independently justified.
- Created by Gemini Argon, human supervised.

# bls12-381-fixed-g1-invert-g2: Special Tower Number Field Sieve (STNFS) via $G_T$ DLP Reduction (Revision v3)

## 0. Revision Record (Addressing Review of Prior Submission `631ee3e7` and `57d7de1`)

### 0.1 Baseline / Prior Work (`631ee3e7` and `57d7de1`)
In prior submission `57d7de1` (`stnfs-guillevic-singh`, claiming `time_log2 = 115.56`) and revision `631ee3e7` (`stnfs-guillevic-singh-v2`, claiming `time_log2 = 118.86`), we formulated a Special Tower Number Field Sieve (STNFS) attack on `bls12-381-fixed-g1-invert-g2` via the reduction $g = e(P, Q_0) \in G_T$, $g^a = z \text{ in } \mathbb{F}_{p^{12}}^*$, and $Q = [a]Q_0 \in G_2$, using the polynomial tuple $h(Y) = Y^6 - Y^2 + 1$, $g_y(X) = X^2 - u y$, and $f_y(X) = X^{12} - 2yX^{10} + 2y^3X^6 + y^5X^2 + y^2 - 1$ listed in Table 11 of Guillevic–Singh (2021).

While the judge confirmed in both reviews that the pairing-to-$G_T$-DLP reduction is sound and that supported STNFS heuristics can qualify under `exploratory-v1` without a production-scale run, submission `631ee3e7` was rejected for four material reasons:
1. **Fatal Algebraic Flaw in the Guillevic–Singh Table 11 Polynomial Tuple ($X^2 - u y$ Splits over $\mathbb{F}_{p^6}$):**
   For $h(Y) = Y^6 - Y^2 + 1$, $\deg h = 6$ is even and $h(0) = +1$, so if $h(Y)$ is irreducible modulo $p$, $\text{Norm}_{\mathbb{F}_{p^6}/\mathbb{F}_p}(y) = (-1)^6 h(0) = +1$ and $\text{Norm}_{\mathbb{F}_{p^6}/\mathbb{F}_p}(u y) = u^6$. Consequently, $(u y)^{(p^6 - 1)/2} = (u^6)^{(p - 1)/2} = (u^3)^{p - 1} \equiv +1 \pmod p$. Thus $u y$ is a square in $\mathbb{F}_{p^6}$, so $g_y(X) = X^2 - u y$ splits into two linear factors over $\mathbb{F}_{p^6}$ rather than defining $\mathbb{F}_{p^{12}}$! (More generally, over any prime $p \equiv 3 \pmod 4$, no even polynomial $h(Y) = g_3(Y^2)$ can simultaneously be irreducible over $\mathbb{F}_p$ and have $y$ be a non-square in $\mathbb{F}_{p^6}$, because $\text{Norm}_{\mathbb{F}_{p^3}/\mathbb{F}_p}(y^2) = -g_3(0)$ and $\text{Norm}_{\mathbb{F}_{p^6}/\mathbb{F}_p}(y) = +g_3(0)$ have opposite signs.)
2. **Runtime `ValueError` at `log2(avg_p_fg)` and Experimental Gaps:**
   - In `DickmanRho.__init__`, integrating $u \rho'(u) = -\rho(u - 1)$ via subtractive forward quadrature from $\rho(1) = 1.0$ suffered from IEEE-754 `float64` catastrophic cancellation around $u \approx 11.2$ ($\rho(u) \approx 10^{-11}$), clamping `rho(uf)` to `0.0` and causing `math.log2(avg_p_fg)` to exit with `ValueError`.
   - `observe_small_norm_smooth_relations()` tested smoothness only of $N_f$ rather than simultaneous smoothness of $N_f$ and $N_g$; Algorithm 6.1 checked only integer coefficient $\gcd = 1$ rather than ideal coprimality $a\mathcal{O}_h + b\mathcal{O}_h = \mathcal{O}_h$; and the $25\%$ degradation check was hypothetical rather than backed by a simultaneous-smoothness calibration experiment.
3. **Initial Descent (Step 4a) Joint-Probability Defect:**
   Section 4.3 required both numerator and denominator LLL norms $N_g(a,b)$ and $N_g(c,d)$ ($2472$ bits each) to be $B_{\text{init}}$-smooth ($B_{\text{init}} = 2^{165}$), but used the single-norm smoothness probability $\rho(14.98) \approx 2^{-60.4}$ instead of the joint probability $\rho(14.98)^2 \approx 2^{-120.8}$.
4. **Unestablished 4-Fold Filtering Reduction ($\text{flt} = 4$) and Random Row Combinations:**
   Section 4.2 assumed a $4\times$ filtering reduction ($\text{flt} = 4, \text{wt} \le 64$) without proving occupancy/rank, and proposed taking dense random linear combinations of excess rows (which destroys row sparsity).

### 0.2 Corrections & Implementation Changes in Revision v3
In this revision (`stnfs-guillevic-singh-v3`), we resolve all four defects from first principles:
1. **Provably Valid Algebraic Tower over $h(Y) = Y^6 - 2Y^3 + 2$ (Fixes Defect 1):**
   - From `target-profile.json`, the BLS12-381 target field $\mathbb{F}_{p^{12}}$ is defined by the irreducible polynomial $W^{12} - 2W^6 + 2 \in \mathbb{F}_p[W]$. Setting $Y = W^2$ gives the degree-$6$ tower polynomial:
     $$h(Y) = Y^6 - 2Y^3 + 2 \in \mathbb{Z}[Y].$$
   - Because $h(W^2) = W^{12} - 2W^6 + 2$ is irreducible of degree $12$ over $\mathbb{F}_p$, $h(Y)$ is **provably irreducible of degree $6$ over $\mathbb{F}_p$**.
   - Furthermore, $\text{Norm}_{\mathbb{F}_{p^6}/\mathbb{F}_p}(y) = (-1)^6 h(0) = 2$, so $\text{Norm}_{\mathbb{F}_{p^6}/\mathbb{F}_p}(u y) = 2 u^6$. Because the BLS12-381 prime $p$ satisfies $p \equiv 3 \pmod 8$ (its hex representation ends in `...aaab`, and $\text{0xab} = 171 \equiv 3 \pmod 8$), $2$ is a quadratic non-residue modulo $p$ ($\left(\frac{2}{p}\right) = (-1)^{(p^2-1)/8} = -1$). Therefore:
     $$(u y)^{(p^6 - 1)/2} = (2 u^6)^{(p - 1)/2} = 2^{(p - 1)/2} (u^3)^{p - 1} \equiv (-1) \cdot 1 = -1 \pmod p.$$
     Thus $u y$ is **provably a quadratic non-residue in $\mathbb{F}_{p^6}$**, $g_y(X) = X^2 - u y$ is **provably irreducible of degree $2$ over $\mathbb{F}_{p^6}$**, and the $p^6$-Frobenius acts on $X$ as $X^{p^6} = X(u y)^{(p^6-1)/2} = -X = \tau(X)$.
   - Using $y^6 = 2y^3 - 2$ in $\mathbb{Z}[y]/(h(y))$ and $3p = P_6(u) = u^6 - 2u^5 + 2u^3 + u + 1$, clearing denominators in $y^6 P_6(X^2/y)$ gives:
     $$f_y(X) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + 2 y^3 - 2 \in \mathbb{Z}[y][X],$$
     which has $\deg_X f_y = 12$, $\|f_y\|_\infty = 2$, only even powers of $X$, and satisfies $f_y(X) \equiv (2y^3 - 2)\cdot 3p \equiv 0 \pmod{p, h(y), g_y(X)}$.
   - In `experiment.py`, `verify_algebraic_and_tower_identities()` machine-checks Rabin's irreducibility conditions $\gcd(Y^{p^d} - Y, h(Y)) \pmod p = 1$ ($d \in \{1,2,3\}$), $Y^{p^6} \equiv Y \pmod{(p, h(Y))}$, and $(u y)^{(p^6-1)/2} \equiv -1 \pmod{(p, h(y))}$ on the exact 381-bit prime $p$.
2. **Strictly Positive Dickman $\rho$ Solver, Ideal Coprimality, and Simultaneous $(N_f, N_g)$ Calibration (Fixes Defect 2):**
   - We replace subtractive quadrature in `DickmanRho` with the strictly positive integral identity $\rho(u) = \frac{1}{u}\int_{u-1}^u \rho(t)\,dt \iff \rho(u_i) = \frac{h(\frac{1}{2}\rho(u_{i-M}) + \sum_{j=1}^{M-1}\rho(u_{i-j}))}{u_i - h/2}$. Because every term is positive with zero subtractions, $\rho(u)$ is strictly positive for all $u > 0$, eliminating catastrophic cancellation and the `ValueError`.
   - In Algorithm 6.1, we enforce **ideal coprimality** by computing the exact integer norms $N_a = N_{K_h/\mathbb{Q}}(a(y))$ and $N_b = N_{K_h/\mathbb{Q}}(b(y))$ and checking $\gcd(N_a, N_b)$ (verifying $\gcd_{\mathbb{F}_\ell[y]}(a, b, h) = 1$ at any common prime divisor $\ell \mid \gcd(N_a, N_b)$).
   - In `observe_simultaneous_smooth_relations_and_calibration()`, we test ideal-coprime pairs $(a(y), b(y))$ on a scaled BLS12 instance ($u_{\text{scaled}} = -6$, $B_{\text{scaled}} = 300$), factor **both** $N_f(a,b)$ and $N_g(a,b)$ completely to **observe simultaneous $(N_f, N_g)$ smooth relations**, and compare the observed simultaneous relation count directly against the Murphy-$E$ prediction on the same sample set.
3. **True Joint Smoothness Probability $\rho(u_{\text{init}})^2$ in Initial Descent (Fixes Defect 3):**
   - In Section 4.3 (Step 4a), we set the initial bootstrapping smoothness bound to $B_{\text{init}} = 2^{210}$. For each $2484$-bit LLL-reduced norm $N_g(a,b)$ and $N_g(c,d)$, $u_{\text{init}} = \frac{2484}{210} = 11.8286$ gives single-norm smoothness probability $\rho(11.8286) \approx 2^{-41.13}$, and we use the **true joint probability** that **both** numerator and denominator norms are simultaneously $B_{\text{init}}$-smooth:
     $$P_{\text{joint}} = \rho(11.8286)^2 \approx \left(2^{-41.13}\right)^2 = 2^{-82.26}.$$
   - Testing $N_{\text{init}} = 2^{88}$ LLL linear combinations yields $\lambda_{\text{init}} = 2^{88 - 82.26} \approx 53.4$ expected joint smooth pairs per target (failure probability $< 10^{-23}$) and costs $2 \times 2^{88} \times 2^{32} = 2^{121}$ instructions ($\ll 2^{138}$ main sieve work), followed by a 5-level special-$\mathfrak{q}$ descent tree from $2^{210}$ to $B = 2^{65.316}$.
4. **Unfiltered ($\text{flt} = 1$) Block-Wiedemann Linear Algebra, No Row Combinations, and $A = 1024$ Rebalancing (Fixes Defect 4):**
   - We **completely eliminate merge filtering ($\text{flt} = 1$)** from our cost and memory bounds: the matrix dimension is bounded directly by the full unmerged $\mathcal{A}=2$ factor-base size $N_{\text{mat}} = \text{SP}_{\text{orb}} = \lceil \text{Li}(B) \rceil + 16 \approx 1.0379 \times 10^{18} \approx 2^{59.8648}$.
   - We **never combine rows**: each row of the sparse matrix $M$ is a single distinct raw smooth relation $(a(y), b(y))$, whose prime-ideal multiplicity sum is deterministically bounded by $\log_2 N_f + \log_2 N_g \le 1320$ and has mean weight $\approx 28.2$ (we charge $\text{wt} = 128$ across all $N_{\text{mat}}$ rows).
   - Because $\text{flt} = 1$ increases linear-algebra work to $R_{\text{LA}} \approx 1.39 \times 10^{41} \approx 2^{136.68}$ instructions, we set $A = 1024$ (with $B = 2^{65.316}$) so that 12D sieving ($R_{\text{sieve}} \approx 4.03 \times 10^{41} \approx 2^{138.21}$ instructions) and unfiltered $\text{flt} = 1$ linear algebra are naturally balanced, while boosting the expected yield to $\lambda_{\text{orb}} \approx 2^{63.16}$ against $\text{SP}_{\text{orb}} \approx 2^{59.86}$—an **oversampling ratio of $\lambda_{\text{orb}} / \text{SP}_{\text{orb}} \approx 9.8\times$** that preserves $P(\text{RF} \ge \text{SP}_{\text{orb}}) > 0.9999 > 0.40$ even under an **$80\%$ ($5\times$) systematic downward error** in the Dickman smoothness model.

### 0.3 Revised Experimental Evidence Summary
Running `experiment.py` with the default initial inputs (`alpha_prime_bound = 1000`, `mc_samples = 4000`, `seed = 42`, `A = 1024`, `B_log2_milli = 65316`) completes in $< 3$ seconds and outputs:

| Metric / Audit Item | Prior (`631ee3e7`) | Revised (`v3`) | Verification in `experiment.py` |
| --- | --- | --- | --- |
| Tower polynomial $h(Y)$ | $Y^6 - Y^2 + 1$ ($X^2-uy$ splits!) | **$Y^6 - 2Y^3 + 2$ ($h(w^2)=w^{12}-2w^6+2$)** | Rabin irreducibility mod $p$ verified |
| Quadratic non-residuosity $(u y)^{(p^6-1)/2} \bmod p$ | $+1$ (reducible $g_y$) | **$-1 \equiv p - 1 \pmod p$ (irreducible $g_y$)** | Verified via $\text{Norm}(u y) = 2 u^6$ & $p \equiv 3 \pmod 8$ |
| Dickman $\rho(u)$ solver at $u \in [1, 16]$ | Subtractive (`ValueError` at $0.0$) | **Positive integral average ($\rho(u) > 0$)** | Verified strictly positive & accurate |
| Simultaneous $(N_f, N_g)$ smooth calibration | $N_f$ only | **Observed joint $(N_f, N_g)$ smooth pairs** | Calibrated against Murphy-$E$ prediction |
| Coprimality check in Algorithm 6.1 | Integer $\gcd(a,b)=1$ only | **Ideal coprimality $\gcd(N(a),N(b), \dots)$** | Empirical coprime ratio matches $1/\zeta_{K_h}(2)$ |
| Oversampling ratio $\lambda_{\text{orb}} / \text{SP}_{\text{orb}}$ | $1.41\times$ ($A = 686$) | **$\approx 9.8\times$ ($A = 1024, B = 2^{65.316}$)** | Survives $80\%$ ($5\times$) systematic model error |
| Initial descent Step 4a joint probability | Used single-norm $2^{-60.4}$ | **True joint $\rho(11.83)^2 \approx 2^{-82.26}$ ($B_{\text{init}}=2^{210}$)** | $2^{88}$ trials $\implies \lambda_{\text{init}} \approx 53.4$ ($2^{121}$ inst) |
| Linear algebra filtering assumption | Assumed $\text{flt} = 4$ & row combos | **Unfiltered $\text{flt} = 1$, raw sparse rows** | $N_{\text{mat}} = \lceil\text{Li}(B)\rceil + 16$, $\text{wt} = 128$ charged |
| Submitted `time_log2` ($R_{\text{total}} / 20{,}000$) | $118.86$ (`5943/50`) | **$124.51$ (`12451/100`)** | Strictly bounds $\log_2(W) = 124.501227$ |
| Submitted `memory_log2_bytes` | $73.0$ | **$76.0$ ($2^{76}$ bytes)** | Covers $2^{64}$ raw relations & $\text{flt}=1$ BW vectors |

---

## 1. Problem Reduction and Complete Attack Architecture

The challenge supplies a uniformly sampled nonidentity point $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and an independently uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^*$ on BLS12-381, where:
- $u = -\text{0xd201000000010000} = -15132376222941642752 = -(2^{63} + 2^{62} + 2^{60} + 2^{57} + 2^{48} + 2^{16})$,
- $p = \frac{(u-1)^2(u^4 - u^2 + 1)}{3} + u = \frac{u^6 - 2u^5 + 2u^3 + u + 1}{3}$ ($381$ bits, $\log_2 p \approx 380.0640$, $p \equiv 3 \pmod 8$),
- $r = u^4 - u^2 + 1$ ($255$-bit prime subgroup order, $\log_2 r \approx 254.8571$),
- $Q_0 \in G_2$ is the fixed public generator of the prime-order-$r$ subgroup $G_2$.

As permitted by `rules.md`, we reduce fixed-$G_1$ pairing inversion into $G_2$ to the Discrete Logarithm Problem (DLP) in the order-$r$ cyclotomic subgroup $G_T \subset \mathbb{F}_{p^{12}}^*$:
1. **Phase 0 (Input Validation & Pairing Evaluation):** Decode $P \in G_1$ and $z \in G_T$ canonically and verify curve/subgroup membership. If $z = 1$, return $Q = \mathcal{O}_{G_2}$ (`null`). Otherwise, compute $g = e(P, Q_0) \in G_T \subset \mathbb{F}_{p^{12}}^*$ and apply the explicit isomorphism $\iota : \mathbb{F}_p[w]/(w^{12}-2w^6+2) \xrightarrow{\sim} \mathbb{F}_p[y,X]/(y^6-2y^3+2, X^2-uy)$ of Section 2.1. Since $P \ne \mathcal{O}_{G_1}$ and $Q_0$ generates $G_2$, $g$ generates $G_T$ and there is a unique $a \in [0, r)$ with $g^a = z$.
2. **Phase 1 (STNFS Tower & Factor Base Enumeration):** Construct the STNFS tower over $K_h = \mathbb{Q}[Y]/(Y^6 - 2Y^3 + 2)$ ($\eta = 6, \kappa = 2$), enumerate the factor base $\mathcal{B} = \mathcal{B}_f \cup \mathcal{B}_g$ of degree-$1$ prime-ideal orbits up to norm bound $B = 2^{65.316}$, and initialize the $n_{\text{Sch}} = 16$ Schirokauer unit maps modulo $r$.
3. **Phase 2 (Relation Collection — 12D Bucket/Line Sieving & Cofactorization):** Sieve over degree-$5$ polynomial pairs $(a(y), b(y)) \in \mathbb{Z}_y^2$ with $\|a\|_\infty, \|b\|_\infty \le A = 1024$ and $\text{lc}(b) \ge 0$, factor smooth norms into prime ideals of $\mathcal{B}_f \cup \mathcal{B}_g$, fold conjugate prime ideals under the Galois involution $\tau : X \mapsto -X$ via $\log(\tau(\mathfrak{q})) \equiv -\log(\mathfrak{q}) \pmod r$, and evaluate the Schirokauer maps modulo $r$.
4. **Phase 3 (Unfiltered $\text{flt} = 1$ Block-Wiedemann Linear Algebra Modulo $r$):** Select $N_{\text{mat}} = \text{SP}_{\text{orb}} = \lceil\text{Li}(B)\rceil + 16$ distinct raw smooth relation rows directly (with zero merge filtering and zero row combinations, preserving raw row sparsity), and solve for a non-zero right kernel vector modulo $r$ via Kaltofen's Block-Wiedemann algorithm ($n=2, m=4$, charging Krylov `smvp`, dense Schirokauer columns, quadratic Matrix Berlekamp–Massey generator synthesis, and `Mksol` solution reconstruction in full).
5. **Phase 4 (Individual Discrete Logarithm — LLL Smoothing at $B_{\text{init}} = 2^{210}$ & 5-Level Special-$\mathfrak{q}$ Descent):** Apply 24-dimensional weighted LLL smoothing ($B_{\text{init}} = 2^{210}$, accounting for the joint probability $\rho(u_{\text{init}})^2 \approx 2^{-82.26}$ across $2^{88}$ trials) and a 5-level special-$\mathfrak{q}$ descent tree in $\mathbb{F}_{p^{12}}^*$ for $g$ and $z$ down to $B = 2^{65.316}$, recovering $a \equiv \log(z)\cdot\log(g)^{-1} \pmod r$.
6. **Phase 5 (Output Construction & Verification):** Compute $Q = [a]Q_0 \in G_2$, verify $e(P, Q) = z$ by one full pairing evaluation, and output $Q$.

---

## 2. Provably Valid STNFS Tower over $h(Y) = Y^6 - 2Y^3 + 2$ and Exact Norms

### 2.1 Base Field $K_h$, Irreducibility Modulo $p$, and Explicit Isomorphism to `field_Fp12`
In `target-profile.json`, the BLS12-381 target field is defined as:
$$\mathbb{F}_{p^{12}} = \mathbb{F}_p[w] / (w^{12} - 2w^6 + 2), \qquad \text{with } i = w^6 - 1 \text{ satisfying } i^2 = -1.$$
Setting $Y = w^2$ yields the degree-$6$ tower polynomial:
$$h(Y) = Y^6 - 2Y^3 + 2 = (Y^3 - 1)^2 + 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 2.$$
Let $y = Y \bmod h(Y)$ and $K_h = \mathbb{Q}(y)$.
1. **Irreducibility over $\mathbb{Q}$ and over $\mathbb{F}_p$:**
   - Over $\mathbb{Q}$, $h(Y) = Y^6 - 2Y^3 + 2$ is $2$-Eisenstein ($\ell = 2$ divides $-2$ and $2$, and $2^2 \nmid 2$), hence irreducible over $\mathbb{Q}$ and $\mathbb{Q}_2$.
   - Over $\mathbb{F}_p$, because $h(w^2) = w^{12} - 2w^6 + 2$ is irreducible of degree $12$ over $\mathbb{F}_p$, $h(Y)$ cannot have any non-trivial factor over $\mathbb{F}_p$; thus $h(Y)$ is **provably irreducible of degree $6$ over $\mathbb{F}_p$**, and $K = \mathbb{F}_p[y]/(h(y)) \cong \mathbb{F}_{p^6}$. (This is also verified directly on the 381-bit prime $p$ via Rabin's test in `experiment.py`.)
2. **Proof that $g_y(X) = X^2 - u y$ is Irreducible over $\mathbb{F}_{p^6}$:**
   In $K = \mathbb{F}_{p^6} = \mathbb{F}_p[y]/(y^6 - 2y^3 + 2)$, since $\deg h = 6$ is even and $h(0) = 2$, the norm of $y$ from $\mathbb{F}_{p^6}$ to $\mathbb{F}_p$ is:
   $$\text{Norm}_{K/\mathbb{F}_p}(y) = (-1)^6 h(0) = 2, \qquad \text{Norm}_{K/\mathbb{F}_p}(u y) = 2 u^6.$$
   The BLS12-381 prime $p$ ends in hex bytes `...aaab` ($\text{0xab} = 171 = 8 \times 21 + 3$), so $p \equiv 3 \pmod 8$. By the Second Supplement to Quadratic Reciprocity, $2$ is a quadratic non-residue modulo $p$: $2^{(p-1)/2} \equiv (-1)^{(p^2-1)/8} = -1 \pmod p$. Since $u \not\equiv 0 \pmod p$, Fermat's Little Theorem gives $(u^3)^{p-1} \equiv 1 \pmod p$, and therefore in $\mathbb{F}_{p^6}$:
   $$(u y)^{(p^6 - 1)/2} = \big(\text{Norm}_{K/\mathbb{F}_p}(u y)\big)^{(p - 1)/2} = (2 u^6)^{(p - 1)/2} = 2^{(p - 1)/2} (u^3)^{p - 1} \equiv -1 \cdot 1 = -1 \pmod p.$$
   Consequently:
   - $u y$ is **provably a quadratic non-residue** in $\mathbb{F}_{p^6}$,
   - $g_y(X) = X^2 - u y$ is **provably irreducible of degree $\kappa = 2$ over $\mathbb{F}_{p^6}$**, defining $\mathbb{F}_{p^6}[X]/(X^2 - u y) \cong \mathbb{F}_{p^{12}}$,
   - The $p^6$-Frobenius automorphism of $\mathbb{F}_{p^{12}} / \mathbb{F}_{p^6}$ acts on $X$ as $X^{p^6} = X(u y)^{(p^6-1)/2} = -X = \tau(X)$.
3. **Explicit Isomorphism with `target-profile.json`'s $\mathbb{F}_p[w]/(w^{12} - 2w^6 + 2)$:**
   In $\mathbb{F}_{p^6} = \mathbb{F}_p[y]/(y^6 - 2y^3 + 2)$, the element $i = y^3 - 1$ satisfies $i^2 = y^6 - 2y^3 + 1 = -1$. Because $p \equiv 3 \pmod 4$, for $u \in \mathbb{F}_p^*$, either $u^{(p-1)/2} \equiv 1 \pmod p$ (in which case $s = u^{(p+1)/4} \in \mathbb{F}_p$ satisfies $s^2 = u$) or $u^{(p-1)/2} \equiv -1 \pmod p$ (in which case $s = (-u)^{(p+1)/4}(y^3 - 1) \in \mathbb{F}_{p^6}$ satisfies $s^2 = (-u)(-1) = u$). In fact, as verified in `experiment.py`, $u^{(p-1)/2} \equiv -1 \pmod p$, so $s(y) = (-u)^{(p+1)/4}(y^3 - 1)$ satisfies $s(y)^2 = u$, and the map:
   $$\iota : \mathbb{F}_p[w]/(w^{12} - 2w^6 + 2) \xrightarrow{\sim} \mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2,\; X^2 - u y), \qquad w^2 \mapsto y,\quad w \mapsto s(y)^{-1} X$$
   is an explicit $O(1)$-cost $\mathbb{F}_p$-algebra isomorphism.
4. **Top Polynomial $f_y(X)$ and Exact Automorphism Reduction ($\mathcal{A} = 2$):**
   Since $3p = P_6(u) = u^6 - 2u^5 + 2u^3 + u + 1$ and $y^6 = 2y^3 - 2$ in $\mathbb{Z}_y = \mathbb{Z}[y]/(y^6 - 2y^3 + 2)$, substituting $u = X^2/y$ into $y^6 P_6(u)$ gives:
   $$f_y(X) = y^6 P_6\!\left(\frac{X^2}{y}\right) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + 2 y^3 - 2 \in \mathbb{Z}_y[X],$$
   which has $\deg_X f_y = 12$, $\|f_y\|_\infty = 2$, and satisfies $f_y(X) \equiv (2y^3 - 2)\cdot 3p \equiv 0 \pmod{p, h(y), g_y(X)}$.
   Because $X^{p^6} = -X = \tau(X)$ and $r \mid (p^4 - p^2 + 1) \mid (p^6 + 1) \implies p^6 \equiv -1 \pmod r$, every degree-$1$ prime ideal $\mathfrak{q} \in \mathcal{B}_f \cup \mathcal{B}_g$ and its conjugate $\tau(\mathfrak{q})$ satisfy:
   $$\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r,$$
   reducing the number of unknown factor-base virtual logarithms by an exact factor of $\mathcal{A} = 2$ to $\text{SP}_{\text{orb}} = \lceil\text{Li}(B)\rceil + n_{\text{Sch}}$.

### 2.2 Exact Integer Resultant Norms via $K_2 = \mathbb{Q}(i)$ and Complex Conjugate Pairs
Because $h(Y) = g_2(Y^3)$ with $g_2(Z) = Z^2 - 2Z + 2 = (Z - 1)^2 + 1$, $K_h = \mathbb{Q}(y)$ contains the quadratic subfield $K_2 = \mathbb{Q}[z]/(z^2 - 2z + 2) \cong \mathbb{Q}(i)$ where $z = y^3 = 1 + i$.
- **Exact Integer Norm Formula:** Any $\gamma(y) = \sum_{j=0}^5 \gamma_j y^j \in \mathbb{Z}_y$ writes as $\gamma(y) = A_0(z) + A_1(z) y + A_2(z) y^2$ with $A_0(z) = \gamma_0 + \gamma_3 z$, $A_1(z) = \gamma_1 + \gamma_4 z$, $A_2(z) = \gamma_2 + \gamma_5 z$. Its relative norm to $K_2$ is the $3 \times 3$ cyclic determinant in $K_2[y]/(y^3 - z)$:
  $$N_{K_h/K_2}(\gamma(y)) = A_0^3 + z A_1^3 + z^2 A_2^3 - 3 z A_0 A_1 A_2 \equiv c_0 + c_1 z \pmod{z^2 - 2z + 2},$$
  and its exact integer norm to $\mathbb{Q}$ is the sum of two squares in $\mathbb{Z}$:
  $$N_{K_h/\mathbb{Q}}(\gamma(y)) = N_{K_2/\mathbb{Q}}(c_0 + c_1 z) = c_0^2 + 2 c_0 c_1 + 2 c_1^2 = (c_0 + c_1)^2 + c_1^2 \in \mathbb{Z}_{\ge 0}.$$
- **Three Distinct Complex Conjugate Pairs:** Since $z_1 = 1 + i = \sqrt{2}\,e^{i\pi/4}$ and $\bar{z}_1 = 1 - i = \sqrt{2}\,e^{-i\pi/4}$, the three cube roots of $z_1$ are:
  $$\theta_k = 2^{1/6} \exp\!\left(i\,\frac{\pi/4 + 2\pi k}{3}\right), \qquad k \in \{0, 1, 2\},$$
  with arguments $\pi/12 = 15^\circ$, $9\pi/12 = 135^\circ$, and $17\pi/12 = 255^\circ$. Their complex conjugates $\bar{\theta}_k$ are the three cube roots of $\bar{z}_1$ (with arguments $-15^\circ, -135^\circ, +105^\circ$). Thus $\{\theta_0, \theta_1, \theta_2\}$ are representatives of the three distinct conjugate pairs of $h(Y) = 0$, and $\prod_{k=0}^2 |\gamma(\theta_k)|^2 = N_{K_h/\mathbb{Q}}(\gamma(y))$ holds identically. In particular, for $a(y) = 1 + y$ and $b(y) = 0$, $N_{K_h/\mathbb{Q}}(1+y) = h(-1) = 1 + 2 + 2 = 5$, so $N_g(1+y, 0) = 5^2 = 25$ and $N_f(1+y, 0) = 5^{12} = 244{,}140{,}625$, verified by both methods in `experiment.py`.

---

## 3. Computation of $1/\zeta_{K_h}(2)$, $\alpha(f_y, h, 10^3)$, $\alpha(g_y, h, 10^3)$, and Simultaneous Smoothness Calibration

### 3.1 Prime-Ideal Factorization of $h(Y) = Y^6 - 2Y^3 + 2$ and $1/\zeta_{K_h}(2)$
The discriminant of $h(Y) = Y^6 - 2Y^3 + 2$ is $\text{Disc}(h) = -2^8 \cdot 3^6 = -186{,}624$.
- At $\ell = 2$: $h(Y)$ is $2$-Eisenstein, so $\mathbb{Z}_y$ is maximal at $2$ and $2\mathcal{O}_h = \mathfrak{l}_2^6$ with a single prime ideal $\mathfrak{l}_2 = (2, y) = (y)$ of inertia degree $\mathfrak{f}(\mathfrak{l}_2, 2) = 1$ and norm $N(\mathfrak{l}_2) = 2$ (contributing $1 - 2^{-2} = 3/4 = 0.75$ to $1/\zeta_{K_h}(2)$).
- At $\ell = 3$: $h(Y) \equiv (Y^2 + Y - 1)^3 \pmod 3$, so above $\ell = 3$ there is a single prime ideal $\mathfrak{l}_3$ of inertia degree $\mathfrak{f}(\mathfrak{l}_3, 3) = 2$ and norm $N(\mathfrak{l}_3) = 3^2 = 9$ (contributing $1 - 9^{-2} = 80/81 \approx 0.987654$).
- At all primes $\ell \ge 5$: $\ell \nmid \text{Disc}(h)$, so $\mathbb{Z}_y$ is locally maximal at $\ell$ and $h(Y) \pmod \ell$ is square-free. By distinct-degree factorization $G_d(Y) = \gcd(Y^{\ell^d} - Y, h(Y)) \pmod \ell$, the number of prime ideals of inertia degree $d$ above $\ell$ is $\deg(G_d)/d$.

Evaluating the Euler product over all primes $\ell \le 1000$ in `experiment.py` yields:
$$\frac{1}{\zeta_{K_h}(2)} = \prod_{\ell \le 1000} \prod_{\mathfrak{l} \mid \ell} \left(1 - \frac{1}{N(\mathfrak{l})^2}\right) \approx 0.7208,$$
which matches the empirical fraction of ideal-coprime pairs $(a(y), b(y))$ observed in `experiment.py`.

### 3.2 TNFS Murphy $\alpha$-Values and Simultaneous Smoothness Calibration
For all unramified primes $\ell \ge 5$ ($\ell \nmid \text{Disc}_X(f_y)\cdot\text{Disc}_X(g_y)$), $v_{\mathfrak{l}}(F_y) = \frac{n_{\mathfrak{l}}^{\text{aff}} q}{q^2 - 1}$ is computed by exact finite-field root counting over $\mathbb{F}_\ell[y]/(h_i(y))$, while at the bad primes ($\ell \in \{2, 3, 7, 13\}$) $v_\ell(F_y)$ is evaluated via exact local root lifting / empirical resultant valuation, giving:
$$\alpha(f_y, h, 10^3) \approx 0.89, \qquad \alpha(g_y, h, 10^3) \approx 0.24 \quad (\text{in base } e).$$
Notice that because $\mathfrak{l}_2 = (y)$ has norm $N(\mathfrak{l}_2) = 2$ and $r = 0$ is a root of both $f_y(X) \equiv X^{12} \pmod y$ and $g_y(X) \equiv X^2 \pmod y$, both $N_f(a,b)$ and $N_g(a,b)$ enjoy higher $2$-adic valuations than under $Y^6 - Y^2 + 1$, lowering $\alpha(f_y, h, 10^3)$ from $2.18$ to $0.89$ and $\alpha(g_y, h, 10^3)$ from $0.96$ to $0.24$.

**Simultaneous $(N_f, N_g)$ Smoothness Calibration in `experiment.py`:**
To validate the Murphy-$E$ / Dickman joint smoothness model against actual simultaneous smooth relations, `observe_simultaneous_smooth_relations_and_calibration()` in `experiment.py` tests $1{,}500$ ideal-coprime pairs $(a(y), b(y))$ on a scaled BLS12 instance ($u_{\text{scaled}} = -6$, $B_{\text{scaled}} = 300$), factors **both** exact integer norms $N_f(a,b)$ and $N_g(a,b)$ by trial division up to $B_{\text{scaled}} = 300$, counts the pairs where **both $N_f(a,b)$ and $N_g(a,b)$ are simultaneously $B_{\text{scaled}}$-smooth**, and compares the observed joint count against the Murphy-$E$ Dickman integral prediction over the same pairs, confirming that the empirical-to-modeled joint yield ratio lies in $[0.75, 1.35]$.

---

## 4. Relation Yield at $(A = 1024, B = 2^{65.316})$, Unfiltered Linear Algebra ($\text{flt} = 1$), and Descent

### 4.1 Expected Yield (Algorithm 6.1) with $\approx 9.8\times$ Relation Oversampling
To eliminate any merge-filtering assumption ($\text{flt} = 1$) while keeping sieving and linear algebra balanced, we set $A = 1024$ and $B = 2^{65.316}$:
- Unreduced sieving volume: $V_0 = \frac{(2A+1)^{12}}{2} = \frac{2049^{12}}{2} \approx 2.751880 \times 10^{39} \approx 2^{131.008451}$.
- $\mathcal{A}=2$ orbit sieving volume ($\text{lc}(b) > 0, \text{lc}(a) \ge 0$): $V_{0,\text{orb}} = \frac{V_0}{2} \approx 1.375940 \times 10^{39} \approx 2^{130.008451}$.
- Valid ideal-coprime orbit space (at $1/\zeta_{K_h}(2) \approx 0.7208$):
  $$V_{1,\text{orb}} = V_{0,\text{orb}} \cdot \frac{1}{\zeta_{K_h}(2)} \approx 9.918 \times 10^{38} \approx 2^{129.536}.$$
- Across $N = 4{,}000$ ideal-coprime samples in `experiment.py` (using the strictly positive integral-average Dickman $\rho$ solver), the mean norm bit-lengths are $\mathbb{E}[\log_2 N_f] \approx 773.5\text{ bits}$ and $\mathbb{E}[\log_2 N_g] \approx 506.0\text{ bits}$, and the mean joint $B$-smoothness probability is $\bar{p} \approx 2^{-66.38} \approx 1.03 \times 10^{-20}$.
- The expected number of distinct ideal-coprime orbit relations is:
  $$\lambda_{\text{orb}} = V_{1,\text{orb}} \cdot \bar{p} \approx 1.02 \times 10^{19} \approx 2^{63.16},$$
  compared against the full unmerged $\mathcal{A}=2$ factor-base orbit threshold:
  $$\text{SP}_{\text{orb}} = \lceil\text{Li}(2^{65.316})\rceil + 16 \approx 1.037897 \times 10^{18} \approx 2^{59.8648}.$$
- **Massive Oversampling Margin ($\lambda_{\text{orb}} / \text{SP}_{\text{orb}} \approx 9.8\times$):**
  Because $\lambda_{\text{orb}}$ exceeds $\text{SP}_{\text{orb}}$ by nearly an order of magnitude ($+3.29\text{ bits}$), even if the true mean smoothness probability $\bar{p}_{\text{true}}$ is systematically **5 times lower ($80\%$ downward modeling error, $\bar{p}_{\text{true}} = 0.20\,\bar{p}$)** than Murphy's Dickman estimate, the degraded expected yield is still $\lambda_{\text{degraded}} = 0.20\,\lambda_{\text{orb}} \approx 2.04 \times 10^{18} > 1.96 \cdot \text{SP}_{\text{orb}}$, with Poisson–Normal z-score:
  $$Z_{\text{degraded}} = \frac{\text{SP}_{\text{orb}} - \lambda_{\text{degraded}}}{\sqrt{\lambda_{\text{degraded}}}} \approx \frac{-1.00 \times 10^{18}}{1.43 \times 10^9} \approx -7.0 \times 10^8 \ll -10,$$
  giving $P(\text{RF} \ge \text{SP}_{\text{orb}}) = 1 - \Phi(Z_{\text{degraded}}) > 0.9999 > 0.40$.

### 4.2 Unfiltered ($\text{flt} = 1$) Sparse Matrix Construction and Useful Rank
We do **not** assume any $20$-fold or $4$-fold merge-filtering reduction ($\text{flt} = 1$), and we do **not** take linear combinations of rows:
1. **Raw Sparse Rows & Deterministic Weight Bound:**
   Let $K_{\text{active}} \le \text{SP}_{\text{orb}} \approx 1.0379 \times 10^{18}$ be the number of factor-base prime-ideal orbits (plus the $16$ Schirokauer unit columns) that appear with non-zero valuation among the collected relations. Because $\text{RF} \ge 2\,\text{SP}_{\text{orb}} \ge 2 K_{\text{active}}$, we select $N_{\text{mat}} = K_{\text{active}} \le \text{SP}_{\text{orb}}$ distinct raw relation rows directly (discarding the remaining excess rows). Every row of the resulting $N_{\text{mat}} \times N_{\text{mat}}$ matrix $M$ is a single unmerged relation $(a(y), b(y))$. Since $N_f(a,b) < 2^{800}$ and $N_g(a,b) < 2^{520}$ and every prime ideal has norm $\ge 2$, the sum of prime-ideal multiplicities in any row is **deterministically bounded by $800 + 520 = 1320$**, while the mean number of distinct prime-ideal factors per row is $\approx \ln\ln N_f + \ln\ln N_g + 16 \approx 28.2$. Charging $\text{wt} = 128$ non-zero entries per row across all $N_{\text{mat}}$ rows covers the unmerged matrix with a $>4.5\times$ row-weight safety factor.
2. **Useful Rank and Block-Wiedemann Kernel Recovery:**
   By standard index-calculus sparse-matrix rank heuristics (Pomerance–Smith 1992; Kaltofen 1995 [4]), selecting $K_{\text{active}}$ independent random smooth relations over the $K_{\text{active}}$ active columns (where each relation is generated by an independent coprime polynomial pair $(a(y), b(y))$ and the true virtual-logarithm vector $v^* \pmod r$ lies in the right kernel of every valid relation row by construction of the commutative TNFS diagram) produces a system whose right kernel recovers $v^*$ on the active factor base. In Kaltofen's Block-Wiedemann algorithm [4, Theorem 7] with blocking parameters $n = 2, m = 4$ over $\mathbb{F}_r$ ($r > 2^{254}$), randomization is performed via the left and right projection blocks $X \in \mathbb{F}_r^{n \times N_{\text{mat}}}, Y \in \mathbb{F}_r^{N_{\text{mat}} \times m}$ without altering the sparsity of $M$, failing to recover a kernel vector with probability at most $2 N_{\text{mat}} / r < 2^{-193}$.

### 4.3 Parameter-Level Derivation of Phase 4 (Initial Smoothing with Joint Probability & 5-Level Descent)
We derive the exact parameters and work of Phase 4 for $\tau \in \{g, z\} \subset G_T \subset \mathbb{F}_{p^{12}}^*$:
1. **Step 4a — 24-Dimensional LLL Smoothing at $B_{\text{init}} = 2^{210}$ (True Joint Probability $\rho(u_{\text{init}})^2$):**
   For a random exponent $e \in [1, r-1]$, set $\tilde{\tau} = \tau^e \in \mathbb{F}_{p^{12}}^*$. In the 24-dimensional lattice $\Lambda_{\tilde{\tau}} = \{(a(y)+b(y)X,\; c(y)+d(y)X) \in \mathbb{Z}_y^2 : (a+bX) \equiv \tilde{\tau}(c+dX) \pmod{p, h(y), X^2-uy}\}$ of determinant $p^{12}$, weighting $b, d$ by $|u|^{1/2} \approx 2^{31.855}$ and applying LLL reduction produces short vectors with $\|a\|_\infty, \|c\|_\infty \approx p^{1/2}|u|^{1/4} \approx 2^{206}$ and $\|b\|_\infty, \|d\|_\infty \approx p^{1/2}|u|^{-1/4} \approx 2^{174}$.
   Thus both $a(y)^2 - u y b(y)^2$ and $c(y)^2 - u y d(y)^2$ have coefficients of size $\approx 2^{413}$, and their degree-$6$ norms $N_g(a,b)$ and $N_g(c,d)$ in $K_h$ have bit-length:
   $$\log_2 N_g(a,b),\; \log_2 N_g(c,d) \le 6 \times 414 = 2484\text{ bits}.$$
   We set the initial smoothness bound to:
   $$B_{\text{init}} = 2^{210}.$$
   Each $2484$-bit norm has Dickman parameter $u_{\text{init}} = \frac{2484}{210} = 11.8286$, giving single-norm smoothness probability $\rho(11.8286) \approx 2^{-41.13}$. Because **both** the numerator norm $N_g(a,b)$ and the denominator norm $N_g(c,d)$ must be simultaneously $B_{\text{init}}$-smooth, the **joint smoothness probability** is:
   $$P_{\text{joint}} = \rho(11.8286)^2 \approx \left(2^{-41.13}\right)^2 = 2^{-82.26}.$$
   Testing $N_{\text{init}} = 2^{88}$ LLL-reduced linear combinations per target gives an expected yield of:
   $$\lambda_{\text{init}} = N_{\text{init}} \cdot P_{\text{joint}} = 2^{88} \cdot 2^{-82.26} = 2^{5.74} \approx 53.4 \text{ simultaneously } B_{\text{init}}\text{-smooth pairs per target},$$
   with failure probability $\exp(-53.4) < 10^{-23}$. Charging $2^{32}$ scalar AArch64 instructions per candidate trial (for LLL linear combination norm evaluation and early-abort trial division / ECM smoothness testing), Step 4a for both $g$ and $z$ costs at most:
   $$R_{\text{descent, 4a}} = 2 \times 2^{88} \times 2^{32} = 2^{121}\text{ scalar AArch64 instructions}.$$
2. **Step 4b — 5-Level Special-$\mathfrak{q}$ Descent Tree from $B_{\text{init}} = 2^{210}$ to $B = 2^{65.316}$:**
   Each prime ideal $\mathfrak{q}$ dividing $N_g(a,b)$ or $N_g(c,d)$ with norm $Q = N(\mathfrak{q}) \in (2^{65.316}, 2^{210}]$ is eliminated by a 5-level special-$\mathfrak{q}$ lattice sieve with exponent contraction factor $0.78$ per level:
   - Level 0 ($Q_0 \le 2^{210}$) $\to$ Level 1 ($Q_1 \le 2^{163.8}$) $\to$ Level 2 ($Q_2 \le 2^{127.8}$) $\to$ Level 3 ($Q_3 \le 2^{99.7}$) $\to$ Level 4 ($Q_4 \le 2^{77.7}$) $\to$ Level 5 ($Q_5 \le 2^{60.6} < B = 2^{65.316}$).
   - At each node, sieving a sublattice $L_{\mathfrak{q}} \subset \mathbb{Z}^{12}$ of volume $V_{\text{desc}} = 2^{60}$ produces residual norms with $\log_2 N_f' \le 6\log_2 Q + 360$ and $\log_2 N_g' \le 448$, for which $u_f + u_g \le 19.2$ and $\rho(u_f)\rho(u_g) \ge 2^{-52}$, yielding on average $\ge 2^8 = 256$ smooth splits per node (allowing us to select a split whose factor-base ideals all lie in the active set $K_{\text{active}}$ and have branching $\le 6$).
   - Across at most $38 \times (1 + 6 + 6^2 + 6^3 + 6^4) = 59{,}090 < 2^{16}$ nodes at $\le 2^{68}$ instructions per node, Step 4b costs $\le 2^{85}$ instructions.
   - Total Phase 4 work is therefore strictly bounded by $R_{\text{descent}} = 2^{122}$ scalar AArch64 instructions ($\ll R_{\text{sieve}} \approx 2^{138.3}$).

---

## 5. Complete Scalar AArch64 Instruction Ledger and Memory Accounting

Per `COST_MODEL.md` and `classical-group-work-v2.json`, total computation is $W = H + R / C$ where $C = 20,000$ scalar AArch64 instructions per reference dense $\mathbb{F}_{p^{12}}$ multiplication ($H = 0$, all operations charged in $R$).

### 5.1 Phase 2: 12D Bucket/Line Sieving at $A = 1024, B = 2^{65.316}$
Over $V_0 = \frac{1}{2}(2A+1)^{12} = \frac{1}{2}(2049)^{12} \approx 2.751880 \times 10^{39} \approx 2^{131.008451}$, by Mertens' Second Theorem the total number of prime-ideal hits across both sides ($\mathcal{B}_f \cup \mathcal{B}_g$) over $V_{0,\text{orb}} = V_0 / 2$ is $N_{\text{hits}} = V_0(\ln\ln B + 0.2615) \approx 4.0742 \cdot V_0$. Charging $32$ scalar AArch64 instructions per bucket/line hit and $16$ instructions per location for array initialization and 64-bit word-masked threshold scanning gives:
$$R_{\text{sieve, derived}} = \big(32(\ln\ln B + 0.2615) + 16\big) V_0 \approx 146.375 \cdot V_0 \approx 4.0281 \times 10^{41} \approx 2^{138.214021}\text{ instructions}.$$
We allocate the explicit upper bound:
$$R_{\text{sieve}} = 4.4 \times 10^{41} \approx 2^{138.341398}\text{ scalar AArch64 instructions}.$$

### 5.2 Phase 3: Unfiltered ($\text{flt} = 1$) Block-Wiedemann Linear Algebra + Reconstruction
For the unmerged sparse matrix of dimension $N_{\text{mat}} \le \text{SP}_{\text{orb}} = \lceil\text{Li}(2^{65.316})\rceil + 16 \approx 1.037897 \times 10^{18} \approx 2^{59.864815}$ and charged row weight $\text{wt} = 128$ (including the $16$ Schirokauer unit columns) over $\mathbb{F}_r$ ($4$ limbs of $64$ bits):
- **Krylov Stage ($n = 2, m = 4$):** $1.5 N_{\text{mat}}$ `smvp` calls $\times 128 N_{\text{mat}} = 192 N_{\text{mat}}^2$ multiply-accumulates in $\mathbb{F}_r$.
- **Matrix Berlekamp–Massey / Generator Synthesis:** $(m+n)^3 N_{\text{mat}}^2 = 216 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$.
- **`Mksol` Solution Reconstruction & Verification:** $0.5 N_{\text{mat}}$ `smvp` calls ($64 N_{\text{mat}}^2$) plus $N_{\text{mat}}$ dense vector updates ($32 N_{\text{mat}}^2$) $= 96 N_{\text{mat}}^2$ operations in $\mathbb{F}_r$.
- **Total $\mathbb{F}_r$ Operations:** $(192 + 216 + 96) N_{\text{mat}}^2 = 504 N_{\text{mat}}^2 \approx 5.4292 \times 10^{38}$ operations in $\mathbb{F}_r$.
- Charging $256$ scalar AArch64 instructions per $4$-limb $\mathbb{F}_r$ multiply-accumulate + sparse index/vector load gives:
  $$R_{\text{LA, derived}} = 256 \times 504 \times N_{\text{mat}}^2 \approx 1.3899 \times 10^{41} \approx 2^{136.678610}\text{ instructions}.$$
  We allocate the explicit upper bound:
  $$R_{\text{LA}} = 1.6 \times 10^{41} \approx 2^{136.881966}\text{ scalar AArch64 instructions}.$$

### 5.3 Summary Instruction Ledger and Score Normalization ($C = 20,000$)

| Phase / Component | Scalar AArch64 Instructions ($R_i$) | $\log_2(R_i)$ | Target-Group-Op Equivalents ($R_i / 20{,}000$) |
| --- | ---: | ---: | ---: |
| Phase 0: Input check, $g = e(P, Q_0)$ & $\iota$ | $16,777,216$ | $24.000000$ | $838.8608$ |
| Phase 1: Poly selection & factor base | $75,557,863,725,914,323,419,136$ | $76.000000$ | $3.777893 \times 10^{18}$ |
| Phase 2a: 12D Bucket/Line Sieving ($A=1024$) | $440,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $138.341398$ | $2.200000 \times 10^{37}$ |
| Phase 2b: ECM cofactorization & Schirokauer | $20,282,409,603,651,670,423,947,251,286,016$ | $104.000000$ | $1.014120 \times 10^{27}$ |
| Phase 3: Unfiltered ($\text{flt}=1$) BW, BM & `Mksol` | $160,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $136.881966$ | $8.000000 \times 10^{36}$ |
| Phase 4: Joint LLL smoothing ($2^{88}$) & descent | $5,316,911,983,139,663,491,615,228,241,121,378,304$ | $122.000000$ | $2.658456 \times 10^{32}$ |
| Phase 5: $Q = [a]Q_0$ & $e(P, Q) = z$ check | $16,777,216$ | $24.000000$ | $838.8608$ |
| **Total Attack Work** | **$600,005,316,932,265,672,729,820,664,539,814,689,767,424$** | **$138.788940$** | **$3.000027 \times 10^{37}$** |

Dividing $R_{\text{total}}$ by $C = 20,000$ gives:
$$W = \frac{R_{\text{total}}}{20,000} \approx 3.00002658 \times 10^{37} \approx 2^{124.501228} < 2^{124.510000}.$$
We therefore submit the exact rational upper bound:
$$\text{time\_log2} = 124.510000 = \frac{12451}{100}, \qquad \text{preprocessing\_log2} = 124.510000.$$

### 5.4 Complete Peak Memory Ledger (`memory_log2_bytes = 76.000000`)
We account for all storage simultaneously resident at peak under unfiltered ($\text{flt} = 1$) operation:
1. **Factor Base Table ($\mathcal{B}_f \cup \mathcal{B}_g$):** $2.0758 \times 10^{18}$ prime-ideal records at $32$ bytes each $\le 6.643 \times 10^{19}\text{ bytes} \approx 2^{65.86}\text{ bytes}$.
2. **Collected Raw Relation Storage ($\le 2^{64}$ relations):** Up to $2^{64}$ raw smooth relations stored at $2{,}048 = 2^{11}$ bytes per relation ($12$ `int16` coefficients, up to $112$ prime-ideal orbit indices and `int16` multiplicities, and $16$ precomputed $32$-byte Schirokauer map values modulo $r$):
   $$M_{\text{relations}} = 2^{64} \times 2^{11}\text{ bytes} = 2^{75.00}\text{ bytes}.$$
3. **Unfiltered Matrix & Block-Wiedemann Working Storage ($N_{\text{mat}} \le 2^{60}$):**
   - $N_{\text{mat}} \le 2^{60}$ selected raw relation rows at $2{,}048$ bytes/row ($2^{71.00}$ bytes).
   - $32$ Block-Wiedemann Krylov and `Mksol` vectors of length $N_{\text{mat}} \le 2^{60}$ over $\mathbb{F}_r$ at $32$ bytes/entry ($32 \times 2^{60} \times 32 = 2^{70.00}$ bytes), plus the virtual logarithm output table ($\text{SP}_{\text{orb}} \times 32 \approx 2^{64.86}$ bytes) and sieve bucket buffers ($< 2^{40}$ bytes).
4. **Total Peak Memory:**
   $$M_{\text{total}} \le 2^{65.86} + 2^{75.00} + 2^{71.00} + 2^{70.05} < 2^{75.14}\text{ bytes} < 2^{76.000000}\text{ bytes}.$$
   Nonuniform advice is $0$ bytes ($\text{nonuniform\_advice\_log2\_bytes} = 0$).

---

## 6. Summary of Disclosed Heuristics and References
1. **Murphy $\alpha$/$E$ Smoothness Model & $\approx 9.8\times$ Oversampling:** Validated locally in `experiment.py` via simultaneous $(N_f, N_g)$ smooth relation observation on a scaled BLS12 instance; the $\approx 9.8\times$ oversampling margin ($\lambda_{\text{orb}} / \text{SP}_{\text{orb}} \ge 8.5$) ensures $P(\text{RF} \ge \text{SP}_{\text{orb}}) > 0.40$ even under an $80\%$ ($5\times$) systematic downward error in $\bar{p}$.
2. **Unfiltered ($\text{flt} = 1$) Sparse Matrix Rank:** Avoids all heuristic merge-filtering reductions ($\text{flt} = 1$) and preserves raw row sparsity; relies only on standard index-calculus sparse kernel sufficiency over active columns.
3. **Joint-Smoothness LLL Bootstrapping & 5-Level Special-$\mathfrak{q}$ Descent:** Explicitly accounts for the joint numerator/denominator smoothness probability $\rho(u_{\text{init}})^2 \approx 2^{-82.26}$ at $B_{\text{init}} = 2^{210}$ across $2^{88}$ trials.

### References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. E. Kaltofen, "Analysis of Coppersmith's Block Wiedemann Algorithm for the Parallel Solution of Sparse Linear Systems," *Mathematics of Computation*, 64(210):777–806, 1995.