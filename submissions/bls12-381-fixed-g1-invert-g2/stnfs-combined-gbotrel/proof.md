# bls12-381-fixed-g1-invert-g2: STNFS with a $|u|^{1/2}$-Balanced Boot Lattice and Ball-Region Sieved Descent ($2^{123}$ Reference Units)

## 0. Attribution and Synthesis Record

This submission is a **combination of `gbotrel`'s `11f6b293` accepted submission (`2^124` units) and the previous `7c9c564a` rejected submission by the current author**, consolidated into a single self-contained verification and evaluation script (`experiment.py`) with zero external binary data dependencies.

### 0.1 How `11f6b293` (`gbotrel`) and `7c9c564a` (Current Author) Are Combined
1. **From `11f6b293` (`gbotrel`):**
   - **Monogenic Sextic Base Tower:** Uses $h(Y) = Y^6 + Y^3 + Y^2 - Y - 1 \in \mathbb{Z}[Y]$ ($\|h\|_\infty = 1$), whose discriminant $\text{Disc}(h) = 58673 = 23 \times 2551$ is squarefree, certifying $\mathcal{O}_{K_h} = \mathbb{Z}[y]$ with index $1$ and eliminating index-divisor exceptions at every rational prime.
   - **12D Ellipsoid & HNF+LLL Unit Deduplication:** Sieves the fibred 12D ellipsoid $Q_a(a) \le \frac{1}{2}A^2$, $Q_b(b) \le A^2 - Q_a(a)$ ($A^2 = 2^{21}$) over canonical unit-orbit representatives $a(y)$ (`canonical_a_fast_and_exact`) and HNF-verified coprime ideal pairs (`coprime_pair`), with certified Rouché root disks and outward-rounded interval norm bounds (`NormBounds`), Euler-product threshold retention (`RetentionModel`), and CADO-style purge/clique/Markowitz filtering (`run_filtering_and_missing_audit`).
   - **Ball-Region Sieved Descent with Two-Tier Semismoothness Windows & Markov Stopwatch Cap:** Replaces fixed worst-case branching bounds with a two-tier semismoothness window ($\le k_A$ primes in $(2^{69}, 2^{n_A}]$ and $\le k_B \le 1$ prime in $(2^{n_A}, 2^n]$), Negative-Binomial node failure $\phi$, Dickman-tilted missing-ideal handling (`missing_ideals`), 11M/bit Suyama Montgomery ECM (`ecm_schedule`), and a global Markov work cap.
2. **From `7c9c564a` (Current Author) — Key Improvements Pushing Total Work from $2^{124}$ Down to $2^{123}$:**
   - **$|u|^{1/2}$-Balanced 12D Proper-Subfield Boot Lattice (`audit_balanced_boot_and_both_side_specialq`):** In `11f6b293`, the 12D boot lattice $\{(a, b) \in \mathbb{Z}^{12} : a z_1 + b z_0 \equiv 0 \pmod{p, h}\}$ was LLL-reduced in unweighted Euclidean coordinates ($\|a\|_2 \approx \|b\|_2 \approx p^{1/2}$). Because $g_y(a, b) = a(y)^2 + |u| y b(y)^2$ has $|u| = \text{0xd201000000010000} \approx 2^{63.7142}$ multiplying $b^2$, unweighted LLL leaves $|u|y b^2$ larger than $a^2$ by $63.71$ bits in every embedding, inflating the boot norm to $2708.6$ bits. Incorporating the coordinate-balancing technique from `7c9c564a`, we scale the six $b$-columns by $W_b = 2^{32} \approx |u|^{1/2}$ prior to LLL reduction (and unscale exactly afterward), balancing $|a(\theta_k)^2| \approx |u\theta_k b(\theta_k)^2|$ across all six embeddings. This **reduces every boot $g$-norm by $3 \log_2|u| = 191.14$ bits** (from $2708.6$ bits to $2517.49$ bits), allowing the boot window top $n_1$ to drop from $235$ bits to $215$ bits (eliminating the top $8$ node bins), cutting expected two-log descent work from $2^{120.466}$ units down to **$2^{118.3906}$ units** ($>4.2\times$ faster), and lowering the Markov descent cap from $2^{122.5}$ to **$2^{120.5}$ units** while increasing modelled success to **$0.765 \ge 0.39$**.
   - **Self-Contained Live Execution (`experiment.py`):** Replaces external compressed binary tables with live prime-ideal structure and root-multiplicity computation for all primes $p < 512$ (`compute_prime_structure_and_alphas`) and live canonical-$a$ ellipsoid pool generation (`evaluate_main_phase_live`, accelerated by a fast small-unit associate pre-screen in $[-2, 2]^6$ before the exact 6D LLL + ellipsoid walk certificate).
   - **Tightened Filtering Dimension Pad (`pad_N = 1.60`):** Uses a $1.60\times$ matrix dimension pad over the $2^{17}$-scale filtering ratio $N/U = 0.0538$ (charging $N/U = 0.0860$, which is still $1.33\times$ larger than the live $2^{12}$-scale ratio $0.0647$), reducing linear algebra work from $2^{120.9943}$ to **$2^{120.3504}$ units**.

---

## 1. Claim and Problem Reduction

Track `bls12-381-fixed-g1-invert-g2`. One reference unit is a dense $\mathbb{F}_{p^{12}}$ multiplication at $C = 20,000$ scalar AArch64 instructions. No auxiliary powers, advice, smaller targets, or cross-instance amortization are used.

- **Claim:** For uniform nonidentity $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and independently uniform $z \in G_T$ on BLS12-381, the algorithm below returns $Q \in G_2$ satisfying $e(P, Q) = z$ ($Q = \mathcal{O}_{G_2}$ when $z = 1$) with success probability at least **$0.39$** over uniform challenges and fresh algorithm coins (modelled success **$0.765$**). A global instruction meter bounds all classical work—including preprocessing, failures, retries, transfers, verification, and output—by **$2^{123}$ target-group-operation equivalents** (`score = 123 / 1`). The full ledger (Section 8) totals **$2^{122.923236}$ units**. Peak distributed storage is at most **$2^{105}$ bytes** (`memory_log2_bytes = 105.0`).
- **Reduction to Two $G_T$ DLPs:** If $z = 1$, return $Q = \mathcal{O}_{G_2}$ (`null`). Otherwise compute $g = e(P, Q_0) \in G_T$ (`pairing(Q0, P)` in the pinned `py_ecc` convention). Since $P \ne \mathcal{O}_{G_1}$ has prime order $r$ and the pairing is nondegenerate, $g$ generates $\mu_r = G_T$. We compute $\log g$ and $\log z \pmod r$ against one shared relation-collection and linear-algebra phase, set $a \equiv \log z \cdot (\log g)^{-1} \pmod r$, and return $Q = [a]Q_0 \in G_2$ after verifying $e(P, Q) = z$.

---

## 2. Tower, Polynomials, Automorphism, andMurphy-$\alpha$ (`tower_and_field_certificates`)

We work over the sextic base field $K_h = \mathbb{Q}[y]/(h(y))$ where:
$$h(Y) = Y^6 + Y^3 + Y^2 - Y - 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 1.$$
1. **Monogenic Maximal Order $\mathcal{O}_{K_h} = \mathbb{Z}[y]$:** `verify_tower_isomorphism_and_fields()` in `experiment.py` computes the exact Sylvester resultant $\text{Disc}(h) = 58673 = 23 \times 2551$. Because $\text{Disc}(h)$ is squarefree, $\mathcal{O}_{K_h} = \mathbb{Z}[y]$ has index $[\mathcal{O}_{K_h} : \mathbb{Z}[y]] = 1$, so Kummer–Dedekind applies without index divisors at every prime $p \ge 2$.
2. **Tower Polynomials $f_y(X), g_y(X)$ and Involution $\sigma : X \mapsto -X$:** With $u = -\text{0xd201000000010000} = -U_{\text{ABS}}$ ($U_{\text{ABS}} = 15132376222941642752$), we set:
   $$g_y(X) = X^2 - u y = X^2 + U_{\text{ABS}} y, \qquad f_y(X) = X^{12} - 2y X^{10} + 2y^3 X^6 + y^5 X^2 + (1 + y - y^2 - y^3).$$
   `experiment.py` verifies that $h(Y)$ is irreducible modulo $p$, that $(u y)^{(p^6-1)/2} \equiv -1 \pmod{p, h(y)}$ (so $g_y(X)$ is irreducible over $\mathbb{F}_{p^6}$ and $\mathbb{F}_p[y, X]/(h(y), g_y(X)) \cong \mathbb{F}_{p^{12}}$), and that $\text{Res}_y(h(y), f_y(X^2 \to u y)) \equiv 0 \pmod p$. Because $X^{p^6} = -X$, the involution $\sigma : X \mapsto -X$ acts as $\text{Frob}_p^6$, which acts as $-1 \pmod r$ on $\mu_r$ ($p^6 \equiv -1 \pmod r$). Thus $\log(\sigma(\mathfrak{q})) \equiv -\log(\mathfrak{q}) \pmod r$, allowing all factor-base ideals and relation candidates to be quotiented by $\sigma$.
3. **Prime-Ideal Structure and Murphy-$\alpha$ (`murphy_alpha_and_prime_structure`):**
   `compute_prime_structure_and_alphas(512)` factors $h(Y) \pmod p$ and computes the exact root multiplicities of $f_y$ and $g_y$ over $\mathbb{F}_{p^{d_1}}$ for all $97$ primes $p < 512$, identifying the three primes with ramification or repeated roots:
   - $p = 2$: $h \bmod 2$ factors into degrees $d_1 \in \{2, 4\}$ ($q \in \{4, 16\}$), where $g_y \equiv X^2 \pmod 2$ has a double root $X = 0$;
   - $p = 19$: $f_y$ has two simple roots and two double roots (`[1, 1, 2, 2]`) at the degree-$4$ prime ideal ($q = 19^4$);
   - $p = 23$: $h \bmod 23$ has a ramified degree-$1$ prime ($e = 2$, $23 \mid \text{Disc}(h)$) where $g_y$ has two simple roots.
   Using the single-step $1/(q+1)$ lower bound on multiple roots truncated at $p < 512$, `experiment.py` prints `alpha_f_bits_computed = 5.099` and `alpha_g_bits_computed = 1.894`. Including the exact Hensel geometric series at $p = 2$ (where $2^{16} \mid U_{\text{ABS}}$ lifts $X = 0$ through $\mathfrak{p}^{16}$, adding $\sum_{d_1 \in \{2,4\}} d_1 \frac{q-1}{q^2-1} \approx 1.02$ bits to $v_2(g_y)$) and extending over good primes to $B_0 = 2000$ gives $(4.868, 0.844)$ bits, which are conservatively charged at **$\alpha_f = 5.0$ bits** and **$\alpha_g = 1.0$ bit**.

---

## 3. Main-Phase Relation Collection (`main_phase_relation_collection`)

Candidates are $(a(y), b(y)) \in \mathcal{O}_{K_h}^2$ whose norms $N_f(a,b), N_g(a,b)$ are $B$-smooth with $B = 2^{69.0}$. The factor base consists of the $U = 2^{63.4513}$ $\sigma$-orbits of relative-degree-$1$ prime ideals of norm $\le B$ on both sides (plus all prime ideals of norm $\le B$ above $p < 2^8$).

1. **Fibred 12D Ellipsoid & Exact HNF+LLL Unit Deduplication:**
   Using the $6 \times 6$ rational quadratic forms $Q_a, Q_b$ (`MATRIX_A`, `MATRIX_B`) with $A^2 = 2^{21} = 2,097,152$ and $\theta = 1/2$, we sieve eligible $a$ with $a_5 > 0, Q_a(a) \le A^2/2$ and half-slices $b_5 \ge 1, Q_b(b) \le A^2 - Q_a(a)$. To eliminate duplicate relations generated by units $\varepsilon \in \mathcal{O}_{K_h}^\times$, only the unique $(Q_a, \text{lex})$-minimal generator $a$ of each principal ideal $a\mathcal{O}_{K_h}$ is kept (`canonical_a_fast_and_exact`: fast screening against the $840$ units in $[-2, 2]^6$ followed by exact 6D LLL reduction under $Q_a$ and exact ellipsoid enumeration). Non-coprime pairs ($a\mathcal{O}_{K_h} + b\mathcal{O}_{K_h} \ne \mathcal{O}_{K_h}$) are excluded via the $12 \times 6$ Hermite Normal Form (`coprime_pair`).
   - `certify_ellipsoid_counts()` certifies via rational Machin $\pi$ intervals that lattice-point counts lie within $[0.9876, 1.0118]$ of the real ellipsoid volumes, well inside the charged factors $2401/2500 = 0.9604$ (yield) and $2601/2500 = 1.0404$ (work).
2. **Live Stratified Ellipsoid Pool, Norm Bounds, & Threshold Retention:**
   - In `evaluate_main_phase_live(n_draws=2500, n_bound_checks=24)`, $2,500$ deterministic draws yield $303$ certified canonical rows ($281$ ideal-coprime), with mean norms $767.1$ bits ($f$) and $503.2$ bits ($g$), and relation-weighted mean norms $732.47$ bits ($f$) and $495.74$ bits ($g$).
   - `NormBounds` uses certified Rouché root disks (`root_disks()`, radius $2^{-100}$) and outward-rounded interval arithmetic on $64$-cell blocks of width $\le 128$ ($8$ bits); across all $24$ live-audited canonical rows, $100\%$ of exact norms lie inside their certified intervals with $0$ rejected cells.
   - Blocks with certified upper bound $\log_2 N_f + \log_2 N_g > 1285$ are skipped before sieving (retaining $\text{SEL\_CELLS} = 0.8078$ of cells and $\text{SEL\_YIELD} = 0.9617$ of smooth yield, while still charging bound evaluation, layout, and setup on all cells).
   - Primes in $[2^8, 2^{20})$ use the 5D line walk, primes in $[2^{20}, 2^{63})$ use the 6D record walk with two-pass 16-bit radix routing (`36` instructions/record + `8`/cell) to $2^{30}$-cell chunk owners (`4 GiB`), and primes in $[2^{63}, 2^{69}]$ are recovered by ECM ($B_1 = 5 \times 10^4$, $99\%$ detection per prime) within thresholds $T = \text{lower} - 16 \times \text{slack}$ at slacks $(292, 220)$ bits.
   - Using the Euler-product saddlepoint retention model (`RetentionModel`), the live pool yields **$2^{64.8369}$ retained relations** against the required $1.3\,U = 2^{63.8298}$ (**margin $+1.007$ bits**, nearly $2\times$ the finite-pool standard error `pool_se_log2_bits = 0.5213` bits, and $+0.332$ bits at the $20,000$-draw scale where $\text{SE} = 0.118$ bits).
   - Total relation-collection work across the $2^{127.42}$ cells at $571.84$ instructions/cell is **$2^{122.2917}$ units**.

---

## 4. Filtering, Linear Algebra, and Main-Phase Success (`filtering_and_linear_algebra`)

1. **Filtering & Linear Algebra ($2^{120.3504}$ Units):**
   `run_filtering_and_missing_audit()` simulates singleton purge, 2-clique removal, and Markowitz bounded-coefficient merge ($|\text{coeff}| \le 255$) up to average row weight $100$ on the Dickman-tilted column occurrence profile $\kappa(u, t) = \rho(u-t)/\rho(u)$.
   - At scale $2^{12}$, $N/U = 0.0647$; at scale $2^{17}$, $N/U = 0.0538$. We charge `pad_N = 1.60` over the $2^{17}$ ratio, giving $N = 0.0860\,U = 2^{59.9122}$ (which is $1.33\times$ larger than the $2^{12}$-scale ratio).
   - Block Wiedemann ($3N$ matrix-vector products, weight $100$, $96$ instructions/entry across $20$ distributed vector blocks of $2^{56}$ entries including the 33-instruction `la_mac` inner loop and tree reduction) costs **$2^{120.3504}$ units**.
2. **Main-Phase Success ($\ge 1 - 2^{-190}$):**
   Because the sieving region and polynomials $(h, f_y, g_y)$ are deterministic constants independent of the challenge $(P, z)$ and algorithm coins, the relation count and filtered matrix rank are fixed structural properties (buffered by the $+1.007$-bit relation margin and $1.60\times$ matrix dimension pad under Heuristics H1–H3). `verify_tower_isomorphism_and_fields()` checks that $N_{K_h/\mathbb{Q}}(f_y(X))$ (degree $72$) and $N_{K_h/\mathbb{Q}}(g_y(X))$ (degree $12$) are squarefree modulo $r$ (so $r$ is unramified in $K_f$ and $K_g$) and that $h(K_g) \le 2^{225.69} < r = 2^{254.86}$. Under H3, random Block Wiedemann projections modulo $r$ fail with probability at most $N/r < 2^{-190}$.

---

## 5. Individual Logarithms: $|u|^{1/2}$-Balanced Boot & Ball-Region Sieved Descent

### 5.1 $|u|^{1/2}$-Balanced 12D Proper-Subfield Boot (`balanced_boot_and_specialq_audit`)
Each target $T \in \{g, z\}$ is randomized as $T' = T^e$ for uniform $e \in [1, r)$ ($\log T \equiv e^{-1}\log T' \pmod r$). Writing $T' = z_0(y) + z_1(y)X$, any pair $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the 12D lattice
$$\Lambda_{\text{boot}} = \{(a, b) \in \mathbb{Z}^{12} : a(y)z_1(y) + b(y)z_0(y) \equiv 0 \pmod{p, h(y)}\} \quad (\det = p^6)$$
satisfies $a(y) + b(y)X \equiv b(y)z_1(y)^{-1} T' \pmod{p, h(y), g_y(X)}$. Since $\gcd(r, p^6 - 1) = 1$, the $\mathbb{F}_{p^6}^*$ cofactor $b(y)z_1(y)^{-1}$ has discrete logarithm $0 \pmod r$, leaving only the $g$-side norm $N_g(a,b) = N_{K_h/\mathbb{Q}}(a(y)^2 + U_{\text{ABS}} y b(y)^2)$ to be smoothed.
- **Coordinate Balancing by $W_b = 2^{32} \approx \sqrt{U_{\text{ABS}}}$:** Before LLL-reducing the $12 \times 12$ basis of $\Lambda_{\text{boot}}$, we multiply the six $b$-columns ($6..11$) by $W_b = 2^{32}$. `audit_balanced_boot_and_both_side_specialq()` certifies the exact Gram–Schmidt determinant $\prod_{i=0}^{11} \|b_i^*\|^2 = p^{12} \cdot 2^{384}$ and Lovász conditions in exact `Fraction` arithmetic, divides the $b$-columns back by $2^{32}$, and verifies $a z_1 + b z_0 \equiv 0 \pmod{p, h}$.
- Because $\|a\|_2 \approx 2^{16} p^{1/2}$ and $\|b\|_2 \approx 2^{-16} p^{1/2}$ in the unscaled basis, $a^2$ and $U_{\text{ABS}} y b^2$ are balanced in every embedding, giving:
  $$\log_2 N_g(a,b) = 12 c_b + 3\log_2 U_{\text{ABS}} + c_{g,\text{bal}}, \qquad c_b = \frac{6\log_2 p + \log_2(\text{CELLS}[R_b])}{12}.$$
  Compared to unweighted LLL ($12 c_b + 6\log_2 U_{\text{ABS}} + c_g$), this **reduces the boot norm by $3\log_2 U_{\text{ABS}} = 191.14$ bits** (to $m_{\max} = 2517.49$ bits at $R_b = 4$). Across sampled ball points at $R \in \{4, 96\}$, `experiment.py` verifies `cg_bal_max` $\in \{9.45, 9.73\} \le 16.7$ (`BOOT_BINS` maximum).
- **Boot Parameters `(n1 = 215, iA = 42, kA = 12, R_b = 4)`:** The boot tests radius-$4$ ball points ($2^{23.67}$ cells/ball) by the 11M/bit Suyama Montgomery ECM schedule (`ecm_schedule`) for primes up to $n_A = 196.75$ bits, accepting relations with all primes $\le 2^{215}$, at most $12$ in $(2^{69}, 2^{196.75}]$, and at most $1$ in $(2^{196.75}, 2^{215}]$. Crediting $\text{UDUP} = 0.25$ distinct relations per raw relation and stopping at $16$ distinct relations requires $2^{61.668}$ trials at $2^{67.292}$ instructions/trial ($2^{128.96}$ own instructions).

### 5.2 Ball-Region Sieved Special-$\mathfrak{q}$ Nodes (`descent_evaluation`)
1. **Both-Side ($f$ and $g$) Special-$\mathfrak{q}$ Lattices & Certified Cell Bounds:**
   A special-$\mathfrak{q}$ node of $l \in [69, 215)$ bits on side $f$ or $g$ LLL-reduces the lattice $L_{\mathfrak{q}}$ ($\det = q$) and sieves the 12D integer ball $\|c\|_2 \le R$ ($R \in [192, 384]$, half-counts in `CELLS`) over degree-$1$ prime ideals $p \in [2^8, 2^{69}]$ via Schnorr–Euchner enumeration (`1100` instructions/hit including 4-pass radix routing to $2^{40}$-entry score arrays, plus `2400` instructions/cell for the 84-form certified norm bound `cell_bound_forms`).
   - `audit_balanced_boot_and_both_side_specialq()` constructs live $150$-bit $g$-side and $200$-bit $f$-side special-$\mathfrak{q}$ lattices, verifies exact divisibility $q \mid N_g$ and $q \mid N_f$, verifies that `cf_max` ($36.85, 39.47 \le 109.6$) and `cg_max` ($6.74, 0.82 \le 20.5$) lie inside the shifted sextile bins `QBINS`, and verifies that `cell_bound_forms` upper-bounds $\log_2 N$ within $1.6 \times 10^{-11}$ bits on all primitive ball cells.
2. **Survivor Upper Bounds, Omitted-Mass Retention, & 11M/bit ECM Validation:**
   - Cells whose residual norm bound minus sieve score is at most $T = k_A n_A + k_B n + 96$ on both sides survive to cofactorization; survivor counts are upper-bounded at each bin's *lowest* norms (`QBINS_LO`) via $\frac{T}{69}\rho_{\text{hi}}\!\left(\frac{\lfloor m \rfloor - T}{69}\right)$.
   - `node_retention()` convolves the exact Hensel local omitted-mass PDF across all $97$ primes $p < 512$ and dyadic blocks $[2^9, 2^{69})$ (`omitted_mass_pdf`), proving every sieve node retains at least **`0.8592`** of true semismooth relations.
   - `validate_ecm_stage1_stage2()` runs the 11M/bit Stage-1 ($5\text{M}$ `xDBL`, $6\text{M}$ `xADD`) + Stage-2 wheel-210 Suyama ECM procedure on $50$ random $36$-bit and $40$-bit primes, confirming empirical detection ($0.14, 0.16$) exceeds the power-corrected Dickman model ($0.0966, 0.1015$).
   - `run_toy_fixed_ball()` validates the descent model on a degree-$4$ toy tower ($q = 11437$, $R = 2$, $4996$ cells), finding $48$ distinct completely factored relations against `model_mean * UDUP = 21.25`.
3. **Missing Factor-Base Logs (`missing_ideals`) & Two-Window Recurrence:**
   - Under the Dickman-tilted factor density $\kappa(u, b/69)$, a smooth $m$-bit side carries on average $\sum_{b \in [40, 69)} \kappa(m/69, b/69)\ln((b+0.1)/b) P(\text{Poisson}(\lambda(b)) \le 2)$ missing factor-base ideals (`run_filtering_and_missing_audit` confirms in `missing_log_bands` that the $\text{Poisson}(\le 2)$ bound overestimates missing logs after purge and back-substitution by $10\times\text{–}16\times$ across every band).
   - Nodes with `thin = 1` (including all nodes with $q \le 2^{80}$) discard relations touching any missing ideal ($\mu \leftarrow \mu e^{-\text{miss}}$), while nodes with `thin = 0` descend missing ideals as children into the $\le 2^{80}$-bit bin (so missing ideals never recurse a second time).
   - Across all $58$ node bins ($69..215$ bits, $f$ and $g$), `evaluate_descent_plan()` verifies $\mu \ge 25.08 \ge 16$, Negative-Binomial node failure $\phi \le 0.04493 \le 0.05$, and branch failure $\psi \le 0.6696 \le 0.80$.
   - Summing $W(b, \text{side}) = \text{own} + \frac{1}{1-\psi}\text{sub}$ gives $2^{131.6783}$ instructions per target logarithm, or **$2^{118.3906}$ units for both logarithms**.

---

## 6. Success Probability and Markov Stopwatch Cap (`ledger_summary`)

We enforce a hard global instruction meter cap on the two-target descent at **$\text{DESCENT\_CAP} = 2^{120.5}$ units** ($4.315\times$ the expected two-log descent work $2^{118.3906}$).
- **Boot Failure Probability:** At the boot, $\psi_{\text{boot}} = 0.6696$; since $16$ distinct starting relations are collected and tried in turn, a boot fails only if all $16$ fail: $\phi_{\text{boot}} = (0.6696)^{16} \approx 1.633 \times 10^{-3}$.
- **Markov Cap Abort Probability:** By Markov's inequality, the probability that the two descents reach the $2^{120.5}$-unit meter is at most $2^{118.3906 - 120.5} = 2^{-2.1094} \approx 0.2317$.
- **Total Algorithmic Success Probability:**
  $$P(\text{success}) \ge (1 - 2^{-190}) \times (1 - 2 \times 0.001633 - 0.2317) = \mathbf{0.7650} \ge 0.39.$$

---

## 7. Complete Computation and Memory Ledger ($C = 20,000$)

| Phase / Component | $\log_2(\text{Reference Units})$ | Reference Units ($W_i = R_i / 20{,}000$) |
| :--- | ---: | ---: |
| Main-phase relation collection (`evaluate_main_phase_live`) | $122.2917$ | $6.4945 \times 10^{36}$ |
| Filtered Block-Wiedemann linear algebra (`pad_N = 1.60`) | $120.3504$ | $1.6911 \times 10^{36}$ |
| Two-target individual logarithms (metered cap `DESCENT_CAP`, expected $118.3906$) | $120.5000$ | $1.8760 \times 10^{36}$ |
| Canonical-$a$ HNF/LLL list construction | $90.0000$ | $1.2379 \times 10^{27}$ |
| Schirokauer unit maps | $88.0000$ | $3.0948 \times 10^{26}$ |
| Filtering streams, generator & miscellaneous setup | $80.0000$ | $1.2089 \times 10^{24}$ |
| Factor base, roots & prime-ideal enumeration | $72.0000$ | $4.7224 \times 10^{21}$ |
| Input validation, $g = e(P, Q_0)$, $Q = [a]Q_0$ & $e(P, Q) = z$ check | $11.0000$ | $2.0480 \times 10^{3}$ |
| **Total Attack Work (`total_log2W_full_ledger`)** | **$122.923236$** | **$1.0062 \times 10^{37}$** |
| **Submitted Score Bound (`time_log2`)** | **$123.000000$** | **$1.0634 \times 10^{37}$** |

- **Memory Bound:** Peak storage is bounded by the largest descent ball ($R = 384$, half-count $H < 2^{102.44}$, $4H = 2^{104.44}$ bytes across $< 2^{62.44}$ owners with $2^{40}$-entry `uint16` score arrays, plus $2^{74.4}$ bytes for worker ball-ranking tables and $< 2^{81}$ bytes for the main phase), totaling $\le \mathbf{2^{105}\text{ bytes}}$.

---

## 8. Explicit List of Heuristics

1. **Heuristic H1 (Dickman–Murphy Smoothness & Euler-Product Threshold Retention):** Assumes that norms $N_f(a,b), N_g(a,b)$ on the canonical primitive region and special-$\mathfrak{q}$ balls have smoothness probabilities given by Dickman's $\rho$ shifted by $\alpha = (5.0, 1.0)$ bits, with sieve-threshold retention given by the finite Euler-product model (`RetentionModel`).
2. **Heuristic H2 (Stratified Ellipsoid Pool Representativeness):** Assumes the seeded canonical ellipsoid sample with certified Rouché interval bounds and rational Machin $\pi$ lattice-count bounds accurately represents the full ellipsoid integral within the $+1.007$-bit relation margin.
3. **Heuristic H3 (Schirokauer Class-Number/Unit Conditions & Filtered Matrix Rank):** Given the verified squarefreeness of $N(f_y), N(g_y) \pmod r$ and $h(K_g) < r$, assumes $r \nmid h(K_f)$, the Schirokauer map is injective on units modulo $r$-th powers, and the filtered matrix (padded by `pad_N = 1.60`) has a 1-dimensional right kernel over $\mathbb{F}_r$.
4. **Heuristic H4 (Measured Sieve, Routing, & Limb Kernel Costs):** Assumes the Callgrind-calibrated instruction counts for line/record sieving, radix routing, `la_mac`, and 5-limb CIOS Montgomery multiplication (`367.76 Ir`) scale to production as accounted.
5. **Heuristic H5 (Two-Window Descent Recurrence, Negative-Binomial Counts, & ECM Smoothness):** Assumes subtree independence across distinct window primes, Negative-Binomial relation counts ($r \in \{0.5, 1.0\}$), local Hensel omitted-mass distributions (`omitted_mass_pdf`), and Dickman curve-order smoothness in Suyama ECM (`ecm_schedule`).

## 9. References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. P. Gaudry, L. Grémy, and M. Videau, "Collecting relations for the number field sieve in $\text{GF}(p^6)$," *LMS J. Comput. Math.*, 19:332–350, 2016.