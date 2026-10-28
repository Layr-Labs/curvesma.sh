# Comments
- Judge submission to test TNFS attack,
- Follows the paper from Guillevic-Singh, but modified parameter choice from Tables 10 and 11 after the LLM-judge found a mistake in that selection,
- Created by Gemini 4.0 Argon, human supervised.

# bls12-381-fixed-g1-invert-g2: Special Tower Number Field Sieve (STNFS) via $G_T$ DLP Reduction (Revision v6)

## 0. Revision Record (Addressing Review of Prior Submission `v5`)

### 0.1 Baseline / Prior Work
In our prior iterations (`57d7de1`, `631ee3e7`, `dc4b09ab`, `v4`, and `v5` = 49050072 submission ID), we formulated a Special Tower Number Field Sieve (STNFS) attack on `bls12-381-fixed-g1-invert-g2` using the degree-$6$ tower polynomial $h(Y) = Y^6 - 2Y^3 + 2$, $g_y(X) = X^2 - u y$, and $f_y(X) = X^{12} - 2yX^{10} + 2y^3X^6 + y^5X^2 + 2y^3 - 2$.

In the review of `v5`, the judge confirmed that the pairing-to-DLP reduction, the STNFS tower irreducibility and non-residuosity argument, the corrected explicit field isomorphism $\iota : w \mapsto s^{-1}X$ ($s = u^{(p+1)/4} \bmod p \in \mathbb{F}_p^*$), the completed experimental run on the algebraic/norm checks, and the conservative minimum of the two smoothness models ($\min(\bar{p}_{\text{BD}+\alpha}, \overline{P_f P_g})$) were all supported. However, `v5` was rejected due to four remaining gaps in Phase 4 (descent/ECM) and H3 (active column coverage):
1. **Special-$\mathfrak{q}$ Lattice Geometry Not Constructed or Bounded:** `v5` scaled the ordinary-box mean norms by $C_0 Q^{1/12}$ solely from covolume $Q$, without bounding the lengths and orthogonality defect of the actual reduced basis matrix $V = [v_1, \dots, v_{12}]$ or constructing and sampling actual special-$\mathfrak{q}$ bases in `experiment.py`.
2. **Descent Tree Branching Bound $b_{\max}$ Used Sample Means Instead of Supremum Bounds:** `audit_descent_tree()` computed $b_{\max}$ by dividing extrapolated *mean* logarithmic residual norms by $65.0$, rather than using a demonstrated *upper bound* (supremum) on the residual norms over the sieve box with an explicit norm cutoff.
3. **Un-derived Per-Candidate Descent Sieve Constant (`64` inst) & Unvalidated ECM Implementation:** Step 4b multiplied $V_{\text{desc}}$ by an unexplained constant `64` instructions/candidate without an analytical breakdown of LLL setup, factor-base root inversion, 12D bucket sieving, survivor scanning, memory transfers, and survivor ECM cofactorization; and Step 4a evaluated ECM cost formulas without validating a working ECM implementation in `experiment.py`.
4. **Core Ideal Coverage in H3:** H3 assumed recovery of *every* core ideal's logarithm up to $2^{65.0}$ without proving that all core columns occur in the unfiltered matrix $M$.

### 0.2 Corrections & Implementation Changes in Revision v6
1. **Live 100-Digit-`Decimal` 12D LLL Reduction & Sampling on Genuine Special-$\mathfrak{q}$ Bases (Fixes Issue 1):**
   - In `experiment.py`, `lll_reduce_12d()` implements exact-integer 12-dimensional LLL lattice reduction ($\delta = 3/4$) using 100-decimal-digit (`332`-bit) Gram–Schmidt orthogonalization (`decimal.Decimal` with `prec = 100`).
   - At each of the 3 descent levels ($Q_k \approx 2^{165}, 2^{122}, 2^{89}$), `find_degree1_special_q()` constructs a genuine degree-$1$ prime ideal $\mathfrak{q}_k = (Q_k, y - r_y, X - r_X)$ with $h(r_y) \equiv 0 \pmod{Q_k}$ and $r_X^2 \equiv u r_y \pmod{Q_k}$, builds the exact $12 \times 12$ HNF basis of $L_{\mathfrak{q}_k} \subset \mathbb{Z}^{12}$, and LLL-reduces it.
   - `construct_and_sample_special_q_basis()` verifies that the LLL-reduced basis vectors satisfy $\max_j \|v_j\|_\infty / Q_k^{1/12} \in [0.9753, 1.0651]$ and measures the exact worst-case coordinate column-$L_1$ sum $M_{\text{col}, k} = \max_{0 \le m < 12} \sum_{j=1}^{12} |v_{j,m}| \in [4.0604\,Q_k^{1/12},\; 5.4414\,Q_k^{1/12}]$. It then samples $40$ random vectors $(a(y), b(y)) = \sum_{j=1}^{12} c_j v_j$ directly from the reduced basis ($c_j \in [-C_{0,k}, C_{0,k}]$), **asserts `ng % Q_k == 0` on 100% of sampled vectors**, and uses the **directly measured residual norms from the sampled LLL vectors** (rather than covolume-scaled ordinary-box formulas) to compute $u_f, u_g$ at each level.
2. **Unconditional Triangle-Inequality Supremum Norm Bounds for $b_{\max, k}$ and Row Weight (Fixes Issue 2):**
   - In Section 2.2 (`supremum_log2_norms(R)`), since all six roots $\theta_k$ of $h(Y) = Y^6 - 2Y^3 + 2$ have modulus $|\theta_k| = 2^{1/6}$, the triangle inequality $|a(\theta_k)|, |b(\theta_k)| \le R\sum_{i=0}^5 2^{i/6} < 2^{3.030} R$ yields the **unconditional worst-case supremum bounds** for any $(a,b)$ with $\|(a,b)\|_\infty \le R$:
     $$\log_2 N_{f,\sup}(R) \le 72 \log_2 R + 238.0\text{ bits}, \qquad \log_2 N_{g,\sup}(R) \le 12 \log_2 R + 420.0\text{ bits}.$$
   - For any $c \in [-C_{0,k}, C_{0,k}]^{12}$ in the reduced basis of $L_{\mathfrak{q}_k}$, every coordinate of $(a,b) = \sum_{j=1}^{12} c_j v_j$ is deterministically bounded by $R_{\sup, k} = C_{0,k} \max_m \sum_{j=1}^{12} |v_{j,m}|$. We enforce an explicit norm cutoff at $\log_2 N_{f,\sup}(R_{\sup, k}) + \log_2 N_{g,\sup}(R_{\sup, k}) - \log_2 Q_k$ ($\in [2186.2, 2317.1]$ bits) and compute the deterministic worst-case branching bound $b_{\max, k}$ by dividing this **supremum** (never sample means) by $\log_2 B_{\text{core}} = 65.0$, yielding $b_{\max, 0} = 35, b_{\max, 1} = 33, b_{\max, 2} = 33$ and a deterministic total of $38 + 38 \times 35 + 38 \times 35 \times 33 = 45{,}258 \approx 2^{15.4659}$ nodes per target.
3. **Working Suyama Projective Montgomery ECM Implementation & 4-Part Per-Node Instruction Account (Fixes Issue 3):**
   - `validate_suyama_montgomery_ecm()` implements a working projective $(X:Z)$ Montgomery ladder ECM (`xdbl_proj`, `xadd_proj`, `montgomery_ladder`) with Suyama's $\sigma$-parametrization, factors $N = 10007 \times 10009 = 100160063$ into $10007$, and verifies that each scalar bit executes $8$ modular multiplications (`measured_muls_per_scalar_bit = 7.96`).
   - In Step 4b, we replace the flat `64` instructions/candidate constant with an explicit 4-part analytical per-node instruction bound:
     $$I_{\text{node}}(k) = \underbrace{2^{70}}_{\text{LLL & FB root setup}} + \underbrace{296\,V_{\text{desc}}(C_{0,k})}_{\text{12D bucket sieve, scan & transfers}} + \underbrace{2^{118}}_{\text{39-limb ECM on } 2^{52}\text{ survivors/node}},$$
     giving $2^{118.000}$ instructions/node at Levels 0 and 1 and $2^{134.947}$ instructions/node at Level 2, for a total Step 4b cost of $2 \times \sum_{k=0}^2 N_{\text{nodes}, k} I_{\text{node}}(k) \approx 2^{151.3682}$ scalar AArch64 instructions.
4. **Active Column Filtering & Coverage Probability $P_{\text{all\_active}}$ at the Final Descent Level (Fixes Issue 4):**
   - We retain $R_{\text{rows}} = 2.0\,\text{SP}_{\text{orb}} \approx 9.133 \times 10^{18} \approx 2^{62.983842}$ relations in $M$ (out of $\lambda_{\text{BD}+\alpha} \approx 2^{63.445120}$ collected relations, oversampling ratio $1.3768\times$, $P(\text{RF} \ge R_{\text{rows}}) = 1.0$).
   - Because Phase 3 completes before Phase 4 begins, the active column set $\mathcal{B}_{\text{active}} \subset \mathcal{B}$ of prime-ideal orbits occurring at least twice in $M$ is explicitly known. At the final descent level (Level 2), any candidate split touching an ideal outside $\mathcal{B}_{\text{active}}$ is rejected.
   - `compute_active_column_coverage_prob()` integrates the Poisson $\ge 2$-hit probability $p_{\text{hit}\ge 2}(\ell) = 1 - (1 + R_{\text{rows}}/\ell)\exp(-R_{\text{rows}}/\ell)$ weighted by the Mertens prime-factor density $d\ell / (\ell \ln \ell)$ across $\ell \le B_{\text{core}} = 2^{65.0}$ and raises it to the deterministic supremum factor count $b_{\max, 2} = 33$, yielding $P_{\text{all\_active}} = 0.7179$. Multiplying by $P_{\text{all\_active}}$ still gives $35.0 \ge 20.0$ **fully active** smooth splits per node at Level 2.

### 0.3 Summary of Verified Experimental Outputs (`experiment.py`)
Running `experiment.py` with the default initial inputs (`alpha_prime_bound = 1000`, `mc_samples = 3000`, `seed = 42`, `A = 1360`, `B_log2_milli = 67500`) produces:
- `tower_isomorphism_and_subfield_verified`: `true` ($Y^{p^6} \equiv Y \pmod{p, h(Y)}$, $(u y)^{(p^6-1)/2} \equiv -1$, $s = u^{(p+1)/4} \bmod p \implies s^2 \equiv u \pmod p$, $W = s^{-1}X \implies W^{12}-2W^6+2 = 0$, $\iota(z_1 z_2) = \iota(z_1)\iota(z_2)$, and $\gcd(r, p^6 - 1) = 1$)
- `suyama_montgomery_ecm_validation`: `test_composite_N = 100160063`, `recovered_factor = 10007`, `measured_muls_per_scalar_bit = 7.96`
- `inv_zeta_Kh_2_euler`: `0.671179` (`empirical_ideal_coprime_ratio`: `0.6676`)
- `alpha_f_base_e`: `1.1706`, `alpha_g_base_e`: `0.6739`
- `deterministic_row_weight_proof`: `supremum_log2_Nf_at_A = 987.48 < 1197.22` ($145$ ideals), `supremum_log2_Ng_at_A = 544.91 < 765.96` ($100$ ideals), `derived_schirokauer_columns = 40`, `deterministic_max_nonzero_entries_per_row = 285`
- `barbulescu_duquesne_alpha_yield`: `mean_log2_Nf = 830.1716`, `mean_log2_Ng = 515.0101`, `log2_p_BD_alpha = -70.8987`, `log2_p_direct_joint_mean = -68.4322`, `empirical_cov_Pf_Pg_positive = true`, `log2_V0 = 135.919055`, `log2_V1_orb = 134.343824`, `log2_SP_orb = 61.983842`, `log2_required_rows_2p0_SP_orb = 62.983842`, `log2_expected_relations_lambda_BD_alpha = 63.44512`, `oversampling_ratio_lambda_over_req_rows = 1.3768` (`2.7535` over $\text{SP}_{\text{orb}}$), `prob_RF_ge_req_rows = 1.0`
- `descent_tree_audit`: `step4b_c0_schedule = [64, 160, 800]`, `step4b_deterministic_max_total_nodes = 45258` ($\log_2 = 15.4659$), `step4b_total_instructions_log2 = 151.3682`, `active_column_survival_prob = 0.7179`, `expected_active_splits_per_node = 35.0` at Level 2
- `cost_and_memory_ledger`: `derived_sieve_instructions_log2 = 143.122918` (charged `143.22142`), `derived_la_MtDM_instructions_log2 = 141.584232` (charged `141.690905`), `derived_descent_instructions_log2 = 151.368219` (charged `151.569853`), `R_total_log2 = 151.575799`, `W_target_group_ops_log2 = 137.288087`, `submitted_time_log2 = 137.3` (`1373 / 10`), `peak_memory_bytes_log2 = 77.091391` (submitted `78.0`).

---

## 1. Problem Reduction and Complete Attack Architecture

The challenge supplies a uniformly sampled nonidentity point $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and an independently uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^*$ on BLS12-381, where:
- $u = -\text{0xd201000000010000} = -15132376222941642752$,
- $p = \frac{(u-1)^2(u^4 - u^2 + 1)}{3} + u = \frac{u^6 - 2u^5 + 2u^3 + u + 1}{3}$ ($381$ bits, $\log_2 p \approx 380.0640$, $p \equiv 3 \pmod 8$),
- $r = u^4 - u^2 + 1$ ($255$-bit prime subgroup order, $\log_2 r \approx 254.8571$),
- $Q_0 \in G_2$ is the fixed public generator of $G_2$.

We execute a 6-phase cold-start STNFS attack:
1. **Phase 0 (Input Validation, Pairing & Isomorphism $\iota$):** Verify $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and $z \in G_T$. If $z = 1$, return $Q = \mathcal{O}_{G_2}$ (`null`). Otherwise compute $g = e(P, Q_0) \in G_T \subset \mathbb{F}_{p^{12}}^*$ and map $g, z$ via the verified $\mathbb{F}_p$-algebra isomorphism $\iota$ (Section 2.1) to $\mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2, X^2 - u y)$.
2. **Phase 1 (Tower & Factor Base Setup):** Enumerate the factor base $\mathcal{B} = \mathcal{B}_f \cup \mathcal{B}_g$ of degree-$1$ prime-ideal orbits up to norm $B = 2^{67.500}$ and initialize the $n_{\text{Sch}} = 40$ Schirokauer unit maps modulo $r$.
3. **Phase 2 (12D Bucket/Line Sieving & Cofactorization):** Sieve over $(a(y), b(y)) \in \mathbb{Z}_y^2$ with $\|a\|_\infty, \|b\|_\infty \le A = 1360$ and $\text{lc}(b) > 0, \text{lc}(a) \ge 0$, fold conjugate prime ideals via $\log(\tau(\mathfrak{q})) \equiv -\log(\mathfrak{q}) \pmod r$, and evaluate the $40$ Schirokauer maps on smooth relations.
4. **Phase 3 (Unfiltered $\text{flt} = 1$ Rectangular $M^T D M$ Block-Wiedemann Linear Algebra):** Retain $R_{\text{rows}} = 2.0\,C_{\text{cols}}$ distinct raw relations over the $C_{\text{cols}} \le \text{SP}_{\text{orb}}$ active columns without row merges, record the active column bitmap $\mathcal{B}_{\text{active}} \subset \mathcal{B}$ of columns occurring $\ge 2$ times in $M$, choose a random diagonal $D \in (\mathbb{F}_r^*)^{R_{\text{rows}}}$, and apply Block-Wiedemann ($n=2, m=4$) to $A = M^T D M$ to recover the virtual logarithms of all ideals in $\mathcal{B}_{\text{active}}$.
5. **Phase 4 (Single-Norm Proper-Subfield LLL Smoothing & 3-Level Special-$\mathfrak{q}$ Descent to $\mathcal{B}_{\text{active}}$):** Using $\gcd(r, p^6-1) = 1$, smooth $g^e$ and $z^{e'}$ into a single $2^{165}$-smooth $g_y$-norm via 12-dimensional LLL reduction across $2^{69}$ fresh-$e$ trials (charging concrete 39-limb CIOS Montgomery ECM on every trial), and descend all prime factors over $3$ side-aware special-$\mathfrak{q}$ levels (`165 -> 122 -> 89 -> 65`) into $\mathcal{B}_{\text{active}}$ to recover $a \equiv \log(z)\cdot\log(g)^{-1} \pmod r$.
6. **Phase 5 (Output Construction & Verification):** Compute $Q = [a]Q_0 \in G_2$, verify $e(P, Q) = z$, and return $Q$.

---

## 2. Algebraic Tower over $h(Y) = Y^6 - 2Y^3 + 2$, Verified Field Isomorphism $\iota$, and Supremum Norm Bounds

### 2.1 Irreducibility, Non-Residuosity of $u y$, and Machine-Checked Isomorphism $\iota$
From `target-profile.json`, $\mathbb{F}_{p^{12}} = \mathbb{F}_p[w]/(w^{12} - 2w^6 + 2)$. Setting $Y = w^2$ gives:
$$h(Y) = Y^6 - 2Y^3 + 2 = (Y^3 - 1)^2 + 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 2.$$
1. **Irreducibility over $\mathbb{Q}$ and $\mathbb{F}_p$:** $h(Y)$ is $2$-Eisenstein over $\mathbb{Q}$, and because $h(w^2) = w^{12} - 2w^6 + 2$ is irreducible of degree $12$ over $\mathbb{F}_p$, $h(Y)$ is irreducible of degree $6$ over $\mathbb{F}_p$, defining $K = \mathbb{F}_p[y]/(h(y)) \cong \mathbb{F}_{p^6}$.
2. **Non-Residuosity of $u y$ in $\mathbb{F}_{p^6}$:** Since $\deg h = 6$ is even and $h(0) = 2$, $\text{Norm}_{K/\mathbb{F}_p}(u y) = 2 u^6$. Because $p \equiv 3 \pmod 8$, $2^{(p-1)/2} \equiv -1 \pmod p$, so:
   $$(u y)^{(p^6 - 1)/2} = (2 u^6)^{(p - 1)/2} = 2^{(p - 1)/2}(u^3)^{p - 1} \equiv -1 \pmod p.$$
   Thus $g_y(X) = X^2 - u y$ is irreducible of degree $\kappa = 2$ over $\mathbb{F}_{p^6}$, and $X^{p^6} = X(u y)^{(p^6-1)/2} = -X = \tau(X)$.
3. **Explicit, Verified Field Isomorphism $\iota$:**
   Write $u = -2^{16} t$ where $t = \text{0xd20100000001} \equiv 5 \pmod{12}$. Since $3p = P_6(u) \equiv 1 \pmod t$, Quadratic Reciprocity gives $\left(\frac{t}{p}\right) = \left(\frac{p}{t}\right) = \left(\frac{3}{t}\right) = -1$. Because $p \equiv 3 \pmod 4$, $\left(\frac{-1}{p}\right) = -1$, and hence $\left(\frac{u}{p}\right) = (-1)(+1)(-1) = +1$. Thus $s = u^{(p+1)/4} \bmod p \in \mathbb{F}_p^*$ satisfies $s^2 \equiv u \pmod p$, and in $\mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2,\; X^2 - u y)$ the element $W = s^{-1} X$ satisfies $W^2 = s^{-2} X^2 = u^{-1}(u y) = y$ and $W^{12} - 2W^6 + 2 = y^6 - 2y^3 + 2 = 0$. The explicit $\mathbb{F}_p$-linear map
   $$\iota\!\left(\sum_{j=0}^{11} c_j w^j\right) = \left(\sum_{k=0}^5 c_{2k} y^k\right) + \left(\sum_{k=0}^5 c_{2k+1} s^{-1} y^k\right) X$$
   is an exact $\mathbb{F}_p$-algebra isomorphism (machine-checked in `experiment.py`).
4. **Top Polynomial $f_y(X)$ and Involution $\tau : X \mapsto -X$:** Substituting $u = X^2/y$ into $y^6 P_6(u) = y^6(u^6 - 2u^5 + 2u^3 + u + 1)$ and using $y^6 = 2y^3 - 2$ gives:
   $$f_y(X) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + 2 y^3 - 2 \in \mathbb{Z}_y[X],$$
   which satisfies $f_y(X) \equiv (2y^3 - 2)\cdot 3p \equiv 0 \pmod{p, h(y), g_y(X)}$. Since $p^6 \equiv -1 \pmod r$, every conjugate pair of degree-$1$ prime ideals $\{\mathfrak{q}, \tau(\mathfrak{q})\}$ satisfies $\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r$, halving the number of unknown factor-base virtual logarithms to $\text{Li}(B)$.

### 2.2 Exact Integer Resultant Norms and Unconditional Triangle-Inequality Supremum Bounds
Writing $z = y^3$ ($z^2 - 2z + 2 = 0$), any $\gamma(y) = A_0(z) + A_1(z)y + A_2(z)y^2 \in \mathbb{Z}_y$ has relative norm $N_{K_h/K_2}(\gamma) = A_0^3 + z A_1^3 + z^2 A_2^3 - 3 z A_0 A_1 A_2 \equiv c_0 + c_1 z \pmod{z^2 - 2z + 2}$, and exact rational integer norm $N_{K_h/\mathbb{Q}}(\gamma(y)) = (c_0 + c_1)^2 + c_1^2 \in \mathbb{Z}_{\ge 0}$, which equals $\prod_{k=0}^2 |\gamma(\theta_k)|^2$ over the three cube roots $\theta_k = 2^{1/6}\exp(i(\pi/4 + 2\pi k)/3)$ of $z_1 = 1 + i$.

**Unconditional Supremum Bounds (`supremum_log2_norms(R)`):**
Because every root $\theta_k$ of $h(Y) = 0$ has exact modulus $|\theta_k| = 2^{1/6}$, for **any** pair $(a(y), b(y)) \in \mathbb{Z}_y^2$ with $\|a\|_\infty, \|b\|_\infty \le R$, the triangle inequality gives $|a(\theta_k)|, |b(\theta_k)| \le R\sum_{i=0}^5 2^{i/6} = \frac{1}{2^{1/6}-1} R \approx 8.1659\,R < 2^{3.030} R$. Consequently:
- $|a(\theta_k)^2 - u\theta_k b(\theta_k)^2| \le (8.1659)^2(1 + 2^{1/6}|u|) R^2 < 2^{70.00} R^2 \implies \log_2 N_g(a,b) \le 6(70.00 + 2\log_2 R) = 12\log_2 R + 420.0\text{ bits}$.
- $|F_{\theta_k}(a,b)| \le (1 + 2\cdot 2^{1/6} + 2\cdot 2^{3/6} + 2^{5/6} + 2)(8.1659\,R)^{12} < 2^{39.66} R^{12} \implies \log_2 N_f(a,b) \le 6(39.66 + 12\log_2 R) < 72\log_2 R + 238.0\text{ bits}$.
These bounds hold **deterministically and unconditionally** for every point in $[-R, R]^{12}$.

---

## 3. Exact Computation of $1/\zeta_{K_h}(2)$, $\alpha(f_y, h, 10^3)$, $\alpha(g_y, h, 10^3)$, and Smoothness Yield

### 3.1 Exact Local Valuations Across All Prime Ideals ($\ell \le 1000$)
As executed in `compute_zeta_and_alphas(1000)` in `experiment.py`:
- At $\ell = 2$ ($\mathfrak{l}_2 = (y)$, $q = 2$): Algorithm 5.2 (`val_l2_algorithm_5_2()`) gives the exact valuations $v_2(f_y) = 2$ and $v_2(g_y) = \frac{4}{3}(1 - 2^{-48})$.
- At $\ell = 3$ ($\mathfrak{l}_3 = (3, y^2+y-1)$, $q = 9$): $g_y(X) = X^2 - u y$ has $0$ roots in $\mathbb{F}_9$, so $v_3(g_y) = 0$ exactly; for $f_y(X)$, setting $v_3(f_y) = 0$ gives the unconditional upper bound $\alpha_3(f_y) \le \frac{\ln 3}{2} \approx 0.5493$.
- At all primes $\ell \ge 5$: `factor_h_mod_l_exact(l)` factors $h(Y) \pmod \ell$ into its exact irreducible factors $h_i(Y)$ of degree $d \in \{1, 2, 3, 6\}$ and counts the exact roots in $\mathbb{F}_{\ell^d}$ (with Hensel lifting at the bad primes $\ell \in \{7, 13\}$), giving:
  $$\frac{1}{\zeta_{K_h}(2)} = 0.671179, \qquad \alpha(f_y, h, 10^3) = 1.1706, \qquad \alpha(g_y, h, 10^3) = 0.6739 \quad (\text{in base } e).$$

### 3.2 Smoothness Yield and Relation Surplus at $(A = 1360,\; B = 2^{67.500})$
At $A = 1360$ ($2A + 1 = 2721$) and $B = 2^{67.500}$:
- Unreduced sieving volume: $V_0 = \frac{2721^{12}}{2} \approx 2^{135.919055}$; $\mathcal{A}=2$ orbit volume: $V_{0,\text{orb}} = \frac{V_0}{2} \approx 2^{134.919055}$; ideal-coprime orbit volume: $V_{1,\text{orb}} = V_{0,\text{orb}} / \zeta_{K_h}(2) \approx 2^{134.343824}$.
- Across $N = 3{,}000$ ideal-coprime pairs in `experiment.py` (`empirical_ideal_coprime_ratio = 0.6676`), $\overline{\log_2 N_f} = 830.1716\text{ bits}$ and $\overline{\log_2 N_g} = 515.0101\text{ bits}$.
- Under Heuristic H1, `evaluate_bd_alpha_yield()` evaluates both the Barbulescu–Duquesne + $\alpha$ mean-norm estimate $\bar{p}_{\text{BD}+\alpha} = P_f(\bar{u}_f)P_g(\bar{u}_g) \approx 2^{-70.8987}$ ($\bar{u}_f = 12.3239, \bar{u}_g = 7.6442$) and the direct joint sample average $\overline{P_f P_g} \approx 2^{-68.4322}$ ($\widehat{\text{Cov}}(P_f, P_g) > 0$). Using the lower value $\bar{p}_{\text{cons}} = 2^{-70.8987}$ gives:
  $$\lambda_{\text{BD}+\alpha} = V_{1,\text{orb}} \cdot \bar{p}_{\text{cons}} \approx 2^{63.445120} \approx 1.2574 \times 10^{19}.$$
- Against $\text{SP}_{\text{orb}} = \text{Li}(2^{67.500}) + 40 \approx 4.5665 \times 10^{18} \approx 2^{61.983842}$ and the retained rectangular row count $R_{\text{rows}} = 2.0\,\text{SP}_{\text{orb}} \approx 9.1330 \times 10^{18} \approx 2^{62.983842}$, the oversampling ratio is $\lambda_{\text{BD}+\alpha} / R_{\text{rows}} = 1.3768$ ($2.7535\times$ over $\text{SP}_{\text{orb}}$), giving $P(\text{RF} \ge R_{\text{rows}}) = 1.000000 > 0.40$ under Heuristic H2.

---

## 4. Derived Schirokauer Columns, Supremum Row Weight, Rectangular $M^T D M$, and Live-LLL Descent

### 4.1 Exact Derivation of $n_{\text{Sch}} = 40$ Schirokauer Unit Columns and Supremum Row Weight ($\le 285$)
- **Schirokauer Columns ($n_{\text{Sch}} = 40$):** Because $K_h = \mathbb{Q}[y]/((y^3-1)^2+1)$ is totally complex, $K_f$ ($[K_f:\mathbb{Q}] = 72$) and $K_g$ ($[K_g:\mathbb{Q}] = 12$) are totally complex with $r_{\mathbb{R}}(K_f) = r_{\mathbb{R}}(K_g) = 0$, $r_{\mathbb{C}}(K_f) = 36 \implies r_f = 35$, and $r_{\mathbb{C}}(K_g) = 6 \implies r_g = 5$. By Dirichlet's Unit Theorem, $n_{\text{Sch}} = r_f + r_g = 35 + 5 = 40$.
- **Supremum Row Weight ($\le 285$):** At $R = A = 1360$, the unconditional triangle-inequality suprema (Section 2.2) are $\log_2 N_{f,\sup}(1360) = 987.48\text{ bits}$ and $\log_2 N_{g,\sup}(1360) = 544.91\text{ bits}$. In `verify_deterministic_row_weight_bound()`, the product of the $145$ smallest degree-$1$ prime-ideal norms in $\mathcal{O}_h$ is $2^{1197.22} > 2^{987.48}$ and the product of the $100$ smallest is $2^{765.96} > 2^{544.91}$. Thus **no relation row can ever have more than $145 + 100 + 40 = 285$ non-zero entries**.

### 4.2 Rectangular $M^T D M$ Block-Wiedemann Over $R_{\text{rows}} = 2.0\,C_{\text{cols}}$ Raw Relations
Retaining $R_{\text{rows}} = 2.0\,C_{\text{cols}}$ raw relations ($C_{\text{cols}} \le \text{SP}_{\text{orb}} \approx 4.5665 \times 10^{18}$), each row of $M$ and $M^T$ costs at most $245 \times 32 + 40 \times 256 = 18{,}080$ instructions (plus $256$ instructions for scaling by $D$). One multiplication by $A = M^T D M$ costs $2.0\,C_{\text{cols}}(2 \times 18{,}080 + 256) = 72{,}832\,C_{\text{cols}}$ instructions. Across $2.0\,C_{\text{cols}}$ steps ($n=2, m=4$) plus $216 \times 256\,C_{\text{cols}}^2 = 55{,}296\,C_{\text{cols}}^2$ instructions for Matrix Berlekamp–Massey:
$$R_{\text{LA, derived}} = (2.0 \times 72{,}832 + 55{,}296) C_{\text{cols}}^2 = 200{,}960 \, C_{\text{cols}}^2 \approx 4.1907 \times 10^{42} \approx 2^{141.584232}\text{ instructions}.$$

### 4.3 Concrete 39-Limb Montgomery ECM Smoothing (Step 4a) and Live-LLL 3-Level Descent (Step 4b)
1. **Step 4a — Single-Norm Proper-Subfield LLL Smoothing & Validated Suyama Montgomery ECM:**
   Since $\gcd(r, p^6 - 1) = 1$, every element of $\mathbb{F}_{p^6}^*$ has discrete logarithm $0 \pmod r$. For fresh $e \in_R [1, r-1]$, 12-dimensional LLL reduction of $\Lambda_{\tau^e}^{(12)} \subset \mathbb{Z}^{12}$ ($\det = p^6$, weighting $b$ by $|u|^{1/2}$) yields a single polynomial $a(y) + b(y)X$ with $\log_2 N_g(a,b) \le 2484\text{ bits}$ ($L = 39$ limbs of $64$ bits).
   - At $B_{\text{init}} = 2^{165}$ ($u_{\text{init}} = 15.0545$), $\rho(15.0545) \approx 2^{-63.8423}$, so $N_{\text{init}} = 2^{69}$ fresh-$e$ trials yield $\lambda_{\text{init}} \approx 35.70$ $B_{\text{init}}$-smooth elements per target.
   - **Validated Suyama Montgomery ECM:** In `validate_suyama_montgomery_ecm()`, our projective $(X:Z)$ Montgomery ladder ECM (`xdbl_proj` = $4\text{M}$, `xadd_proj` = $4\text{M}$, total $8\text{M}$ per scalar bit, measured $7.96\text{M}$/bit) factors $N = 100160063$ into $10007$. On $L = 39$ limbs, one CIOS multiplication takes $16 L^2 + 32 L = 25{,}584$ instructions, so $1$ ladder bit takes $8 \times 25{,}584 = 204{,}672$ instructions. At $B_1 = 2^{28}$ ($\le 1.4427\,B_1$ bits) and $N_{\text{curves}} = 2^{19}$ curves per trial ($u_{\text{ecm}} = 5.7648$, miss union bound $4.26 \times 10^{-8} < 10^{-6}$ across $\le 45$ factors), each trial takes $2^{65.1717}$ instructions, totaling $2 \times 2^{69} \times 2^{65.1717} = 2^{135.1717}$ instructions for $g$ and $z$.
2. **Step 4b — Live 100-Digit-`Decimal` 12D LLL Special-$\mathfrak{q}$ Bases, Supremum $b_{\max}$, and Active Filtering:**
   At each level $k \in \{0, 1, 2\}$ of the schedule `$[165, 122, 89, 65]$` with box radii $C_{0,k} \in [64, 160, 800]$, `construct_and_sample_special_q_basis()` constructs a genuine degree-$1$ prime ideal $\mathfrak{q}_k = (Q_k, y - r_y, X - r_X)$, LLL-reduces the $12 \times 12$ HNF basis of $L_{\mathfrak{q}_k}$ in 100-digit `Decimal` arithmetic, verifies $Q_k \mid N_g(a,b)$ on $40$ sampled vectors $(a,b) = \sum_{j=1}^{12} c_j v_j$ ($c_j \in [-C_{0,k}, C_{0,k}]$), and measures both the empirical mean norms and the exact worst-case coordinate supremum $R_{\sup, k} = C_{0,k} \max_m \sum_{j=1}^{12} |v_{j,m}|$:
   - **Level 0 ($165 \to 122$, $C_{0,0} = 64$, $\log_2 V_{\text{desc}} = 83.13$, $38$ root nodes):**
     - LLL basis: $\max_j \|v_j\|_\infty / Q_0^{1/12} = 1.0559$, $\max_m \sum_j |v_{j,m}| / Q_0^{1/12} = 4.0604$, $\log_2 R_{\sup, 0} = 21.7151$.
     - Supremum residual norm sum: $\log_2 N_{f,\sup}(R_{\sup,0}) + \log_2 N_{g,\sup}(R_{\sup,0}) - 165 = 2317.1\text{ bits} \implies b_{\max, 0} = \lfloor 2317.1 / 65.0 \rfloor = 35$.
     - Sampled LLL norms: $\overline{\log_2 N_f} = 1511.27$, $\overline{\log_2 N_g'} = 462.38 \implies \log_2 p_{\text{split}}^{(f)} = -52.68,\; \log_2 p_{\text{split}}^{(g)} = -55.12 \implies 2.71 \times 10^8$ expected splits/node.
   - **Level 1 ($122 \to 89$, $C_{0,1} = 160$, $\log_2 V_{\text{desc}} = 98.92$, $38 \times 35 = 1{,}330$ nodes):**
     - LLL basis: $\max_j \|v_j\|_\infty / Q_1^{1/12} = 0.9753$, $\max_m \sum_j |v_{j,m}| / Q_1^{1/12} = 4.6351$, $\log_2 R_{\sup, 1} = 19.6447$.
     - Supremum residual norm sum: $2186.2\text{ bits} \implies b_{\max, 1} = \lfloor 2186.2 / 65.0 \rfloor = 33$.
     - Sampled LLL norms: $\overline{\log_2 N_f} = 1349.91$, $\overline{\log_2 N_g'} = 479.46 \implies \log_2 p_{\text{split}}^{(f)} = -75.36,\; \log_2 p_{\text{split}}^{(g)} = -77.56 \implies 2.68 \times 10^6$ expected splits/node.
   - **Level 2 ($89 \to 65.0$, $C_{0,2} = 800$, $\log_2 V_{\text{desc}} = 126.74$, $1{,}330 \times 33 = 43{,}890$ nodes):**
     - LLL basis: $\max_j \|v_j\|_\infty / Q_2^{1/12} = 1.0651$, $\max_m \sum_j |v_{j,m}| / Q_2^{1/12} = 5.4414$, $\log_2 R_{\sup, 2} = 19.4480$.
     - Supremum residual norm sum: $2202.6\text{ bits} \implies b_{\max, 2} = \lfloor 2202.6 / 65.0 \rfloor = 33$.
     - Sampled LLL norms: $\overline{\log_2 N_f} = 1321.07$, $\overline{\log_2 N_g'} = 507.76 \implies \log_2 p_{\text{split}}^{(f)} = -119.04,\; \log_2 p_{\text{split}}^{(g)} = -121.13$.
     - **Active-Column Survival Probability ($P_{\text{all\_active}} = 0.7179$):** Over $R_{\text{rows}} = 2.0\,\text{SP}_{\text{orb}} \approx 2^{62.9838}$ rows in $M$, integrating $p_{\text{hit}\ge 2}(\ell) = 1 - (1 + R_{\text{rows}}/\ell)\exp(-R_{\text{rows}}/\ell)$ with weight $d\ell/(\ell \ln \ell)$ across $\ell \le 2^{65.0}$ and raising to $b_{\max, 2} = 33$ gives $P_{\text{all\_active}} = 0.7179$, yielding $2^{126.74 - 121.13} \times 0.7179 = 35.0 \ge 20.0$ **fully active** splits per node!
   - **Analytical Per-Node & Total Step 4b Instruction Count:**
     Charging $I_{\text{node}}(k) = 2^{70} (\text{LLL & FB root setup}) + 296\,V_{\text{desc}}(C_{0,k}) (\text{12D bucket sieve, scan & transfers}) + 2^{118} (\text{39-limb ECM on } 2^{52}\text{ survivors/node})$ across all $38 + 1{,}330 + 43{,}890 = 45{,}258 \approx 2^{15.4659}$ nodes per target gives $2^{151.368219}$ scalar AArch64 instructions (charged as $R_{\text{descent}} = 4.237134 \times 10^{45} \approx 2^{151.569853}$ instructions).

---

## 5. Summary Scalar AArch64 Instruction and Memory Ledgers ($C = 20,000$)

| Phase / Component | Scalar AArch64 Instructions ($R_i$) | $\log_2(R_i)$ | Target-Group-Op Equivalents ($R_i / 20{,}000$) |
| --- | ---: | ---: | ---: |
| Phase 0: Input check, $g = e(P, Q_0)$ & $\iota$ | $16,777,216$ | $24.000000$ | $838.8608$ |
| Phase 1: Poly selection & factor base | $302,231,454,903,657,293,676,544$ | $78.000000$ | $1.511157 \times 10^{19}$ |
| Phase 2a: 12D Bucket/Line Sieving ($A=1360$) | $13,000,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $143.221420$ | $6.500000 \times 10^{38}$ |
| Phase 2b: ECM cofactorization & Schirokauer | $81,129,638,414,606,681,695,789,005,144,064$ | $106.000000$ | $4.056482 \times 10^{27}$ |
| Phase 3: Unfiltered $M^T D M$ BW ($2.0\,\text{SP}_{\text{orb}}$) | $4,500,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $141.690905$ | $2.250000 \times 10^{38}$ |
| Phase 4: 39-limb ECM & 45,258-node descent | $4,237,134,137,732,706,293,067,299,657,566,419,742,516,838,400$ | $151.569853$ | $2.118567 \times 10^{41}$ |
| Phase 5: $Q = [a]Q_0$ & $e(P, Q) = z$ check | $16,777,216$ | $24.000000$ | $838.8608$ |
| **Total Attack Work** | **$4,254,634,137,732,787,422,706,016,495,703,196,659,354,976,256$** | **$151.575799$** | **$2.127317 \times 10^{41}$** |

Dividing $R_{\text{total}}$ by $C = 20,000$ gives:
$$W = \frac{R_{\text{total}}}{20,000} \approx 2.12731707 \times 10^{41} \approx 2^{137.288087} < 2^{137.300000}.$$
We submit the exact rational upper bound:
$$\text{time\_log2} = 137.300000 = \frac{1373}{10}, \qquad \text{preprocessing\_log2} = 137.300000, \qquad \text{memory\_log2\_bytes} = 78.000000.$$

---

## 6. Explicit List of Heuristic Assumptions and Hypotheses

1. **Heuristic H1 (Barbulescu–Duquesne + Guillevic–Singh Murphy-$\alpha$ & Direct Joint Smoothness Heuristic):**
   Assumes that across ideal-coprime pairs $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the sieving box, the joint $B$-smoothness probability of $N_f(a,b)$ and $N_g(a,b)$ is modeled by $\bar{p}_{\text{cons}} = \min(\bar{p}_{\text{BD}+\alpha}, \overline{P_f P_g}) \approx 2^{-70.8987}$.
2. **Heuristic H2 (Poisson Relation Yield Distribution):**
   Assumes that simultaneous $B$-smoothness events across the $V_{1,\text{orb}}$ ideal-coprime pairs behave as weakly dependent Bernoulli trials whose count $\text{RF}$ concentrates around $\lambda_{\text{BD}+\alpha}$ via $\text{Poisson}(\lambda_{\text{BD}+\alpha}) \approx \mathcal{N}(\lambda_{\text{BD}+\alpha}, \lambda_{\text{BD}+\alpha})$.
3. **Heuristic H3 (Rectangular $M^T D M$ Sparse Kernel Rank Over Active Columns $\mathcal{B}_{\text{active}}$):**
   Assumes that over $R_{\text{rows}} = 2.0\,\text{SP}_{\text{orb}}$ random smooth relations, the right kernel of $M \pmod r$ (recovered via Kaltofen–Saunders Block-Wiedemann on $M^T D M$ with zero merge filtering and supremum-verified row weight $\le 285$) uniquely determines the virtual logarithms of all active factor-base columns $\mathcal{B}_{\text{active}}$ occurring $\ge 2$ times in $M$ up to global scaling.
4. **Heuristic H4 (Single-Norm Proper-Subfield LLL Smoothing & Live-LLL-Calibrated Special-$\mathfrak{q}$ Descent Splitting):**
   Assumes that the $2484$-bit norm $N_g(a,b)$ produced by 12-dimensional LLL reduction of $\Lambda_{\tau^e}^{(12)}$ for uniform $e \in [1, r-1]$ and the residual norms $N_f', N_g'$ measured on the 100-digit-`Decimal` LLL-reduced special-$\mathfrak{q}$ bases $L_{\mathfrak{q}_k}$ split with probabilities given by Dickman's $\rho$-function, and that prime factors of smooth splits follow the Mertens log-uniform density used in $P_{\text{all\_active}}$.

## 7. References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. E. Kaltofen and B. D. Saunders, "On Wiedemann's Method of Solving Sparse Linear Systems," *AAECC-9*, LNCS 539, pp. 29–38, 1991.