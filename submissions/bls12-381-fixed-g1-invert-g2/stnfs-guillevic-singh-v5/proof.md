# Comments
- Judge submission to test TNFS attack,
- Follows the paper from Guillevic-Singh, but modified parameter choice from Tables 10 and 11 after the LLM-judge found a mistake in that selection,
- Created by Gemini 4.0 Argon, human supervised.

# bls12-381-fixed-g1-invert-g2: Special Tower Number Field Sieve (STNFS) via $G_T$ DLP Reduction (Revision v5)

## 0. Revision Record (Addressing Review of Prior Submission `v4` / `dc4b09ab`)

### 0.1 Baseline / Prior Work
In our prior iterations (`57d7de1`, `631ee3e7`, `dc4b09ab`, and `v4`), we formulated a Special Tower Number Field Sieve (STNFS) attack on `bls12-381-fixed-g1-invert-g2` using the degree-$6$ tower polynomial $h(Y) = Y^6 - 2Y^3 + 2$, $g_y(X) = X^2 - u y$, and $f_y(X) = X^{12} - 2yX^{10} + 2y^3X^6 + y^5X^2 + 2y^3 - 2$.

In the review of `v4`, the judge confirmed that the pairing-to-DLP reduction, the tower irreducibility argument ($h(w^2) = w^{12} - 2w^6 + 2$ and $\text{Norm}(u y) = 2u^6$ with $p \equiv 3 \pmod 8$), and the completed experimental execution on the checked tower and ledger identities were supported, but rejected the submission due to three remaining defects:
1. **Incorrect Sign of $\left(\frac{u}{p}\right)$ in the Explicit Field Isomorphism $\iota$ (`proof.md` §2.1):**
   For $u = -\text{0xd201000000010000} = -2^{16} t$ ($t = \text{0xd20100000001} \equiv 5 \pmod{12}$) and $3p = P_6(u) \equiv 1 \pmod t$, Quadratic Reciprocity gives $\left(\frac{t}{p}\right) = \left(\frac{p}{t}\right) = \left(\frac{3}{t}\right) = -1$. Since $p \equiv 3 \pmod 4$, $\left(\frac{-1}{p}\right) = -1$, so $\left(\frac{u}{p}\right) = (-1)(+1)(-1) = +1$ ($u$ is a **quadratic residue** in $\mathbb{F}_p^*$, not a non-square). Consequently, `v4`'s formula $s = (-u)^{(p+1)/4}(y^3 - 1)$ had $s^2 = -u$, making $W = s^{-1}X$ satisfy $W^2 = -y$ and $W^{12} - 2W^6 + 2 = 4y^3 \ne 0$.
2. **Descent Region Construction, $g$-Side vs. $f$-Side Residual Norms, Node Count, and Concrete ECM Cost (Section 4.4):**
   - Step 4a relied on the leading asymptotic expression $L_{2^{165}}(1/2, \sqrt{2})$ without a concrete multi-limb arithmetic implementation and per-factor failure-budget analysis on $2484$-bit trials.
   - Step 4b (`audit_descent_tree()`) used hard-coded residual norm bounds without tying them to the explicit sieve box radius $C_0$, applied the $f$-side residual formula ($5\log_2 Q + \text{const}$) uniformly instead of tracking $g$-side ideals ($\mathfrak{q} \in \mathcal{B}_g \implies Q \mid N_g$, leaving $N_f$ with $6\log_2 Q + \text{const}$), and asserted $< 2^{18}$ nodes without a deterministic derivation.
3. **1D Jensen Claim on a Product of Dependent Variables:**
   `v4` claimed that 1-variable convexity of $\rho(u)$ implies $\mathbb{E}[P_f P_g] \ge P_f(\bar{u}_f) P_g(\bar{u}_g)$ via separate Jensen inequalities, whereas separate Jensen inequalities bound $\mathbb{E}[P_f]\mathbb{E}[P_g]$, not $\mathbb{E}[P_f P_g] = \mathbb{E}[P_f]\mathbb{E}[P_g] + \text{Cov}(P_f, P_g)$.

### 0.2 Corrections & Implementation Changes in Revision v5
1. **Corrected & Machine-Checked Field Isomorphism $\iota$ with $s = u^{(p+1)/4} \bmod p \in \mathbb{F}_p^*$ (Fixes Issue 1):**
   - Since $\left(\frac{u}{p}\right) = +1$ and $p \equiv 3 \pmod 4$, the scalar $s = u^{(p+1)/4} \bmod p \in \mathbb{F}_p^*$ lies in the base field $\mathbb{F}_p$ and satisfies $s^2 \equiv u^{(p+1)/2} = u \cdot u^{(p-1)/2} \equiv u \pmod p$.
   - Setting $W = s^{-1} X \in \mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2,\; X^2 - u y)$ gives:
     $$W^2 = s^{-2} X^2 = u^{-1}(u y) = y, \qquad W^{12} - 2W^6 + 2 = (W^2)^6 - 2(W^2)^3 + 2 = y^6 - 2y^3 + 2 = 0.$$
   - In `verify_algebraic_and_tower_identities()` in `experiment.py`, we compute `s_scalar = pow(u_mod_p, (P_MOD + 1) // 4, P_MOD)`, verify `(s_scalar * s_scalar) % P_MOD == u_mod_p`, verify $W^2 = y$ and $W^{12} - 2W^6 + 2 = 0$ in the tower, and verify $\iota(z_1 \cdot z_2) = \iota(z_1) \cdot \iota(z_2)$ on elements of $\mathbb{F}_{p^{12}}$.
2. **Concrete 39-Limb CIOS Montgomery ECM Implementation & Side-Aware 3-Level Special-$\mathfrak{q}$ Descent (Fixes Issue 2):**
   - **Step 4a (Concrete 39-Limb Montgomery ECM):** Each $2484$-bit norm $N_g(a,b)$ fits in $L = 39$ limbs of $64$ bits. One Coarsely Integrated Operand Scanning (CIOS) Montgomery multiplication takes $16 L^2 + 32 L = 25{,}584 < 2^{15}$ scalar AArch64 instructions, and one projective Montgomery ladder step (`xDBLADD`) takes $8 \times 25{,}584 = 204{,}672 < 2^{18}$ instructions per scalar bit. With Stage-$1$ bound $B_1 = 2^{28}$ ($\le 1.4427\,B_1$ scalar bits), each curve costs $7.92 \times 10^{13} \approx 2^{46.1717}$ instructions. For any prime factor $q \le B_{\text{init}} = 2^{165}$, Suyama-curve smoothness of $|E(\mathbb{F}_q)|/12$ has $u_{\text{ecm}} = (165 - \log_2 12)/28 = 5.7648$ and success probability $\rho(5.7648) \approx 2^{-14.66}$. Running $N_{\text{curves}} = 2^{19}$ curves per trial gives per-factor miss probability $\exp(-2^{19}\rho(5.7648)) \approx 9.46 \times 10^{-10}$ (union bound over $\le 45$ prime factors: $4.26 \times 10^{-8} < 10^{-6}$) and costs $2^{19 + 46.1717} = 2^{65.1717} < 2^{66}$ instructions per trial ($2 \times 2^{69} \times 2^{65.1717} \approx 2^{135.17}$ instructions total for $g$ and $z$).
   - **Step 4b (Side-Aware 3-Level Special-$\mathfrak{q}$ Descent):** We use a 3-level schedule `$[165.0, 122.0, 89.0, 65.0]$` down to $B_{\text{core}} = 2^{65.0}$ ($2.5$ bits below $B = 2^{67.5}$) with explicit per-level box radii $C_{0,k} \in [64, 128, 550]$ in the 12-dimensional index-$Q$ sublattice $L_{\mathfrak{q}}$. A box $[-C_{0,k}, C_{0,k}]^{12}$ in $L_{\mathfrak{q}}$ (covolume $Q$) has sieve volume $V_{\text{desc}}(C_{0,k}) = \frac{1}{2}(2C_{0,k}+1)^{12}$ and coordinate half-width $C_{0,k} Q^{1/12}$ (scaled by $S = C_{0,k} Q^{1/12} / A$ relative to $A = 1360$). `audit_descent_tree()` computes **both** the $f$-side ($\mathfrak{q} \in \mathcal{B}_f \implies Q \mid N_f$) and $g$-side ($\mathfrak{q} \in \mathcal{B}_g \implies Q \mid N_g$) residual norm bounds at every level, verifies $\min(p_{\text{split}}^{(f)}, p_{\text{split}}^{(g)}) \cdot V_{\text{desc}}(C_{0,k}) \ge 29.4 \ge 20.0$ at every level, and bounds the total tree size deterministically via $b_{\max, k} = \lfloor(\log_2 N_f' + \log_2 N_g')/65.0\rfloor \in \{30, 27\}$ ($33{,}640 \approx 2^{15.0379}$ total nodes across all 3 levels, costing $2^{142.2389}$ instructions).
3. **Both BD+$\alpha$ and Direct Joint Smoothness Evaluation Without 1D-Jensen-on-Product Claim (Fixes Issue 3):**
   In `evaluate_bd_alpha_yield()`, we compute both the Barbulescu–Duquesne + $\alpha$ mean-norm smoothness model $\bar{p}_{\text{BD}+\alpha} = P_f(\bar{u}_f)P_g(\bar{u}_g) \approx 2^{-70.8987}$ and the direct joint sample average $\overline{P_f P_g} \approx 2^{-68.4322}$ (verifying empirically that $\text{Cov}(P_f, P_g) > 0$), and use the smaller value $\min(\bar{p}_{\text{BD}+\alpha}, \overline{P_f P_g}) = 2^{-70.8987}$ as Heuristic H1 without claiming that 1D Jensen alone bounds a product.

### 0.3 Summary of Verified Experimental Outputs (`experiment.py`)
Running `experiment.py` with the default initial inputs (`alpha_prime_bound = 1000`, `mc_samples = 3000`, `seed = 42`, `A = 1360`, `B_log2_milli = 67500`) produces:
- `tower_isomorphism_and_subfield_verified`: `true` ($Y^{p^6} \equiv Y \pmod{p, h(Y)}$, $(u y)^{(p^6-1)/2} \equiv -1 \pmod{p, h(y)}$, $u^{(p-1)/2} \equiv +1 \pmod p$, $s = u^{(p+1)/4} \bmod p \implies s^2 \equiv u \pmod p$, $W = s^{-1}X \implies W^{12}-2W^6+2 = 0$, $\iota(z_1 z_2) = \iota(z_1)\iota(z_2)$, and $\gcd(r, p^6 - 1) = 1$)
- `inv_zeta_Kh_2_euler`: `0.671179` (`empirical_ideal_coprime_ratio`: `0.6676`)
- `alpha_f_base_e`: `1.1706`, `alpha_g_base_e`: `0.6739`
- `deterministic_row_weight_proof`: `max_distinct_f = 145` ($\log_2 \text{prod} = 1197.22$), `max_distinct_g = 100` ($\log_2 \text{prod} = 765.96$), `derived_schirokauer_columns = 40`, `deterministic_max_nonzero_entries_per_row = 285`
- `barbulescu_duquesne_alpha_yield`: `mean_log2_Nf = 830.1716`, `mean_log2_Ng = 515.0101`, `u_f_BD_alpha = 12.3239`, `u_g_BD_alpha = 7.6442`, `log2_p_BD_alpha = -70.8987`, `log2_p_direct_joint_mean = -68.4322`, `empirical_cov_Pf_Pg_positive = true`, `log2_V0 = 135.919055`, `log2_V1_orb = 134.343824`, `log2_SP_orb = 61.983842`, `log2_required_rows_1p25_SP_orb = 62.30577`, `log2_expected_relations_lambda_BD_alpha = 63.44512`, `oversampling_ratio_lambda_over_req_rows = 2.2028` (`2.7535` over $\text{SP}_{\text{orb}}$), `prob_RF_ge_req_rows = 1.0`
- `descent_tree_audit`: `step4a_single_norm_u_init = 15.0545`, `step4a_single_norm_prob_log2 = -63.8423`, `step4a_log2_trials = 69.0`, `step4a_expected_smooth_elements = 35.7`, `step4a_39limb_cios_mul_instructions = 25584`, `step4a_ecm_B1_log2 = 28.0`, `step4a_ecm_curves_per_trial_log2 = 19.0`, `step4a_ecm_miss_prob_union_bound = 4.2553e-8`, `step4a_instructions_per_trial_log2 = 65.1717`, `step4b_c0_schedule = [64, 128, 550]`, `step4b_deterministic_max_total_nodes_log2 = 15.0379` ($33{,}640$ nodes), `step4b_total_instructions_log2 = 142.2389`
- `cost_and_memory_ledger`: `derived_sieve_instructions_log2 = 143.122918` (charged `143.22142`), `derived_la_MtDM_instructions_log2 = 141.126608` (charged `141.328335`), `charged_descent_instructions_log2 = 142.52098` (`8e42`), `R_total_log2 = 144.13569`, `W_target_group_ops_log2 = 129.847977`, `submitted_time_log2 = 129.9` (`1299 / 10`), `peak_memory_bytes_log2 = 77.091391` (submitted `78.0`).

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
4. **Phase 3 (Unfiltered $\text{flt} = 1$ Rectangular $M^T D M$ Block-Wiedemann Linear Algebra):** Retain $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ distinct raw relations over the $C_{\text{cols}} \le \text{SP}_{\text{orb}}$ active columns without row merges, choose a random diagonal $D \in (\mathbb{F}_r^*)^{R_{\text{rows}}}$, and apply Block-Wiedemann ($n=2, m=4$) to $A = M^T D M$ to recover the virtual logarithms of all active factor-base ideals (including all core ideals of norm $\le B_{\text{core}} = 2^{65.0}$).
5. **Phase 4 (Single-Norm Proper-Subfield LLL Smoothing & 3-Level Special-$\mathfrak{q}$ Descent):** Using $\gcd(r, p^6-1) = 1$, smooth $g^e$ and $z^{e'}$ into a single $2^{165}$-smooth $g_y$-norm via 12-dimensional LLL reduction across $2^{69}$ fresh-$e$ trials (charging concrete 39-limb CIOS Montgomery ECM on every trial), and descend all prime factors down to $B_{\text{core}} = 2^{65.0}$ over $3$ side-aware special-$\mathfrak{q}$ levels to recover $a \equiv \log(z)\cdot\log(g)^{-1} \pmod r$.
6. **Phase 5 (Output Construction & Verification):** Compute $Q = [a]Q_0 \in G_2$, verify $e(P, Q) = z$, and return $Q$.

---

## 2. Algebraic Tower over $h(Y) = Y^6 - 2Y^3 + 2$, Verified Field Isomorphism $\iota$, and Exact Norms

### 2.1 Irreducibility, Non-Residuosity of $u y$, and Machine-Checked Isomorphism $\iota$
From `target-profile.json`, $\mathbb{F}_{p^{12}} = \mathbb{F}_p[w]/(w^{12} - 2w^6 + 2)$. Setting $Y = w^2$ gives:
$$h(Y) = Y^6 - 2Y^3 + 2 = (Y^3 - 1)^2 + 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 2.$$
1. **Irreducibility over $\mathbb{Q}$ and $\mathbb{F}_p$:** $h(Y)$ is $2$-Eisenstein over $\mathbb{Q}$, and because $h(w^2) = w^{12} - 2w^6 + 2$ is irreducible of degree $12$ over $\mathbb{F}_p$, $h(Y)$ is irreducible of degree $6$ over $\mathbb{F}_p$, defining $K = \mathbb{F}_p[y]/(h(y)) \cong \mathbb{F}_{p^6}$.
2. **Non-Residuosity of $u y$ in $\mathbb{F}_{p^6}$:** Since $\deg h = 6$ is even and $h(0) = 2$, $\text{Norm}_{K/\mathbb{F}_p}(u y) = 2 u^6$. Because $p \equiv 3 \pmod 8$, $2^{(p-1)/2} \equiv -1 \pmod p$, so:
   $$(u y)^{(p^6 - 1)/2} = (2 u^6)^{(p - 1)/2} = 2^{(p - 1)/2}(u^3)^{p - 1} \equiv -1 \pmod p.$$
   Thus $g_y(X) = X^2 - u y$ is irreducible of degree $\kappa = 2$ over $\mathbb{F}_{p^6}$, and $X^{p^6} = X(u y)^{(p^6-1)/2} = -X = \tau(X)$.
3. **Explicit, Verified Field Isomorphism $\iota$:**
   Write $u = -2^{16} t$ where $t = \text{0xd20100000001} \equiv 5 \pmod{12}$. Since $3p = P_6(u) \equiv 1 \pmod t$, Quadratic Reciprocity gives $\left(\frac{t}{p}\right) = \left(\frac{p}{t}\right) = \left(\frac{3}{t}\right) = -1$. Because $p \equiv 3 \pmod 4$, $\left(\frac{-1}{p}\right) = -1$, and hence:
   $$\left(\frac{u}{p}\right) = \left(\frac{-1}{p}\right)\left(\frac{2^{16}}{p}\right)\left(\frac{t}{p}\right) = (-1)(+1)(-1) = +1.$$
   Because $u$ is a quadratic residue in $\mathbb{F}_p^*$ and $p \equiv 3 \pmod 4$, the scalar
   $$s = u^{(p+1)/4} \bmod p \in \mathbb{F}_p^*$$
   satisfies $s^2 \equiv u^{(p+1)/2} = u \cdot u^{(p-1)/2} \equiv u \pmod p$. In the STNFS tower $\mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2,\; X^2 - u y)$, the element $W = s^{-1} X$ therefore satisfies:
   $$W^2 = s^{-2} X^2 = u^{-1}(u y) = y, \qquad W^{12} - 2W^6 + 2 = (W^2)^6 - 2(W^2)^3 + 2 = y^6 - 2y^3 + 2 = 0.$$
   Thus the explicit $\mathbb{F}_p$-linear map
   $$\iota : \mathbb{F}_p[w]/(w^{12} - 2w^6 + 2) \xrightarrow{\sim} \mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2,\; X^2 - u y), \qquad \iota\!\left(\sum_{j=0}^{11} c_j w^j\right) = \left(\sum_{k=0}^5 c_{2k} y^k\right) + \left(\sum_{k=0}^5 c_{2k+1} s^{-1} y^k\right) X$$
   is an exact $\mathbb{F}_p$-algebra isomorphism (machine-checked on $W^2 = y$, $W^{12}-2W^6+2=0$, and $\iota(z_1 z_2) = \iota(z_1)\iota(z_2)$ in `experiment.py`), preserving all multiplicative and discrete-logarithm relations between $g = e(P, Q_0)$ and $z$.
4. **Top Polynomial $f_y(X)$ and Involution $\tau : X \mapsto -X$:** Substituting $u = X^2/y$ into $y^6 P_6(u) = y^6(u^6 - 2u^5 + 2u^3 + u + 1)$ and using $y^6 = 2y^3 - 2$ gives:
   $$f_y(X) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + 2 y^3 - 2 \in \mathbb{Z}_y[X],$$
   which satisfies $f_y(X) \equiv (2y^3 - 2)\cdot 3p \equiv 0 \pmod{p, h(y), g_y(X)}$. Since $p^6 \equiv -1 \pmod r$, every conjugate pair of degree-$1$ prime ideals $\{\mathfrak{q}, \tau(\mathfrak{q})\}$ satisfies $\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r$, halving the number of unknown factor-base virtual logarithms to $\text{Li}(B)$.

### 2.2 Exact Integer Resultant Norms via $K_2 = \mathbb{Q}(i)$
Writing $z = y^3$ ($z^2 - 2z + 2 = 0$), any $\gamma(y) = A_0(z) + A_1(z)y + A_2(z)y^2 \in \mathbb{Z}_y$ has relative norm $N_{K_h/K_2}(\gamma) = A_0^3 + z A_1^3 + z^2 A_2^3 - 3 z A_0 A_1 A_2 \equiv c_0 + c_1 z \pmod{z^2 - 2z + 2}$, and exact rational integer norm:
$$N_{K_h/\mathbb{Q}}(\gamma(y)) = (c_0 + c_1)^2 + c_1^2 \in \mathbb{Z}_{\ge 0}.$$
This agrees to machine precision ($\le 2.01 \times 10^{-16}$ relative error) with $\prod_{k=0}^2 |\gamma(\theta_k)|^2$ over the three cube roots $\theta_k = 2^{1/6}\exp(i(\pi/4 + 2\pi k)/3)$ ($k \in \{0, 1, 2\}$) of $z_1 = 1 + i$.

---

## 3. Exact Computation of $1/\zeta_{K_h}(2)$, $\alpha(f_y, h, 10^3)$, $\alpha(g_y, h, 10^3)$, and Smoothness Yield

### 3.1 Exact Local Valuations Across All Prime Ideals ($\ell \le 1000$)
The discriminant of $h(Y) = Y^6 - 2Y^3 + 2$ is $\text{Disc}(h) = -2^8 \cdot 3^6$:
- **At $\ell = 2$ (Exact Algorithm 5.2 Lifting):** $h(Y)$ is $2$-Eisenstein, so $\mathbb{Z}_y$ is locally maximal at $2$ with a single principal prime ideal $\mathfrak{l}_2 = (y)$ of norm $q = 2$ and uniformizer $\gamma = y$.
  - For $f_y(X)$: modulo $y$, $f_y(X) \equiv X^{12}$ has a single multiple root $\bar{s} = 0$. Substituting $X \gets y X$ gives $f_y(yX) = y^{12}X^{12} - 2y^{11}X^{10} + 2y^9 X^6 + y^7 X^2 + (2y^3 - 2)$. Because $2y^3 - 2 = y^6$, dividing by content $y^6$ leaves $y^6 X^{12} - 2y^5 X^{10} + 2y^3 X^6 + y X^2 + 1 \equiv 1 \pmod y$, which has $0$ roots mod $y$. Thus Algorithm 5.2 terminates after $1$ step with $\text{val\_aff} = 6/2 = 3$ and $v_2(f_y) = \frac{2}{2+1}\cdot 3 = 2$.
  - For $g_y(X) = X^2 - u y$: since $2 = y^3(2 - y^3)$ has $\text{val}_{\mathfrak{l}_2}(2) = 6$ and $2^{16} \parallel u$, $\text{val}_{\mathfrak{l}_2}(u y) = 97$. Algorithm 5.2 lifts $\bar{s} = 0$ for $48$ steps of valuation $2$ before terminating at step $49$ ($y X_{49}^2 - \text{unit} \equiv 1 \pmod y$), giving the exact value $v_2(g_y) = \frac{2}{3}\sum_{k=1}^{48} \frac{2}{2^k} = \frac{4}{3}(1 - 2^{-48}) \approx 1.333333$.
- **At $\ell = 3$:** $h(Y) \equiv (Y^2 + Y - 1)^3 \pmod 3$, so $\mathfrak{l}_3 = (3, y^2 + y - 1)$ has $\mathfrak{f}(\mathfrak{l}_3, 3) = 2$ and norm $q = 9$. Since $u \equiv 1 \pmod 3$ and $y^4 \equiv -1 \pmod{\mathfrak{l}_3}$, $g_y(X) = X^2 - u y$ has $0$ roots in $\mathbb{F}_9$, so $v_3(g_y) = 0$ exactly. For $f_y(X)$, since $v_3(f_y) \ge 0$ unconditionally, we use the conservative upper bound $\alpha_3(f_y) \le \frac{\ln 3}{2} \approx 0.5493$ (setting $v_3(f_y) = 0$, which only increases $\alpha(f_y, h)$ and lowers the modeled smoothness probability).
- **At All Primes $\ell \ge 5$:** $\ell \nmid \text{Disc}(h)$, so $\mathbb{Z}_y$ is locally maximal at $\ell$ and $h(Y) \pmod \ell$ is square-free. In `experiment.py`, `factor_h_mod_l_exact(l)` factors $h(Y) \pmod \ell$ into its exact irreducible factors $h_i(Y)$ of degree $d \in \{1, 2, 3, 6\}$ and counts the exact roots $n_{\mathfrak{l}}^{\text{aff}}$ of $f_y(X)$ and $g_y(X)$ in $\mathbb{F}_{\ell^d} = \mathbb{F}_\ell[y]/(h_i(y))$ for every factor $h_i(Y)$ (including Hensel lifting at the bad primes $\ell \in \{7, 13\}$ where $\text{Disc}_U(P_6) = 3^3 \cdot 7 \cdot 13 \cdot 50933$).
- Summing over all primes $\ell \le 1000$ in `experiment.py` yields:
  $$\frac{1}{\zeta_{K_h}(2)} = 0.671179, \qquad \alpha(f_y, h, 10^3) = 1.1706, \qquad \alpha(g_y, h, 10^3) = 0.6739 \quad (\text{in base } e).$$

### 3.2 Barbulescu–Duquesne + Murphy-$\alpha$ and Direct Joint Smoothness Evaluation at $(A = 1360,\; B = 2^{67.500})$
At $A = 1360$ ($2A + 1 = 2721$) and $B = 2^{67.500}$:
- Unreduced sieving volume: $V_0 = \frac{2721^{12}}{2} \approx 8.2738 \times 10^{40} \approx 2^{135.919055}$; $\mathcal{A}=2$ orbit volume: $V_{0,\text{orb}} = \frac{V_0}{2} \approx 2^{134.919055}$.
- Valid ideal-coprime orbit volume (with $1/\zeta_{K_h}(2) = 0.671179 \approx 2^{-0.575231}$):
  $$V_{1,\text{orb}} = V_{0,\text{orb}} \cdot \frac{1}{\zeta_{K_h}(2)} \approx 2.7759 \times 10^{40} \approx 2^{134.343824}.$$
- Across $N = 3{,}000$ ideal-coprime pairs $(a(y), b(y))$ in `experiment.py` (empirical ideal-coprime ratio $0.6676$), the mean logarithmic norms are $\overline{\log_2 N_f} = 830.1716\text{ bits}$ and $\overline{\log_2 N_g} = 515.0101\text{ bits}$.
- Under Heuristic H1, `evaluate_bd_alpha_yield()` in `experiment.py` evaluates **both**:
  1. The **Barbulescu–Duquesne + Murphy-$\alpha$ mean-norm formula** (Section 7.4.2 of Guillevic–Singh, p. 26):
     $$\bar{u}_f = \frac{\overline{\ln N_f} + \alpha(f_y, h, 10^3)}{\ln B} = 12.3239, \qquad \bar{u}_g = \frac{\overline{\ln N_g} + \alpha(g_y, h, 10^3)}{\ln B} = 7.6442,$$
     $$P_f(\bar{u}_f) = \rho(\bar{u}_f) + (1-\gamma)\frac{\rho(\bar{u}_f-1)}{\overline{\ln N_f}}, \qquad P_g(\bar{u}_g) = \rho(\bar{u}_g) + (1-\gamma)\frac{\rho(\bar{u}_g-1)}{\overline{\ln N_g}}, \qquad \bar{p}_{\text{BD}+\alpha} = P_f(\bar{u}_f)P_g(\bar{u}_g) \approx 2^{-70.8987}.$$
  2. The **direct joint sample average** over the $3{,}000$ coprime pairs:
     $$\overline{P_f P_g} = \frac{1}{N}\sum_{i=1}^N P_f(u_{f,i}) P_g(u_{g,i}) \approx 2^{-68.4322},$$
     for which the empirical sample covariance $\widehat{\text{Cov}}(P_f, P_g) = \overline{P_f P_g} - \overline{P_f}\cdot\overline{P_g} > 0$ is strictly positive (because small-coefficient pairs $(a,b)$ reduce both $N_f$ and $N_g$ simultaneously).
- Taking the more conservative of the two estimates ($\bar{p}_{\text{cons}} = \min(\bar{p}_{\text{BD}+\alpha}, \overline{P_f P_g}) = 2^{-70.8987}$) gives the modeled expected orbit relation count:
  $$\lambda_{\text{BD}+\alpha} = V_{1,\text{orb}} \cdot \bar{p}_{\text{cons}} \approx 2^{134.343824 - 70.898704} = 2^{63.445120} \approx 1.2574 \times 10^{19}.$$
- Compared against the $\mathcal{A}=2$ factor-base size $\text{SP}_{\text{orb}} = \text{Li}(2^{67.500}) + 40 \approx 4.5665 \times 10^{18} \approx 2^{61.983842}$ and the required rectangular matrix row count $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}} \approx 5.7081 \times 10^{18} \approx 2^{62.305770}$:
  $$\frac{\lambda_{\text{BD}+\alpha}}{R_{\text{rows}}} = 2.2028, \qquad \frac{\lambda_{\text{BD}+\alpha}}{\text{SP}_{\text{orb}}} = 2.7535.$$
  Conditional on Heuristics H1 and H2 ($\text{RF} \sim \text{Poisson}(\lambda_{\text{BD}+\alpha})$), the survival probability to collect at least $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ relations has z-score $Z = (R_{\text{rows}} - \lambda_{\text{BD}+\alpha})/\sqrt{\lambda_{\text{BD}+\alpha}} \approx -1.94 \times 10^9 \ll -10$, so $P(\text{RF} \ge R_{\text{rows}}) = 1 - \Phi(Z) = 1.000000 > 0.40$.

---

## 4. Derived Schirokauer Columns, Deterministic Row Weight, Rectangular $M^T D M$, and Concrete ECM + Descent

### 4.1 Exact Derivation of $n_{\text{Sch}} = 40$ Schirokauer Unit Columns
Since $h(Y) = (Y^3 - 1)^2 + 1 > 0$ for all $Y \in \mathbb{R}$, $K_h = \mathbb{Q}(y)$ has $r_{\mathbb{R}}(K_h) = 0$ real embeddings and $r_{\mathbb{C}}(K_h) = 3$ complex conjugate pairs. Because $K_f = K_h[X]/(f_y(X))$ ($[K_f:\mathbb{Q}] = 72$) and $K_g = K_h[X]/(g_y(X))$ ($[K_g:\mathbb{Q}] = 12$) are extensions of the totally complex field $K_h$, neither $K_f$ nor $K_g$ has any real embeddings:
$$r_{\mathbb{R}}(K_f) = 0,\quad r_{\mathbb{C}}(K_f) = \frac{72}{2} = 36, \qquad r_{\mathbb{R}}(K_g) = 0,\quad r_{\mathbb{C}}(K_g) = \frac{12}{2} = 6.$$
By Dirichlet's Unit Theorem, the free unit ranks of $\mathcal{O}_{K_f}^*$ and $\mathcal{O}_{K_g}^*$ are $r_f = 36 - 1 = 35$ and $r_g = 6 - 1 = 5$. Since finite torsion roots of unity have order coprime to the $255$-bit prime $r$, they lie in the kernel of the $r$-adic Schirokauer map modulo $r$. Thus the total number of Schirokauer unit-map columns modulo $r$ is **exactly $n_{\text{Sch}} = r_f + r_g = 35 + 5 = 40$**.

### 4.2 Machine-Verified Deterministic Worst-Case Row Weight ($\text{wt}_{\max} \le 285$)
In an unfiltered ($\text{flt} = 1$) relation row $(a(y), b(y))$ with $\|a\|_\infty, \|b\|_\infty \le A = 1360$:
- Each non-zero prime-ideal entry in row $(a,b)$ corresponds to a **distinct** degree-$1$ prime-ideal orbit $\mathfrak{q}$ dividing $N_f(a,b) < 2^{900}$ (in $\mathcal{B}_f$) or $N_g(a,b) < 2^{540}$ (in $\mathcal{B}_g$), with its integer valuation $e(\mathfrak{q}) \ge 1$ stored as the matrix entry value.
- Ordering all degree-$1$ prime ideals $\mathfrak{q}_1, \mathfrak{q}_2, \dots$ of $\mathcal{O}_h$ by non-decreasing norm $N(\mathfrak{q}_1) \le N(\mathfrak{q}_2) \le \dots$, `verify_deterministic_row_weight_bound()` in `experiment.py` verifies that:
  $$\prod_{i=1}^{145} N(\mathfrak{q}_i) = 2^{1197.22} > 2^{900} > \max_{\|a\|_\infty,\|b\|_\infty \le 1360} N_f(a,b), \qquad \prod_{i=1}^{100} N(\mathfrak{q}_i) = 2^{765.96} > 2^{540} > \max_{\|a\|_\infty,\|b\|_\infty \le 1360} N_g(a,b).$$
- Consequently, **no relation row can ever be divisible by more than $145$ distinct prime ideals in $\mathcal{B}_f$ or more than $100$ distinct prime ideals in $\mathcal{B}_g$**. Including the $n_{\text{Sch}} = 40$ Schirokauer columns, the number of non-zero entries in any row of $M$ is **deterministically and unconditionally bounded by $145 + 100 + 40 = 285$**.

### 4.3 Rectangular $M^T D M$ Block-Wiedemann Over $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ Raw Relations
We retain $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ distinct raw relations over the $C_{\text{cols}} \le \text{SP}_{\text{orb}} = \text{Li}(2^{67.500}) + 40 \approx 4.5665 \times 10^{18} \approx 2^{61.983842}$ active columns, and apply Kaltofen–Saunders rectangular Block-Wiedemann ($n=2, m=4$) to $A = M^T D M : \mathbb{F}_r^{C_{\text{cols}}} \to \mathbb{F}_r^{C_{\text{cols}}}$, where $D = \text{diag}(d_1, \dots, d_{R_{\text{rows}}})$ has entries sampled uniformly from $\mathbb{F}_r^*$.
- For each row of $M$ (and $M^T$), multiplying the $\le 245$ integer-valuation prime-ideal entries takes at most $245 \times 32 = 7{,}840$ scalar AArch64 instructions, and multiplying the $40$ full $\mathbb{F}_r$ Schirokauer entries takes $40 \times 256 = 10{,}240$ instructions, totaling **$18{,}080$ scalar AArch64 instructions per row** (plus $256$ instructions per row for diagonal scaling by $D$).
- One application of $A = M^T D M$ costs $1.25\,C_{\text{cols}}(2 \times 18{,}080 + 256) = 45{,}520 \, C_{\text{cols}}$ scalar AArch64 instructions.
- Across all $2.0 C_{\text{cols}}$ applications of $A$ ($1.5 C_{\text{cols}}$ Krylov steps + $0.5 C_{\text{cols}}$ `Mksol` steps) plus $216 \times 256 \, C_{\text{cols}}^2 = 55{,}296 \, C_{\text{cols}}^2$ instructions for quadratic Matrix Berlekamp–Massey generator synthesis:
  $$R_{\text{LA, derived}} = (2.0 \times 45{,}520 + 55{,}296) C_{\text{cols}}^2 = 146{,}336 \, C_{\text{cols}}^2 \approx 3.0518 \times 10^{42} \approx 2^{141.126608}\text{ instructions}.$$

### 4.4 Concrete 39-Limb Montgomery ECM Smoothing (Step 4a) and Side-Aware 3-Level Descent (Step 4b)
1. **Step 4a — Single-Norm Proper-Subfield LLL Smoothing with Concrete 39-Limb CIOS Montgomery ECM:**
   Because $r = u^4 - u^2 + 1 \mid (p^4 - p^2 + 1)$ and $\gcd(r, p^6 - 1) = 1$, every element of $\mathbb{F}_{p^6}^* = (\mathbb{F}_p[y]/(h(y)))^*$ has discrete logarithm $0 \pmod r$ (Guillevic [3, Section 3.1]). For a fresh random exponent $e \in [1, r-1]$, let $\tilde{\tau} = \tau^e = \tau_0(y) + \tau_1(y)X \in \mathbb{F}_{p^{12}}^*$ and $s_\tau(y) = \tau_0(y)\tau_1(y)^{-1} \in \mathbb{F}_{p^6}$. Any non-zero $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the **12-dimensional lattice**:
   $$\Lambda_{\tilde{\tau}}^{(12)} = \big\{(a(y), b(y)) \in \mathbb{Z}_y^2 : a(y) \equiv b(y)s_\tau(y) \pmod{p, h(y)}\big\} \subset \mathbb{Z}^{12}$$
   satisfies $e \cdot \log(\tau) \equiv \log(a(y) + b(y)X) \pmod r$. Because $\Lambda_{\tilde{\tau}}^{(12)}$ has dimension $12$ and determinant $p^6$, weighting $b$ by $|u|^{1/2} \approx 2^{31.855}$ and applying 12-dimensional LLL reduction yields a shortest vector $(a, b)$ with $\|a\|_\infty \le 2^{206}$ and $\|b\|_\infty \le 2^{174}$, producing a **single polynomial $a(y) + b(y)X$** whose $g_y$-norm $N_g(a,b)$ has bit-length $\le 2484$ bits ($L = \lceil 2484/64 \rceil = 39$ limbs of $64$ bits).
   - Setting $B_{\text{init}} = 2^{165}$ gives $u_{\text{init}} = \frac{2484}{165} = 15.0545$ and single-norm smoothness probability $\rho(15.0545) \approx 2^{-63.8423}$. Drawing $N_{\text{init}} = 2^{69}$ **fresh random exponents $e \in [1, r-1]$** per target yields $\lambda_{\text{init}} = 2^{69 - 63.8423} \approx 35.70$ $B_{\text{init}}$-smooth elements per target.
   - **Concrete 39-Limb Montgomery ECM Cost & Failure Analysis:**
     - One $39$-limb CIOS Montgomery multiplication uses at most $16 L^2 + 32 L = 25{,}584$ scalar AArch64 instructions, and one projective Montgomery ladder step (`xDBLADD`, $6\text{M}+5\text{A}$) uses $8 \times 25{,}584 = 204{,}672$ instructions per scalar bit.
     - With Stage-$1$ prime bound $B_1 = 2^{28}$, the Stage-$1$ multiplier $E = \prod_{\ell \le B_1} \ell^{\lfloor \log_\ell B_1 \rfloor}$ has $\le 1.4427\,B_1 \approx 3.873 \times 10^8$ bits, so one ECM curve executes $3.873 \times 10^8 \times 204{,}672 \approx 7.926 \times 10^{13} \approx 2^{46.1717}$ scalar AArch64 instructions.
     - For any prime factor $q \le 2^{165}$, a random Suyama curve has group order $|E(\mathbb{F}_q)|$ divisible by $12$, so $|E(\mathbb{F}_q)|/12$ is $B_1$-smooth with Dickman parameter $u_{\text{ecm}} \le \frac{165 - \log_2 12}{28} = 5.7648$ and probability $\ge \rho(5.7648) \approx 2^{-14.66}$.
     - Running $N_{\text{curves}} = 2^{19}$ independent curves per trial misses a prime factor $q \le 2^{165}$ with probability at most $\exp(-2^{19}\rho(5.7648)) \approx 9.456 \times 10^{-10}$; by the union bound over at most $45$ prime factors of $N_g(a,b)$, ECM fails to factor a $2^{165}$-smooth integer with probability $\le 45 \times 9.456 \times 10^{-10} \approx 4.255 \times 10^{-8} < 10^{-6}$.
     - Total instruction cost per trial (including 12-dim LLL reduction): $2^{19} \times 2^{46.1717} = 2^{65.1717} < 2^{66}$ instructions. Across $2 \times 2^{69}$ trials (for $g$ and $z$), Step 4a costs $2^{70 + 65.1717} = 2^{135.1717}$ instructions.
2. **Step 4b — Side-Aware 3-Level Special-$\mathfrak{q}$ Descent from $B_{\text{init}} = 2^{165}$ to $B_{\text{core}} = 2^{65.0}$:**
   For a prime ideal $\mathfrak{q}$ of norm $Q \in (2^{65.0}, 2^{165}]$, the condition $\mathfrak{q} \mid (a'(y) - b'(y)X)$ defines an index-$Q$ sublattice $L_{\mathfrak{q}} \subset \mathbb{Z}^{12}$ of covolume $Q$. A box $[-C_{0,k}, C_{0,k}]^{12}$ in the reduced basis of $L_{\mathfrak{q}}$ contains $V_{\text{desc}}(C_{0,k}) = \frac{1}{2}(2C_{0,k}+1)^{12}$ candidates (up to sign) with coordinate half-width $C_{0,k} Q^{1/12}$. Scaling relative to the main sieve box $[-A, A]^{12}$ ($A = 1360$, where $\overline{\log_2 N_f} = 830.1716$ and $\overline{\log_2 N_g} = 515.0101$) by $S(Q, C_{0,k}) = C_{0,k} Q^{1/12} / 1360$, the pre-division mean norms are:
   $$\overline{\log_2 N_f}(Q, C_{0,k}) = 830.1716 + 72\!\left(\frac{\log_2 Q}{12} + \log_2\!\left(\frac{C_{0,k}}{1360}\right)\right),$$
   $$\overline{\log_2 N_g}(Q, C_{0,k}) = 515.0101 + 12\!\left(\frac{\log_2 Q}{12} + \log_2\!\left(\frac{C_{0,k}}{1360}\right)\right).$$
   We track the side of $\mathfrak{q}$ explicitly:
   - **If $\mathfrak{q} \in \mathcal{B}_f$ ($f$-side ideal):** $Q \mid N_f$, so $\overline{\log_2 N_f'}^{(f)} = \overline{\log_2 N_f}(Q, C_{0,k}) - \log_2 Q$ and $\overline{\log_2 N_g'}^{(f)} = \overline{\log_2 N_g}(Q, C_{0,k})$.
   - **If $\mathfrak{q} \in \mathcal{B}_g$ ($g$-side ideal, including all root ideals from Step 4a):** $Q \mid N_g$, so $\overline{\log_2 N_f'}^{(g)} = \overline{\log_2 N_f}(Q, C_{0,k})$ and $\overline{\log_2 N_g'}^{(g)} = \overline{\log_2 N_g}(Q, C_{0,k}) - \log_2 Q$.
   Auditing the 3-level schedule `$[165.0, 122.0, 89.0, 65.0]$` with box radii $C_{0,k} \in [64, 128, 550]$ in `audit_descent_tree()` gives:
   - **Level 0 ($165.0 \to 122.0$, $C_{0,0} = 64$, $\log_2 V_{\text{desc}} = 83.13$, $\le 40$ root nodes):**
     - $f$-side: $u_f = 10.965, u_g = 5.140 \implies \log_2 p_{\text{split}}^{(f)} = -52.28$
     - $g$-side: $u_f = 12.317, u_g = 3.788 \implies \log_2 p_{\text{split}}^{(g)} = -54.71$
     - Expected smooth splits per node: $2^{83.13 - 54.71} = 3.59 \times 10^8 \ge 20.0$; deterministic max children $> 2^{65.0}$: $b_{\max, 0} = \lfloor(\overline{\log_2 N_f'} + \overline{\log_2 N_g'})/65.0\rfloor = 30$.
   - **Level 1 ($122.0 \to 89.0$, $C_{0,1} = 128$, $\log_2 V_{\text{desc}} = 95.07$, $\le 40 \times 30 = 1{,}200$ nodes):**
     - $f$-side: $u_f = 13.424, u_g = 6.698 \implies \log_2 p_{\text{split}}^{(f)} = -72.90$
     - $g$-side: $u_f = 14.794, u_g = 5.327 \implies \log_2 p_{\text{split}}^{(g)} = -75.06$
     - Expected smooth splits per node: $2^{95.07 - 75.06} = 1.05 \times 10^6 \ge 20.0$; deterministic max children $> 2^{65.0}$: $b_{\max, 1} = 27$.
   - **Level 2 ($89.0 \to 65.0$, $C_{0,2} = 550$, $\log_2 V_{\text{desc}} = 120.26$, $\le 1{,}200 \times 27 = 32{,}400$ nodes):**
     - $f$-side: $u_f = 18.171, u_g = 9.051 \implies \log_2 p_{\text{split}}^{(f)} = -113.35$
     - $g$-side: $u_f = 19.540, u_g = 7.682 \implies \log_2 p_{\text{split}}^{(g)} = -115.38$
     - Expected smooth splits per node: $2^{120.26 - 115.38} = 29.4 \ge 20.0$.
   - **Deterministic Total Node Count & Instruction Bound:** Across all 3 levels, the total number of nodes per target is at most $40 + 1{,}200 + 32{,}400 = 33{,}640 \approx 2^{15.0379}$ nodes. Charging $2^6 = 64$ instructions per candidate in $V_{\text{desc}}(C_{0,k})$ across both targets ($g$ and $z$) gives total Step 4b work $2^{142.2389}$ instructions.
   - Combining Step 4a ($2^{135.1717}$ instructions) and Step 4b ($2^{142.2389}$ instructions), total Phase 4 work is bounded by $R_{\text{descent}} = 8 \times 10^{42} \approx 2^{142.520980}$ scalar AArch64 instructions.

---

## 5. Summary Scalar AArch64 Instruction and Memory Ledgers ($C = 20,000$)

| Phase / Component | Scalar AArch64 Instructions ($R_i$) | $\log_2(R_i)$ | Target-Group-Op Equivalents ($R_i / 20{,}000$) |
| --- | ---: | ---: | ---: |
| Phase 0: Input check, $g = e(P, Q_0)$ & $\iota$ | $16,777,216$ | $24.000000$ | $838.8608$ |
| Phase 1: Poly selection & factor base | $302,231,454,903,657,293,676,544$ | $78.000000$ | $1.511157 \times 10^{19}$ |
| Phase 2a: 12D Bucket/Line Sieving ($A=1360$) | $13,000,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $143.221420$ | $6.500000 \times 10^{38}$ |
| Phase 2b: ECM cofactorization & Schirokauer | $81,129,638,414,606,681,695,789,005,144,064$ | $106.000000$ | $4.056482 \times 10^{27}$ |
| Phase 3: Unfiltered $M^T D M$ BW, BM & `Mksol` | $3,500,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $141.328335$ | $1.750000 \times 10^{38}$ |
| Phase 4: 39-limb ECM ($2^{69}$) & 3-level descent | $8,000,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $142.520980$ | $4.000000 \times 10^{38}$ |
| Phase 5: $Q = [a]Q_0$ & $e(P, Q) = z$ check | $16,777,216$ | $24.000000$ | $838.8608$ |
| **Total Attack Work** | **$24,500,000,000,081,129,638,716,838,136,599,446,332,375,040$** | **$144.135690$** | **$1.225000 \times 10^{39}$** |

Dividing $R_{\text{total}}$ by $C = 20,000$ gives:
$$W = \frac{R_{\text{total}}}{20,000} \approx 1.22500000 \times 10^{39} \approx 2^{129.847977} < 2^{129.900000}.$$
We submit the exact rational upper bound:
$$\text{time\_log2} = 129.900000 = \frac{1299}{10}, \qquad \text{preprocessing\_log2} = 129.900000, \qquad \text{memory\_log2\_bytes} = 78.000000.$$
Peak memory ($2^{77.091391} < 2^{78.000000}$ bytes) covers the factor base ($2^{67.98}$ bytes), up to $2^{65}$ raw relation records at $4{,}096$ bytes each ($245$ prime-ideal `(uint64, int16)` pairs + $40$ $32$-byte $\mathbb{F}_r$ Schirokauer values $\le 3{,}754 < 4{,}096$ bytes/row, totaling $2^{77}$ bytes), and $32$ Block-Wiedemann vectors of length $R_{\text{rows}} \le 2^{63}$ ($2^{73}$ bytes).

---

## 6. Explicit List of Heuristic Assumptions and Hypotheses

1. **Heuristic H1 (Barbulescu–Duquesne + Guillevic–Singh Murphy-$\alpha$ & Direct Joint Smoothness Heuristic):**
   Assumes that across ideal-coprime pairs $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the sieving box, the joint $B$-smoothness probability of $N_f(a,b)$ and $N_g(a,b)$ is accurately modeled by $\bar{p}_{\text{cons}} = \min(\bar{p}_{\text{BD}+\alpha}, \overline{P_f P_g}) \approx 2^{-70.8987}$. *Support & Limitations:* Standard in TNFS security analyses [1, 2]; buffered here by a $2.2028\times$ oversampling ratio over $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ ($2.7535\times$ over $\text{SP}_{\text{orb}}$).
2. **Heuristic H2 (Poisson Relation Yield Distribution):**
   Assumes that simultaneous $B$-smoothness events across the $V_{1,\text{orb}}$ ideal-coprime pairs behave as weakly dependent Bernoulli trials whose count $\text{RF}$ concentrates around $\lambda_{\text{BD}+\alpha}$ via $\text{Poisson}(\lambda_{\text{BD}+\alpha}) \approx \mathcal{N}(\lambda_{\text{BD}+\alpha}, \lambda_{\text{BD}+\alpha})$.
3. **Heuristic H3 (Rectangular $M^T D M$ Sparse Kernel Rank Over Core Ideals):**
   Assumes that over $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ random smooth relations, the right kernel of $M \pmod r$ (recovered via Kaltofen–Saunders Block-Wiedemann on $M^T D M$ with zero merge filtering and deterministic row weight $\le 285$) uniquely determines the virtual logarithms of all core factor-base ideals of norm $\le B_{\text{core}} = 2^{65.0}$ up to global scaling.
4. **Heuristic H4 (Single-Norm Proper-Subfield LLL Smoothing & Side-Aware Special-$\mathfrak{q}$ Descent Splitting):**
   Assumes that the $2484$-bit norm $N_g(a,b)$ produced by 12-dimensional LLL reduction of $\Lambda_{\tau^e}^{(12)}$ for uniform $e \in [1, r-1]$ and the side-aware residual norms $N_f', N_g'$ in the 3-level special-$\mathfrak{q}$ descent tree split with probabilities given by Dickman's $\rho$-function.

## 7. References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. E. Kaltofen and B. D. Saunders, "On Wiedemann's Method of Solving Sparse Linear Systems," *AAECC-9*, LNCS 539, pp. 29–38, 1991.