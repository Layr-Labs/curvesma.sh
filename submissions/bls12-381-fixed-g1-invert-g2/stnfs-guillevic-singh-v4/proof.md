# Comments
- Judge submission to test TNFS attack,
- Follows the paper from Guillevic-Singh, but modified parameter choice from Tables 10 and 11 after the LLM-judge found a mistake in that selection,
- Created by Gemini 4.0 Argon, human supervised.

# bls12-381-fixed-g1-invert-g2: Special Tower Number Field Sieve (STNFS) via $G_T$ DLP Reduction (Revision v4)

## 0. Revision Record (Addressing Review of Prior Submission `dc4b09ab`)

### 0.1 Baseline / Prior Work (`dc4b09ab`, `631ee3e7`, and `57d7de1`)
In prior submission `dc4b09ab` (`stnfs-guillevic-singh-v3`, following earlier iterations `631ee3e7` and `57d7de1`), we formulated a Special Tower Number Field Sieve (STNFS) attack on `bls12-381-fixed-g1-invert-g2` using the degree-$6$ tower polynomial $h(Y) = Y^6 - 2Y^3 + 2$, $g_y(X) = X^2 - u y$, and $f_y(X) = X^{12} - 2yX^{10} + 2y^3X^6 + y^5X^2 + 2y^3 - 2$.

In the review of `dc4b09ab`, the judge **confirmed and accepted the algebraic tower construction**:
> *"The pairing-to-DLP reduction is sound, and the revised tower's irreducibility/nonresiduosity argument follows from the protected field definition and $\text{Norm}(uy)=2u^6$. Rejection does not require an unconditional proof or a production-scale attack."*

However, `dc4b09ab` was rejected due to a test-script bug and three remaining methodological gaps:
0. **Test Implementation Defect in `verify_algebraic_and_tower_identities()`:** Line 176 assigned `y_p6 = poly_pow_mod(y_p3, P_MOD, ...)`, which computed $(Y^{p^3})^p = Y^{p^4}$ rather than $Y^{p^6}$, causing `assert y_p6 == [0, 1]` to fail before the rest of `experiment.py` executed.
1. **Approximate Local Valuations in `compute_zeta_and_alphas` & Overclaiming a Proved Model-Error Bound:** `compute_zeta_and_alphas` hardcoded valuations at $\ell = 2, 3$ and substituted $1/(q - 1)$ for higher-degree ideals ($d \ge 2$) instead of performing exact root counting across all ideal degrees; and `submission.json` phrased the oversampling margin as if it proved production Dickman model error is bounded by $5\times$.
2. **Special-$\mathfrak{q}$ Descent Residual Norms and Initial Smoothing Cost:**
   - In Step 4a, testing a 24-dimensional rational reconstruction $(a+bX)/(c+dX)$ required two norms to be smooth simultaneously, did not justify charging only $2^{32}$ instructions per trial against ECM smoothness-testing costs, and varied linear combinations of a single LLL basis without re-randomizing the exponent $e$.
   - In Step 4b, the residual norm bound $\log_2 N_f' \le 6\log_2 Q + 360$ gave $u_f + u_g = (6 \cdot 77.7 + 360 + 448)/60.6 \approx 21.0 > 19.2$.
3. **Matrix Rank, Row Weight, and Schirokauer Column Derivation:**
   - Selecting $K_{\text{active}}$ rows and discarding the remaining relations can omit columns or leave a larger kernel.
   - Charging $128$ entries per row relied on an average-weight estimate rather than a deterministic bound on distinct prime-ideal factors per row.
   - The number of Schirokauer unit columns was asserted rather than derived from Dirichlet's Unit Theorem for $K_f$ and $K_g$.

### 0.2 Corrections & Implementation Changes in Revision v4
In this revision (`stnfs-guillevic-singh-v4`), we fix every issue above, verify `experiment.py` execution end-to-end, and replace Algorithm 6.1 with the deterministic **Barbulescu–Duquesne mean-norm formula shifted by Guillevic–Singh's exact Murphy $\alpha$** at $(A = 1360,\; B = 2^{67.500})$:
1. **Fixed Sequential Frobenius Powers in `verify_algebraic_and_tower_identities()` (Fixes Problem 0):**
   `experiment.py` computes `y_p1`, `y_p2`, `y_p3`, `y_p4`, `y_p5`, `y_p6` sequentially by raising each previous polynomial to `P_MOD`, verifying $\gcd(Y^{p^d} - Y, h(Y)) \pmod p = 1$ for $d \in \{1, 2, 3\}$ and $Y^{p^6} \equiv Y \pmod{(p, h(Y))}$, and completes cleanly with `"status": "ok"`.
2. **Exact Local Root Counting Across All Degrees $d \in \{1, 2, 3, 6\}$ and Bad Primes + Barbulescu–Duquesne + $\alpha$ Yield Formula (Fixes Problem 1):**
   - For every prime $\ell \ge 5$, `factor_h_mod_l_exact(l)` factors $h(Y) = (Y^3 - 1)^2 + 1 \pmod \ell$ into its **exact monic irreducible factors $h_i(Y) \in \mathbb{F}_\ell[Y]$** of degree $d \in \{1, 2, 3, 6\}$. For every prime ideal $\mathfrak{l} = (\ell, h_i(y))$ of every degree $d$, `count_roots_fy_good_ideal` and `count_roots_gy_good_ideal` count the **exact roots** in $\mathbb{F}_{\ell^d} = \mathbb{F}_\ell[y]/(h_i(y))$, with exact Hensel root lifting at the bad primes $\ell \in \{7, 13\}$.
   - At $\ell = 2$, $h(Y) = Y^6 - 2Y^3 + 2$ is $2$-Eisenstein, so $\mathbb{Z}[y]$ is locally maximal at $2$ with a single principal prime ideal $\mathfrak{l}_2 = (y)$ of norm $q = 2$ and uniformizer $\gamma = y$. Because $2y^3 - 2 = y^6$ in $\mathbb{Z}[y]/(h(y))$, $f_y(yX)/y^6 = y^6 X^{12} - 2y^5 X^{10} + 2y^3 X^6 + y X^2 + 1 \equiv 1 \pmod y$ has $0$ roots mod $y$, so Algorithm 5.2 (`val_l2_algorithm_5_2()`) terminates in $1$ step and gives the **exact rational valuation** $v_2(f_y) = \frac{2}{3}\cdot 3 = 2$. Similarly, for $g_y(X) = X^2 - u y$ ($\text{val}_{\mathfrak{l}_2}(u y) = 97$), Algorithm 5.2 gives $v_2(g_y) = \frac{2}{3}\sum_{k=1}^{48} \frac{2}{2^k} = \frac{4}{3}(1 - 2^{-48})$.
   - At $\ell = 3$, since $u \equiv 1 \pmod 3$ and $y$ is a non-square in $\mathbb{F}_9 = \mathbb{F}_3[y]/(y^2 + y - 1)$, $g_y(X) = X^2 - u y$ has $0$ roots mod $\mathfrak{l}_3$, giving $v_3(g_y) = 0$ exactly; for $f_y(X)$, since $v_3(f_y) \ge 0$ unconditionally, setting $v_3(f_y) = 0$ gives the unconditional upper bound $\alpha_3(f_y) \le \frac{\ln 3}{2} \approx 0.5493$ (which conservatively lowers the modeled smoothness probability). Across all $\ell \le 1000$, this yields $1/\zeta_{K_h}(2) = 0.671179$, $\alpha(f_y, h, 10^3) = 1.1706$, and $\alpha(g_y, h, 10^3) = 0.6739$.
   - Instead of Algorithm 6.1, we evaluate the expected relation yield via the **Barbulescu–Duquesne [2] mean-norm formula shifted by Guillevic–Singh's $\alpha$** (Section 7.4.2 of Guillevic–Singh, p. 26):
     $$\bar{u}_F = \frac{\overline{\ln N_F} + \alpha(F_y, h, 10^3)}{\ln B}, \qquad P_F = \rho(\bar{u}_F) + (1 - \gamma)\frac{\rho(\bar{u}_F - 1)}{\overline{\ln N_F}}, \qquad \bar{p}_{\text{BD}+\alpha} = P_f \cdot P_g.$$
     By Jensen's inequality (since $u \mapsto \rho(u)$ is strictly convex on $u \ge 1$), evaluating at $\overline{\ln N_F}$ is a deterministic lower bound on the expected Dickman smoothness probability. We explicitly list the Dickman smoothness model as Heuristic H1 in Section 6 rather than claiming a proved bound on production model error.
3. **Proper-Subfield Single-Norm LLL Smoothing ($\gcd(r, p^6-1)=1$) + Full ECM Budget & Verified Descent Table (Fixes Problem 2):**
   - Because $r = u^4 - u^2 + 1 \mid (p^4 - p^2 + 1)$, we have $\gcd(r, p^6 - 1) = 1$ (verified in `experiment.py`). Therefore, **every non-zero element of the intermediate subfield $\mathbb{F}_{p^6} = \mathbb{F}_p(y)$ has discrete logarithm $0 \pmod r$** (Guillevic [3, Section 3.1]).
   - For a fresh random exponent $e \in [1, r-1]$, writing $\tilde{\tau} = \tau^e = \tau_0(y) + \tau_1(y)X$ and $s(y) = \tau_0(y)\tau_1(y)^{-1} \in \mathbb{F}_{p^6}$, any pair $(a(y), b(y)) \in \mathbb{Z}_y^2$ satisfying the $6$ $\mathbb{F}_p$-linear congruences $a(y) \equiv b(y)s(y) \pmod{p, h(y)}$ satisfies $a(y) + b(y)X \equiv b(y)\tau_1(y)^{-1}\tilde{\tau} \pmod{p, h(y), X^2-uy}$ with $b(y)\tau_1(y)^{-1} \in \mathbb{F}_{p^6}^*$, and hence:
     $$e \cdot \log(\tau) \equiv \log\big(a(y) + b(y)X\big) \pmod r.$$
   - Thus $\Lambda_{\tilde{\tau}}^{(12)} \subset \mathbb{Z}^{12}$ is a **12-dimensional lattice of determinant $p^6$** producing a **single polynomial $a(y) + b(y)X$** (no denominator!) with a **single $2484$-bit norm $N_g(a,b)$**.
   - At $B_{\text{init}} = 2^{165}$ ($u_{\text{init}} = \frac{2484}{165} = 15.0545$), the single-norm smoothness probability is $\rho(15.0545) \approx 2^{-63.8423}$. Testing $N_{\text{init}} = 2^{69}$ **fresh random exponents $e \in [1, r-1]$** (one fresh 12-dim LLL reduction per trial, so zero coefficient growth) gives $\lambda_{\text{init}} = 2^{69 - 63.8423} \approx 35.70$ smooth elements per target, and we charge the **full ECM smoothness-testing cost** $L_{2^{165}}(1/2, \sqrt{2}) = \exp(\sqrt{2\ln(2^{165})\ln\ln(2^{165})}) \approx 2^{47.5015}$ ops ($\le 2^{57}$ scalar AArch64 instructions) on every single trial, totaling $2 \times 2^{69} \times 2^{57} = 2^{127}$ instructions.
   - In Step 4b, `audit_descent_tree()` in `experiment.py` audits all $8$ transitions ($165 \to 142 \to 123 \to 107 \to 93 \to 81 \to 71 \to 62.5 \to 55.0 = \log_2 B_{\text{core}}$), verifying that $\log_2 p_{\text{split}} \ge -75.45$ at every level.
4. **Derived $n_{\text{Sch}} = 40$ Schirokauer Columns, Machine-Verified Deterministic Row-Weight Bound ($\le 285$), and Rectangular $M^T D M$ Block-Wiedemann (Fixes Problem 3):**
   - In Section 4.1, we derive the exact unit ranks of $K_f$ ($r_f = 35$) and $K_g$ ($r_g = 5$) from Dirichlet's Unit Theorem, giving **exactly $n_{\text{Sch}} = 35 + 5 = 40$ Schirokauer unit columns**.
   - In `verify_deterministic_row_weight_bound()`, we prove and machine-check that the product of the $145$ smallest degree-$1$ prime-ideal norms in $\mathcal{O}_h$ is $2^{1197.22} > 2^{900} > \max N_f(a,b)$ and the product of the $100$ smallest is $2^{765.96} > 2^{540} > \max N_g(a,b)$. Thus **no relation row can ever contain more than $145 + 100 = 245$ distinct prime-ideal entries** alongside the $40$ Schirokauer columns ($\le 285$ non-zero entries total per row).
   - Instead of discarding rows down to a square submatrix, we retain $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ raw relations ($C_{\text{cols}} \le \text{SP}_{\text{orb}}$) and run Kaltofen–Saunders rectangular Block-Wiedemann on $A = M^T D M : \mathbb{F}_r^{C_{\text{cols}}} \to \mathbb{F}_r^{C_{\text{cols}}}$ for a random diagonal $D$, charging the deterministic worst-case row cost ($18{,}080$ instructions/row) across all $R_{\text{rows}}$ rows.

### 0.3 Summary of Verified Experimental Outputs (`experiment.py`)
Running `experiment.py` with the default initial inputs (`alpha_prime_bound = 1000`, `mc_samples = 3000`, `seed = 42`, `A = 1360`, `B_log2_milli = 67500`) produces:
- `tower_and_subfield_identities_verified`: `true` ($Y^{p^6} \equiv Y \pmod{p, h(Y)}$, $(u y)^{(p^6-1)/2} \equiv -1 \pmod{p, h(y)}$, $\gcd(r, p^6 - 1) = 1$)
- `inv_zeta_Kh_2_euler`: `0.671179` (empirical ideal-coprime ratio: `0.6676`)
- `alpha_f_base_e`: `1.1706`, `alpha_g_base_e`: `0.6739`
- `deterministic_row_weight_proof`: `max_distinct_f = 145` ($\log_2 \text{prod} = 1197.22$), `max_distinct_g = 100` ($\log_2 \text{prod} = 765.96$), `derived_schirokauer_columns = 40`, `deterministic_max_nonzero_entries_per_row = 285`
- `barbulescu_duquesne_alpha_yield`: `mean_log2_Nf = 830.1716`, `mean_log2_Ng = 515.0101`, `u_f_BD_alpha = 12.3239`, `u_g_BD_alpha = 7.6442`, `log2_p_BD_alpha = -70.8987`, `log2_V0 = 135.919055`, `log2_V1_orb = 134.343824`, `log2_SP_orb = 61.983842`, `log2_required_rows_1p25_SP_orb = 62.30577`, `log2_expected_relations_lambda_BD_alpha = 63.44512`, `oversampling_ratio_lambda_over_req_rows = 2.2028` (`2.7535` over $\text{SP}_{\text{orb}}$), `prob_RF_ge_req_rows = 1.0`
- `descent_tree_audit`: `step4a_single_norm_u_init = 15.0545`, `step4a_single_norm_prob_log2 = -63.8423`, `step4a_log2_trials = 69.0`, `step4a_expected_smooth_elements = 35.7`, `step4a_ecm_ops_per_trial_log2 = 47.5015`
- `cost_and_memory_ledger`: `derived_sieve_instructions_log2 = 143.122918` (charged `143.22142`), `derived_la_MtDM_instructions_log2 = 141.126608` (charged `141.328335`), `R_total_log2 = 143.565389`, `W_target_group_ops_log2 = 129.277677`, `submitted_time_log2 = 129.3` (`1293 / 10`), `peak_memory_bytes_log2 = 77.091391` (submitted `78.0`).

---

## 1. Problem Reduction and Complete Attack Architecture

The challenge supplies a uniformly sampled nonidentity point $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and an independently uniform target $z \in G_T \subset \mathbb{F}_{p^{12}}^*$ on BLS12-381, where:
- $u = -\text{0xd201000000010000} = -15132376222941642752$,
- $p = \frac{(u-1)^2(u^4 - u^2 + 1)}{3} + u = \frac{u^6 - 2u^5 + 2u^3 + u + 1}{3}$ ($381$ bits, $\log_2 p \approx 380.0640$, $p \equiv 3 \pmod 8$),
- $r = u^4 - u^2 + 1$ ($255$-bit prime subgroup order, $\log_2 r \approx 254.8571$),
- $Q_0 \in G_2$ is the fixed public generator of $G_2$.

We execute a 6-phase cold-start STNFS attack:
1. **Phase 0 (Input Validation & Pairing):** Verify $P \in G_1 \setminus \{\mathcal{O}_{G_1}\}$ and $z \in G_T$. If $z = 1$, return $Q = \mathcal{O}_{G_2}$ (`null`). Otherwise compute $g = e(P, Q_0) \in G_T \subset \mathbb{F}_{p^{12}}^*$ and map $g, z$ via the explicit isomorphism $\iota$ (Section 2.1) to $\mathbb{F}_p[y, X]/(y^6 - 2y^3 + 2, X^2 - u y)$.
2. **Phase 1 (Tower & Factor Base Setup):** Enumerate the factor base $\mathcal{B} = \mathcal{B}_f \cup \mathcal{B}_g$ of degree-$1$ prime-ideal orbits up to norm $B = 2^{67.500}$ and initialize the $n_{\text{Sch}} = 40$ Schirokauer unit maps modulo $r$.
3. **Phase 2 (12D Bucket/Line Sieving & Cofactorization):** Sieve over $(a(y), b(y)) \in \mathbb{Z}_y^2$ with $\|a\|_\infty, \|b\|_\infty \le A = 1360$ and $\text{lc}(b) > 0, \text{lc}(a) \ge 0$, fold conjugate prime ideals via $\log(\tau(\mathfrak{q})) \equiv -\log(\mathfrak{q}) \pmod r$, and evaluate the $40$ Schirokauer maps on smooth relations.
4. **Phase 3 (Unfiltered $\text{flt} = 1$ Rectangular $M^T D M$ Block-Wiedemann Linear Algebra):** Retain $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ distinct raw relations over the $C_{\text{cols}} \le \text{SP}_{\text{orb}}$ active columns without row merges, choose a random diagonal $D \in (\mathbb{F}_r^*)^{R_{\text{rows}}}$, and apply Block-Wiedemann ($n=2, m=4$) to $A = M^T D M$ to recover the virtual logarithms of all active factor-base ideals (including all core ideals of norm $\le B_{\text{core}} = 2^{55}$).
5. **Phase 4 (Single-Norm Proper-Subfield LLL Smoothing & Special-$\mathfrak{q}$ Descent):** Using $\gcd(r, p^6-1) = 1$, smooth $g^e$ and $z^{e'}$ into a single $2^{165}$-smooth $g_y$-norm via 12-dimensional LLL reduction across $2^{69}$ fresh-$e$ trials (charging full ECM on every trial), and descend all prime factors down to $B_{\text{core}} = 2^{55}$ over $8$ special-$\mathfrak{q}$ steps to recover $a \equiv \log(z)\cdot\log(g)^{-1} \pmod r$.
6. **Phase 5 (Output Construction & Verification):** Compute $Q = [a]Q_0 \in G_2$, verify $e(P, Q) = z$, and return $Q$.

---

## 2. Algebraic Tower over $h(Y) = Y^6 - 2Y^3 + 2$ and Exact Integer Norms

### 2.1 Irreducibility, Non-Residuosity of $u y$, and Isomorphism to `field_Fp12`
From `target-profile.json`, $\mathbb{F}_{p^{12}} = \mathbb{F}_p[w]/(w^{12} - 2w^6 + 2)$. Setting $Y = w^2$ gives:
$$h(Y) = Y^6 - 2Y^3 + 2 = (Y^3 - 1)^2 + 1 \in \mathbb{Z}[Y], \qquad \|h\|_\infty = 2.$$
1. **Irreducibility over $\mathbb{Q}$ and $\mathbb{F}_p$:** $h(Y)$ is $2$-Eisenstein over $\mathbb{Q}$, and because $h(w^2) = w^{12} - 2w^6 + 2$ is irreducible of degree $12$ over $\mathbb{F}_p$, $h(Y)$ is irreducible of degree $6$ over $\mathbb{F}_p$, defining $K = \mathbb{F}_p[y]/(h(y)) \cong \mathbb{F}_{p^6}$.
2. **Non-Residuosity of $u y$ in $\mathbb{F}_{p^6}$:** Since $\deg h = 6$ is even and $h(0) = 2$, $\text{Norm}_{K/\mathbb{F}_p}(u y) = 2 u^6$. Because $p \equiv 3 \pmod 8$, $2^{(p-1)/2} \equiv -1 \pmod p$, so:
   $$(u y)^{(p^6 - 1)/2} = (2 u^6)^{(p - 1)/2} = 2^{(p - 1)/2}(u^3)^{p - 1} \equiv -1 \pmod p.$$
   Thus $g_y(X) = X^2 - u y$ is irreducible of degree $\kappa = 2$ over $\mathbb{F}_{p^6}$, and $X^{p^6} = X(u y)^{(p^6-1)/2} = -X = \tau(X)$.
3. **Explicit Isomorphism $\iota$:** In $\mathbb{F}_{p^6}$, $(y^3 - 1)^2 = -1$ and $u^{(p-1)/2} \equiv -1 \pmod p$, so $s(y) = (-u)^{(p+1)/4}(y^3 - 1) \in \mathbb{F}_{p^6}$ satisfies $s(y)^2 = u$, giving the explicit isomorphism $\iota : w^2 \mapsto y,\; w \mapsto s(y)^{-1} X$.
4. **Top Polynomial $f_y(X)$ and Involution $\tau : X \mapsto -X$:** Substituting $u = X^2/y$ into $y^6 P_6(u) = y^6(u^6 - 2u^5 + 2u^3 + u + 1)$ and using $y^6 = 2y^3 - 2$ gives:
   $$f_y(X) = X^{12} - 2 y X^{10} + 2 y^3 X^6 + y^5 X^2 + 2 y^3 - 2 \in \mathbb{Z}_y[X],$$
   which satisfies $f_y(X) \equiv (2y^3 - 2)\cdot 3p \equiv 0 \pmod{p, h(y), g_y(X)}$. Since $p^6 \equiv -1 \pmod r$, every conjugate pair of degree-$1$ prime ideals $\{\mathfrak{q}, \tau(\mathfrak{q})\}$ satisfies $\log(\tau(\mathfrak{q})) \equiv p^6 \log(\mathfrak{q}) \equiv -\log(\mathfrak{q}) \pmod r$, halving the number of unknown factor-base virtual logarithms to $\text{Li}(B)$.

### 2.2 Exact Integer Resultant Norms via $K_2 = \mathbb{Q}(i)$
Writing $z = y^3$ ($z^2 - 2z + 2 = 0$), any $\gamma(y) = A_0(z) + A_1(z)y + A_2(z)y^2 \in \mathbb{Z}_y$ has relative norm $N_{K_h/K_2}(\gamma) = A_0^3 + z A_1^3 + z^2 A_2^3 - 3 z A_0 A_1 A_2 \equiv c_0 + c_1 z \pmod{z^2 - 2z + 2}$, and exact rational integer norm:
$$N_{K_h/\mathbb{Q}}(\gamma(y)) = (c_0 + c_1)^2 + c_1^2 \in \mathbb{Z}_{\ge 0}.$$
This agrees to machine precision ($\le 2.01 \times 10^{-16}$ relative error) with $\prod_{k=0}^2 |\gamma(\theta_k)|^2$ over the three cube roots $\theta_k = 2^{1/6}\exp(i(\pi/4 + 2\pi k)/3)$ ($k \in \{0, 1, 2\}$) of $z_1 = 1 + i$.

---

## 3. Exact Computation of $1/\zeta_{K_h}(2)$, $\alpha(f_y, h, 10^3)$, $\alpha(g_y, h, 10^3)$, and BD + $\alpha$ Yield

### 3.1 Exact Local Valuations Across All Prime Ideals ($\ell \le 1000$)
The discriminant of $h(Y) = Y^6 - 2Y^3 + 2$ is $\text{Disc}(h) = -2^8 \cdot 3^6$:
- **At $\ell = 2$ (Exact Algorithm 5.2 Lifting):** $h(Y)$ is $2$-Eisenstein, so $\mathbb{Z}_y$ is locally maximal at $2$ with a single principal prime ideal $\mathfrak{l}_2 = (y)$ of norm $q = 2$ and uniformizer $\gamma = y$.
  - For $f_y(X)$: modulo $y$, $f_y(X) \equiv X^{12}$ has a single multiple root $\bar{s} = 0$. Substituting $X \gets y X$ gives $f_y(yX) = y^{12}X^{12} - 2y^{11}X^{10} + 2y^9 X^6 + y^7 X^2 + (2y^3 - 2)$. Because $2y^3 - 2 = y^6$, dividing by content $y^6$ leaves $y^6 X^{12} - 2y^5 X^{10} + 2y^3 X^6 + y X^2 + 1 \equiv 1 \pmod y$, which has $0$ roots mod $y$. Thus Algorithm 5.2 terminates after $1$ step with $\text{val\_aff} = 6/2 = 3$ and $v_2(f_y) = \frac{2}{2+1}\cdot 3 = 2$.
  - For $g_y(X) = X^2 - u y$: since $2 = y^3(2 - y^3)$ has $\text{val}_{\mathfrak{l}_2}(2) = 6$ and $2^{16} \parallel u$, $\text{val}_{\mathfrak{l}_2}(u y) = 97$. Algorithm 5.2 lifts $\bar{s} = 0$ for $48$ steps of valuation $2$ before terminating at step $49$ ($y X_{49}^2 - \text{unit} \equiv 1 \pmod y$), giving the exact value $v_2(g_y) = \frac{2}{3}\sum_{k=1}^{48} \frac{2}{2^k} = \frac{4}{3}(1 - 2^{-48}) \approx 1.333333$.
- **At $\ell = 3$:** $h(Y) \equiv (Y^2 + Y - 1)^3 \pmod 3$, so $\mathfrak{l}_3 = (3, y^2 + y - 1)$ has $\mathfrak{f}(\mathfrak{l}_3, 3) = 2$ and norm $q = 9$. Since $u \equiv 1 \pmod 3$ and $y^4 \equiv -1 \pmod{\mathfrak{l}_3}$, $g_y(X) = X^2 - u y$ has $0$ roots in $\mathbb{F}_9$, so $v_3(g_y) = 0$ exactly. For $f_y(X)$, since $v_3(f_y) \ge 0$ unconditionally, we use the conservative upper bound $\alpha_3(f_y) \le \frac{\ln 3}{2} \approx 0.5493$ (setting $v_3(f_y) = 0$, which only increases $\alpha(f_y, h)$ and lowers the modeled smoothness probability).
- **At All Primes $\ell \ge 5$:** $\ell \nmid \text{Disc}(h)$, so $\mathbb{Z}_y$ is locally maximal at $\ell$ and $h(Y) \pmod \ell$ is square-free. In `experiment.py`, `factor_h_mod_l_exact(l)` factors $h(Y) \pmod \ell$ into its exact irreducible factors $h_i(Y)$ of degree $d \in \{1, 2, 3, 6\}$ and counts the exact roots $n_{\mathfrak{l}}^{\text{aff}}$ of $f_y(X)$ and $g_y(X)$ in $\mathbb{F}_{\ell^d} = \mathbb{F}_\ell[y]/(h_i(y))$ for every factor $h_i(Y)$ (including Hensel lifting at the bad primes $\ell \in \{7, 13\}$ where $\text{Disc}_U(P_6) = 3^3 \cdot 7 \cdot 13 \cdot 50933$).
- Summing over all primes $\ell \le 1000$ in `experiment.py` yields:
  $$\frac{1}{\zeta_{K_h}(2)} = 0.671179, \qquad \alpha(f_y, h, 10^3) = 1.1706, \qquad \alpha(g_y, h, 10^3) = 0.6739 \quad (\text{in base } e).$$

### 3.2 Barbulescu–Duquesne + Murphy-$\alpha$ Relation Yield at $(A = 1360,\; B = 2^{67.500})$
At $A = 1360$ ($2A + 1 = 2721$) and $B = 2^{67.500}$:
- Unreduced sieving volume: $V_0 = \frac{2721^{12}}{2} \approx 8.2738 \times 10^{40} \approx 2^{135.919055}$; $\mathcal{A}=2$ orbit volume: $V_{0,\text{orb}} = \frac{V_0}{2} \approx 2^{134.919055}$.
- Valid ideal-coprime orbit volume (with $1/\zeta_{K_h}(2) = 0.671179 \approx 2^{-0.575231}$):
  $$V_{1,\text{orb}} = V_{0,\text{orb}} \cdot \frac{1}{\zeta_{K_h}(2)} \approx 2.7759 \times 10^{40} \approx 2^{134.343824}.$$
- Across $N = 3{,}000$ ideal-coprime pairs $(a(y), b(y))$ in `experiment.py` (empirical ideal-coprime ratio $0.6676$), the mean logarithmic norms are $\overline{\log_2 N_f} = 830.1716\text{ bits}$ and $\overline{\log_2 N_g} = 515.0101\text{ bits}$.
- Evaluating the **Barbulescu–Duquesne + Murphy-$\alpha$ formula** (Heuristic H1):
  $$\bar{u}_f = \frac{\overline{\ln N_f} + \alpha(f_y, h, 10^3)}{\ln B} = 12.3239, \qquad \bar{u}_g = \frac{\overline{\ln N_g} + \alpha(g_y, h, 10^3)}{\ln B} = 7.6442,$$
  $$P_f = \rho(\bar{u}_f) + (1-\gamma)\frac{\rho(\bar{u}_f-1)}{\overline{\ln N_f}}, \qquad P_g = \rho(\bar{u}_g) + (1-\gamma)\frac{\rho(\bar{u}_g-1)}{\overline{\ln N_g}}, \qquad \bar{p}_{\text{BD}+\alpha} = P_f \cdot P_g \approx 2^{-70.8987},$$
  yields the expected orbit relation count:
  $$\lambda_{\text{BD}+\alpha} = V_{1,\text{orb}} \cdot \bar{p}_{\text{BD}+\alpha} \approx 2^{134.343824 - 70.898704} = 2^{63.445120} \approx 1.2574 \times 10^{19}.$$
- Compared against the $\mathcal{A}=2$ factor-base size $\text{SP}_{\text{orb}} = \text{Li}(2^{67.500}) + 40 \approx 4.5665 \times 10^{18} \approx 2^{61.983842}$ and the required rectangular matrix row count $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}} \approx 5.7081 \times 10^{18} \approx 2^{62.305770}$:
  $$\frac{\lambda_{\text{BD}+\alpha}}{R_{\text{rows}}} = 2.2028, \qquad \frac{\lambda_{\text{BD}+\alpha}}{\text{SP}_{\text{orb}}} = 2.7535.$$
  Under Heuristic H2 ($\text{RF} \sim \text{Poisson}(\lambda_{\text{BD}+\alpha})$), the survival probability to collect at least $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ relations has z-score $Z = (R_{\text{rows}} - \lambda_{\text{BD}+\alpha})/\sqrt{\lambda_{\text{BD}+\alpha}} \approx -1.94 \times 10^9 \ll -10$, so:
  $$P(\text{RF} \ge R_{\text{rows}}) = 1 - \Phi(Z) = 1.000000 > 0.40.$$

---

## 4. Derived Schirokauer Columns, Deterministic Row Weight, Rectangular $M^T D M$, and Single-Norm Descent

### 4.1 Exact Derivation of $n_{\text{Sch}} = 40$ Schirokauer Unit Columns
Since $h(Y) = (Y^3 - 1)^2 + 1 > 0$ for all $Y \in \mathbb{R}$, $K_h = \mathbb{Q}(y)$ has $r_{\mathbb{R}}(K_h) = 0$ real embeddings and $r_{\mathbb{C}}(K_h) = 3$ complex conjugate pairs. Because $K_f = K_h[X]/(f_y(X))$ ($[K_f:\mathbb{Q}] = 72$) and $K_g = K_h[X]/(g_y(X))$ ($[K_g:\mathbb{Q}] = 12$) are extensions of the totally complex field $K_h$, neither $K_f$ nor $K_g$ has any real embeddings:
$$r_{\mathbb{R}}(K_f) = 0,\quad r_{\mathbb{C}}(K_f) = \frac{72}{2} = 36, \qquad r_{\mathbb{R}}(K_g) = 0,\quad r_{\mathbb{C}}(K_g) = \frac{12}{2} = 6.$$
By Dirichlet's Unit Theorem, the free unit ranks of $\mathcal{O}_{K_f}^*$ and $\mathcal{O}_{K_g}^*$ are:
$$r_f = r_{\mathbb{R}}(K_f) + r_{\mathbb{C}}(K_f) - 1 = 35, \qquad r_g = r_{\mathbb{R}}(K_g) + r_{\mathbb{C}}(K_g) - 1 = 5.$$
Since finite torsion roots of unity have order coprime to the $255$-bit prime $r$, they lie in the kernel of the $r$-adic Schirokauer map modulo $r$. Thus the total number of Schirokauer unit-map columns modulo $r$ is **exactly $n_{\text{Sch}} = r_f + r_g = 35 + 5 = 40$**.

### 4.2 Machine-Verified Deterministic Worst-Case Row Weight ($\text{wt}_{\max} \le 285$)
In an unfiltered ($\text{flt} = 1$) relation row $(a(y), b(y))$ with $\|a\|_\infty, \|b\|_\infty \le A = 1360$:
- Each non-zero prime-ideal entry in row $(a,b)$ corresponds to a **distinct** degree-$1$ prime-ideal orbit $\mathfrak{q}$ dividing $N_f(a,b) < 2^{900}$ (in $\mathcal{B}_f$) or $N_g(a,b) < 2^{540}$ (in $\mathcal{B}_g$), with its integer valuation $e(\mathfrak{q}) \ge 1$ stored as the matrix entry value.
- Ordering all degree-$1$ prime ideals $\mathfrak{q}_1, \mathfrak{q}_2, \dots$ of $\mathcal{O}_h$ by non-decreasing norm $N(\mathfrak{q}_1) \le N(\mathfrak{q}_2) \le \dots$, `verify_deterministic_row_weight_bound()` in `experiment.py` computes the cumulative norm product and verifies that:
  $$\prod_{i=1}^{145} N(\mathfrak{q}_i) = 2^{1197.22} > 2^{900} > \max_{\|a\|_\infty,\|b\|_\infty \le 1360} N_f(a,b),$$
  $$\prod_{i=1}^{100} N(\mathfrak{q}_i) = 2^{765.96} > 2^{540} > \max_{\|a\|_\infty,\|b\|_\infty \le 1360} N_g(a,b).$$
- Consequently, **no relation row can ever be divisible by more than $145$ distinct prime ideals in $\mathcal{B}_f$ or more than $100$ distinct prime ideals in $\mathcal{B}_g$**. Including the $n_{\text{Sch}} = 40$ Schirokauer columns, the number of non-zero entries in any row of $M$ is **deterministically and unconditionally bounded by $145 + 100 + 40 = 285$** (and the sum of integer valuations across the $\le 245$ prime-ideal entries is deterministically $\le 900 + 540 = 1440$).

### 4.3 Rectangular $M^T D M$ Block-Wiedemann Over $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ Raw Relations
Instead of discarding excess relations down to a square submatrix or combining rows, we retain $R_{\text{rows}} = 1.25\,C_{\text{cols}}$ distinct raw relations over the $C_{\text{cols}} \le \text{SP}_{\text{orb}} = \text{Li}(2^{67.500}) + 40 \approx 4.5665 \times 10^{18} \approx 2^{61.983842}$ active columns, and apply Kaltofen–Saunders rectangular Block-Wiedemann ($n=2, m=4$) to the symmetric operator:
$$A = M^T D M : \mathbb{F}_r^{C_{\text{cols}}} \to \mathbb{F}_r^{C_{\text{cols}}},$$
where $D = \text{diag}(d_1, \dots, d_{R_{\text{rows}}})$ has entries sampled uniformly from $\mathbb{F}_r^*$.
- Every retained relation is preserved intact with zero fill-in.
- For each row of $M$ (and $M^T$), multiplying the $\le 245$ integer-valuation prime-ideal entries takes at most $245 \times 32 = 7{,}840$ scalar AArch64 instructions, and multiplying the $40$ full $\mathbb{F}_r$ Schirokauer entries takes $40 \times 256 = 10{,}240$ instructions, totaling **$18{,}080$ scalar AArch64 instructions per row** (plus $256$ instructions per row for diagonal scaling by $D$).
- One application of $A = M^T D M$ traverses $M$ ($R_{\text{rows}}$ rows), $D$ ($R_{\text{rows}}$ entries), and $M^T$ ($R_{\text{rows}}$ rows), costing $1.25\,C_{\text{cols}}(2 \times 18{,}080 + 256) = 45{,}520 \, C_{\text{cols}}$ scalar AArch64 instructions.
- Across all $2.0 C_{\text{cols}}$ applications of $A$ ($1.5 C_{\text{cols}}$ Krylov steps + $0.5 C_{\text{cols}}$ `Mksol` steps) plus $216 \times 256 \, C_{\text{cols}}^2 = 55{,}296 \, C_{\text{cols}}^2$ instructions for quadratic Matrix Berlekamp–Massey generator synthesis:
  $$R_{\text{LA, derived}} = (2.0 \times 45{,}520 + 55{,}296) C_{\text{cols}}^2 = 146{,}336 \, C_{\text{cols}}^2 \approx 3.0518 \times 10^{42} \approx 2^{141.126608}\text{ instructions}.$$

### 4.4 Single-Norm Proper-Subfield Smoothing (Step 4a) and Special-$\mathfrak{q}$ Descent Tree to $B_{\text{core}} = 2^{55}$ (Step 4b)
1. **Step 4a — Single-Norm Proper-Subfield LLL Smoothing ($\gcd(r, p^6 - 1) = 1$):**
   Because $r = u^4 - u^2 + 1 \mid (p^4 - p^2 + 1)$ and $\gcd(r, p^6 - 1) = 1$, every element of $\mathbb{F}_{p^6}^* = (\mathbb{F}_p[y]/(h(y)))^*$ has discrete logarithm $0 \pmod r$ (Guillevic [3, Section 3.1]). For a fresh random exponent $e \in [1, r-1]$, let $\tilde{\tau} = \tau^e = \tau_0(y) + \tau_1(y)X \in \mathbb{F}_{p^{12}}^*$ and $s(y) = \tau_0(y)\tau_1(y)^{-1} \in \mathbb{F}_{p^6}$. Any non-zero $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the **12-dimensional lattice**:
   $$\Lambda_{\tilde{\tau}}^{(12)} = \big\{(a(y), b(y)) \in \mathbb{Z}_y^2 : a(y) \equiv b(y)s(y) \pmod{p, h(y)}\big\} \subset \mathbb{Z}^{12}$$
   satisfies $a(y) + b(y)X \equiv b(y)\tau_1(y)^{-1}\tilde{\tau} \pmod{p, h(y), X^2 - u y}$, where $b(y)\tau_1(y)^{-1} \in \mathbb{F}_{p^6}^*$ has logarithm $0 \pmod r$, and therefore:
   $$e \cdot \log(\tau) \equiv \log\big(a(y) + b(y)X\big) \pmod r.$$
   Because $\Lambda_{\tilde{\tau}}^{(12)}$ has dimension $12$ and determinant $p^6$, weighting $b$ by $|u|^{1/2} \approx 2^{31.855}$ and applying 12-dimensional LLL reduction yields a shortest vector $(a, b)$ with $\|a\|_\infty \le 2^{206}$ and $\|b\|_\infty \le 2^{174}$, producing a **single polynomial $a(y) + b(y)X$** whose $g_y$-norm $N_g(a,b)$ has bit-length $\le 2484$ bits.
   - Setting $B_{\text{init}} = 2^{165}$ gives $u_{\text{init}} = \frac{2484}{165} = 15.0545$ and single-norm smoothness probability $\rho(15.0545) \approx 2^{-63.8423}$.
   - Drawing $N_{\text{init}} = 2^{69}$ **fresh random exponents $e \in [1, r-1]$** per target (one fresh 12-dim LLL reduction per trial, avoiding any linear-combination norm growth) yields $\lambda_{\text{init}} = 2^{69 - 63.8423} \approx 35.70$ $B_{\text{init}}$-smooth elements per target (failure probability $e^{-35.70} < 10^{-15}$).
   - Charging the **full ECM smoothness-testing cost** $L_{2^{165}}(1/2, \sqrt{2}) = \exp(\sqrt{2\ln(2^{165})\ln\ln(2^{165})}) \approx 2^{47.5015}$ ops ($\le 2^{57}$ scalar AArch64 instructions, including 12-dim LLL and trial division) on every single trial gives:
     $$R_{\text{descent, 4a}} = 2 \times 2^{69} \times 2^{57} = 2^{127}\text{ scalar AArch64 instructions}.$$
2. **Step 4b — 8-Step Special-$\mathfrak{q}$ Descent from $B_{\text{init}} = 2^{165}$ to $B_{\text{core}} = 2^{55}$:**
   For a prime ideal $\mathfrak{q}$ of norm $Q \in (2^{55}, 2^{165}]$, LLL reduction of the index-$Q$ sublattice $L_{\mathfrak{q}} \subset \mathbb{Z}^{12}$ with residual mean norm bounds $\overline{\log_2 N_f'} \le 5\log_2 Q + 342\text{ bits}$ and $\overline{\log_2 N_g'} \le \log_2 Q + 434\text{ bits}$ is audited across all $8$ steps ($165 \to 142 \to 123 \to 107 \to 93 \to 81 \to 71 \to 62.5 \to 55.0 = \log_2 B_{\text{core}}$) by `audit_descent_tree()` in `experiment.py`:
   - Level 0 ($165 \to 142$): $u_f = 8.218, u_g = 4.218, \log_2 p_{\text{split}} = -34.41$
   - Level 1 ($142 \to 123$): $u_f = 8.553, u_g = 4.683, \log_2 p_{\text{split}} = -37.83$
   - Level 2 ($123 \to 107$): $u_f = 8.944, u_g = 5.206, \log_2 p_{\text{split}} = -41.88$
   - Level 3 ($107 \to 93$): $u_f = 9.430, u_g = 5.817, \log_2 p_{\text{split}} = -46.93$
   - Level 4 ($93 \to 81$): $u_f = 9.963, u_g = 6.506, \log_2 p_{\text{split}} = -52.75$
   - Level 5 ($81 \to 71$): $u_f = 10.521, u_g = 7.254, \log_2 p_{\text{split}} = -59.16$
   - Level 6 ($71 \to 62.5$): $u_f = 11.152, u_g = 8.080, \log_2 p_{\text{split}} = -66.56$ (already $< B = 2^{67.5}$)
   - Level 7 ($62.5 \to 55.0$): $u_f = 11.900, u_g = 9.027, \log_2 p_{\text{split}} = -75.45$.
   Sieving $V_{\text{desc}} = 2^{82}$ per node yields $\ge 2^{6.55} \approx 94$ smooth splits per node; across $< 2^{18}$ total nodes, Step 4b costs $\le 2^{108} \ll 2^{127}$ instructions. Because all terminal prime ideals have norm $\le B_{\text{core}} = 2^{55} \ll B = 2^{67.5}$, each terminal ideal appears on average $R_{\text{rows}} / 2^{55} \approx 2^{62.306 - 55} \approx 158.2$ times in the $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ retained relations of $M$ (probability of $0$ hits in $M$ is $e^{-158.2} < 10^{-68}$).

---

## 5. Summary Scalar AArch64 Instruction and Memory Ledgers ($C = 20,000$)

| Phase / Component | Scalar AArch64 Instructions ($R_i$) | $\log_2(R_i)$ | Target-Group-Op Equivalents ($R_i / 20{,}000$) |
| --- | ---: | ---: | ---: |
| Phase 0: Input check, $g = e(P, Q_0)$ & $\iota$ | $16,777,216$ | $24.000000$ | $838.8608$ |
| Phase 1: Poly selection & factor base | $302,231,454,903,657,293,676,544$ | $78.000000$ | $1.511157 \times 10^{19}$ |
| Phase 2a: 12D Bucket/Line Sieving ($A=1360$) | $13,000,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $143.221420$ | $6.500000 \times 10^{38}$ |
| Phase 2b: ECM cofactorization & Schirokauer | $81,129,638,414,606,681,695,789,005,144,064$ | $106.000000$ | $4.056482 \times 10^{27}$ |
| Phase 3: Unfiltered $M^T D M$ BW, BM & `Mksol` | $3,500,000,000,000,000,000,000,000,000,000,000,000,000,000$ | $141.328335$ | $1.750000 \times 10^{38}$ |
| Phase 4: Single-norm LLL ($2^{69}$ ECM) & descent | $170,141,183,460,469,231,731,687,303,715,884,105,728$ | $127.000000$ | $8.507059 \times 10^{33}$ |
| Phase 5: $Q = [a]Q_0$ & $e(P, Q) = z$ check | $16,777,216$ | $24.000000$ | $838.8608$ |
| **Total Attack Work** | **$16,500,170,141,264,590,107,948,569,823,903,162,216,480,768$** | **$143.565389$** | **$8.250085 \times 10^{38}$** |

Dividing $R_{\text{total}}$ by $C = 20,000$ gives:
$$W = \frac{R_{\text{total}}}{20,000} \approx 8.25008507 \times 10^{38} \approx 2^{129.277677} < 2^{129.300000}.$$
We submit the exact rational upper bound:
$$\text{time\_log2} = 129.300000 = \frac{1293}{10}, \qquad \text{preprocessing\_log2} = 129.300000, \qquad \text{memory\_log2\_bytes} = 78.000000.$$
Peak memory ($2^{77.091391} < 2^{78.000000}$ bytes) covers the factor base ($2^{67.98}$ bytes), up to $2^{65}$ raw relation records at $4{,}096$ bytes each ($245$ prime-ideal `(uint64, int16)` pairs + $40$ $32$-byte $\mathbb{F}_r$ Schirokauer values $\le 3{,}754 < 4{,}096$ bytes/row, totaling $2^{77}$ bytes), and $32$ Block-Wiedemann vectors of length $R_{\text{rows}} \le 2^{63}$ ($2^{73}$ bytes).

---

## 6. Explicit List of Heuristic Assumptions and Hypotheses

1. **Heuristic H1 (Barbulescu–Duquesne + Guillevic–Singh Murphy-$\alpha$ Smoothness Heuristic):**
   Assumes that across ideal-coprime pairs $(a(y), b(y)) \in \mathbb{Z}_y^2$ in the sieving box, the joint $B$-smoothness probability of $N_f(a,b)$ and $N_g(a,b)$ is accurately modeled by $\bar{p}_{\text{BD}+\alpha} = P_f \cdot P_g$ evaluated at the mean logarithmic norms $\overline{\ln N_F}$ shifted by $\alpha(F_y, h, 10^3)$. *Support & Limitations:* Standard in TNFS security analyses [1, 2]; buffered here by a $2.2028\times$ oversampling ratio over $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ ($2.7535\times$ over $\text{SP}_{\text{orb}}$).
2. **Heuristic H2 (Poisson Relation Yield Distribution):**
   Assumes that simultaneous $B$-smoothness events across the $V_{1,\text{orb}}$ ideal-coprime pairs behave as weakly dependent Bernoulli trials whose count $\text{RF}$ concentrates around $\lambda_{\text{BD}+\alpha}$ via $\text{Poisson}(\lambda_{\text{BD}+\alpha}) \approx \mathcal{N}(\lambda_{\text{BD}+\alpha}, \lambda_{\text{BD}+\alpha})$.
3. **Heuristic H3 (Rectangular $M^T D M$ Sparse Kernel Rank Over Core Ideals):**
   Assumes that over $R_{\text{rows}} = 1.25\,\text{SP}_{\text{orb}}$ random smooth relations, the right kernel of $M \pmod r$ (recovered via Kaltofen–Saunders Block-Wiedemann on $M^T D M$ with zero merge filtering and deterministic row weight $\le 285$) uniquely determines the virtual logarithms of all core factor-base ideals of norm $\le B_{\text{core}} = 2^{55}$ up to global scaling.
4. **Heuristic H4 (Single-Norm Proper-Subfield LLL Smoothing & Special-$\mathfrak{q}$ Descent Splitting):**
   Assumes that the $2484$-bit norm $N_g(a,b)$ produced by 12-dimensional LLL reduction of $\Lambda_{\tau^e}^{(12)}$ for uniform $e \in [1, r-1]$ and the residual norms $N_f', N_g'$ in the 8-step special-$\mathfrak{q}$ descent tree split with probabilities given by Dickman's $\rho$-function.

## 7. References
1. A. Guillevic and S. Singh, "On the Alpha Value of Polynomials in the Tower Number Field Sieve Algorithm," *Mathematical Cryptology*, 1(1):1–39, 2021.
2. R. Barbulescu and S. Duquesne, "Updating Key Size Estimations for Pairings," *Journal of Cryptology*, 32(4):1298–1336, 2019.
3. A. Guillevic, "Faster individual discrete logarithms in finite fields of composite extension degree," *Mathematics of Computation*, 88(317):1273–1301, 2019.
4. E. Kaltofen and B. D. Saunders, "On Wiedemann's Method of Solving Sparse Linear Systems," *AAECC-9*, LNCS 539, pp. 29–38, 1991.