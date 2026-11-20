# Zero-Retry Deterministic-DAG Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-`q` nodes**, **`boot_rels = 1`**, **`k = 1` (ZERO internal retries inside the DAG!)**, **root-only exact conditional i.i.d. fresh-`e_(j,t)` boot attempts (`B_max = 10`, `e_(j,t) in {0..r-1}`)**, **heavy-tailed `NB(r = 1.0, alpha = 1.0)` reached-ideal `H5` (`p0(mu) = (1 + mu)^(-1)`)**, **explicit deterministic per-node special-`q` norm cutoff `|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, **expected survivor/ECM work framed under `H-EA`**, **zero-loss degree-preserving excess-row trim + Markowitz deg-0 guard + LIFO `repurge_stack` (`post_peel_cols_lost = post_merge_cols_lost = n_reconstruct_unit = 0`!)**, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`, explicit `D_skew` boot LLL (`S_boot = 2^21`, **run LIVE from scratch on ALL 8/8 genuine `G_T` targets with ZERO supplied LLL vectors!**), and an explicit pointwise aborting combined logarithm phase meter **`t_desc = 2^124.2`** core units (`P(D >= t_desc) <= E[D] / t_desc` by Markov).

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, empirical `q`-lattice offsets, the 24-base production `(4, 3)` `Q(zeta_15)` special-`q` witness set (`prod_sq24`), the 224-ideal **same-field `Q(zeta_15)` real-2-core-FB literal-`k = 1`** depth-1 + depth-2 campaign (`q15_v5k_campaign`), and boot config (`G0`, `m_sv = 2560.0`b, `s_boot = 2^21`, ZERO supplied LLL vectors!). `audit_cm.py` calls `cm_models.eval_descent_live` to compute seeded peeling closure `(h, exc_crit)` and every node's `pr`, `sv`, `mu`, `miss`, `sieve`, expected survivors `cells * (sv + pr)`, `surv_cost`, `own`, `p0 = 1/(1+mu)`, `Phi_dag`, `W_det`, `C_det`, single-attempt `Phi_1`, root-retry `2*Q_1^10`, and Wald expected work 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`E[D] = 2^111.525535 <= 2^111.526`, `Phi_1 = 0.34198`, `2*Q_1^10 = 4.39e-05 <= 0.000044`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`E[D] = 2^111.509266`, `Phi_1 = 0.33468`).

## 1. Fresh-Coin Zero-Retry i.i.d. Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = diag(1..1, 2^21..2^21)`)

For each discrete logarithm (`log g` with `g = e(P0, Q0)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy, Live `F_(p^12) -> F_(p^6)` Subfield Elimination, & Exact `D_skew` 16D LLL (ALL 8/8 LIVE!)
In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). Any `(a, b) in L_Guillevic` with `b != 0 (mod P)` gives **`Psi_g(a - b X) = -(b mod P) w_1^(-1) c_6 z'`** with exact relative norm **`N_(K_g/K_H)(a - b X) = a^3 - omega^2 u b^3`**.
- **Exact `D_skew`-Weighted 16D Integer LLL (`S_boot = 2^21`, `delta = 0.99`):** Multiply the 8 `b`-coordinates by **`S_boot = 2^21 = 2,097,152`**, forming `D_skew * L_Guillevic` (`det = p^4 S_boot^8`), run exact de Weger integer LLL (`cm_models.lll_16d_int`, `delta = 0.99`), and unscale `b = b_scaled / S_boot`. Every reduced basis satisfies the exact `D_skew`-weighted Lovasz bound **`||D_skew b_1||_2 <= (0.99 - 0.25)^(-15/4) (p^4 S_boot^8)^(1/16) = 2^107.30`**.
- **Single-Relation i.i.d. Boot (`boot_rels = 1`, `e_(j,t) in {0, ..., r - 1}`, ZERO excluded `G_T` element, ZERO `BALL_PAIR`):**
  1. For root attempt `j in {1..10}` and boot trial `t = 1..T_boot_max`, draw fresh uniform **`e_(j,t) in {0, ..., r - 1}`** (inclusive of `0`, so for ANY fixed nonidentity `(g', z')`, **`z_(j,t) = z' * g'^(e_(j,t))`** — and `g_(j,t) = g' * g'^(e_(j,t))` for `log g'` — is **EXACTLY uniform over ALL of `G_T`**, excluding zero elements!; if `z_(j,t) == 1`, prob exact `1/r < 2^(-254)`, that trial fails cleanly and is charged in `Q_1`), run `F_(p^12) -> F_(p^6)` elimination to obtain `gamma_t`, run `lll_16d_int` on `D_skew * L_Guillevic(gamma_t)`, and take the **single shortest vector `b_1 = (a_1, b_1)`** (`cells_b = 1`).
  2. Apply the deterministic norm gate **`|N_(K_g/Q)(a_1^3 - omega^2 u b_1^3)| <= 2^2560`** (`m_sv = 2560.0`b). **`RET_BOOT = 0.50`** is declared as a production-law model inequality supported (not proved) by 8 genuine production `G_T` challenges executed LIVE below (`8/8 <= 2543.21`b, Clopper-Pearson 95% LB **`0.6877 >> 0.50`**, plus `32/32 <= 2553.16`b in `[11]`).
  3. Run worst-case multi-factor ECM + limb-scaled 100-round MR on `b_1`: `k_lx = 129` primes in `(2^20, 2^70]`, `kA = 13` in `(2^70, 2^234.53125]`, `1` in `(n_A, 2^265]` (`mr_instr(2560) = 1.21e10` Ir, `det_b = 0.999^143 = 0.8667`).
  4. Take the **FIRST** qualifying boot relation (`boot_rels = 1`), cap trials at **`T_boot_max = ceil(ln(10^4) / p_boot_trial)`** (`<= 10^(-4)`/log, `<= 2e-4` for 2 logs), expected two-log boot work **`2^102.8155`** core units.
- **ALL-8/8 LIVE Genuine Production `G_T` Boot LLL Verification (`audit_cm.py`: `voronoi_and_boot_live`, ZERO supplied LLL vectors!):** (a) Verifies genuine production BLS12-381 `G_T` generator `G_0 != 1, G_0^r == 1` in `F_(p^12)`, generates 8 deterministic seeded exponents `x_t in [1, r-1]` LIVE, and executes exact **`F_(p^12) -> F_(p^6)` subfield elimination LIVE** on ALL 8 challenges `z_t = G_0^(x_t)` to recover `gamma_t`; (b) runs exact de Weger integer **`lll_16d_int(D_skew * B_hnf, delta=0.99)` LIVE FROM SCRATCH ON ALL 8/8 TARGETS** (`boot_supplied_lll_vectors = 0` — ZERO precomputed `v1`/`v2` vectors stored anywhere in the package!), certifying on EVERY SINGLE TARGET `det^2 == p^8 S_boot^16`, all 15 exact Lovasz `99/100` inequalities, all 120 exact size-reduction checks, `b_scaled == 0 mod S_boot`, `a_1 == gamma_t b_1 (mod P)`, `||D_skew b_1||_2 <= 2^106.09 <= 2^107.30`, and **`8/8` first-vector norms `<= 2543.21 <= m_sv = 2560.0`b** (Clopper-Pearson 95% LB **`0.6877 >> RET_BOOT = 0.50`**; on targets `0..3`, finds `v_2` live with `|N(det)| == p^4` and certifies `4/4` 17-cell deciles!), plus constructive 16D LLL on `p_0 = 7207` (`k = 2, index = 1`).

## 2. Exact Seeded Peeling Closure (`h = 0.256845`), Zero-Loss `repurge_stack` Back-Substitution, & `<= LX = 2^70` Discard

- **Seeded Peeling Fixed Point (`cm_models._peel_fixed_point_bls`, live in `audit_cm.py`):** On the two-sided BLS12-381 factor base (`B_f = 2^71.29`, `B_g = 2^71.69`, `lrels = log2(1.09 U') = 63.441701`, `u_f = 548.9 / 71.29`, `u_g = 677.8 / 71.69`), ideals with `lam >= 100.0` form the protected seed. Unprotected ideals satisfy `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))`, giving **`h = 0.256845`** and **`exc_crit = 0.7625 << exc = 1.09`** (**`+43.0%` margin** above peeling collapse).
- **Zero-Loss Degree-Preserving Excess-Row Trim + Markowitz Guard + LIFO `repurge_stack` (`RANK-BRIDGE.md` §3):** During excess-row trimming from the post-singleton-peeling 2-core down to `D = N + 4`, CADO-NFS degree-2 clique removal drops a candidate row ONLY when every incident column has degree `>= 2` and every cascading degree-1 peel has a unique degree-1 column in its row (pushed onto LIFO `repurge_stack`; rollback if any row would acquire `>= 2` degree-1 columns), and Markowitz merging rejects any pivot-row cancellation that would drop a nonpivot column to degree `0`. Consequently **`post_peel_cols_lost == post_merge_cols_lost == 0` by deterministic invariant** (`0` unrecovered post-peeling or post-merge columns, `0` back-substitution errors, `n_reconstruct_unit == 0` across `12/1`, `12/6`, `13/1`, `13/2`, `14/1`). After BW solves the `N` active core columns, strictly 3-stack triangular back-substitution (`L_rec` forward -> `merge_stack` reverse -> `repurge_stack` reverse, `< 2^64` u in `filtering_and_misc = 2^80`; **ZERO 4th unit-prop scan and ZERO residual-cycle solver**) recovers **100% of all post-peeling 2-core columns**, making `miss_side`'s post-peeling formula exact with zero post-trim/post-merge loss!
- **Exact `<= LX = 2^70` Leaf Discard & Conservative `> 2^70` Children:** Leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = kB = 0`). Any prime `<= 2^70` absent from the post-peeling 2-core discards the candidate cleanly (`exp(-miss)` integrated over `[40, 70)`, charging `exp(-x)(1 + x)` at `x = lam(1 - h)` for `lam < 100` and `exp(-lam)` for `lam >= 100`). Every prime **`> 2^70`** (including `(2^70, B_s]`) is conservatively counted as an unknown Window-A/B child in `ss_eval` and descended strictly downward.

## 3. Root-Only i.i.d. Boot Retries (`B_max = 10`), Zero-Internal-Retry DAG (`k = 1`), `NB(r = 1.0)` `H5`, Expected Work under `H-EA`, & Same-Field `Q(zeta_15)` Real-2-Core-FB Evidence

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`, and enumerates `alpha v_1 + beta v_2` (`UDUP_NODE = 0.95`; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b; `[7e]` shift `(+3.40, +2.35)`b inside `(+10, +4)`b).

1. **Zero-Internal-Retry DAG (`k = 1`) + Explicit Norm Cutoff + Expected Work under `H-EA` + Ideal-Log Cache (`P(repeat) = 0`):** Inside a single root attempt, every reached special-`q` node `v` enforces the deterministic norm cutoff **`|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, sieves its fixed ball ONCE, and cofactors qualifying survivors (with expected survivor count `cells * (sv + pr)` and expected sieve-hit / ECM / 100-round MR work priced under declared expected-work model inequality **`H-EA`**, wrapped globally inside the hard aborting meter `t_desc = 2^124.2`). Each node takes the **FIRST** qualifying usable relation (`M_usable(v) = 0 =>` abort this attempt; zero internal retries). Strictly decreasing child norms + a per-attempt ideal-log hash table give `P(repeat) = 0` and deterministic **`C_det_1log = 4,554` nodes/log** (`P(node cap) = 0`).
2. **Unconditional Additive Boole Union Bound under Heavy-Tailed Fixed-`FB` `NB(r = 1.0, alpha = 1.0)` `H5`:** Because collection uses a fixed region, fixed polynomials, and a fixed public ECM seed (independent of the challenge and fresh run coins), the post-peeling 2-core factor base **`FB` is a single fixed deterministic set**. Under declared fixed-`FB` reached-ideal ensemble production-law model inequality **`H5`** (**geometric / `NB(r = 1.0, alpha = 1.0)`**: `Pr_(fresh e)(M_usable(q_sigma) = 0 | slot sigma reached; fixed FB) <= p0(mu_model) = (1 + mu_model)^(-1)`, part of the PGF `(1 + mu_model (1 - t))^(-1)` for `t in [0, 0.8]`; `p0_max = 0.00648` at `min_mu = 153.32`, allowing full geometric overdispersion `CV = 1, Var = mu + mu^2` — **`100x+` wider than Poisson**!), Boole's inequality holds unconditionally over the fresh boot coins `e` for arbitrary slot dependence within a single attempt: **`Phi_1 = 0.341982 <= 0.35`** (`0.334684` at `384/384`), giving single-attempt conditional failure **`Q_1 = Phi_1 + 10^(-4) + 2^(-254) = 0.342082 < 0.35`** uniformly for EVERY `(g', z')` given the fixed `FB`.
3. **Root-Only Exact Conditional i.i.d. Boot Retries (`B_max = 10`, `e_(j,t) in {0..r-1}`) + Exact Conditional Wald Work (`H-EA`) + `t_desc = 2^124.2` Meter (Markov):** For each discrete logarithm, run up to **`B_max = 10` independent root attempts** using fresh uniform **`e_(j,t) in {0, ..., r - 1}`** (clearing the DAG cache between attempts). Conditional on the fixed deterministic `FB` and ANY fixed `(g', z')`, because `z_(j,t) = z' * g'^(e_(j,t))` (resp. `g' * g'^(e_(j,t))`) is **EXACTLY uniform on ALL of `G_T`** (zero excluded element!), attempts `1..B_max` are **100% EXACTLY i.i.d. conditional on `(FB, g', z')`**:
   - **Exact conditional i.i.d. amplification:** `Pr(all B_max fail | fixed FB, g', z') = Pr(attempt 1 fails | fixed FB)^10 <= Q_1^10` per log -> **`2 * Q_1^10 = 4.39e-05 <= 0.000044`** for BOTH logs uniformly over `(g', z')`!
   - **Exact conditional Wald's identity (`E[W_j | N_att >= j; fixed FB, g', z'] = E[W_1 | fixed FB] <= boot_own_exp + W_dag_1`):** Expected two-log work is `E[D | fixed FB, g', z'] <= 2 * ((1 - Q_1^10) / (1 - Q_1)) * (boot_own_exp + W_dag_1) = 1.5199x` single-attempt work = **`2^111.525535 <= 2^111.526`** core units (`2^111.509266` at `384/384`) uniformly over `(g', z')`.
   - **Hard Combined-Log Meter `t_desc = 2^124.2` (Markov):** Pointwise `min(D, t_desc) <= D` preserves `E[min(D, t_desc) | fixed FB, g', z'] <= E[D | fixed FB, g', z']` with **Markov meter-abort `P(D >= t_desc | fixed FB, g', z') <= E[D | fixed FB, g', z'] / t_desc = 0.00015302 <= 0.000157`**.
4. **Two Scoped Same-Field `Q(zeta_15)` Empirical Compatibility Checks (`q15_v5k_descent_live` + `prod_sq24_live` in `audit_cm.py`):**
   - **(a) SAME-FIELD `Q(zeta_15)` Real-2-Core-FB Literal-`k = 1` First-Usable Depth-1 + Depth-2 Campaign (`224` ideals; Stage 1 REAL 2-core FB rebuilt LIVE + `1` complete 2-hop chain + `1` opposite-side control = `4` full 16D `2,048`-pt balls replayed LIVE!):** Exact same cyclotomic base field **`K_H = Q(zeta_15)` (`Phi_15`, deg 8) & order-12 Frobenius `tau`** as production, with `tau`-invariant `g_t(X) = X^3 + 2 omega^2` and `f_t(X) = X^3 - (2 + omega)`:
     - **Stage 1 (REAL collection + REAL 2-core singleton peeling LIVE):** Collects `6,636` genuine `B_FB = 2^14 = 16,384` smooth `Q(zeta_15)` relations from `18,000` short `(a, b)` pairs and executes real 2-core singleton peeling (`6,611` 2-core rels, `|FB_g| = 108, |FB_f| = 100`) — **ZERO prescribed missing-prime list and ZERO synthetic SplitMix64 ECM thinning!**
     - **Stage 2 (16D LLL + LITERAL `k = 1` first-usable against REAL `FB_2core`):** Full 16D LLL on `L_q` (`det = q`), canonical sparse `2,048`-pt ball (`S_b = 11.0`), `D1 = [17.5, 19.2)` (`nc = 17.2`), `D2 = [15.5, 17.2)` (`nc = 15.3`), usable iff BOTH norms strip to `<= 2^nc` AND every prime `<= 2^14` actually belongs to the real peeled 2-core `(FB_g, FB_f)`. Inspects **ONLY literal `M_use == 1`** (ZERO skipping, ZERO side preference) against pre-calibrated `qbins` envelopes (`mu_u_d1 = 16.58 g / 14.47 f`, `mu_u_d2 = 12.74 g / 11.37 f`):

     | Cohort | Ideals | `M_u` range | Empty | `NB(1.0) p0` vs `CP95` | `M_u/mu_u` | `var` vs `NB(1.0)` | `sv2/sv_up` (max) | PGF `emp/nb10` (`0/.2/.5/.8`) |
     |---|---:|---|---:|---|---:|---|---:|---:|
     | `ALL_224` | 224 | `3..165` | **`0/224`** | `0.0625` vs **`0.0133`** (`4.7x`) | `1.7790` | `2.649` vs `3.283` | **`0.0671`** (`0.3980`) | `0 / .00065 / .0181 / .2015` |
     | `D1_100` | 100 | `4..148` | **`0/100`** | `0.0608` vs `0.0295` (`2.1x`) | `1.8409` | `2.277` vs `3.508` | **`0.0744`** (`0.3968`) | `0 / .00031 / .0152 / .1726` |
     | `D2_24` | 24 | `3..65` | **`0/24`** | `0.0768` | `1.4773` | `0.912` vs `2.305` | **`0.0391`** (`0.1228`) | `0 / .00354 / .0403 / .2392` |
     | `CTRL_100` | 100 | `5..165` | **`0/100`** | `0.0608` vs `0.0295` (`2.1x`) | `1.7895` | `3.451` vs `3.318` | **`0.0663`** (`0.3980`) | `0 / .00011 / .0144 / .2194` |

     **Zero empties (`0/224`, `CP95 UB = 0.01328` is `4.7x` LOWER than `NB(1.0)` `p0 = 0.06248`!)**, **`D1/CTRL = 1.0287`** (`>= 1.0`, zero depth-1 degradation!), ALL four cohorts satisfy `PGF_emp(t) / PGF_NB(1.0)(t) <= 0.24 < 1.0` (and `<= 0.17` vs `NB(1.5)`!) across `t in {0.0, 0.2, 0.5, 0.8}`, and actual two-sided **`sv2 / sv_upper = 0.0671`** (`max = 0.3980 < 1.0`)!
   - **(b) 24-Base Actual Production BLS12-381 `(4, 3)` `Q(zeta_15)` `(f_9, g_3)` Special-`q` 16D Lattice Check (`audit_cm.py`: `prod_sq24_live`):** Across `12 g-side + 12 f-side` actual production primes spanning `73..258` bits (`384` exact production deg-72/deg-24 norms): `det(v1, v2) % q == 0` and exact `q | N` verified LIVE; all 6 empirical sextiles `p16..max` lie strictly below `qbins` `c_hi + cf, c_hi + cg` (`g-side` headroom `>= 39.68`b `f` / `>= 23.87`b `g`; `f-side` `>= 37.39`b `f` / `>= 27.31`b `g`); smoothness-model survivor rates at the actual empirical production norms vs the `qbins_lo` plan envelope give mean **`0.0503` (`g`, max `0.1226`)** and **`0.2602` (`f`, max `1.0000`)**.
   - **Scale-Transfer & Explicit Exploratory Scope:** Scaling within the SAME base field `Q(zeta_15)` and `16D` special-`q` lattice family from `(B_FB = 2^14, mu ~ 11..17)` to production `(LX = 2^70, mu >= 153)` increases the number of primes in `FB` by `> 10^16x` and replaces `2,048`-pt balls by `> 2^90`-pt balls (where relative norm variation across 16D balls concentrates, directly confirmed at `73..258`-bit production `q` by `prod_sq24_live`), while `NB(1.0)` preserves full geometric overdispersion (`p0 = 1/(1+mu)`, `CV = 1`) rather than Poisson concentration. Together, `q15_v5k_descent_live` (`0/224`, `CP95 0.0133 < 0.0625`, `PGF_emp/nb10 <= 0.24`) and `prod_sq24_live` provide scoped same-field finite-scale + production-norm compatibility evidence (NOT statistical proof of production tails); production tail and expected-work transfer rest explicitly on declared exploratory fixed-`FB` model inequalities `H5` and `H-EA`.

## 4. Complete Live-Computed Deterministic-DAG Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|p0|Phi_dag|Nf/Ng|own|W_det|C_det|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|76|70|70.0|0|0|16983.8|5.9e-5|0.00006|841/816|110.55|110.55|1|
|[70,76)|g|sieve|80|70|70.0|0|0|25377.3|3.9e-5|0.00004|916/748|111.74|111.74|1|
|[76,81)|f|sieve|68|76|76.0|1|0|45631.2|2.2e-5|0.00012|846/820|108.01|112.34|3|
|[76,81)|g|sieve|68|76|76.0|1|0|43766.8|2.3e-5|0.00012|922/744|108.01|112.34|3|
|[81,86)|f|sieve|72|76|76.0|1|0|19258.2|5.2e-5|0.00015|870/829|109.31|112.44|3|
|[81,86)|g|sieve|76|76|76.0|1|0|29242.7|3.4e-5|0.00013|956/750|110.56|112.65|3|
|[86,91)|f|sieve|68|76|76.0|2|0|20379.9|4.9e-5|0.00025|881/835|110.45|113.46|5|
|[86,91)|g|sieve|72|76|76.0|2|0|29705.0|3.4e-5|0.00023|973/751|111.25|113.58|5|
|[91,96)|f|sieve|76|76|76.0|2|0|14563.3|6.9e-5|0.00027|910/846|111.45|113.63|5|
|[91,96)|g|sieve|80|76|76.0|2|0|19078.8|5.2e-5|0.00025|1006/757|112.28|113.86|5|
|[96,101)|f|sieve|84|76|76.0|2|0|9030.6|1.1e-4|0.00031|938/857|113.07|114.17|5|
|[96,101)|g|sieve|84|76|76.0|2|0|7228.9|1.4e-4|0.00033|1034/761|113.02|114.15|5|
|[101,106)|f|sieve|84|76|76.0|2|0|2321.3|4.3e-4|0.00063|956/864|112.93|114.11|5|
|[101,106)|g|sieve|80|81|75.5|1|1|4670.7|2.1e-4|0.00055|1051/762|111.84|114.21|9|
|[106,111)|f|sieve|88|81|75.5|1|1|3432.0|2.9e-4|0.00063|978/874|113.97|114.93|9|
|[106,111)|g|sieve|84|81|81.0|2|0|4850.1|2.1e-4|0.00069|1079/766|112.91|114.79|13|
|[111,116)|f|sieve|76|81|75.5|2|1|1512.4|6.6e-4|0.00110|980/876|115.22|115.84|11|
|[111,116)|g|sieve|80|81|75.5|2|1|1636.1|6.1e-4|0.00105|1096/767|115.18|115.81|11|
|[116,121)|f|sieve|84|86|86.0|2|0|2627.4|3.8e-4|0.00095|1008/887|112.88|114.94|13|
|[116,121)|g|sieve|80|96|83.0|1|1|4259.2|2.3e-4|0.00103|1119/769|111.77|115.39|17|
|[121,126)|f|sieve|84|96|83.0|1|1|2612.8|3.8e-4|0.00118|1026/894|112.88|115.52|17|
|[121,126)|g|sieve|84|101|85.5|1|1|4700.9|2.1e-4|0.00114|1147/773|112.88|115.78|17|
|[126,131)|f|sieve|84|96|96.0|2|0|2715.4|3.7e-4|0.00140|1043/902|112.88|115.93|21|
|[126,131)|g|sieve|88|116|100.2|2|0|5846.9|1.7e-4|0.00146|1174/778|113.95|116.44|21|
|[131,136)|f|sieve|80|96|83.0|2|1|3662.5|2.7e-4|0.00135|1055/908|114.13|116.08|23|
|[131,136)|g|sieve|84|96|83.0|2|1|2866.4|3.5e-4|0.00143|1192/778|113.97|116.04|23|
|[136,141)|f|sieve|84|101|85.5|2|1|6429.5|1.6e-4|0.00136|1078/917|114.78|116.44|23|
|[136,141)|g|sieve|88|101|85.5|2|1|4439.8|2.3e-4|0.00143|1219/783|114.71|116.42|23|
|[141,146)|f|sieve|72|96|96.0|3|0|1119.6|8.9e-4|0.00243|1080/919|113.16|116.48|31|
|[141,146)|g|sieve|88|101|85.5|2|1|1085.0|9.2e-4|0.00213|1241/785|114.21|116.28|23|
|[146,151)|f|sieve|68|106|93.6|2|1|767.3|1.3e-3|0.00351|1091/925|111.46|116.54|35|
|[146,151)|g|sieve|68|101|101.0|3|0|434.4|2.3e-3|0.00422|1237/779|111.40|116.78|31|
|[151,156)|f|sieve|80|106|93.6|2|1|820.1|1.2e-3|0.00343|1125/938|112.86|116.61|35|
|[151,156)|g|sieve|84|106|98.1|2|1|691.4|1.4e-3|0.00391|1282/788|113.48|116.90|35|
|[156,161)|f|sieve|76|116|100.2|2|1|2086.2|4.8e-4|0.00391|1138/943|113.35|117.61|45|
|[156,161)|g|sieve|92|106|106.0|3|0|1573.2|6.4e-4|0.00418|1313/794|115.70|117.33|43|
|[161,166)|f|sieve|84|116|100.2|2|1|1294.9|7.7e-4|0.00420|1166/954|113.91|117.64|45|
|[161,166)|g|sieve|92|116|105.9|2|1|1548.1|6.5e-4|0.00515|1336/797|115.42|117.84|53|
|[166,171)|f|sieve|88|126|98.0|2|1|1354.8|7.4e-4|0.00434|1188/964|114.43|117.69|55|
|[166,171)|g|sieve|84|131|100.5|2|1|806.1|1.2e-3|0.00538|1349/796|113.16|117.83|63|
|[171,176)|f|sieve|88|126|106.8|2|1|1607.1|6.2e-4|0.00558|1205/971|115.00|118.03|79|
|[171,176)|g|sieve|88|131|110.0|2|1|1390.7|7.2e-4|0.00621|1376/800|114.47|118.16|87|
|[176,181)|f|sieve|76|141|105.5|2|1|649.2|1.5e-3|0.00675|1208/973|112.39|117.98|75|
|[176,181)|g|sieve|64|131|100.5|3|1|500.5|2.0e-3|0.00677|1366/792|114.61|118.13|73|
|[181,186)|f|sieve|68|146|108.0|2|1|179.1|5.6e-3|0.01276|1214/977|110.81|118.20|99|
|[181,186)|g|sieve|68|131|100.5|3|1|204.0|4.9e-3|0.00966|1395/796|113.92|118.08|73|
|[186,191)|f|sieve|72|151|110.5|2|1|211.7|4.7e-3|0.01507|1237/987|111.72|118.33|111|
|[186,191)|g|sieve|60|141|105.5|3|1|219.3|4.5e-3|0.01093|1404/794|112.82|118.18|89|
|[191,196)|f|sieve|72|161|115.5|2|1|396.0|2.5e-3|0.01495|1255/994|113.28|119.21|137|
|[191,196)|g|sieve|72|141|105.5|3|1|252.6|3.9e-3|0.01034|1446/803|114.03|118.23|89|
|[196,201)|f|sieve|80|166|118.0|2|1|644.4|1.5e-3|0.01520|1283/1005|114.96|119.42|159|
|[196,201)|g|sieve|76|151|110.5|3|1|1084.9|9.2e-4|0.01261|1474/807|116.02|118.79|133|
|[201,206)|f|sieve|84|171|120.5|2|1|679.2|1.5e-3|0.01548|1306/1014|115.71|119.48|179|
|[201,206)|g|sieve|68|146|119.9|3|1|627.0|1.6e-3|0.01260|1485/806|116.00|119.18|145|
|[206,211)|f|sieve|68|181|125.5|2|1|202.6|4.9e-3|0.02307|1301/1015|114.79|119.65|235|
|[206,211)|g|sieve|64|161|115.5|3|1|153.3|6.5e-3|0.02106|1501/807|112.76|119.46|161|
|[211,216)|f|sieve|80|181|125.5|2|1|204.5|4.9e-3|0.02303|1335/1028|115.43|119.68|235|
|[211,216)|g|sieve|68|166|118.0|3|1|191.8|5.2e-3|0.02099|1530/811|113.40|119.60|189|
|[216,221)|f|ecm|20|171|136.3|5|1|128915.8|7.8e-6|0.02398|1209/987|124.24|124.33|349|
|[216,221)|g|sieve|72|171|120.5|3|1|237.5|4.2e-3|0.02035|1558/816|114.06|119.63|209|
|[221,226)|f|ecm|20|176|139.6|5|1|203540.6|4.9e-6|0.02605|1227/995|125.07|125.13|397|
|[221,226)|g|sieve|76|181|125.5|3|1|865.7|1.2e-3|0.02163|1586/820|116.22|119.91|269|
|[226,231)|f|ecm|20|171|148.9|5|1|595310.3|1.7e-6|0.04838|1244/1002|126.59|126.61|449|
|[226,231)|g|sieve|88|181|125.5|3|1|664.2|1.5e-3|0.02198|1624/828|116.86|119.97|269|
|[231,236)|f|ecm|20|176|139.6|6|1|866282.7|1.2e-6|0.02890|1262/1010|125.67|125.72|443|
|[231,236)|g|sieve|72|176|139.6|3|1|595.0|1.7e-3|0.02202|1626/823|117.41|120.27|305|
|[236,241)|f|ecm|20|176|139.6|6|1|350349.0|2.9e-6|0.02890|1279/1017|125.67|125.72|443|
|[236,241)|g|sieve|48|191|149.4|3|1|185.7|5.4e-3|0.05456|1606/812|114.83|120.29|399|
|[241,246)|f|ecm|20|171|148.9|6|1|1372964.7|7.3e-7|0.05612|1297/1025|127.35|127.36|515|
|[241,246)|g|sieve|56|191|149.4|3|1|180.4|5.5e-3|0.05472|1644/819|115.48|120.31|399|
|[246,251)|f|ecm|20|191|149.4|6|1|3991958.1|2.5e-7|0.07241|1314/1032|127.43|127.45|597|
|[246,251)|g|sieve|68|191|149.4|3|1|227.1|4.4e-3|0.05359|1687/829|116.76|120.38|399|
|[251,256)|f|ecm|20|191|149.4|6|1|1698034.8|5.9e-7|0.07241|1332/1040|127.43|127.45|597|
|[251,256)|g|sieve|68|196|152.7|3|1|183.6|5.4e-3|0.05463|1710/831|116.76|120.69|437|
|[256,261)|f|ecm|20|191|149.4|6|1|696640.2|1.4e-6|0.07241|1349/1047|127.43|127.45|597|
|[256,261)|g|sieve|72|201|156.0|3|1|243.9|4.1e-3|0.05510|1738/836|117.60|120.95|503|
|[261,266)|f|ecm|20|186|160.6|6|1|4323408.0|2.3e-7|0.07122|1367/1055|129.57|129.57|715|
|[261,266)|g|sieve|76|206|159.2|3|1|310.5|3.2e-3|0.05570|1766/840|118.69|121.44|589|
