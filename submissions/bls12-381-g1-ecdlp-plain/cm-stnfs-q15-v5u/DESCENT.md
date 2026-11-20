# Zero-Retry Deterministic-DAG Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn `rho_dB` on boot & all nodes**, **`boot_rels = 1`**, **`k = 1` (zero internal retries)**, **root-only exact conditional i.i.d. fresh-`e_(j,t)` boot retries (`B_max = 10`, `e_(j,t) in {0..r-1}`)**, **fixed-`FB` `NB(r = 1.0, alpha = 1.0)` `H5` (`p0 = (1 + mu)^(-1)`)**, **per-node cutoff `|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, **expected work under `H-EA`**, **zero-loss `repurge_stack` (`post_peel_cols_lost = post_merge_cols_lost = 0`)**, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`, `D_skew` `S_boot = 2^21` (**ALL 8/8 LIVE from scratch, `boot_supplied_lll_vectors = 0`**), and pointwise meter **`t_desc = 2^124.2`** (Markov `P(D >= t_desc) <= E[D]/t_desc`).

`cm_data.json` stores ONLY the 78 structural tuples `plan_struct`, `q`-lattice offsets, `prod_sq24` + `prod_enc12` (`36` actual production `(f_9, g_3)` lattices), and boot config (`G0`, `m_sv = 2560.0`b, `s_boot = 2^21`). `audit_cm.py` calls `cm_models.eval_descent_live` to compute `(h, exc_crit)` and all 78 nodes LIVE (`MESH = 96/192`: `E[D] = 2^111.525535 <= 2^111.526`, `Phi_1 = 0.34198`, `2*Q_1^10 = 4.39e-05 <= 0.000044`; strictly tighter at `{"descent_mesh": 384}`: `E[D] = 2^111.509266`, `Phi_1 = 0.33468`).

## 1. Fresh-Coin Zero-Retry i.i.d. Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = diag(1..1, 2^21..2^21)`)

For each discrete logarithm (`log g` with `g = e(P0, Q0)` and `log z`), each boot trial draws fresh uniform `e_(j,t) in {0, ..., r - 1}` (inclusive of `0`, fresh coins) and forms `z_(j,t) = z' * g'^(e_(j,t))` (resp. `g_(j,t) = g' * g'^(e_(j,t))`), exact uniform on all of `G_T` (`P(identity) = 1/r < 2^(-254)` charged in `Q_1`).

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
- **ALL-8/8 LIVE Genuine Production `G_T` Boot LLL (`voronoi_and_boot_live`, `boot_supplied_lll_vectors = 0`):** Verifies `G_0 != 1, G_0^r == 1`, runs `F_(p^12) -> F_(p^6)` elimination AND from-scratch 381-bit pure-integer de Weger `lll_16d_int(99/100)` LIVE on ALL `8/8` targets (`det^2 == p^8 S_boot^16`, 15 Lovasz + 120 size-red. each, `||D_skew b_1||_2 <= 2^106.09 <= 2^107.30`, **`8/8 <= 2544.69 <= 2560.0`b**, CP95 LB **`0.6877 >> RET_BOOT = 0.50`**; `|N(det)| == p^4` & `4/4` deciles on `0..3`), plus constructive 16D pure-integer `lll_16d_int` on `p_0 = 7207` (`k = 2, index = 1`; ONE unified exact integer LLL reducer everywhere, ZERO float LLL).

## 2. Exact Seeded Peeling Closure (`h = 0.256845`), Zero-Loss `repurge_stack` Back-Substitution, & `<= LX = 2^70` Discard

- **Seeded Peeling Fixed Point (`cm_models._peel_fixed_point_bls`, live in `audit_cm.py`):** On the two-sided BLS12-381 factor base (`B_f = 2^71.29`, `B_g = 2^71.69`, `lrels = log2(1.09 U') = 63.441701`, `u_f = 548.9 / 71.29`, `u_g = 677.8 / 71.69`), ideals with `lam >= 100.0` form the protected seed. Unprotected ideals satisfy `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))`, giving **`h = 0.256845`** and **`exc_crit = 0.7625 << exc = 1.09`** (**`+43.0%` margin** above peeling collapse).
- **Zero-Loss Degree-Preserving Excess-Row Trim + Markowitz Guard + LIFO `repurge_stack` (`RANK-BRIDGE.md` §3):** During excess-row trim down to `D = N + 4`, degree-2 clique removal drops a row ONLY when all incident columns have degree `>= 2` and every cascading degree-1 peel touches a row with a unique degree-1 column (pushed onto LIFO `repurge_stack`; rollback otherwise), and Markowitz merging rejects any elimination dropping a nonpivot column to degree `0`. Hence **`post_peel_cols_lost == post_merge_cols_lost == 0` by deterministic invariant** (`n_reconstruct_unit == 0`). Post-BW 3-stack triangular back-substitution (`L_rec` -> `merge_stack` -> `repurge_stack`, `< 2^64` u in `filtering_and_misc = 2^80`) recovers **100%** of post-peeling 2-core columns, making `miss_side` exact!
- **Exact `<= LX = 2^70` Leaf Discard & Conservative `> 2^70` Children:** Leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = kB = 0`). Any prime `<= 2^70` absent from the 2-core discards the candidate (`exp(-miss)` over `[40, 70)`, charging `exp(-x)(1 + x)` at `x = lam(1 - h)` for `lam < 100` and `exp(-lam)` for `lam >= 100`); every prime **`> 2^70`** is conservatively descended.

## 3. Root-Only i.i.d. Boot Retries (`B_max = 10`), Zero-Internal-Retry DAG (`k = 1`), `NB(r = 1.0)` `H5`, Expected Work under `H-EA`, & Direct Production `(f_9, g_3)` `Q(zeta_15)` Audits

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D pure-integer `lll_16d_int` extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`, and enumerates `alpha v_1 + beta v_2` (`UDUP_NODE = 0.95`; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b; live `[7e]` production encountered-vs-control shift `g: (-2.36, -3.35)`b, `f: (+0.35, -2.27)`b — well inside `(+10, +4)`b, alongside offline GP `[7e]` `(+3.40, +2.35)`b).

1. **Zero-Internal-Retry DAG (`k = 1`) + Explicit Norm Cutoff + Expected Work under `H-EA` + Ideal-Log Cache (`P(repeat) = 0`):** Inside a single root attempt, every reached special-`q` node `v` enforces the deterministic norm cutoff **`|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, sieves its fixed ball ONCE, and cofactors qualifying survivors (with expected survivor count `cells * (sv + pr)` and expected sieve-hit / ECM / 100-round MR work priced under declared expected-work model inequality **`H-EA`**, wrapped globally inside the hard aborting meter `t_desc = 2^124.2`). Each node takes the **FIRST** qualifying usable relation (`M_usable(v) = 0 =>` abort this attempt; zero internal retries). Strictly decreasing child norms + a per-attempt ideal-log hash table give `P(repeat) = 0` and deterministic **`C_det_1log = 4,554` nodes/log** (`P(node cap) = 0`).
2. **Unconditional Additive Boole Union Bound under Heavy-Tailed Fixed-`FB` `NB(r = 1.0, alpha = 1.0)` `H5`:** Because collection uses a fixed region, fixed polynomials, and a fixed public ECM seed (independent of the challenge and fresh run coins), the post-peeling 2-core factor base **`FB` is a single fixed deterministic set**. Under declared fixed-`FB` reached-ideal ensemble production-law model inequality **`H5`** (**geometric / `NB(r = 1.0, alpha = 1.0)`**: `Pr_(fresh e)(M_usable(q_sigma) = 0 | slot sigma reached; fixed FB) <= p0(mu_model) = (1 + mu_model)^(-1)`, part of the PGF `(1 + mu_model (1 - t))^(-1)` for `t in [0, 0.8]`; `p0_max = 0.00648` at `min_mu = 153.32`, allowing full geometric overdispersion `CV = 1, Var = mu + mu^2` — **`100x+` wider than Poisson**!), Boole's inequality holds unconditionally over the fresh boot coins `e` for arbitrary slot dependence within a single attempt: **`Phi_1 = 0.341982 <= 0.35`** (`0.334684` at `384/384`), giving single-attempt conditional failure **`Q_1 = Phi_1 + 10^(-4) + 2^(-254) = 0.342082 < 0.35`** uniformly for EVERY `(g', z')` given the fixed `FB`.
3. **Root-Only Exact Conditional i.i.d. Boot Retries (`B_max = 10`, `e_(j,t) in {0..r-1}`) + Exact Conditional Wald Work (`H-EA`) + `t_desc = 2^124.2` Meter (Markov):** For each discrete logarithm, run up to **`B_max = 10` independent root attempts** using fresh uniform **`e_(j,t) in {0, ..., r - 1}`** (clearing the DAG cache between attempts). Conditional on the fixed deterministic `FB` and ANY fixed `(g', z')`, because `z_(j,t) = z' * g'^(e_(j,t))` (resp. `g' * g'^(e_(j,t))`) is **EXACTLY uniform on ALL of `G_T`** (zero excluded element!), attempts `1..B_max` are **100% EXACTLY i.i.d. conditional on `(FB, g', z')`**:
   - **Exact conditional i.i.d. amplification:** `Pr(all B_max fail | fixed FB, g', z') = Pr(attempt 1 fails | fixed FB)^10 <= Q_1^10` per log -> **`2 * Q_1^10 = 4.39e-05 <= 0.000044`** for BOTH logs uniformly over `(g', z')`!
   - **Exact conditional Wald's identity (`E[W_j | N_att >= j; fixed FB, g', z'] = E[W_1 | fixed FB] <= boot_own_exp + W_dag_1`):** Expected two-log work is `E[D | fixed FB, g', z'] <= 2 * ((1 - Q_1^10) / (1 - Q_1)) * (boot_own_exp + W_dag_1) = 1.5199x` single-attempt work = **`2^111.525535 <= 2^111.526`** core units (`2^111.509266` at `384/384`) uniformly over `(g', z')`.
   - **Hard Combined-Log Meter `t_desc = 2^124.2` (Markov):** Pointwise `min(D, t_desc) <= D` preserves `E[min(D, t_desc) | fixed FB, g', z'] <= E[D | fixed FB, g', z']` with **Markov meter-abort `P(D >= t_desc | fixed FB, g', z') <= E[D | fixed FB, g', z'] / t_desc = 0.00015302 <= 0.000154`**.
4. **Direct Production `(f_9, g_3)` `Q(zeta_15)` Audits (`prod_sq36_live` in `audit_cm.py`, ZERO Cubic Surrogate!):**
   - **36 ACTUAL Production BLS12-381 `(f_9, g_3)` `Q(zeta_15)` Special-`q` 16D Lattices (`24` control `73..258`b + `12` genuine parent-verified encountered children `71..255`b):** Every lattice is reduced with pure exact integer de Weger `lll_16d_int(99/100)` and verified LIVE on the ACTUAL production degree-72 `f_9` and degree-24 `g_3` polynomials:
     - **(a) Algebraic Containment, Sextile Envelope, & Modeled Survivor Ratio (`576` exact production norms):** Certifies `det % q == 0`, `q | N_side`, AND for all `12/12` encountered children LIVE parent divisibility **`q | N_side(parent_a, parent_b)`** with non-degenerate `||parent||_2 > q^(1/16)`. Every sextile `p16..max` sits **`23.87–42.83`b** below `qbins` on BOTH control and encountered cohorts; modeled single-ball survivor ratio vs `eval_descent_live` is **`0.0501` (`ctrl_g`), `0.2598` (`ctrl_f`), `0.0321` (`enc_g`), `0.0193` (`enc_f`)** (`max <= 1.0`).
     - **(b) LIVE `[7e]` Production Encountered-vs-Control Norm Check:** Comparing the `12` parent-verified encountered children directly against the `24` control lattices on the ACTUAL production `(f_9, g_3)` polynomials gives **`g`-side `(delta_cf, delta_cg) = (-2.36, -3.35)`b** and **`f`-side `(+0.35, -2.27)`b** — ZERO adverse norm inflation on encountered production ideals!
     - **(c) Direct Production `(f_9, g_3)` 4-Window `NB(1.0)` Smoothness-Window & Encountered/Control Check (`4,608` exact production norms = `36` lattices x `128` provably `P^1(Q(zeta_15))`-unique rank-2 points `a_1 + beta a_2, b_1 + beta b_2`):** Strips primes `<= 4096` from actual `N_f(a, b) / (q if f else 1)` and `N_g(a, b) / (q if g else 1)` across 4 windows (`12b..16b`):

     | Window | `n_zero` | `M_min` | `mean` | `enc/ctrl` | `var` / `NB(1.0)` | PGF `emp/nb10` (`0/.2/.5/.8`) |
     |---|---:|---:|---:|---:|---|---:|
     | `thr12b` | **`0/36`** | `4` | `13.00` | `0.8024` | `17.60` / `182.00` (`10.3x`) | `0 / .0006 / .0241 / .2957` |
     | `thr14b` | **`0/36`** | `3` | `8.25` | `0.7887` | `7.11` / `76.31` (`10.7x`) | `0 / .0026 / .0674 / .4965` |
     | `thr15b` | **`0/36`** | `3` | `7.03` | `0.9249` | `5.00` / `56.42` (`11.3x`) | `0 / .0050 / .0953 / .5633` |
     | `thr16b` | **`0/36`** | `2` | `5.44` | `0.8613` | `3.91` / `35.09` (`9.0x`) | `0 / .0236 / .1915 / .6794` |

   - **Exploratory Scope:** `prod_sq36_live` and `voronoi_and_boot_live` provide direct production `(f_9, g_3)` compatibility support (NOT statistical proof of `4e-4` production descent tails); production tail and expected-work transfer across `[2^70, 2^265)` rest explicitly on declared exploratory fixed-`FB` hypotheses `H5` and `H-EA`.

## 4. Complete Live-Computed Deterministic-DAG Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|p0|Phi_dag|Nf/Ng|own|W_det|C_det|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|76|70|70|0|0|16984|5.9e-5|0.000|841/816|110.6|110.6|1|
|[70,76)|g|sieve|80|70|70|0|0|25377|3.9e-5|0.000|916/748|111.7|111.7|1|
|[76,81)|f|sieve|68|76|76|1|0|45631|2.2e-5|0.000|846/820|108.0|112.3|3|
|[76,81)|g|sieve|68|76|76|1|0|43767|2.3e-5|0.000|922/744|108.0|112.3|3|
|[81,86)|f|sieve|72|76|76|1|0|19258|5.2e-5|0.000|870/829|109.3|112.4|3|
|[81,86)|g|sieve|76|76|76|1|0|29243|3.4e-5|0.000|956/750|110.6|112.6|3|
|[86,91)|f|sieve|68|76|76|2|0|20380|4.9e-5|0.000|881/835|110.5|113.5|5|
|[86,91)|g|sieve|72|76|76|2|0|29705|3.4e-5|0.000|973/751|111.2|113.6|5|
|[91,96)|f|sieve|76|76|76|2|0|14563|6.9e-5|0.000|910/846|111.4|113.6|5|
|[91,96)|g|sieve|80|76|76|2|0|19079|5.2e-5|0.000|1006/757|112.3|113.9|5|
|[96,101)|f|sieve|84|76|76|2|0|9031|1.1e-4|0.000|938/857|113.1|114.2|5|
|[96,101)|g|sieve|84|76|76|2|0|7229|1.4e-4|0.000|1034/761|113.0|114.1|5|
|[101,106)|f|sieve|84|76|76|2|0|2321|4.3e-4|0.001|956/864|112.9|114.1|5|
|[101,106)|g|sieve|80|81|76|1|1|4671|2.1e-4|0.001|1051/762|111.8|114.2|9|
|[106,111)|f|sieve|88|81|76|1|1|3432|2.9e-4|0.001|978/874|114.0|114.9|9|
|[106,111)|g|sieve|84|81|81|2|0|4850|2.1e-4|0.001|1079/766|112.9|114.8|13|
|[111,116)|f|sieve|76|81|76|2|1|1512|6.6e-4|0.001|980/876|115.2|115.8|11|
|[111,116)|g|sieve|80|81|76|2|1|1636|6.1e-4|0.001|1096/767|115.2|115.8|11|
|[116,121)|f|sieve|84|86|86|2|0|2627|3.8e-4|0.001|1008/887|112.9|114.9|13|
|[116,121)|g|sieve|80|96|83|1|1|4259|2.3e-4|0.001|1119/769|111.8|115.4|17|
|[121,126)|f|sieve|84|96|83|1|1|2613|3.8e-4|0.001|1026/894|112.9|115.5|17|
|[121,126)|g|sieve|84|101|86|1|1|4701|2.1e-4|0.001|1147/773|112.9|115.8|17|
|[126,131)|f|sieve|84|96|96|2|0|2715|3.7e-4|0.001|1043/902|112.9|115.9|21|
|[126,131)|g|sieve|88|116|100|2|0|5847|1.7e-4|0.001|1174/778|114.0|116.4|21|
|[131,136)|f|sieve|80|96|83|2|1|3662|2.7e-4|0.001|1055/908|114.1|116.1|23|
|[131,136)|g|sieve|84|96|83|2|1|2866|3.5e-4|0.001|1192/778|114.0|116.0|23|
|[136,141)|f|sieve|84|101|86|2|1|6430|1.6e-4|0.001|1078/917|114.8|116.4|23|
|[136,141)|g|sieve|88|101|86|2|1|4440|2.3e-4|0.001|1219/783|114.7|116.4|23|
|[141,146)|f|sieve|72|96|96|3|0|1120|8.9e-4|0.002|1080/919|113.2|116.5|31|
|[141,146)|g|sieve|88|101|86|2|1|1085|9.2e-4|0.002|1241/785|114.2|116.3|23|
|[146,151)|f|sieve|68|106|94|2|1|767|1.3e-3|0.004|1091/925|111.5|116.5|35|
|[146,151)|g|sieve|68|101|101|3|0|434|2.3e-3|0.004|1237/779|111.4|116.8|31|
|[151,156)|f|sieve|80|106|94|2|1|820|1.2e-3|0.003|1125/938|112.9|116.6|35|
|[151,156)|g|sieve|84|106|98|2|1|691|1.4e-3|0.004|1282/788|113.5|116.9|35|
|[156,161)|f|sieve|76|116|100|2|1|2086|4.8e-4|0.004|1138/943|113.4|117.6|45|
|[156,161)|g|sieve|92|106|106|3|0|1573|6.4e-4|0.004|1313/794|115.7|117.3|43|
|[161,166)|f|sieve|84|116|100|2|1|1295|7.7e-4|0.004|1166/954|113.9|117.6|45|
|[161,166)|g|sieve|92|116|106|2|1|1548|6.5e-4|0.005|1336/797|115.4|117.8|53|
|[166,171)|f|sieve|88|126|98|2|1|1355|7.4e-4|0.004|1188/964|114.4|117.7|55|
|[166,171)|g|sieve|84|131|100|2|1|806|1.2e-3|0.005|1349/796|113.2|117.8|63|
|[171,176)|f|sieve|88|126|107|2|1|1607|6.2e-4|0.006|1205/971|115.0|118.0|79|
|[171,176)|g|sieve|88|131|110|2|1|1391|7.2e-4|0.006|1376/800|114.5|118.2|87|
|[176,181)|f|sieve|76|141|106|2|1|649|1.5e-3|0.007|1208/973|112.4|118.0|75|
|[176,181)|g|sieve|64|131|100|3|1|500|2.0e-3|0.007|1366/792|114.6|118.1|73|
|[181,186)|f|sieve|68|146|108|2|1|179|5.6e-3|0.013|1214/977|110.8|118.2|99|
|[181,186)|g|sieve|68|131|100|3|1|204|4.9e-3|0.010|1395/796|113.9|118.1|73|
|[186,191)|f|sieve|72|151|110|2|1|212|4.7e-3|0.015|1237/987|111.7|118.3|111|
|[186,191)|g|sieve|60|141|106|3|1|219|4.5e-3|0.011|1404/794|112.8|118.2|89|
|[191,196)|f|sieve|72|161|116|2|1|396|2.5e-3|0.015|1255/994|113.3|119.2|137|
|[191,196)|g|sieve|72|141|106|3|1|253|3.9e-3|0.010|1446/803|114.0|118.2|89|
|[196,201)|f|sieve|80|166|118|2|1|644|1.5e-3|0.015|1283/1005|115.0|119.4|159|
|[196,201)|g|sieve|76|151|110|3|1|1085|9.2e-4|0.013|1474/807|116.0|118.8|133|
|[201,206)|f|sieve|84|171|120|2|1|679|1.5e-3|0.015|1306/1014|115.7|119.5|179|
|[201,206)|g|sieve|68|146|120|3|1|627|1.6e-3|0.013|1485/806|116.0|119.2|145|
|[206,211)|f|sieve|68|181|126|2|1|203|4.9e-3|0.023|1301/1015|114.8|119.6|235|
|[206,211)|g|sieve|64|161|116|3|1|153|6.5e-3|0.021|1501/807|112.8|119.5|161|
|[211,216)|f|sieve|80|181|126|2|1|204|4.9e-3|0.023|1335/1028|115.4|119.7|235|
|[211,216)|g|sieve|68|166|118|3|1|192|5.2e-3|0.021|1530/811|113.4|119.6|189|
|[216,221)|f|ecm|20|171|136|5|1|128916|7.8e-6|0.024|1209/987|124.2|124.3|349|
|[216,221)|g|sieve|72|171|120|3|1|237|4.2e-3|0.020|1558/816|114.1|119.6|209|
|[221,226)|f|ecm|20|176|140|5|1|203541|4.9e-6|0.026|1227/995|125.1|125.1|397|
|[221,226)|g|sieve|76|181|126|3|1|866|1.2e-3|0.022|1586/820|116.2|119.9|269|
|[226,231)|f|ecm|20|171|149|5|1|595310|1.7e-6|0.048|1244/1002|126.6|126.6|449|
|[226,231)|g|sieve|88|181|126|3|1|664|1.5e-3|0.022|1624/828|116.9|120.0|269|
|[231,236)|f|ecm|20|176|140|6|1|866283|1.2e-6|0.029|1262/1010|125.7|125.7|443|
|[231,236)|g|sieve|72|176|140|3|1|595|1.7e-3|0.022|1626/823|117.4|120.3|305|
|[236,241)|f|ecm|20|176|140|6|1|350349|2.9e-6|0.029|1279/1017|125.7|125.7|443|
|[236,241)|g|sieve|48|191|149|3|1|186|5.4e-3|0.055|1606/812|114.8|120.3|399|
|[241,246)|f|ecm|20|171|149|6|1|1372965|7.3e-7|0.056|1297/1025|127.3|127.4|515|
|[241,246)|g|sieve|56|191|149|3|1|180|5.5e-3|0.055|1644/819|115.5|120.3|399|
|[246,251)|f|ecm|20|191|149|6|1|3991958|2.5e-7|0.072|1314/1032|127.4|127.4|597|
|[246,251)|g|sieve|68|191|149|3|1|227|4.4e-3|0.054|1687/829|116.8|120.4|399|
|[251,256)|f|ecm|20|191|149|6|1|1698035|5.9e-7|0.072|1332/1040|127.4|127.4|597|
|[251,256)|g|sieve|68|196|153|3|1|184|5.4e-3|0.055|1710/831|116.8|120.7|437|
|[256,261)|f|ecm|20|191|149|6|1|696640|1.4e-6|0.072|1349/1047|127.4|127.4|597|
|[256,261)|g|sieve|72|201|156|3|1|244|4.1e-3|0.055|1738/836|117.6|120.9|503|
|[261,266)|f|ecm|20|186|161|6|1|4323408|2.3e-7|0.071|1367/1055|129.6|129.6|715|
|[261,266)|g|sieve|76|206|159|3|1|310|3.2e-3|0.056|1766/840|118.7|121.4|589|
