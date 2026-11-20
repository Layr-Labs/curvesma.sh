# Zero-Retry Deterministic-DAG Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-`q` nodes**, **`boot_rels = 1`**, **`k = 1` (ZERO retries at every node!)**, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`, explicit `D_skew` boot LLL (`S_boot = 2^21`), and an explicit aborting combined logarithm phase meter **`t = 2^121.4`** core units (`P(D >= t) <= E[D] / t` by Markov).

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, empirical `q`-lattice offsets, the 24-base production `(4, 3)` `Q(zeta_15)` special-`q` witness set (`prod_sq24`), the 200-ideal asymmetric `(4, 3)` `Q(zeta_12)` campaign (`q12_43_campaign`), and boot config (`m_sv = 2560.0`b, `s_boot = 2^21`). `audit_cm.py` calls `cm_models.eval_descent_live` to compute seeded peeling closure `(h, exc_crit)` and every node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `p0`, `Phi_dag`, `W_det`, `C_det`, and boot work/cap 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`E[D] = 2^108.267402 <= 2^108.268`, `2*Phi_dag = 0.095083 <= 0.09509`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`E[D] = 2^108.265147`, `2*Phi_dag = 0.090189`).

## 1. Fresh-Coin Zero-Retry i.i.d. Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = diag(1..1, 2^21..2^21)`)

For each discrete logarithm (`log g` with `g = e(P0, Q)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy, Live `F_(p^12) -> F_(p^6)` Subfield Elimination, & Exact `D_skew` 16D LLL
In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). Any `(a, b) in L_Guillevic` with `b != 0 (mod P)` gives **`Psi_g(a - b X) = -(b mod P) w_1^(-1) c_6 z'`** with exact relative norm **`N_(K_g/K_H)(a - b X) = a^3 - omega^2 u b^3`**.
- **Exact `D_skew`-Weighted 16D Integer LLL (`S_boot = 2^21`, `delta = 0.99`):** Multiply the 8 `b`-coordinates by **`S_boot = 2^21 = 2,097,152`**, forming `D_skew * L_Guillevic` (`det = p^4 S_boot^8`), run exact de Weger integer LLL (`cm_models.lll_16d_int`, `delta = 0.99`), and unscale `b = b_scaled / S_boot`. Every reduced basis satisfies the exact `D_skew`-weighted Lovasz bound **`||D_skew b_1||_2 <= (0.99 - 0.25)^(-15/4) (p^4 S_boot^8)^(1/16) = 2^107.30`**.
- **Single-Relation i.i.d. Boot (`boot_rels = 1`, ZERO boot retries, ZERO `BALL_PAIR`):**
  1. For trial `t = 1..T_boot_max`, draw fresh uniform `e_t in {1..r-1}`, form `z_t = z' * g'^(e_t)` (`< 2^25` Ir, absorbed in `LLL_NODE = 2^38`), run `F_(p^12) -> F_(p^6)` elimination to obtain `gamma_t`, run `lll_16d_int` on `D_skew * L_Guillevic(gamma_t)`, and take the **single shortest vector `b_1 = (a_1, b_1)`** (`cells_b = 1`).
  2. Apply the deterministic norm gate **`|N_(K_g/Q)(a_1^3 - omega^2 u b_1^3)| <= 2^2560`** (`m_sv = 2560.0`b). **`RET_BOOT = 0.50`** is declared as a production-law model inequality supported (not proved) by 20 genuine production `G_T` challenges below (`20/20 <= 2543.21`b, Clopper-Pearson 95% LB `0.8609 >> 0.50`, plus `32/32 <= 2553.16`b in `[11]`).
  3. Run worst-case multi-factor ECM + limb-scaled 100-round MR on `b_1`: `k_lx = 129` primes in `(2^20, 2^70]`, `kA = 13` in `(2^70, 2^234.53125]`, `1` in `(n_A, 2^265]` (`mr_instr(2560) = 1.21e10` Ir, `det_b = 0.999^143 = 0.8667`).
  4. Take the **FIRST** qualifying boot relation (`boot_rels = 1`), cap trials at **`T_boot_max = ceil(ln(10^4) / p_boot_trial)`** (`<= 10^(-4)`/log, `<= 2e-4` for 2 logs), expected two-log boot work **`2^102.8155`** core units.
- **Live 20-Challenge Genuine Production `G_T` Verification (`audit_cm.py`: `voronoi_and_boot_live` + `boot_witness.json`):** (a) Verifies genuine production BLS12-381 `G_T` generator `G_0 != 1, G_0^r == 1` in `F_(p^12)`, generates 20 deterministic seeded exponents `x_t in [1, r-1]` LIVE, and executes exact **`F_(p^12) -> F_(p^6)` subfield elimination LIVE** on all 20 challenges `z_t = G_0^(x_t)` to recover `gamma_t` (zero precomputed `gamma_t` stored!); (b) runs exact de Weger integer `lll_16d_int(D_skew * B_hnf, delta=0.99)` **LIVE from scratch on Target 0**, certifying `det^2 == p^8 S_boot^16`, all 15 Lovasz `0.99` & 120 size-reduction checks, and bit-for-bit match vs `boot_witness.json`; (c) verifies `a_1 == gamma_t b_1 (mod P)`, `||D_skew b_1||_2 <= 2^106.14 <= 2^107.30`, and **`20/20` norms `2532.70..2543.21 <= m_sv = 2560.0`b** (mean `2537.64`b, Clopper-Pearson 95% LB **`0.8609 >> RET_BOOT = 0.50`**; `|N(det)| = p^4` and 17-cell deciles `4/4` on `full4`), plus constructive 16D LLL on `p_0 = 7207` (`k = 2, index = 1`).

## 2. Exact Seeded Peeling Closure (`h = 0.256845`), Zero-Column-Loss 2-Core `reconstructlog`, & `<= LX = 2^70` Discard

- **Seeded Peeling Fixed Point (`cm_models._peel_fixed_point_bls`, live in `audit_cm.py`):** On the two-sided BLS12-381 factor base (`B_f = 2^71.29`, `B_g = 2^71.69`, `lrels = log2(1.09 U') = 63.441701`, `u_f = 548.9 / 71.29`, `u_g = 677.8 / 71.69`), ideals with `lam >= 100.0` form the protected seed. Unprotected ideals satisfy `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))`, giving **`h = 0.256845`** and **`exc_crit = 0.7625 << exc = 1.09`** (**`+43.0%` margin** above peeling collapse).
- **Exact Zero-Column-Loss Post-BW `reconstructlog` (`cols_lost_after_2core == 0`, `RANK-BRIDGE.md` §3):** After BW solves the `N` active core columns, 4-stage post-BW back-substitution (`L_rec` -> `merge_stack` -> unit-propagation `reconstructlog` over the retained post-singleton-peeling 2-core relations -> tiny residual cycle solve, `< 2^64` core units absorbed inside `filtering_and_misc = 2^80`) recovers exact `vlog` for **100%** of post-singleton-peeling 2-core columns (`filtsim_and_hx_live` certifies `cols_lost_after_2core == 0` and `backsub_2core_exact == True` on `6/6` seeds `log2U in 12..14`). Thus `miss_side`'s post-peeling formula has **zero unmodeled post-trim column loss**!
- **Exact `<= LX = 2^70` Leaf Discard & Conservative `> 2^70` Children:** Leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = kB = 0`). Any prime `<= 2^70` missing from the solved 2-core discards the candidate cleanly (`exp(-miss)` over `[40, 70)`, charging `exp(-x)(1 + x)` at `x = lam(1 - h)` for `lam < 100` and `exp(-lam)` for `lam >= 100`). Every prime **`> 2^70`** (including `(2^70, B_s]`) is conservatively counted as an unknown Window-A/B child in `ss_eval` and descended strictly downward.

## 3. Zero-Retry Deterministic-DAG (`k = 1`), Ideal-Log Cache, `t = 2^121.4` Meter, & Exploratory `H5` / `H-EA` Evidence

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`, and enumerates `alpha v_1 + beta v_2` (`UDUP_NODE = 0.95`; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b; `[7e]` shift `(+3.40, +2.35)`b inside `(+10, +4)`b).

1. **Zero-Retry Rule (`k = 1`):** Every reached special-`q` node `v` sieves/ECMs its fixed ball ONCE (charging **100% expected full survivor work** `cells * sv + cells * pr` under `H-EA` + worst-case multi-factor ECM + limb-scaled 100-round MR) and takes the **FIRST** qualifying relation (`M_v = 0 =>` abort; zero retries anywhere).
2. **Strictly Decreasing Norm DAG + Ideal-Log Cache (`P(repeat abort) = 0`):** Every child norm is strictly smaller than its parent; an ideal-log hash table reuses any already-descended `(side, q') -> vlog(q')` (zero cycles, zero repeat-ideal aborts).
3. **Deterministic Max Node Count `C_det`, Expected Work `E[D]`, & Hard Combined-Log Meter `t = 2^121.4` (Markov):** `C_det_1log = 4,554` nodes/log deterministically (`P(node cap) = 0`). Under `H-EA`, expected DAG subtree work is `2^108.2341` u and expected two-log work is **`E[D] <= 2^108.267402 <= 2^108.268`** core units (`2^108.265147` at `384/384`). Wrapping both logs in an aborting meter **`t = 2^121.4`** preserves `E[min(D, t)] <= E[D]` with **Markov meter-abort `P(D >= t) <= E[D] / t = 0.00011140 <= 0.000112`** (zero independence needed).
4. **Unconditional Additive Boole Union Bound under Reached-Ideal Ensemble `H5`:** Under declared reached-ideal ensemble production-law model inequality **`H5`** (`Pr(M(q_sigma) = 0 | slot sigma reached) <= p0(mu_model) = (1 + 0.25 mu_model)^(-2)`, part of `E[t^M | slot sigma reached] <= (1 + 0.25 mu_model (1 - t))^(-2)` for `t in [0, 0.8]`, `NB2` with `r = 2`), Boole's inequality holds unconditionally for arbitrary slot dependence: **`2 * Phi_1log = 0.095083 <= 0.09509`** (`0.090189` at `384/384`, `min_mu = 196.31`).
5. **Two Scoped Empirical Compatibility Checks (`q12_43_descent_live` + `prod_sq24_live` in `audit_cm.py`):**
   - **(a) Structurally Matching Asymmetric Degree-`(4, 3)` `Q(zeta_12)` Two-Window `parent -> k=1 child` Campaign (`200` ideals; `12` chains + `12` controls LIVE):** Irreducible `f_4(X) = X^4 + X^2 + 2` (deg 16), `g_3(X) = X^3 + w^2` (deg 12), `LX = 2^16`, child bin `[18.0, 19.5)`, two-window `(nA = 17.2, nc = 18.0, kA = 1, kB = 1, T = 35.2)`, `ball_pts = 569`, first-principles Murphy `alpha = (+3.55, +3.04)`, pre-calibrated **`qbins`/`qbins_sub`/`qbins_lo` norm-offset envelope** (NEVER per-cell exact norms!: `mu_env = 3.98 g / 7.06 f`, two-sided `sv_upper = 548.69 g / 576.06 f`). `audit_cm.py` (`q12_43_descent_live`) reconstructs **12 full `parent -> k=1 first qualifying relation -> genuine child` chains + 12 controls LIVE** (exact bit-for-bit `M` and two-sided `sv_2` match) and audits all 200 ideals (`100` slot-reached + `100` matched-pair fresh controls):

     | Cohort | Ideals | `M` range | Empty | `nb2_p0` | `M/mu_env` | `var(M/mu)` vs `NB(2)` | `sv2 / sv_upper` (max) | PGF `emp/nb2` (`0 / 0.2 / 0.5 / 0.8`) |
     |---|---:|---|---:|---:|---:|---|---:|---:|
     | `ALL_200` | 200 | `1..29` | **`0/200`** | `0.19085` | `1.7381` | `0.5407` vs `1.852` | **`0.2798`** (`0.5729`) | `0 / 0.0139 / 0.0651 / 0.2950` |
     | `SLOT_100` | 100 | `1..29` | **`0/100`** | `0.19085` | `1.7532` | `0.5621` vs `1.869` | **`0.2822`** (`0.5729`) | `0 / 0.0145 / 0.0654 / 0.2891` |
     | `CTRL_100` | 100 | `1..29` | **`0/100`** | `0.19085` | `1.7230` | `0.5236` vs `1.835` | **`0.2774`** (`0.5381`) | `0 / 0.0134 / 0.0648 / 0.3009` |

     Slot-vs-control ratios: `M/mu` = **`1.0175`** (zero slot degradation), actual two-sided residual survivor upper-bound check `sv2 / sv_upper = 0.2798 < 1.0` (`max = 0.5729 < 1.0`, slot/ctrl = **`1.0173`**, directly supporting `H-EA`!).
   - **(b) 24-Base Actual Production BLS12-381 `(4, 3)` `Q(zeta_15)` `(f_9, g_3)` Special-`q` 16D Lattice Check (`audit_cm.py`: `prod_sq24_live`):** Across `12 g-side + 12 f-side` production primes spanning `73..258` bits (`384` exact production deg-72/deg-24 norms): `det(v1, v2) % q == 0` and exact `q | N` verified LIVE; all 6 empirical sextiles `p16..max` lie strictly below `qbins` `c_hi + cf, c_hi + cg` (`g-side` headroom `>= 39.68`b `f` / `>= 23.87`b `g`; `f-side` `>= 37.39`b `f` / `>= 27.31`b `g`); single-ball modeled survivor ratio `sv_emp / sv_plan` has mean **`0.0503` (`g`, max `0.1226`)** and **`0.2602` (`f`, max `1.0000`)**.
   - **Explicit Exploratory Scope:** Together, `0/200` empties + `var(M/mu) = 0.5407` vs `NB(2) 1.852` + two-sided `sv2/sv_upper = 0.2798 < 1.0` + the 24-base production `(4, 3)` norm-envelope check provide scoped finite-scale compatibility evidence (NOT statistical proof of `4e-4` production tails); production tail and survivor-work transfer rests explicitly on declared exploratory model inequalities `H5` and `H-EA`.

## 4. Complete Live-Computed Deterministic-DAG Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|p0|Phi_dag|Nf/Ng|own|W_det|C_det|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|64|70|70.0|0|0|3769.2|1.1e-6|0.00000|823/810|106.59|106.59|1|
|[70,76)|g|sieve|64|70|70.0|0|0|3755.2|1.1e-6|0.00000|893/740|106.59|106.59|1|
|[76,81)|f|sieve|56|76|76.0|1|0|8460.1|2.2e-7|0.00000|826/813|103.59|107.67|3|
|[76,81)|g|sieve|56|76|76.0|1|0|8316.5|2.3e-7|0.00000|902/737|103.60|107.68|3|
|[81,86)|f|sieve|56|76|76.0|1|0|2252.2|3.1e-6|0.00001|843/820|103.53|107.67|3|
|[81,86)|g|sieve|64|76|76.0|1|0|6523.0|3.8e-7|0.00000|938/744|106.60|108.18|3|
|[86,91)|f|sieve|56|76|76.0|2|0|3739.0|1.1e-6|0.00001|861/828|107.35|109.10|5|
|[86,91)|g|sieve|48|76|76.0|2|0|886.9|2.0e-5|0.00002|931/737|104.92|108.70|5|
|[91,96)|f|sieve|56|76|76.0|2|0|1008.9|1.6e-5|0.00002|878/835|106.01|108.81|5|
|[91,96)|g|sieve|56|76|76.0|2|0|895.8|2.0e-5|0.00002|969/744|105.88|108.79|5|
|[96,101)|f|sieve|64|76|76.0|2|0|840.3|2.2e-5|0.00003|910/848|107.43|109.12|5|
|[96,101)|g|sieve|64|76|76.0|2|0|700.7|3.2e-5|0.00004|1006/752|107.26|109.07|5|
|[101,106)|f|sieve|80|76|76.0|2|0|1551.2|6.6e-6|0.00001|950/863|111.83|111.97|5|
|[101,106)|g|sieve|64|81|75.5|1|1|684.5|3.4e-5|0.00004|1028/754|106.97|109.50|9|
|[106,111)|f|sieve|80|81|75.5|1|1|1505.8|7.0e-6|0.00001|968/870|111.79|112.02|9|
|[106,111)|g|sieve|64|81|81.0|2|0|473.3|7.0e-5|0.00008|1051/757|106.81|109.86|13|
|[111,116)|f|sieve|64|81|75.5|2|1|327.8|1.5e-4|0.00015|962/870|112.64|112.81|11|
|[111,116)|g|sieve|64|81|75.5|2|1|235.3|2.8e-4|0.00029|1073/759|111.80|112.09|11|
|[116,121)|f|sieve|80|86|86.0|2|0|1665.9|5.7e-6|0.00002|1003/885|111.76|112.12|13|
|[116,121)|g|sieve|80|96|83.0|1|1|4259.2|8.8e-7|0.00005|1119/769|111.77|112.28|17|
|[121,126)|f|sieve|80|96|83.0|1|1|1750.5|5.2e-6|0.00006|1020/893|111.75|112.27|17|
|[121,126)|g|sieve|80|101|85.5|1|1|3161.3|1.6e-6|0.00007|1141/772|111.76|112.30|17|
|[126,131)|f|sieve|80|96|96.0|2|0|1726.1|5.3e-6|0.00009|1038/900|111.76|112.41|21|
|[126,131)|g|sieve|80|116|100.2|2|0|2629.5|2.3e-6|0.00013|1164/774|111.77|112.47|21|
|[131,136)|f|sieve|64|96|83.0|2|1|513.3|6.0e-5|0.00012|1032/900|110.73|111.84|23|
|[131,136)|g|sieve|64|96|83.0|2|1|280.5|2.0e-4|0.00026|1163/769|109.50|111.40|23|
|[136,141)|f|sieve|64|101|85.5|2|1|598.8|4.4e-5|0.00012|1050/908|110.53|111.80|23|
|[136,141)|g|sieve|64|101|85.5|2|1|301.2|1.7e-4|0.00025|1186/772|109.18|111.38|23|
|[141,146)|f|sieve|64|96|96.0|3|0|402.1|9.7e-5|0.00023|1067/915|111.42|112.48|31|
|[141,146)|g|sieve|80|101|85.5|2|1|504.0|6.2e-5|0.00014|1231/782|112.21|112.73|23|
|[146,151)|f|sieve|64|106|93.6|2|1|446.8|7.9e-5|0.00024|1085/923|110.47|112.99|35|
|[146,151)|g|sieve|80|101|101.0|3|0|1662.1|5.8e-6|0.00020|1254/784|113.91|114.19|31|
|[151,156)|f|sieve|80|106|93.6|2|1|820.1|2.4e-5|0.00018|1125/938|112.86|113.79|35|
|[151,156)|g|sieve|80|106|98.1|2|1|467.6|7.2e-5|0.00027|1276/787|112.50|113.64|35|
|[156,161)|f|sieve|64|116|100.2|2|1|477.7|6.9e-5|0.00064|1120/938|110.79|113.92|45|
|[156,161)|g|sieve|80|106|106.0|3|0|486.2|6.7e-5|0.00027|1299/789|113.16|114.52|43|
|[161,166)|f|sieve|80|116|100.2|2|1|865.9|2.1e-5|0.00059|1160/953|113.04|114.43|45|
|[161,166)|g|sieve|80|116|105.9|2|1|496.2|6.4e-5|0.00064|1321/792|112.67|114.75|53|
|[166,171)|f|sieve|80|126|98.0|2|1|601.0|4.4e-5|0.00061|1178/960|112.56|114.32|55|
|[166,171)|g|sieve|80|131|100.5|2|1|520.5|5.8e-5|0.00063|1344/794|112.11|114.25|63|
|[171,176)|f|sieve|80|126|106.8|2|1|742.4|2.9e-5|0.00069|1195/968|113.31|114.99|79|
|[171,176)|g|sieve|80|131|110.0|2|1|648.9|3.8e-5|0.00070|1366/797|112.61|114.84|87|
|[176,181)|f|sieve|80|141|105.5|2|1|991.0|1.6e-5|0.00060|1213/975|113.29|114.98|75|
|[176,181)|g|sieve|64|131|100.5|3|1|500.5|6.3e-5|0.00070|1366/792|114.61|115.33|73|
|[181,186)|f|sieve|80|146|108.0|2|1|742.0|2.9e-5|0.00076|1230/983|113.40|115.09|99|
|[181,186)|g|sieve|80|131|100.5|3|1|846.4|2.2e-5|0.00066|1411/802|116.51|116.74|73|
|[186,191)|f|sieve|80|151|110.5|2|1|501.7|6.3e-5|0.00080|1248/990|113.44|115.55|111|
|[186,191)|g|sieve|64|141|105.5|3|1|378.3|1.1e-4|0.00076|1411/797|113.84|115.35|89|
|[191,196)|f|sieve|80|161|115.5|2|1|986.9|1.6e-5|0.00183|1265/998|114.80|116.47|137|
|[191,196)|g|sieve|80|141|105.5|3|1|626.9|4.0e-5|0.00069|1456/807|115.68|116.28|89|
|[196,201)|f|sieve|80|166|118.0|2|1|644.4|3.8e-5|0.00221|1283/1005|114.96|116.69|159|
|[196,201)|g|sieve|64|151|110.5|3|1|251.7|2.4e-4|0.00109|1456/802|113.35|115.68|133|
|[201,206)|f|sieve|80|171|120.5|2|1|456.1|7.6e-5|0.00225|1300/1013|115.00|116.70|179|
|[201,206)|g|sieve|64|146|119.9|3|1|380.7|1.1e-4|0.00196|1478/804|115.05|116.36|145|
|[206,211)|f|sieve|80|181|125.5|2|1|771.9|2.7e-5|0.00230|1318/1020|117.09|117.86|235|
|[206,211)|g|sieve|80|161|115.5|3|1|1002.9|1.6e-5|0.00228|1524/814|116.22|117.20|161|
|[211,216)|f|sieve|80|181|125.5|2|1|204.5|3.7e-4|0.00264|1335/1028|115.43|117.12|235|
|[211,216)|g|sieve|80|166|118.0|3|1|791.3|2.5e-5|0.00264|1546/817|115.92|117.18|189|
|[216,221)|f|ecm|20|171|136.3|5|1|128915.8|9.6e-10|0.00350|1209/987|124.24|124.24|349|
|[216,221)|g|sieve|80|171|120.5|3|1|566.3|4.9e-5|0.00266|1569/819|115.70|117.09|209|
|[221,226)|f|ecm|20|176|139.6|5|1|203540.6|3.9e-10|0.00360|1227/995|125.07|125.08|397|
|[221,226)|g|sieve|64|181|125.5|3|1|196.3|4.0e-4|0.00311|1568/814|113.53|116.90|269|
|[226,231)|f|ecm|20|171|148.9|5|1|595310.3|4.5e-11|0.00391|1244/1002|126.59|126.59|449|
|[226,231)|g|sieve|80|181|125.5|3|1|298.8|1.7e-4|0.00289|1614/824|115.35|117.22|269|
|[231,236)|f|ecm|20|176|139.6|6|1|866282.7|2.1e-11|0.00404|1262/1010|125.67|125.68|443|
|[231,236)|g|sieve|64|176|139.6|3|1|216.5|3.3e-4|0.00304|1613/819|115.60|117.20|305|
|[236,241)|f|ecm|20|176|139.6|6|1|350349.0|1.3e-10|0.00404|1279/1017|125.67|125.68|443|
|[236,241)|g|sieve|56|191|149.4|3|1|722.9|3.0e-5|0.00316|1622/817|117.20|118.57|399|
|[241,246)|f|ecm|20|171|148.9|6|1|1372964.7|8.5e-12|0.00443|1297/1025|127.35|127.35|515|
|[241,246)|g|sieve|64|191|149.4|3|1|564.8|4.9e-5|0.00318|1658/824|117.50|118.69|399|
|[246,251)|f|ecm|20|191|149.4|6|1|3991958.1|1.0e-12|0.00471|1314/1032|127.43|127.43|597|
|[246,251)|g|sieve|80|191|149.4|3|1|872.7|2.1e-5|0.00315|1704/834|119.20|119.68|399|
|[251,256)|f|ecm|20|191|149.4|6|1|1698034.8|5.5e-12|0.00471|1332/1040|127.43|127.43|597|
|[251,256)|g|sieve|80|196|152.7|3|1|750.0|2.8e-5|0.00419|1726/837|119.19|119.78|437|
|[256,261)|f|ecm|20|191|149.4|6|1|696640.2|3.3e-11|0.00471|1349/1047|127.43|127.43|597|
|[256,261)|g|sieve|80|201|156.0|3|1|584.1|4.6e-5|0.00492|1749/839|119.18|119.79|503|
|[261,266)|f|ecm|20|186|160.6|6|1|4323408.0|8.6e-13|0.00703|1367/1055|129.57|129.57|715|
|[261,266)|g|sieve|80|206|159.2|3|1|479.9|6.8e-5|0.00706|1771/842|119.41|119.97|589|
