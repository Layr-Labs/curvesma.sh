# Zero-Retry Deterministic-DAG Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-`q` nodes**, **`boot_rels = 1`**, **`k = 1` (ZERO retries at every node!)**, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`.

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, the empirical `q`-lattice offsets, and the boot config (`m_sv = 2560.0`b). `audit_cm.py` invokes `cm_models.eval_descent_live` to compute seeded peeling closure `(h, exc_crit)` and every node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `p0`, `Phi_dag`, `W_det`, `C_det`, and boot work/cap 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`exp = 2^108.267402 <= 2^108.268`, `2*Phi_dag = 0.095083 <= 0.09509`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`exp = 2^108.265147`, `2*Phi_dag = 0.090189`).

## 1. Fresh-Coin Zero-Retry i.i.d. Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `ball_pair = 0`)

For each discrete logarithm (`log g` with `g = e(P0, Q)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy & i.i.d. Single-Cell Boot Trials
In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). Any `(a, b) in L_Guillevic` with `b != 0 (mod P)` gives **`Psi_g(a - b X) = -(b mod P) w_1^(-1) c_6 z'`** with exact relative norm **`N_(K_g/K_H)(a - b X) = a^3 - omega^2 u b^3`**.
- **Single-Relation i.i.d. Boot (`boot_rels = 1`, ZERO boot retries, ZERO `BALL_PAIR`, ZERO 17-cell dependence):**
  1. For trial `t = 1, 2, ..., T_boot_max`, draw a **fresh independent uniform `e_t in {1, ..., r - 1}`**, compute `z_t = z' * g'^(e_t)` (`GT` exp `< 2^25` Ir, absorbed into `LLL_NODE = 2^38` Ir alongside 16D LLL `< 2^28` Ir), run 16D integer LLL on `L_Guillevic(gamma_t)`, and take the **single shortest vector `b_1 = (a_1, b_1)`** (`cells_b = 1`).
  2. Apply the 1-line deterministic norm gate **`|N_(K_g/Q)(a_1^3 - omega^2 u b_1^3)| <= 2^2560`** (`m_sv = 2560.0`b; passes `32/32` targets in `[11]` with `sv_mean = 2540.82`b, `sv_max = 2553.16`b, and `4/4` live BLS12-381 witnesses in `boot_witness.json` with `b_1` norms `2534.44..2542.66`b; conservatively charged `RET_BOOT = 0.50` = `2x` trials!).
  3. Run the full worst-case multi-factor ECM + limb-scaled 100-round MR extraction loop on `b_1`: `k_lx = 129` primes in `(2^20, 2^70]`, `kA = 13` primes in `(2^70, n_A = 2^234.53125]`, and `1` prime in `(n_A, n_1 = 2^265]`, charging `LLL_NODE + 129*(ecm(70)+mr) + 13*(ecm(nA)+mr) + mr + IDEAL` (`mr_instr(2560) = 1.21e10` Ir, `det_b = 0.999^143 = 0.8667`).
  4. **Take the FIRST qualifying boot relation (`boot_rels = 1`, ZERO boot retries!)** and descend its `<= 14` children once via the deterministic DAG below. Cap boot trials per log at **`T_boot_max = ceil(ln(10^4) / p_boot_trial)`**: because consecutive trials are strictly i.i.d., failure to find a qualifying boot relation within `T_boot_max` trials is **`(1 - p_boot_trial)^T_boot_max <= 10^(-4)`** per log (**`<= 2 x 10^(-4)`** for both logs!). Expected two-log boot work is **`2^102.8155`** core units; deterministic every-run two-log boot cap is **`2^106.0188`** core units!
- **Live verification (`audit_cm.py`: `voronoi_and_boot_live`):** (a) Runs constructive 16D integer LLL + `(v_1 = b_1, v_2 = b_k)` extraction live from unreduced 16 x 16 HNF on split-base `p_0 = 7207 == 7 (mod 15)` (`k = 2, index = 1`); (b) verifies 4 bundled 381-bit BLS12-381 `w'=omega^2` bases (`boot_witness.json`): `a_i == gamma b_i (mod P)`, `|N(det)| = p^4`, **`4/4` single-cell `b_1` norms `<= 2542.66 < m_sv = 2560.0`b**.

## 2. Exact Seeded Peeling Closure (`h = 0.256845`) & Leaf-Domain `<= LX = 2^70` Discard

- **Seeded Peeling Fixed Point (`cm_models._peel_fixed_point_bls`, verified live in `audit_cm.py`):** On the two-sided BLS12-381 factor base (`B_f = 2^71.29`, `B_g = 2^71.69`, `lrels = log2(1.09 U') = 63.441701`, `u_f = 548.9 / 71.29`, `u_g = 677.8 / 71.69`), ideals with initial expected degree `lam >= 100.0` form the protected seed. Unprotected ideals (`lam < 100.0`) satisfy the strictly conservative peeling fixed-point recurrence `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))`, using upper-endpoint divisor density `wt = km(c) ln(c/b)` and lower-endpoint degree `lam = 2^(lrels - c) * 12 * km(b)` on `[b, c] = [b, b + 0.1]`. Live solution: **`h = 0.256845`** and criticality threshold **`exc_crit = 0.7627 << exc = 1.09`** (**`+42.9%` margin** above peeling collapse!).
- **Exact `<= LX = 2^70` Leaf Discard & Conservative `> 2^70` Window-A/B Treatment:** Descent's leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = 0, kB = 0`, requiring every prime strictly `<= 2^70`). Any prime `<= 2^70` missing from the solved FB causes the boot or special-`q` candidate to be **discarded cleanly** (`exp(-miss)` over `[40, 70)`, where unprotected `lam < 100` charges the conservative surviving-degree-`<= 1` probability `exp(-x)(1 + x)` with `x = lam(1 - h)`, strictly larger than `exp(-x)`, and protected `lam >= 100` charges `exp(-lam)`; below `2^40`, `lam >= 8.77e8 => exp(-lam) = 0`). Every prime **`> 2^70`** — **whether in `(2^70, B_s]` or `> B_s`** — is conservatively counted as an unknown Window-A/B child in `ss_eval` and descended strictly downward (zero credit taken for solved FB ideals in `(2^70, B_s]`).

## 3. Zero-Retry Deterministic-DAG Theorem (`k = 1`), Ideal-Log Cache, & Unconditional Boole Union Bound

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`, and enumerates `alpha v_1 + beta v_2` (`UDUP_NODE = 0.95` charged; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b across 56 BLS12-381 lattices; `[7e]` shift `(+3.40, +2.35)`b inside `(+10, +4)`b).

1. **Zero-Retry Rule (`k = 1`):** At every reached special-`q` node `v`, sieve/ECM its fixed lattice ball ONCE (charging **100% full survivor work** + full worst-case multi-factor ECM + limb-scaled 100-round MR) and take the **FIRST** qualifying relation. If `M_v = 0`, abort the run as failure. **Zero retries occur at ANY node or boot.**
2. **Strictly Decreasing Norm DAG + Ideal-Log Cache (`P(repeat abort) = 0`):** Across all 78 rows, every child ideal has strictly smaller norm `q_child < q_parent`, so the visited ideals form a strictly acyclic DAG terminating at `<= 2^70`. Maintain a hash table mapping every already-descended ideal `(side, q') -> vlog(q')`; whenever a child ideal `q'` was already descended earlier in the DAG, simply reuse its cached `vlog(q')`! **Cycles are impossible, repeat-ideal aborts are `0` identically, and NO `COLL_PEAK` collision heuristic exists.**
3. **Deterministic Worst-Case Node Count & DAG Work (ZERO Counter Abort, ZERO DAG Markov Cap Abort):** Because each node expands at most `kA` Window-A children and `kB` Window-B children on each side (at most ONCE!), worst-case node count and DAG work satisfy exact deterministic recurrences:
   - `W_det[bn, side] = own + kA * (wmax_f(nA) + wmax_g(nA)) + kB * (wmax_f(n) + wmax_g(n))`,
   - `C_det[bn, side] = 1 + kA * (cmax_f(nA) + cmax_g(nA)) + kB * (cmax_f(n) + cmax_g(n))`.
   Live evaluation gives **`C_det_1log = 4,554` nodes/log** (deterministic upper bound!), two-log DAG subtree work **`2^108.2341`** core units, expected two-log descent **`2^108.267402 <= 2^108.268`** core units (`2^108.265147` at `384/384`), and deterministic every-run two-log cap **`2^108.515426 <= 2^108.516`** core units!
4. **Unconditional Additive Boole Union Bound (ZERO Subtree Independence, ZERO Sequential Retry Hazard!):**
   Let the deterministic DAG template have slots `sigma` (at most `4,554` slots per log). Under `H5` — which is now **SOLELY** the single-node marginal empty-ball heuristic `Pr(M(q_sigma) = 0 | slot sigma reached) <= p0(mu) = (1 + 0.25 mu)^(-2)` (`NB2`, `r = 2, alpha = 0.5`) — Boole's inequality holds **unconditionally for arbitrary dependence between DAG slots**:
   - `Phi_dag[bn, side] = p0 + kA * (pmax_f(nA) + pmax_g(nA)) + kB * (pmax_f(n) + pmax_g(n))`,
   - `Phi_1log = kA_b * pmax_g(nA_b) + pmax_g(n1_b)`.
   Live evaluation (`min_mu = 196.31`, `232.62` at `384/384`) gives **`2 * Phi_1log = 0.095083 <= 0.09509`** (`0.090189` at `384/384`)!
   **Eliminated completely:** sequential retry hazard `P(F_i | A_i, prehistory)`, conditional child work `E[W_child | A_i, prehistory]`, subtree independence, repeat-ideal collisions, node-counter aborts, and DAG Markov meter caps!
5. **Single-Node `Q(zeta_15)` Empty-Ball Support:** Live `audit_cm.py: single_node_q15_live` runs exact 16D LLL rank-2 zero-retry DAG descent with ideal-log cache on the reduced `tau`-invariant `Q(zeta_15)` tower (`emp_p0 = 0.6000 <= nb2_p0 = 0.7561`, `0` censors). Across our broader 3,200-root `Q(zeta_15)` single-node campaign (`11,866` reached nodes across 4 regimes incl. 5% fixed missing FB), empirical `P(M = 0) / (1 + 0.25 mu)^(-2)` stayed `<= 0.903 < 1.000` in every regime (exploratory reduced-tower support for single-node `H5`, not a proof at `B ~ 2^71.5`).

## 4. Complete Live-Computed Deterministic-DAG Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|p0|Phi_dag|Nf/Ng|own|W_det|C_det|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|64|70|70.0|0|0|3769.2|1.12e-06|0.00000|823/810|106.59|106.59|1|
|[70,76)|g|sieve|64|70|70.0|0|0|3755.2|1.13e-06|0.00000|893/740|106.59|106.59|1|
|[76,81)|f|sieve|56|76|76.0|1|0|8460.1|2.23e-07|0.00000|826/813|103.59|107.67|3|
|[76,81)|g|sieve|56|76|76.0|1|0|8316.5|2.31e-07|0.00000|902/737|103.60|107.68|3|
|[81,86)|f|sieve|56|76|76.0|1|0|2252.2|3.14e-06|0.00001|843/820|103.53|107.67|3|
|[81,86)|g|sieve|64|76|76.0|1|0|6523.0|3.76e-07|0.00000|938/744|106.60|108.18|3|
|[86,91)|f|sieve|56|76|76.0|2|0|3739.0|1.14e-06|0.00001|861/828|107.35|109.10|5|
|[86,91)|g|sieve|48|76|76.0|2|0|886.9|2.02e-05|0.00002|931/737|104.92|108.70|5|
|[91,96)|f|sieve|56|76|76.0|2|0|1008.9|1.56e-05|0.00002|878/835|106.01|108.81|5|
|[91,96)|g|sieve|56|76|76.0|2|0|895.8|1.98e-05|0.00002|969/744|105.88|108.79|5|
|[96,101)|f|sieve|64|76|76.0|2|0|840.3|2.24e-05|0.00003|910/848|107.43|109.12|5|
|[96,101)|g|sieve|64|76|76.0|2|0|700.7|3.22e-05|0.00004|1006/752|107.26|109.07|5|
|[101,106)|f|sieve|80|76|76.0|2|0|1551.2|6.62e-06|0.00001|950/863|111.83|111.97|5|
|[101,106)|g|sieve|64|81|75.5|1|1|684.5|3.37e-05|0.00004|1028/754|106.97|109.50|9|
|[106,111)|f|sieve|80|81|75.5|1|1|1505.8|7.02e-06|0.00001|968/870|111.79|112.02|9|
|[106,111)|g|sieve|64|81|81.0|2|0|473.3|7.02e-05|0.00008|1051/757|106.81|109.86|13|
|[111,116)|f|sieve|64|81|75.5|2|1|327.8|1.45e-04|0.00015|962/870|112.64|112.81|11|
|[111,116)|g|sieve|64|81|75.5|2|1|235.3|2.79e-04|0.00029|1073/759|111.80|112.09|11|
|[116,121)|f|sieve|80|86|86.0|2|0|1665.9|5.74e-06|0.00002|1003/885|111.76|112.12|13|
|[116,121)|g|sieve|80|96|83.0|1|1|4259.2|8.80e-07|0.00005|1119/769|111.77|112.28|17|
|[121,126)|f|sieve|80|96|83.0|1|1|1750.5|5.20e-06|0.00006|1020/893|111.75|112.27|17|
|[121,126)|g|sieve|80|101|85.5|1|1|3161.3|1.60e-06|0.00007|1141/772|111.76|112.30|17|
|[126,131)|f|sieve|80|96|96.0|2|0|1726.1|5.35e-06|0.00009|1038/900|111.76|112.41|21|
|[126,131)|g|sieve|80|116|100.2|2|0|2629.5|2.31e-06|0.00013|1164/774|111.77|112.47|21|
|[131,136)|f|sieve|64|96|83.0|2|1|513.3|5.98e-05|0.00012|1032/900|110.73|111.84|23|
|[131,136)|g|sieve|64|96|83.0|2|1|280.5|1.98e-04|0.00026|1163/769|109.50|111.40|23|
|[136,141)|f|sieve|64|101|85.5|2|1|598.8|4.40e-05|0.00012|1050/908|110.53|111.80|23|
|[136,141)|g|sieve|64|101|85.5|2|1|301.2|1.72e-04|0.00025|1186/772|109.18|111.38|23|
|[141,146)|f|sieve|64|96|96.0|3|0|402.1|9.70e-05|0.00023|1067/915|111.42|112.48|31|
|[141,146)|g|sieve|80|101|85.5|2|1|504.0|6.20e-05|0.00014|1231/782|112.21|112.73|23|
|[146,151)|f|sieve|64|106|93.6|2|1|446.8|7.87e-05|0.00024|1085/923|110.47|112.99|35|
|[146,151)|g|sieve|80|101|101.0|3|0|1662.1|5.76e-06|0.00020|1254/784|113.91|114.19|31|
|[151,156)|f|sieve|80|106|93.6|2|1|820.1|2.36e-05|0.00018|1125/938|112.86|113.79|35|
|[151,156)|g|sieve|80|106|98.1|2|1|467.6|7.19e-05|0.00027|1276/787|112.50|113.64|35|
|[156,161)|f|sieve|64|116|100.2|2|1|477.7|6.89e-05|0.00064|1120/938|110.79|113.92|45|
|[156,161)|g|sieve|80|106|106.0|3|0|486.2|6.66e-05|0.00027|1299/789|113.16|114.52|43|
|[161,166)|f|sieve|80|116|100.2|2|1|865.9|2.11e-05|0.00059|1160/953|113.04|114.43|45|
|[161,166)|g|sieve|80|116|105.9|2|1|496.2|6.40e-05|0.00064|1321/792|112.67|114.75|53|
|[166,171)|f|sieve|80|126|98.0|2|1|601.0|4.37e-05|0.00061|1178/960|112.56|114.32|55|
|[166,171)|g|sieve|80|131|100.5|2|1|520.5|5.82e-05|0.00063|1344/794|112.11|114.25|63|
|[171,176)|f|sieve|80|126|106.8|2|1|742.4|2.87e-05|0.00069|1195/968|113.31|114.99|79|
|[171,176)|g|sieve|80|131|110.0|2|1|648.9|3.75e-05|0.00070|1366/797|112.61|114.84|87|
|[176,181)|f|sieve|80|141|105.5|2|1|991.0|1.62e-05|0.00060|1213/975|113.29|114.98|75|
|[176,181)|g|sieve|64|131|100.5|3|1|500.5|6.29e-05|0.00070|1366/792|114.61|115.33|73|
|[181,186)|f|sieve|80|146|108.0|2|1|742.0|2.87e-05|0.00076|1230/983|113.40|115.09|99|
|[181,186)|g|sieve|80|131|100.5|3|1|846.4|2.21e-05|0.00066|1411/802|116.51|116.74|73|
|[186,191)|f|sieve|80|151|110.5|2|1|501.7|6.26e-05|0.00080|1248/990|113.44|115.55|111|
|[186,191)|g|sieve|64|141|105.5|3|1|378.3|1.09e-04|0.00076|1411/797|113.84|115.35|89|
|[191,196)|f|sieve|80|161|115.5|2|1|986.9|1.63e-05|0.00183|1265/998|114.80|116.47|137|
|[191,196)|g|sieve|80|141|105.5|3|1|626.9|4.02e-05|0.00069|1456/807|115.68|116.28|89|
|[196,201)|f|sieve|80|166|118.0|2|1|644.4|3.81e-05|0.00221|1283/1005|114.96|116.69|159|
|[196,201)|g|sieve|64|151|110.5|3|1|251.7|2.45e-04|0.00109|1456/802|113.35|115.68|133|
|[201,206)|f|sieve|80|171|120.5|2|1|456.1|7.56e-05|0.00225|1300/1013|115.00|116.70|179|
|[201,206)|g|sieve|64|146|119.9|3|1|380.7|1.08e-04|0.00196|1478/804|115.05|116.36|145|
|[206,211)|f|sieve|80|181|125.5|2|1|771.9|2.66e-05|0.00230|1318/1020|117.09|117.86|235|
|[206,211)|g|sieve|80|161|115.5|3|1|1002.9|1.58e-05|0.00228|1524/814|116.22|117.20|161|
|[211,216)|f|sieve|80|181|125.5|2|1|204.5|3.68e-04|0.00264|1335/1028|115.43|117.12|235|
|[211,216)|g|sieve|80|166|118.0|3|1|791.3|2.53e-05|0.00264|1546/817|115.92|117.18|189|
|[216,221)|f|ecm|20|171|136.3|5|1|128915.8|9.63e-10|0.00350|1209/987|124.24|124.24|349|
|[216,221)|g|sieve|80|171|120.5|3|1|566.3|4.92e-05|0.00266|1569/819|115.70|117.09|209|
|[221,226)|f|ecm|20|176|139.6|5|1|203540.6|3.86e-10|0.00360|1227/995|125.07|125.08|397|
|[221,226)|g|sieve|64|181|125.5|3|1|196.3|3.99e-04|0.00311|1568/814|113.53|116.90|269|
|[226,231)|f|ecm|20|171|148.9|5|1|595310.3|4.51e-11|0.00391|1244/1002|126.59|126.59|449|
|[226,231)|g|sieve|80|181|125.5|3|1|298.8|1.75e-04|0.00289|1614/824|115.35|117.22|269|
|[231,236)|f|ecm|20|176|139.6|6|1|866282.7|2.13e-11|0.00404|1262/1010|125.67|125.68|443|
|[231,236)|g|sieve|64|176|139.6|3|1|216.5|3.29e-04|0.00304|1613/819|115.60|117.20|305|
|[236,241)|f|ecm|20|176|139.6|6|1|350349.0|1.30e-10|0.00404|1279/1017|125.67|125.68|443|
|[236,241)|g|sieve|56|191|149.4|3|1|722.9|3.03e-05|0.00316|1622/817|117.20|118.57|399|
|[241,246)|f|ecm|20|171|148.9|6|1|1372964.7|8.49e-12|0.00443|1297/1025|127.35|127.35|515|
|[241,246)|g|sieve|64|191|149.4|3|1|564.8|4.95e-05|0.00318|1658/824|117.50|118.69|399|
|[246,251)|f|ecm|20|191|149.4|6|1|3991958.1|1.00e-12|0.00471|1314/1032|127.43|127.43|597|
|[246,251)|g|sieve|80|191|149.4|3|1|872.7|2.08e-05|0.00315|1704/834|119.20|119.68|399|
|[251,256)|f|ecm|20|191|149.4|6|1|1698034.8|5.55e-12|0.00471|1332/1040|127.43|127.43|597|
|[251,256)|g|sieve|80|196|152.7|3|1|750.0|2.81e-05|0.00419|1726/837|119.19|119.78|437|
|[256,261)|f|ecm|20|191|149.4|6|1|696640.2|3.30e-11|0.00471|1349/1047|127.43|127.43|597|
|[256,261)|g|sieve|80|201|156.0|3|1|584.1|4.63e-05|0.00492|1749/839|119.18|119.79|503|
|[261,266)|f|ecm|20|186|160.6|6|1|4323408.0|8.56e-13|0.00703|1367/1055|129.57|129.57|715|
|[261,266)|g|sieve|80|206|159.2|3|1|479.9|6.83e-05|0.00706|1771/842|119.41|119.97|589|
