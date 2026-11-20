# Probabilistic Sieve Descent with `H-PGF (0.5*mu)`, `THRESHOLD-REPAIR`, & Fresh-Exponent `i.i.d.` Boot (`Q(zeta_15)`)

Target: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, `w' = omega^2` (`alpha = (+2.74, +0.77)`), full de Bruijn `rho_dB`, `boot_rels = 1`, calibrated `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` (`phi >= 0.01`, `phi <= 0.04`, `psi < 0.80`), fresh-exponent `i.i.d.` single-cell boot (`e in {0..r-1}`), `|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`, `H-EA`, zero-loss `repurge_stack`, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`, `D_skew` `S_boot = 2^21` (`boot_supplied_lll_vectors = 0`), node budget `B_node = 2^25.0` + repeated-ideal abort + meter `t_desc = 2^123.5` (Markov).

`cm_data.json` holds the 78 structural tuples `plan_struct`, `q`-offsets, `prod_sq24` + `prod_enc12`, `q15_v5l_campaign`, and boot config (`G0`, `m_sv = 2560.0`b, `s_boot = 2^21`). `cm_models.eval_descent_live` computes `(h, exc_crit)` and all 78 nodes LIVE (`64/192`: `E[D] = 2^110.599603 <= 2^110.600`, `E[N] = 2^14.8147 <= 2^14.82`, `Phi_boot = 0.42625`, `max_phi = 0.0337`, `max_psi = 0.4183`; `384/384`: `E[D] = 2^110.557264`, `Phi_boot = 0.41147`).

## 1. Fresh-Exponent `i.i.d.` Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = diag(1..1, 2^21..2^21)`)

For each log (`log g` with `g = e(P0, Q0)`, and `log z`), each boot trial draws fresh uniform `e in {0, ..., r - 1}` (inclusive of `0`), forming `z_t = z' * g'^e` (resp. `g' * g'^e`), exact `Unif(G_T)` (if `z_t = 1`, redraw).
- **Subfield Elimination & Exact `D_skew` 16D LLL (ALL 4/4 LIVE):** In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` (`|F_(p^4)^*|`, `|F_(p^6)^*|` coprime to `r`), eliminating `X^2` in `c_6 z_t` (`c_6 in F_(p^6)^*`) solves a `2 x 3` homogeneous `F_(p^2)` system giving `c_6 z_t = -w_1 (gamma - X)` with **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** for **100%** of nonidentity targets. Scaling `b` by `S_boot = 2^21` in `L_Guillevic = {(a, b) in O_H^2 : a - gamma b == 0 mod P}` (`det = p^4`) and running pure-integer de Weger `cm_models.lll_16d_int(99/100)` (`||D_skew b_1||_2 <= (0.74)^(-15/4)(p^4 S_boot^8)^(1/16) = 2^107.30`) gives `b_1 = (a_1, b_1)` (`cells_b = 1`) with norm **`a_1^3 - omega^2 u b_1^3`**.
- **Single-Relation Boot (`boot_rels = 1`, `cells_b = 1`, `ball_pair = 0`):** Enforce **`|N_(K_g/Q)(a_1^3 - omega^2 u b_1^3)| <= 2^2560`** (`RET_BOOT = 0.50`, declared model inequality checked on `4/4 <= 2542.56`b deterministic SHA-256 targets live and `32/32 <= 2553.16`b in `[11]`); run worst-case ECM + 100-round MR (`k_lx = 129`, `kA = 13`, `1` in `(n_A, 2^265]`, `mr_instr(2560) = 1.21e10`, `det_b = 0.8667`); descend its children; if any child fails, start a fresh exponent `e`. Expected two-log boot own work **`2^102.8156`** u; Wald retry factor `1/(1 - Phi_boot) = 1.7429x` (`Phi_boot = 13 * phi_max(nA) + phi_max(265) = 0.42625 < 0.50`, NO `Gate E` assumption needed!).
- **ALL-4/4 LIVE Boot LLL (`voronoi_and_boot_live`, `boot_supplied_lll_vectors = 0`):** Verifies `G_0 != 1, G_0^r == 1`, runs `F_(p^12) -> F_(p^6)` elimination AND from-scratch 381-bit pure-integer `lll_16d_int(99/100)` LIVE on ALL `4/4` deterministic SHA-256 targets (`det^2 == p^8 S_boot^16`, 15 Lovasz + 120 size-red. each, `||D_skew b_1||_2 <= 2^106.09 <= 2^107.30`, **`4/4 <= 2542.56 <= 2560.0`b**, deterministic compatibility check for `RET_BOOT = 0.50`, `4/4` deciles) + constructive 16D pure-integer `lll_16d_int` on `p_0 = 7207` (`k = 2, index = 1`, ZERO float LLL).

## 2. Exact Seeded Peeling Closure (`h = 0.256845`), Zero-Loss `repurge_stack`, & `<= LX = 2^70` Discard

- **Seeded Peeling (`_peel_fixed_point_bls`):** Protected seed `lam >= 100.0`; fixed point `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))` gives **`h = 0.256845`**, **`exc_crit = 0.7625 << 1.09`** (`+42.9%` margin).
- **Zero-Loss Trim + LIFO `repurge_stack` (`RANK-BRIDGE.md` §3):** Clique trim + Markowitz deg-0 guard keep **`post_peel_cols_lost == post_merge_cols_lost == 0`** (`n_reconstruct_unit == 0`); 3-stack triangular back-sub recovers 100% of 2-core columns (`< 2^64` u in `filtering_and_misc = 2^80`), making `miss_side` exact.
- **`<= LX = 2^70` Discard:** `[70, 76)` has `n_c = 70, kA = kB = 0`. Missing primes `<= 2^70` discard via `exp(-miss)` over `[40, 70)`; every prime `> 2^70` is conservatively descended.

## 3. `H-PGF (0.5*mu)` + `THRESHOLD-REPAIR` (`phi >= 0.01`), Node Budget & Repeat Loss, & Live `Q(zeta_15)` Audits

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D pure-integer `lll_16d_int` extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`. The deterministic rank-2 search set enumerates **`S_rank2(R) = A(R) x B_8(R)`** (`alpha v_1 + beta v_2`): `A(R)` consists of the first **`N_alpha(R) = floor((B_8(R) - 1) * 0.008700)`** canonical `alpha` vectors (lexicographic order inside `(B_8(R) \ {0}) cap (F_Vor/mu_30)`, `P_CAN_ALPHA = 0.008700 <= p_voronoi/30`), paired with ALL `beta in B_8(R)` (`||beta||_2 <= R`), giving exact cell count **`cells_rank2[R] = N_alpha(R) * B_8(R)`** (`B_8(R)` exact 8D DP for `R <= 32`, conservative lower bound `floor(0.998 (pi^4/24) R^8)` for `R > 32`, **recomputed with `cells_rank2[R] % B_8(R) == 0` asserted LIVE in `audit_cm.py`**; `UDUP_NODE = 0.95`; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b; live `[7e]` production shift `g: (-2.36, -3.35)`b, `f: (+0.35, -2.27)`b inside `(+10, +4)`b, alongside offline GP `[7e]` `(+3.40, +2.35)`b).

1. **Calibrated `H-PGF (0.5*mu)` Reached-Ideal Ensemble Envelope (`H5`, NO exact NB claim!):** Let `M` be the number of distinct qualifying usable relations in a node's ball before child descents, and let `mu` be the raw Dickman/Murphy/retention/detection/norm-envelope scale. **`mu` is NOT asserted to equal the measured mean of `M`.** For `t in [0, 0.8]`, `H-PGF` states the calibrated conservative PGF envelope:
   `E[t^M] <= G_(mu,r)(t) := (1 + 0.5 * mu * (1 - t) / r)^(-r)`,
   with `r = 0.5` when effective window cap `kA + kB <= 1`, and `r = 1.0` otherwise. Here `0.5` is a calibrated conservative envelope parameter; no exact negative-binomial distribution or validated production mean is claimed. `H-PGF` concerns the reached-ideal ensembles used by the recurrence (including parent selection and prior failures, not conditional on a fully fixed `q`) and extrapolates across degree, prime size, and adaptive depth as an explicit exploratory STNFS hypothesis (`H5`).
2. **`THRESHOLD-REPAIR` (ZERO count/solvability independence assumed!):** Every reached node enforces **`|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, charges full ball sieve/bound + full modeled survivor cofactor work (`H-EA`, NO lazy-work discount), and tries qualifying relations sequentially until first child-tree success. Let `psi = kA (phi_max,f(nA) + phi_max,g(nA)) + kB (phi_max,f(n) + phi_max,g(n)) < 0.80` be the additive union bound on a relation's child failures, with sequential hypothesis `P(F_i | M >= i, F_1..F_(i-1), pre-entry history) <= psi`. For every `k >= 1` and `t in (0, 0.8]` (`t = 0` for `k = 1`), Markov's inequality on `t^M >= t^(k-1) 1_(M < k)` + union bound prove:
   `P(node fails) <= P(M < k) + psi^k <= G_(mu,r)(t) / t^(k-1) + psi^k`.
   The plan's `phi` minimizes this certificate over the fixed finite grid `k in 1..64, t in {0.00, 0.01, ..., 0.80}` and **floors at `phi >= 0.01`**, without ever substituting `PGF(psi)` or assuming count/solvability independence. All 78 rows satisfy `0.01 <= phi <= 0.04` (`max_phi = 0.0337`, `max_psi = 0.4183 < 0.80`, `min_mu = 465.0`). Expected attempted relations are `<= 1/(1 - psi)`, giving `W = own + (1/(1 - psi)) * child_W` and `E[N] = 1 + (1/(1 - psi)) * child_N`.
3. **Fresh-Exponent `i.i.d.` Boot + Node Budget (`B_node = 2^25`) + Repeated-Ideal Abort + Meter (`t_desc = 2^123.5`):** Because each boot retry draws fresh uniform `e in {0..r-1}` (`Unif(G_T)`), Wald's identity gives `E[D] = 2 (boot_own + W_dag_1)/(1 - Phi_boot) = 2^110.599603 <= 2^110.600` u (`2^110.557264` at `384/384`) and `E[N_2log] = 2 E[N_1]/(1 - Phi_boot) = 2^14.8147 <= 2^14.82`. Three explicit abort guards bound all tail/conditioning risks via pure Markov & union bounds:
   - **Node invocation budget `B_node = 2^25.0`:** abort if total nodes across both logs exceed `B_node` -> Markov `p_node_abort <= E[N_2log] / B_node <= 8.7e-4`.
   - **Repeated-ideal abort:** abort if any descent ideal is visited twice (preserving the unconditioned reached-ideal ensemble). Every descent ideal has norm `q > 2^70` (`N_ideals(70) = (2^76 - 2^70)/(12 * 76 * ln 2) ~= 2^66.70 > 2^66.6` degree-1 `<tau>`-orbits in the smallest bin `[70, 76)`). Under `H5`'s declared collision premise (conditional on history, a newly generated child ideal hits any specific prior ideal in bin `b` with probability `<= C_rep / N_ideals(b)`, `C_rep = 2^6 = 2^76 / 2^70` exact `1/p`-weighting ratio across the 6-bit bin `[70, 76)` -> `N_eff = N_ideals(70)/C_rep > 2^60.70 > 2^60`), birthday union bound over `<= B_node` nodes gives `p_repeat <= C_rep * B_node^2 / (2 N_ideals(70)) <= B_node^2 / (2 * 2^60) = 2^(-11) <= 4.9e-4`.
   - **Combined-log phase work meter `t_desc = 2^123.5`:** Markov `p_cap <= E[D] / t_desc = 0.00013083 <= 1.4e-4`.
4. **Two Complementary Live Same-Field `Q(zeta_15)` Audits (`prod_sq36_live` + `q15_v5l_descent_live` in `audit_cm.py`, ONE Unified Pure-Int `lll_16d_int`):**
   - **(a) 36 ACTUAL Production BLS12-381 `(f_9, g_3)` `Q(zeta_15)` Special-`q` 16D Lattices (`prod_sq36_live`: `24` control `73..258`b + `12` parent-verified encountered `71..255`b):** All 36 bases were pre-reduced with pure exact integer `lll_16d_int(99/100)` and are verified **LIVE** by `verify_prod_sq36` (`det % q == 0`, `q | N_side`, LIVE `12/12` parent divisibility **`q | N_side(parent_a, parent_b)`** with `||parent||_2 > q^(1/16)`, sextiles `23.87–42.83`b below `qbins`, survivor ratio `0.0193–0.2598 <= 1.0`, LIVE `[7e]` shift `g: -2.36/-3.35`b, `f: +0.35/-2.27`b). On `2,304` exact production norms (`36 x 64` `P^1`-unique pts, primes `<= 4096`), all 4 windows `thr12b..15b` (`mean = 3.28..6.11`) have **`0/36` empty (`M_min = 1..2`), `enc/ctrl = 0.7317..0.8780`, `PGF_emp / H-PGF(0.5*mu, r=1.0) <= 0.6781 < 1.0` & `r=0.5 <= 0.6571 < 1.0`**, with post-`4096` cofactor residuals logged in `prod_cofactor_rem_bits`.
   - **(b) Same-Field `Q(zeta_15)` Exact-`<tau>`-Orbit `<K_H^*, tau>`-Deduplicated Full Usable-Relation Real-2-Core-FB Campaign (`q15_v5l_descent_live`, `224` ideals; Stage 1 REAL `B_FB = 2^14` 2-core FB + `4` live `2,048`-pt balls via `lll_16d_int`):** Exact `Q(zeta_15)` & order-12 `tau` bijection verified LIVE on split `p <= 200`; `393/520` unique 2-core rels (`|FB_orb_g|=112, |FB_orb_f|=107`); across `224` ideals (`100` D1 + `24` D2 + `100` CTRL): **`0/224` empty (`CP95 0.0133 << H-PGF(0.5*mu) p0 = 0.2376`), `M_u/mu_u = 1.4473` (`>= 1.39` all cohorts), `D1/CTRL = 0.9278`, `D2/CTRL = 0.9616`, `PGF_emp / H-PGF(0.5*mu, r=1.0) <= 0.2957` & `r=0.5 <= 0.2720`, actual two-sided `sv2/sv_upper = 0.0235` (`max 0.0579`)!**
   - **Exploratory Scope:** Together, `prod_sq36_live` and `q15_v5l_descent_live` provide scoped empirical compatibility support (NOT statistical proof of production tails); transfer across polynomial degree, prime size, and adaptive depth rests explicitly on declared exploratory hypotheses `H5` (`H-PGF (0.5*mu)`) and `H-EA`.

## 4. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 64/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|psi|phi|Nf/Ng|own|W|EN|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|56|70|70|0|0|1241|0.000|0.028|809/805|103.5|103.5|1.0|
|[70,76)|g|sieve|56|70|70|0|0|1257|0.000|0.028|879/735|103.5|103.5|1.0|
|[76,81)|f|sieve|60|76|76|1|0|15090|0.057|0.024|833/815|105.2|105.9|3.1|
|[76,81)|g|sieve|60|76|76|1|0|14773|0.057|0.025|909/739|105.2|105.9|3.1|
|[81,86)|f|sieve|68|76|76|1|0|11516|0.057|0.027|864/827|108.0|108.1|3.1|
|[81,86)|g|sieve|68|76|76|1|0|10751|0.057|0.028|945/746|108.0|108.1|3.1|
|[86,91)|f|sieve|48|76|76|2|0|942|0.113|0.016|845/823|105.0|106.4|5.5|
|[86,91)|g|sieve|48|76|76|2|0|885|0.113|0.017|931/737|104.9|106.4|5.5|
|[91,96)|f|sieve|52|76|76|2|0|531|0.113|0.027|871/833|104.9|106.3|5.5|
|[91,96)|g|sieve|52|76|76|2|0|475|0.113|0.030|962/742|104.8|106.3|5.5|
|[96,101)|f|sieve|60|76|76|2|0|480|0.113|0.029|903/845|106.3|107.0|5.5|
|[96,101)|g|sieve|64|76|76|2|0|699|0.113|0.021|1006/752|107.3|107.7|5.5|
|[101,106)|f|sieve|72|76|76|2|0|611|0.113|0.023|940/859|109.5|109.6|5.5|
|[101,106)|g|sieve|64|81|76|1|1|681|0.113|0.021|1028/754|107.0|108.2|10.3|
|[106,111)|f|sieve|72|81|76|1|1|625|0.113|0.023|957/867|109.4|109.7|10.3|
|[106,111)|g|sieve|64|81|81|2|0|471|0.113|0.030|1051/757|106.8|108.6|15.1|
|[111,116)|f|sieve|68|81|76|2|1|543|0.170|0.029|969/872|113.6|113.6|13.3|
|[111,116)|g|sieve|72|81|76|2|1|642|0.170|0.026|1086/763|113.5|113.6|13.3|
|[116,121)|f|sieve|72|86|86|2|0|690|0.113|0.021|992/882|109.4|110.9|15.1|
|[116,121)|g|sieve|64|96|83|1|1|646|0.115|0.022|1096/762|106.7|110.4|20.5|
|[121,126)|f|sieve|72|96|83|1|1|686|0.115|0.021|1010/889|109.3|110.9|20.5|
|[121,126)|g|sieve|64|101|86|1|1|465|0.115|0.030|1118/764|106.7|110.4|20.5|
|[126,131)|f|sieve|72|96|96|2|0|714|0.116|0.020|1027/897|109.4|110.9|25.9|
|[126,131)|g|sieve|68|116|100|2|0|680|0.117|0.021|1147/769|108.1|110.6|26.0|
|[131,136)|f|sieve|68|96|83|2|1|838|0.171|0.021|1039/902|111.7|112.4|29.4|
|[131,136)|g|sieve|72|96|83|2|1|752|0.171|0.023|1176/773|111.4|112.2|29.4|
|[136,141)|f|sieve|64|101|86|2|1|585|0.172|0.028|1050/908|110.5|111.8|29.4|
|[136,141)|g|sieve|72|101|86|2|1|810|0.172|0.022|1198/776|111.1|112.0|29.4|
|[141,146)|f|sieve|68|96|96|3|0|654|0.174|0.026|1074/917|112.3|112.8|41.0|
|[141,146)|g|sieve|84|101|86|2|1|733|0.172|0.023|1237/783|113.2|113.5|29.4|
|[146,151)|f|sieve|68|106|94|2|1|746|0.175|0.023|1091/925|111.5|112.4|46.8|
|[146,151)|g|sieve|72|101|101|3|0|680|0.176|0.025|1243/781|112.3|112.8|41.1|
|[151,156)|f|sieve|80|106|94|2|1|797|0.175|0.022|1125/938|112.9|113.3|46.8|
|[151,156)|g|sieve|84|106|98|2|1|672|0.176|0.025|1282/788|113.5|113.8|46.9|
|[156,161)|f|sieve|68|116|100|2|1|788|0.177|0.023|1126/940|111.7|115.1|62.3|
|[156,161)|g|sieve|84|106|106|3|0|736|0.176|0.024|1304/791|114.0|114.3|58.5|
|[161,166)|f|sieve|76|116|100|2|1|541|0.177|0.030|1155/951|112.2|115.1|62.3|
|[161,166)|g|sieve|84|116|106|2|1|712|0.177|0.024|1327/793|113.6|115.4|73.9|
|[166,171)|f|sieve|80|126|98|2|1|573|0.177|0.029|1178/960|112.6|115.2|77.6|
|[166,171)|g|sieve|84|131|100|2|1|767|0.177|0.023|1349/796|113.2|115.3|90.9|
|[171,176)|f|sieve|80|126|107|2|1|709|0.177|0.024|1195/968|113.3|115.4|112.5|
|[171,176)|g|sieve|80|131|110|2|1|618|0.177|0.027|1366/797|112.6|115.2|125.8|
|[176,181)|f|sieve|76|141|106|2|1|612|0.177|0.027|1208/973|112.4|115.2|110.8|
|[176,181)|g|sieve|68|131|100|3|1|807|0.236|0.026|1372/794|115.6|116.3|112.2|
|[181,186)|f|sieve|80|146|108|2|1|696|0.177|0.025|1230/983|113.4|115.4|148.2|
|[181,186)|g|sieve|80|131|100|3|1|795|0.236|0.027|1411/802|116.5|117.0|112.2|
|[186,191)|f|sieve|84|151|110|2|1|735|0.177|0.024|1253/992|114.3|115.7|169.6|
|[186,191)|g|sieve|72|141|106|3|1|970|0.236|0.022|1423/801|115.7|116.4|139.9|
|[191,196)|f|sieve|76|161|116|2|1|594|0.178|0.028|1260/996|114.1|117.1|217.1|
|[191,196)|g|sieve|84|141|106|3|1|870|0.236|0.025|1462/808|116.5|116.9|139.9|
|[196,201)|f|sieve|80|166|118|2|1|595|0.179|0.028|1283/1005|115.0|117.5|253.4|
|[196,201)|g|sieve|76|151|110|3|1|998|0.236|0.022|1474/807|116.0|116.7|215.8|
|[201,206)|f|sieve|84|171|120|2|1|624|0.179|0.027|1306/1014|115.7|117.7|292.7|
|[201,206)|g|sieve|72|146|120|3|1|965|0.237|0.023|1491/808|116.9|117.9|233.1|
|[206,211)|f|sieve|80|181|126|2|1|703|0.180|0.025|1318/1020|117.1|118.4|391.5|
|[206,211)|g|sieve|80|161|116|3|1|914|0.237|0.024|1524/814|116.2|117.9|271.0|
|[211,216)|f|sieve|92|181|126|2|1|578|0.180|0.029|1350/1033|117.5|118.6|391.5|
|[211,216)|g|sieve|80|166|118|3|1|716|0.238|0.029|1546/817|115.9|118.0|319.6|
|[216,221)|f|ecm|14|171|136|5|1|2738|0.359|0.013|1172/975|116.0|118.5|721.8|
|[216,221)|g|sieve|84|171|120|3|1|792|0.238|0.027|1574/821|116.5|118.2|362.0|
|[221,226)|f|ecm|14|176|140|5|1|4268|0.359|0.010|1190/982|116.9|118.8|830.7|
|[221,226)|g|sieve|76|181|126|3|1|772|0.239|0.028|1586/820|116.2|118.4|476.0|
|[226,231)|f|ecm|12|171|149|5|1|2216|0.359|0.016|1191/985|114.8|118.4|949.6|
|[226,231)|g|sieve|92|181|126|3|1|859|0.239|0.025|1628/829|117.6|118.8|476.0|
|[231,236)|f|ecm|12|176|140|6|1|2601|0.418|0.017|1209/992|113.9|118.7|1016.8|
|[231,236)|g|sieve|76|176|140|3|1|859|0.239|0.025|1631/825|118.3|119.0|546.0|
|[236,241)|f|ecm|12|176|140|6|1|1109|0.418|0.033|1226/1000|113.9|118.7|1016.8|
|[236,241)|g|sieve|56|191|149|3|1|644|0.239|0.032|1622/817|117.2|118.9|755.0|
|[241,246)|f|ecm|12|171|149|6|1|3859|0.418|0.012|1244/1007|115.6|118.8|1198.1|
|[241,246)|g|sieve|68|191|149|3|1|824|0.239|0.026|1665/826|118.4|119.4|755.0|
|[246,251)|f|ecm|10|191|149|6|1|1353|0.418|0.028|1242/1008|111.5|119.1|1440.7|
|[246,251)|g|sieve|80|191|149|3|1|771|0.239|0.028|1704/834|119.2|119.8|755.0|
|[251,256)|f|ecm|12|191|149|6|1|4876|0.418|0.010|1279/1022|115.7|119.2|1440.7|
|[251,256)|g|sieve|80|196|153|3|1|659|0.239|0.032|1726/837|119.2|120.0|840.2|
|[256,261)|f|ecm|12|191|149|6|1|2112|0.418|0.020|1296/1030|115.7|119.2|1440.7|
|[256,261)|g|sieve|84|201|156|3|1|776|0.239|0.027|1754/841|119.9|120.5|987.7|
|[261,266)|f|ecm|10|186|161|6|1|1359|0.418|0.028|1295/1031|113.6|119.7|1718.2|
|[261,266)|g|sieve|84|206|159|3|1|615|0.239|0.034|1777/843|120.1|120.9|1168.8|
