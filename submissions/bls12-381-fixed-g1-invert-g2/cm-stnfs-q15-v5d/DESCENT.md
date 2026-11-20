# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `MAIN_LB = 2^71.29`, `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (`<= lr_rethi - 0.009 - 2 SE = 63.4832`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes**, `BOOT_RELS = 48`, `RET_BOOT = 0.50`, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, the empirical `q`-lattice offsets, and the boot deciles (`m_max = 2560.0`b). `audit_cm.py` invokes `cm_models.eval_descent_live` to compute every node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `psi`, `phi`, `W`, `EN`, and boot work/failure 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`109.426407 <= 109.427`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`109.423939`).

## 1. Fresh-Coin Randomization, Constructive 16D LLL Rank-2 Submodule, & Boot

For each discrete logarithm (`log g` with `g = e(P, Q0)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy & Constructive 16D Integer LLL Rank-2 Submodule Extraction

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, neither `w_0` nor `w_1` can vanish for `z' in G_T \ {1}`; setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). For any `(a, b) in L_Guillevic` with `b != 0 (mod P)`, evaluating **`a - b X`** gives **`Psi_g(a - b X) = (b mod P)(gamma - X) = -(b mod P) w_1^(-1) c_6 z'`**, which has exact discrete log `e log_r(z) mod r` and exact relative norm **`N_(K_g/K_H)(a - b X) = Res_X(a - b X, X^3 - omega^2 u) = a^3 - omega^2 u b^3`**.
- **Independent Single-Cell Boot Trials (`cells_b = 1`, ZERO `BALL_PAIR`, ZERO 17-cell dependence!):**
  1. For each boot trial `t = 1, 2, ...`, draw a **fresh independent uniform `e_t in {1, ..., r - 1}`**, compute `z_t = z' * g'^(e_t)` (`GT` exponentiation `< 2^25` Ir, absorbed into `LLL_NODE = 2^38` Ir alongside 16D LLL `< 2^28` Ir), solve the `2 x 3` `F_(p^2)`-linear system for `gamma_t`, and form the 16 x 16 HNF basis of `L_Guillevic(gamma_t)`.
  2. Run ordinary 16D integer LLL under the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)` and take the **single shortest vector `b_1 = (a_1, b_1)`** (`cells_b = 1`).
  3. Apply the 1-line deterministic norm gate **`|N_(K_g/Q)(a_1^3 - omega^2 u b_1^3)| <= 2^2560`** (`m_sv = 2560.0`b; passes `32/32` targets in `[11]` with `sv_mean = 2540.82`b, `sv_max = 2553.16`b, and `4/4` live BLS12-381 witnesses in `boot_witness.json` with `b_1` norms `2534.44..2542.66`b; conservatively charged `RET_BOOT = 0.50` = `2x` trials!).
  4. Because each trial draws an independent uniform `e_t`, **consecutive boot trials are 100% i.i.d.** — requiring **ZERO 17-cell joint-smoothness assumption and ZERO `BALL_PAIR` waiting-time heuristic (`ball_pair = 0`)**! Smoothness is evaluated strictly at the worst-case norm cap `m_sv = 2560.0`b (`ss_eval([2560.0], ...)`), giving exact geometric waiting time `per_rel = 1 / (RET_BOOT * pr_b * det_b)` and exact Wald's identity across up to `BOOT_RELS = 48` independent qualifying boot relations.
- **Live verification (`audit_cm.py`: `voronoi_and_boot_live`):** (a) Runs full constructive 16D integer LLL + `(v_1 = b_1, v_2 = b_k)` rank-2 submodule extraction (used for special-`q` nodes below) live from an unreduced 16 x 16 HNF basis on split-base `p_0 = 7207 == 7 (mod 15)` (`k = 2`, `index = 1`, `a_i == gamma b_i mod P`); (b) verifies 4 bundled reduced 381-bit BLS12-381 `w'=omega^2` bases (`boot_witness.json`): `a_i == gamma b_i (mod P)`, `|N(a_1 b_2 - a_2 b_1)| == p^4` (`index = 1`), and **`4/4` single-cell `b_1` norms `<= 2542.66 < m_sv = 2560.0`b** (`68/68` ball norms `< 2600`b).
- **Explicit multi-factor boot ECM extraction loop & limb-scaled 100-round MR:** After `TDIV` strips primes `<= 2^20`, each `<= 2560`b `b_1` norm has at most `k_lx = 129` primes in `(2^20, 2^70]`, `kA = 13` primes in `(2^70, n_A = 2^234.53125]`, and `1` prime in `(n_A, n_1 = 2^265]`. Every trial is conservatively charged `LLL_NODE +` ALL `129 * (ecm(70, 2560, 0.999) + mr(2560)) + 13 * (ecm(n_A, 2560, 0.999) + mr(2560)) + mr(2560) + IDEAL` instructions (zero early-exit credit; `mr_instr(2560) = 1.21e10` Ir limb-scaled 100-round MR, `<= 4^(-100) = 2^(-200)` error/test) with multiplicative detection **`det_b = 0.999^143 = 0.8667`**.
- **Certified Boot (`audit_cm.py` live):** `n_1 = 265`, `n_A = 234.53125`, `k_A = 13`, `cells_b = 1` (i.i.d.), `m_sv = 2560.0`b, `BOOT_RELS = 48`, `RET_BOOT = 0.50`:
  - Default `MESH = 96/192`: `pr = 2^-38.158758`, `miss(<=2^70) = 0.222155`, `psi = 0.350000`, **`phi_boot = psi^48 = 1.3040e-22 <= 1.0e-12`**, boot own `2^116.6019` + subtree `2^122.6931` -> **`two_log = 2^109.426407 <= 2^109.427`** core units.
  - On-demand `MESH = 384/384` (`{"descent_mesh": 384}`): **`two_log = 2^109.423939`**, `phi_boot = 1.3040e-22`.

## 2. Special-`q` Nodes, Exact `<= 2^70` Leaf Discard, `T_k + psi^k` Union Bound, & Full Survivor Work

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1 = b_1, v_2 = b_k)` (free rank-2 `Z[zeta_15]`-submodule `M <= L_q`, since `rank_Z(L_q) = 16 > 8 = rank_Z(O_H v_1)`) and enumerates `alpha v_1 + beta v_2` (`alpha` unit-canonical in `B_8(R)`, `UDUP_NODE = 0.95` charged). Measured across 56 BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`): worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b (`cm_data.json`). `[7e]` child-conditioning shift across `22` encountered vs `24` random ideals: `(+3.40, +2.35)`b inside `(+10, +4)`b.

- **Exact leaf-domain missing-ideal discard (`<= LX = 2^70`) & conservative `> 2^70` Window-A/B treatment:** Descent's base smoothness / leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = 0, kB = 0`, requiring every prime strictly `<= 2^70`). Any prime `<= 2^70` is a leaf: if it is missing from the solved FB `<= 2^70`, the boot or special-`q` candidate is **discarded cleanly** (`pr_b *= exp(-mb)`, `mu *= exp(-miss_n)`, with `miss_side` integrated exactly over `[40, 70.0]`; below `2^40`, `lam >= 8.77e8` so `exp(-lam) == 0` identically — **zero recursion below `2^70`**, strictly preserving size decrease at every bin `0..38`). Every prime **`> 2^70`** — **whether in `(2^70, B_s]` (`B_f = 2^71.29`, `B_g = 2^71.69`) or `> B_s`** — is conservatively treated as an unknown Window-A/B child, ALREADY counted in `ss_eval`'s `(kA, kB)` budget and descended strictly downward (zero credit taken for solved FB ideals in `(2^70, B_s]`).
- **Half-scale NB-PGF `T_k + psi^k` union bound (`THRESHOLD-REPAIR.md`, NO `PGF(psi)` substitution!), `1.25x` conditional envelope, & full survivor work:** Under `H5` (half-scale Gamma-Poisson / NB PGF envelope `E[t^M] <= G(t) = (1 + beta (1 - t))^(-r)` for `t in [0, 0.8]` with `r = 2`, `alpha = 0.5`, `beta = 0.25 mu`, `DISC = 0`, and explicit `1.25x`-padded conditional child-failure hazard `P(F_i | A_i, prehistory) <= psi = min(0.95, 1.25 * psi_base)` and child work `1.25 * W_base`), each node's failure probability is bounded by **`phi_node = max(phi_floor, min(1, min_{k=1..64}(T_k + psi^k + 1e-10)))`** (evaluated live in `cm_models.phi_nb_pgf_union`; see `THRESHOLD-REPAIR.md` and `ADAPTIVE.md`). Every sieve node is charged **100% full survivor work (`1.0 x survivors * surv_cost`)** (zero lazy-ordering premise!), and every survivor is charged the FULL worst-case multi-factor ECM extraction loop (`k_lx_f`, `k_lx_g` schedules at `<= LX`, `2 kA` schedules at `<= nA`, limb-scaled `mr_instr(nb)` throughout, and `det_n = 0.999^(k_lx_f + k_lx_g + 2*(kA + kB))`).

## 3. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|psi|phi|Nf/Ng|own|W|EN|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|48|70|70.0|0|0|336.9|0.000|0.0200|793/800|99.95|99.95|0.00|
|[70,76)|g|sieve|48|70|70.0|0|0|347.0|0.000|0.0200|863/730|99.95|99.95|0.00|
|[76,81)|f|sieve|40|76|76.0|1|0|448.4|0.050|0.0200|791/801|96.34|101.39|1.86|
|[76,81)|g|sieve|40|76|76.0|1|0|461.7|0.050|0.0200|867/725|96.38|101.39|1.86|
|[81,86)|f|sieve|56|76|76.0|1|0|2335.8|0.050|0.0200|843/820|103.53|103.82|1.86|
|[81,86)|g|sieve|56|76|76.0|1|0|2232.1|0.050|0.0200|924/739|103.53|103.82|1.86|
|[86,91)|f|sieve|40|76|76.0|2|0|192.2|0.100|0.0200|826/816|102.18|103.31|2.71|
|[86,91)|g|sieve|40|76|76.0|2|0|185.5|0.100|0.0200|912/730|102.29|103.36|2.71|
|[91,96)|f|sieve|48|76|76.0|2|0|266.7|0.100|0.0200|862/830|103.57|104.11|2.71|
|[91,96)|g|sieve|48|76|76.0|2|0|242.2|0.100|0.0200|953/739|103.41|104.00|2.71|
|[96,101)|f|sieve|56|76|76.0|2|0|281.6|0.100|0.0200|896/843|104.93|105.16|2.71|
|[96,101)|g|sieve|56|76|76.0|2|0|239.8|0.100|0.0200|992/747|104.76|105.02|2.71|
|[101,106)|f|sieve|64|76|76.0|2|0|230.3|0.100|0.0200|927/855|106.92|106.98|2.71|
|[101,106)|g|sieve|56|81|75.5|1|1|235.7|0.100|0.0200|1014/749|104.22|104.84|3.79|
|[106,111)|f|sieve|64|81|75.5|1|1|235.0|0.100|0.0200|945/863|106.77|106.89|3.79|
|[106,111)|g|sieve|56|81|81.0|2|0|163.4|0.100|0.0200|1037/752|103.98|104.92|4.40|
|[111,116)|f|sieve|64|81|75.5|2|1|341.5|0.150|0.0200|962/870|112.64|112.65|4.13|
|[111,116)|g|sieve|64|81|75.5|2|1|245.1|0.150|0.0200|1073/759|111.80|111.80|4.13|
|[116,121)|f|sieve|64|86|86.0|2|0|262.0|0.100|0.0200|980/878|106.69|107.50|4.40|
|[116,121)|g|sieve|56|96|83.0|1|1|226.9|0.100|0.0200|1082/757|103.75|106.63|4.87|
|[121,126)|f|sieve|64|96|83.0|1|1|264.4|0.100|0.0200|997/885|106.65|107.54|4.87|
|[121,126)|g|sieve|56|101|85.5|1|1|163.9|0.100|0.0200|1104/759|103.67|107.20|4.87|
|[126,131)|f|sieve|56|96|96.0|2|0|89.9|0.100|0.0217|1001/888|103.70|106.72|5.23|
|[126,131)|g|sieve|56|116|100.2|2|0|143.5|0.100|0.0200|1127/762|103.76|107.67|5.23|
|[131,136)|f|sieve|56|96|83.0|2|1|172.6|0.150|0.0200|1018/895|108.55|108.98|5.38|
|[131,136)|g|sieve|64|96|83.0|2|1|292.8|0.150|0.0200|1163/769|109.50|109.75|5.38|
|[136,141)|f|sieve|56|101|85.5|2|1|201.2|0.150|0.0200|1036/903|108.59|109.15|5.38|
|[136,141)|g|sieve|64|101|85.5|2|1|314.7|0.150|0.0200|1186/772|109.18|109.58|5.38|
|[141,146)|f|sieve|64|96|96.0|3|0|420.7|0.150|0.0200|1067/915|111.42|111.49|5.88|
|[141,146)|g|sieve|80|101|85.5|2|1|527.5|0.150|0.0200|1231/782|112.21|112.26|5.38|
|[146,151)|f|sieve|56|106|93.6|2|1|151.1|0.150|0.0200|1071/918|108.53|109.45|6.12|
|[146,151)|g|sieve|64|101|101.0|3|0|266.2|0.150|0.0200|1231/777|110.50|110.77|5.88|
|[151,156)|f|sieve|80|106|93.6|2|1|860.2|0.150|0.0200|1125/938|112.86|112.92|6.12|
|[151,156)|g|sieve|80|106|98.1|2|1|489.9|0.150|0.0200|1276/787|112.50|112.60|6.12|
|[156,161)|f|sieve|56|116|100.2|2|1|162.1|0.150|0.0200|1106/933|108.89|113.91|6.59|
|[156,161)|g|sieve|80|106|106.0|3|0|510.1|0.150|0.0200|1299/789|113.16|113.27|6.51|
|[161,166)|f|sieve|80|116|100.2|2|1|909.7|0.150|0.0200|1160/953|113.04|114.51|6.59|
|[161,166)|g|sieve|80|116|105.9|2|1|520.8|0.150|0.0200|1321/792|112.67|114.40|6.88|
|[166,171)|f|sieve|80|126|98.0|2|1|631.6|0.150|0.0200|1178/960|112.56|114.35|6.97|
|[166,171)|g|sieve|80|131|100.5|2|1|546.6|0.152|0.0200|1344/794|112.11|114.24|7.23|
|[171,176)|f|sieve|80|126|106.8|2|1|781.1|0.150|0.0200|1195/968|113.31|114.63|7.57|
|[171,176)|g|sieve|80|131|110.0|2|1|682.1|0.152|0.0200|1366/797|112.61|114.39|7.75|
|[176,181)|f|sieve|80|141|105.5|2|1|1043.2|0.152|0.0200|1213/975|113.29|114.62|7.52|
|[176,181)|g|sieve|56|131|100.5|3|1|172.2|0.202|0.0200|1352/787|112.53|114.42|7.49|
|[181,186)|f|sieve|80|146|108.0|2|1|781.9|0.152|0.0200|1230/983|113.40|114.77|7.98|
|[181,186)|g|sieve|80|131|100.5|3|1|890.6|0.202|0.0200|1411/802|116.51|116.74|7.49|
|[186,191)|f|sieve|80|151|110.5|2|1|528.9|0.152|0.0200|1248/990|113.44|114.79|8.20|
|[186,191)|g|sieve|64|141|105.5|3|1|398.0|0.202|0.0200|1411/797|113.84|114.92|7.83|
|[191,196)|f|sieve|64|161|115.5|2|1|152.9|0.152|0.0200|1242/990|111.58|116.09|8.61|
|[191,196)|g|sieve|80|141|105.5|3|1|660.9|0.202|0.0200|1456/807|115.68|116.07|7.83|
|[196,201)|f|sieve|80|166|118.0|2|1|680.4|0.152|0.0200|1283/1005|114.96|116.96|8.86|
|[196,201)|g|sieve|64|151|110.5|3|1|265.2|0.202|0.0200|1456/802|113.35|114.82|8.52|
|[201,206)|f|sieve|80|171|120.5|2|1|481.8|0.152|0.0200|1300/1013|115.00|116.97|9.12|
|[201,206)|g|sieve|64|146|119.9|3|1|401.4|0.202|0.0200|1478/804|115.05|116.59|8.63|
|[206,211)|f|sieve|80|181|125.5|2|1|816.3|0.152|0.0200|1318/1020|117.09|117.86|9.59|
|[206,211)|g|sieve|80|161|115.5|3|1|1059.5|0.202|0.0200|1524/814|116.22|117.32|8.90|
|[211,216)|f|sieve|80|181|125.5|2|1|216.3|0.152|0.0200|1335/1028|115.43|117.13|9.59|
|[211,216)|g|sieve|80|166|118.0|3|1|836.4|0.202|0.0200|1546/817|115.92|117.45|9.17|
|[216,221)|f|ecm|20|171|136.3|5|1|135781.4|0.313|0.0200|1209/987|124.24|124.25|10.30|
|[216,221)|g|sieve|80|171|120.5|3|1|598.9|0.202|0.0200|1569/819|115.70|117.38|9.39|
|[221,226)|f|ecm|20|176|139.6|5|1|214482.9|0.313|0.0200|1227/995|125.07|125.08|10.55|
|[221,226)|g|sieve|64|181|125.5|3|1|207.6|0.202|0.0200|1568/814|113.53|117.02|9.83|
|[226,231)|f|ecm|20|171|148.9|5|1|627985.2|0.313|0.0200|1244/1002|126.59|126.59|10.71|
|[226,231)|g|sieve|80|181|125.5|3|1|316.6|0.202|0.0200|1614/824|115.35|117.31|9.83|
|[231,236)|f|ecm|20|176|139.6|6|1|914269.4|0.365|0.0200|1262/1010|125.67|125.68|10.80|
|[231,236)|g|sieve|64|176|139.6|3|1|229.3|0.209|0.0200|1613/819|115.60|117.39|10.02|
|[236,241)|f|ecm|20|176|139.6|6|1|370153.6|0.365|0.0200|1279/1017|125.67|125.68|10.80|
|[236,241)|g|sieve|56|191|149.4|3|1|765.6|0.209|0.0200|1622/817|117.20|118.68|10.48|
|[241,246)|f|ecm|20|171|148.9|6|1|1451277.0|0.365|0.0200|1297/1025|127.35|127.35|11.01|
|[241,246)|g|sieve|64|191|149.4|3|1|599.1|0.209|0.0200|1658/824|117.50|118.80|10.48|
|[246,251)|f|ecm|20|191|149.4|6|1|4221698.8|0.365|0.0200|1314/1032|127.43|127.44|11.31|
|[246,251)|g|sieve|80|191|149.4|3|1|927.1|0.209|0.0200|1704/834|119.20|119.74|10.48|
|[251,256)|f|ecm|20|191|149.4|6|1|1797699.6|0.365|0.0200|1332/1040|127.43|127.44|11.31|
|[251,256)|g|sieve|80|196|152.7|3|1|797.2|0.209|0.0200|1726/837|119.19|119.86|10.68|
|[256,261)|f|ecm|20|191|149.4|6|1|737888.2|0.365|0.0200|1349/1047|127.43|127.44|11.31|
|[256,261)|g|sieve|80|201|156.0|3|1|621.2|0.209|0.0200|1749/839|119.18|119.98|10.95|
|[261,266)|f|ecm|20|186|160.6|6|1|4584369.0|0.365|0.0200|1367/1055|129.57|129.57|11.62|
|[261,266)|g|sieve|80|206|159.2|3|1|511.2|0.209|0.0200|1771/842|119.41|120.19|11.23|

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work:** `2^109.426407` <= **`2^109.427`** core units (`2^109.423939` at `MESH = 384/384`), `+0.563` bits vs `c36` (`108.864`).
- **Metered cap:** `t = 2^121.4` core units; by Markov's inequality, `P(cap) <= 2^(109.427 - 121.4) = 2.487527e-04` (`2.736280e-04` at `+10%` headroom).
- **Node counter (`E[N]` live in `audit_cm.py`):** `E[N] = 2^15.8958`, hard abort `N_MAX = 2^26.2` -> `P(counter abort) <= 2^-10.3042`.
- **Repeat-ideal bound:** Smallest bin `[70, 76)` holds `>= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits; with `2^6` `1/p`-weighting slack, `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-9.26) <= 2^-9.25`.
