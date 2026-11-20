# Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `MAIN_LB = 2^71.29`, `LX = 2^70`, conservative 2-SE lowered relation count `LOG2_RELS = 63.46` (`<= lr_rethi - 0.009 - 2 SE = 63.4846`), twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-q nodes**, `BOOT_RELS = 48`, `RET_BOOT = 0.50`, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`, `BALL_PAIR = 2^24`.

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, the empirical `q`-lattice offsets, and the boot deciles (`m_max = 2600.0`b). `audit_cm.py` invokes `cm_models.eval_descent_live` to compute every node's `pr`, `sv`, `mu`, `miss`, `sieve`, `surv_cost`, `own`, `psi`, `phi`, `W`, `EN`, and boot work/failure 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`107.448555 <= 107.449`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`107.383549`).

## 1. Fresh-Coin Randomization, Constructive 16D LLL Rank-2 Submodule, & Boot

For each discrete logarithm (`log g` with `g = e(P, Q0)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

### Algebraic Non-Degeneracy & Constructive 16D Integer LLL Rank-2 Submodule Extraction

In `F_(p^12) = F_(p^4)[X]/(X^3 - omega^2 u)` with `F_(p^4) = O_H / P = Fix(tau^4)` and `F_(p^6) = F_(p^2)[X]/(X^3 - omega^2 u) = Fix(tau^6)`:
- Both `F_(p^4)^*` and `F_(p^6)^*` have orders coprime to `r`, so every nonzero element in `F_(p^4)^* F_(p^6)^*` has discrete log `0 mod r`.
- Write `z' = z_0 + z_1 X + z_2 X^2` (`z_i in F_(p^4)`). Vanishing of the `X^2`-coefficient of `c_6 z'` (`c_6 in F_(p^6)^*`) is a homogeneous `2 x 3` linear system over `F_(p^2)` and therefore always has a nonzero solution `c_6`. Writing `c_6 z' = w_0 + w_1 X`, neither `w_0` nor `w_1` can vanish for `z' in G_T \ {1}`; setting **`gamma = -w_0 / w_1 in F_(p^4) \ F_(p^2)`** gives `c_6 z' = -w_1 (gamma - X)` for **100% of nonidentity targets**.
- Define the 16D lattice `L_Guillevic = { (a, b) in O_H^2 : a - gamma b == 0 (mod P) } subset Z^16` (`det = p^4`). For any `(a, b) in L_Guillevic` with `b != 0 (mod P)`, evaluating **`a - b X`** gives **`Psi_g(a - b X) = (b mod P)(gamma - X) = -(b mod P) w_1^(-1) c_6 z'`**, which has exact discrete log `e log_r(z) mod r` and exact relative norm **`N_(K_g/K_H)(a - b X) = Res_X(a - b X, X^3 - omega^2 u) = a^3 - omega^2 u b^3`**.
- **Constructive extraction algorithm (NOT assuming LLL generates all of `L_Guillevic`):**
  1. Form the 16 x 16 HNF basis of `L_Guillevic` (4 `F_p`-linear congruences `a == gamma b mod (p, h_pi)`).
  2. Run ordinary 16D integer LLL under the `T_2`-skewed metric `Q_boot = diag(G_t2, |u|^(2/3) G_t2)` (charged `LLL_NODE = 2^38` Ir per ball; `< 2^28` Ir actual).
  3. Set `v_1 = b_1` (shortest LLL vector) and `v_2 = b_k` (the first LLL basis vector `K_H`-linearly independent of `v_1`, i.e., `det_H(v_1, b_k) = a_1 b_k - b_1 a_k != 0` in `O_H`; since `rank_Z(L_Guillevic) = 16 > 8 = rank_Z(O_H v_1)`, such a `k in {2..9}` always exists).
  4. Then **`M = Z[zeta_15] v_1 + Z[zeta_15] v_2 <= L_Guillevic`** is **provably a free rank-2 `Z[zeta_15]`-submodule**: every vector `v_1 + beta v_2` (`beta in B_8(1)`, `17` vectors `0, +-e_i`) lies in `L_Guillevic`, has zero unit duplicates (`UDUP_BOOT = 1.0`), and is base-primitive (`gcd(1, beta) = 1`).
- **Live verification (`audit_cm.py`: `voronoi_and_boot_live`):** (a) Runs the full constructive 16D integer LLL + `(b_1, b_k)` extractor live from an unreduced 16 x 16 HNF basis on split-base `p_0 = 7207 == 7 (mod 15)` (`k = 2`, `index = 1`, `a_i == gamma b_i mod P`); (b) verifies 4 bundled reduced 381-bit BLS12-381 `w'=omega^2` bases (`boot_witness.json`): `a_i == gamma b_i (mod P)`, `|N(a_1 b_2 - a_2 b_1)| == p^4` (`index = 1`), all `68/68` exact relative norms `|N(a^3 - u omega^2 b^3)| < 2600`b (mean `2551.75`b, max `2563.01`b), AND **`4/4` passes of the exact 10-decile discard predicate** (`ball_sorted[idx_j] <= boot_sub_deciles[j]`).
- **Deterministic 50% norm-certificate filter (`RET_BOOT = 0.50`, `m_max = 2600.0`b):** After LLL and `(v_1, v_2)` extraction, the algorithm computes all 17 exact norms of `a - b X` and discards any ball whose sorted deciles exceed `boot_sub_deciles` (charging `RET_BOOT = 0.50` = `2x` balls; supported by `32/32` targets in `[11]` with max `2575.73`b, `kappa_hat = 0.7076 << BALL_PAIR = 2^24`, and `4/4` live witnesses).
- **Explicit multi-factor boot ECM extraction loop & limb-scaled MR:** After `TDIV` strips primes `<= 2^20`, each `<= 2600`b norm has at most `k_lx = 131` primes in `(2^20, 2^70]`, `kA = 13` primes in `(2^70, n_A = 2^234.53125]`, and `1` prime in `(n_A, n_1 = 2^265]`. Every boot candidate is conservatively charged ALL `131 * (ecm(70, 2600, 0.999) + mr(2600)) + 13 * (ecm(n_A, 2600, 0.999) + mr(2600)) + mr(2600) + IDEAL` instructions (zero early-exit credit; `mr_instr(2600) = 5.14e9` Ir limb-scaled 40-round MR) with multiplicative detection **`det_b = 0.999^145 = 0.8650`**.
- **Certified Boot (`audit_cm.py` live):** `n_1 = 265`, `n_A = 234.53125`, `k_A = 13`, `BOOT_RELS = 48`, `RET_BOOT = 0.50`:
  - Default `MESH = 96/192`: `pr = 2^-37.905005`, `miss = 0.225854`, `psi = 0.367379`, **`phi_boot = psi^48 = 1.3350e-21 <= 1.0e-12`**, boot own `2^116.4621` + subtree `2^120.6597` -> **`two_log = 2^107.448555 <= 2^107.449`** core units.
  - On-demand `MESH = 384/384` (`{"descent_mesh": 384}`): **`two_log = 2^107.383549`**, `phi_boot = 4.5566e-23`.

## 2. Special-`q` Nodes, Exact `P(M < k) + psi^k` Union Bound, & Full Survivor Work

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1, v_2)` and enumerates `alpha v_1 + beta v_2` (`alpha` unit-canonical in `B_8(R)`, `UDUP_NODE = 0.95` charged). Measured across 56 BLS12-381 special-`q` lattices (`q in [2^72, 2^260]`): worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b (`cm_data.json`). `[7e]` child-conditioning shift across `22` encountered vs `24` random ideals: `(+3.40, +2.35)`b inside `(+10, +4)`b.

- **Half-scale NB-PGF `T_k + psi^k` union bound (`THRESHOLD-REPAIR.md`, NO `PGF(psi)` substitution!) & full survivor work:** Under `H5` (half-scale Gamma-Poisson / NB PGF envelope `E[t^M] <= G(t) = (1 + beta (1 - t))^(-r)` for `t in [0, 0.8]` with `r = 2`, `alpha = 0.5`, `beta = 0.25 mu`, `DISC = 0`, and explicit conditional child-failure hazard `P(F_i | A_i, prehistory) <= psi`), each node's failure probability is bounded by **`phi_node = max(phi_floor, min(1, min_{k=1..64}(T_k + psi^k + 1e-10)))`** (evaluated live in `cm_models.phi_nb_pgf_union`; see `THRESHOLD-REPAIR.md` and `ADAPTIVE.md`). Every sieve node is charged **100% full survivor work (`1.0 x survivors * surv_cost`)** (zero lazy-ordering premise!), and every survivor is charged the FULL worst-case multi-factor ECM extraction loop (`k_lx_f`, `k_lx_g` schedules at `<= LX`, `2 kA` schedules at `<= nA`, limb-scaled `mr_instr(nb)` throughout, and `det_n = 0.999^(k_lx_f + k_lx_g + 2*(kA + kB))`).

## 3. Complete Live-Computed Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|psi|phi|Nf/Ng|own|W|EN|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|48|70|70.0|0|0|336.9|0.000|0.0200|793/800|99.95|99.95|0.00|
|[70,76)|g|sieve|48|70|70.0|0|0|347.0|0.000|0.0200|863/730|99.95|99.95|0.00|
|[76,81)|f|sieve|40|76|76.0|1|0|448.4|0.040|0.0200|791/801|96.34|101.07|1.62|
|[76,81)|g|sieve|40|76|76.0|1|0|461.7|0.040|0.0200|867/725|96.38|101.07|1.62|
|[81,86)|f|sieve|56|76|76.0|1|0|2335.8|0.040|0.0200|843/820|103.53|103.77|1.62|
|[81,86)|g|sieve|56|76|76.0|1|0|2232.1|0.040|0.0200|924/739|103.53|103.77|1.62|
|[86,91)|f|sieve|40|76|76.0|2|0|192.2|0.080|0.0200|826/816|102.18|103.13|2.42|
|[86,91)|g|sieve|40|76|76.0|2|0|185.5|0.080|0.0200|912/730|102.29|103.19|2.42|
|[91,96)|f|sieve|48|76|76.0|2|0|266.7|0.080|0.0200|862/830|103.57|104.00|2.42|
|[91,96)|g|sieve|48|76|76.0|2|0|242.2|0.080|0.0200|953/739|103.41|103.89|2.42|
|[96,101)|f|sieve|56|76|76.0|2|0|281.6|0.080|0.0200|896/843|104.93|105.11|2.42|
|[96,101)|g|sieve|56|76|76.0|2|0|239.8|0.080|0.0200|992/747|104.76|104.97|2.42|
|[101,106)|f|sieve|56|76|76.0|2|0|91.7|0.084|0.0200|913/850|104.16|104.65|2.60|
|[101,106)|g|sieve|56|81|75.5|1|1|235.7|0.080|0.0200|1014/749|104.22|104.66|3.30|
|[106,111)|f|sieve|56|81|75.5|1|1|94.4|0.084|0.0200|931/858|103.89|104.62|3.41|
|[106,111)|g|sieve|56|81|81.0|2|0|163.4|0.080|0.0200|1037/752|103.98|104.64|3.85|
|[111,116)|f|sieve|56|81|75.5|2|1|109.1|0.120|0.0203|948/865|110.54|110.55|3.65|
|[111,116)|g|sieve|56|81|75.5|2|1|100.9|0.125|0.0234|1059/754|109.84|109.86|3.75|
|[116,121)|f|sieve|56|86|86.0|2|0|106.3|0.084|0.0200|966/873|103.73|106.25|3.93|
|[116,121)|g|sieve|56|96|83.0|1|1|226.9|0.080|0.0200|1082/757|103.75|106.26|4.27|
|[121,126)|f|sieve|56|96|83.0|1|1|108.1|0.084|0.0200|983/880|103.65|106.31|4.33|
|[121,126)|g|sieve|56|101|85.5|1|1|163.9|0.080|0.0200|1104/759|103.67|106.83|4.27|
|[126,131)|f|sieve|56|96|96.0|2|0|112.9|0.085|0.0200|1001/888|103.70|106.39|4.65|
|[126,131)|g|sieve|56|116|100.2|2|0|143.5|0.080|0.0200|1127/762|103.76|107.29|4.60|
|[131,136)|f|sieve|56|96|83.0|2|1|172.6|0.120|0.0200|1018/895|108.55|108.88|4.76|
|[131,136)|g|sieve|56|96|83.0|2|1|123.8|0.125|0.0200|1149/764|107.33|108.04|4.82|
|[136,141)|f|sieve|56|101|85.5|2|1|201.2|0.120|0.0200|1036/903|108.59|109.03|4.76|
|[136,141)|g|sieve|56|101|85.5|2|1|133.9|0.125|0.0200|1172/767|107.14|108.14|4.82|
|[141,146)|f|sieve|56|96|96.0|3|0|135.5|0.120|0.0200|1053/910|109.44|109.65|5.23|
|[141,146)|g|sieve|64|101|85.5|2|1|101.5|0.125|0.0232|1208/774|107.91|108.58|4.82|
|[146,151)|f|sieve|56|106|93.6|2|1|151.1|0.120|0.0200|1071/918|108.52|109.00|5.44|
|[146,151)|g|sieve|56|101|101.0|3|0|115.0|0.125|0.0200|1217/772|108.60|109.27|5.27|
|[151,156)|f|sieve|64|106|93.6|2|1|123.0|0.120|0.0200|1102/930|109.12|109.45|5.44|
|[151,156)|g|sieve|64|106|98.1|2|1|95.2|0.125|0.0259|1253/779|108.56|109.25|5.48|
|[156,161)|f|sieve|56|116|100.2|2|1|162.1|0.124|0.0200|1106/933|108.89|111.73|5.81|
|[156,161)|g|sieve|64|106|106.0|3|0|104.7|0.126|0.0221|1276/782|109.52|109.92|5.82|
|[161,166)|f|sieve|64|116|100.2|2|1|131.0|0.124|0.0200|1137/945|109.33|111.80|5.81|
|[161,166)|g|sieve|64|116|105.9|2|1|103.1|0.129|0.0228|1298/784|108.86|111.73|6.12|
|[166,171)|f|sieve|64|126|98.0|2|1|95.9|0.124|0.0255|1155/953|108.56|111.69|6.14|
|[166,171)|g|sieve|64|131|100.5|2|1|114.0|0.130|0.0200|1321/787|107.73|111.62|6.38|
|[171,176)|f|sieve|64|126|106.8|2|1|113.4|0.124|0.0200|1172/960|109.79|111.89|6.69|
|[171,176)|g|sieve|64|131|110.0|2|1|102.5|0.124|0.0228|1343/789|108.73|111.71|6.84|
|[176,181)|f|sieve|64|141|105.5|2|1|158.7|0.124|0.0200|1190/968|109.71|111.88|6.65|
|[176,181)|g|sieve|56|131|100.5|3|1|231.7|0.170|0.0200|1352/787|112.53|113.15|6.65|
|[181,186)|f|sieve|64|146|108.0|2|1|114.6|0.124|0.0200|1207/975|109.88|111.91|7.06|
|[181,186)|g|sieve|64|131|100.5|3|1|176.7|0.170|0.0200|1388/794|112.95|113.44|6.65|
|[186,191)|f|sieve|64|151|110.5|2|1|107.8|0.129|0.0212|1225/983|109.92|111.93|7.26|
|[186,191)|g|sieve|56|141|105.5|3|1|177.2|0.170|0.0200|1397/792|111.73|112.68|6.98|
|[191,196)|f|sieve|64|161|115.5|2|1|203.8|0.140|0.0200|1242/990|111.58|113.75|7.61|
|[191,196)|g|sieve|64|141|105.5|3|1|133.7|0.170|0.0200|1433/799|112.17|112.92|6.98|
|[196,201)|f|sieve|64|166|118.0|2|1|140.5|0.140|0.0200|1260/998|111.68|114.05|7.84|
|[196,201)|g|sieve|64|151|110.5|3|1|265.2|0.164|0.0200|1456/802|113.35|113.73|7.56|
|[201,206)|f|sieve|64|171|120.5|2|1|95.9|0.145|0.0267|1277/1005|111.79|114.08|8.04|
|[201,206)|g|sieve|56|146|119.9|3|1|181.3|0.182|0.0200|1464/799|112.99|114.29|7.71|
|[206,211)|f|sieve|64|181|125.5|2|1|170.2|0.145|0.0200|1295/1013|113.92|115.15|8.45|
|[206,211)|g|sieve|64|161|115.5|3|1|225.6|0.184|0.0200|1501/807|112.76|114.37|7.91|
|[211,216)|f|sieve|80|181|125.5|2|1|296.1|0.145|0.0200|1335/1028|115.43|115.99|8.45|
|[211,216)|g|sieve|64|166|118.0|3|1|173.2|0.184|0.0200|1523/809|112.45|114.49|8.15|
|[216,221)|f|ecm|20|171|136.3|5|1|180207.7|0.276|0.0200|1209/987|124.24|124.24|9.22|
|[216,221)|g|sieve|80|171|120.5|3|1|598.9|0.183|0.0200|1569/819|115.70|116.11|8.30|
|[221,226)|f|ecm|20|176|139.6|5|1|285781.7|0.276|0.0200|1227/995|125.07|125.07|9.42|
|[221,226)|g|sieve|64|181|125.5|3|1|294.2|0.190|0.0200|1568/814|113.53|115.16|8.70|
|[226,231)|f|ecm|20|171|148.9|5|1|841344.8|0.276|0.0200|1244/1002|126.59|126.59|9.61|
|[226,231)|g|sieve|80|181|125.5|3|1|455.0|0.190|0.0200|1614/824|115.35|116.03|8.71|
|[231,236)|f|ecm|20|176|139.6|6|1|914269.4|0.314|0.0200|1262/1010|125.67|125.67|9.65|
|[231,236)|g|sieve|64|176|139.6|3|1|329.1|0.190|0.0200|1613/819|115.60|116.04|8.90|
|[236,241)|f|ecm|20|176|139.6|6|1|500633.7|0.320|0.0200|1279/1017|125.67|125.67|9.66|
|[236,241)|g|sieve|48|191|149.4|3|1|281.0|0.190|0.0200|1606/812|114.83|115.79|9.35|
|[241,246)|f|ecm|20|171|148.9|6|1|1451277.0|0.314|0.0200|1297/1025|127.35|127.35|9.88|
|[241,246)|g|sieve|56|191|149.4|3|1|275.4|0.190|0.0200|1644/819|115.48|116.16|9.35|
|[246,251)|f|ecm|20|191|149.4|6|1|4221698.8|0.314|0.0200|1314/1032|127.43|127.43|10.14|
|[246,251)|g|sieve|64|191|149.4|3|1|205.8|0.190|0.0200|1681/827|115.82|116.38|9.35|
|[251,256)|f|ecm|20|191|149.4|6|1|1797699.6|0.314|0.0200|1332/1040|127.43|127.43|10.14|
|[251,256)|g|sieve|64|196|152.7|3|1|174.0|0.198|0.0200|1703/829|115.82|116.58|9.52|
|[256,261)|f|ecm|20|191|149.4|6|1|1015645.8|0.321|0.0200|1349/1047|127.43|127.43|10.15|
|[256,261)|g|sieve|64|201|156.0|3|1|140.3|0.198|0.0200|1726/832|115.81|116.68|9.74|
|[261,266)|f|ecm|20|186|160.6|6|1|4584369.0|0.329|0.0200|1367/1055|129.57|129.57|10.43|
|[261,266)|g|sieve|64|206|159.2|3|1|114.4|0.199|0.0250|1748/834|116.05|117.00|9.99|

## 4. Expected Work, Node Counter, Repeat Bound, and Metered Cap

- **Expected two-log work:** `2^107.448555` <= **`2^107.449`** core units (`2^107.383549` at `MESH = 384/384`), `-1.415` bits vs `c36` (`108.864`).
- **Metered cap:** `t = 2^119.5` core units; by Markov's inequality, `P(cap) <= 2^(107.449 - 119.5) = 2.356609e-04` (`2.592270e-04` at `+10%` headroom).
- **Node counter (`E[N]` live in `audit_cm.py`):** `E[N] = 2^14.4757`, hard abort `N_MAX = 2^25.8` -> `P(counter abort) <= 2^-11.3242`.
- **Repeat-ideal bound:** Smallest bin `[70, 76)` holds `>= 2^67.64 > 2^66.66` distinct degree-1 `tau`-orbits; with `2^6` `1/p`-weighting slack, `P(repeat) <= (N_MAX^2 / (2 x 2^66.66)) x 2^6 = 2^(-10.06) <= 2^-10.05`.
