# Zero-Retry Deterministic-DAG Individual-Logarithm Descent (`CM-STNFS` over `Q(zeta_15)`)

Target operating point: `B_f = 2^71.29`, `B_g = 2^71.69`, `LX = 2^70`, twist `w' = omega^2` (`alpha_f = +2.74`, `alpha_g = +0.77`), **de Bruijn local-density correction enabled on BOTH boot and all special-`q` nodes**, **`boot_rels = 1`**, **`k = 1` (ZERO internal retries inside the DAG!)**, **root-only i.i.d. fresh-`e_j` boot attempts (`B_max = 10`)**, **explicit deterministic per-node special-`q` norm cutoff `|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`**, **expected survivor/ECM work framed under `H-EA`**, **zero-loss degree-preserving excess-row trim + Markowitz deg-0 guard + LIFO `repurge_stack` (`post_peel_cols_lost = post_merge_cols_lost = n_reconstruct_unit = 0`!)**, `RET_BOOT = 0.50`, `UDUP_NODE = 0.95`, `cells_b = 1`, `ball_pair = 0`, explicit `D_skew` boot LLL (`S_boot = 2^21`), and an explicit aborting combined logarithm phase meter **`t = 2^121.4`** core units (`P(D >= t) <= E[D] / t` by Markov).

`cm_data.json` stores ONLY the 78 structural node 9-tuples `plan_struct`, empirical `q`-lattice offsets, the 24-base production `(4, 3)` `Q(zeta_15)` special-`q` witness set (`prod_sq24`), the 224-ideal depth-1 + depth-2 high-mean (`mu_u >= 43.7`) usable-relation asymmetric `(4, 3)` `Q(zeta_12)` campaign (`q12_43_campaign`), and boot config (`m_sv = 2560.0`b, `s_boot = 2^21`). `audit_cm.py` calls `cm_models.eval_descent_live` to compute seeded peeling closure `(h, exc_crit)` and every node's `pr`, `sv`, `mu`, `miss`, `sieve`, expected survivors `cells * (sv + pr)`, `surv_cost`, `own`, `p0`, `Phi_dag`, `W_det`, `C_det`, single-attempt `Phi_1`, root-retry `2*Q_1^10`, and Wald expected work 100% live from scratch. Because lower Riemann sums are strictly monotone in mesh density (`pr_96 <= pr_384 <= pr_exact`), pricing at default `MESH = 96/192` (`E[D] = 2^108.841670 <= 2^108.842`, `Phi_1 = 0.32828`, `2*Q_1^10 = 2.92e-05 <= 0.000030`) is strictly more conservative than on-demand `{"descent_mesh": 384}` (`E[D] = 2^108.817423`, `Phi_1 = 0.31796`).

## 1. Fresh-Coin Zero-Retry i.i.d. Single-Cell Boot (`boot_rels = 1`, `cells_b = 1`, `D_skew = diag(1..1, 2^21..2^21)`)

For each discrete logarithm (`log g` with `g = e(P0, Q0)` and `log z`), draw an independent uniform `e in {1, ..., r - 1}` (fresh coins) and boot on `z' = z^e` (resp. `g' = g^e`), making both inputs exact uniform on `G_T \ {1}`.

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

## 2. Exact Seeded Peeling Closure (`h = 0.256845`), Zero-Loss `repurge_stack` Back-Substitution, & `<= LX = 2^70` Discard

- **Seeded Peeling Fixed Point (`cm_models._peel_fixed_point_bls`, live in `audit_cm.py`):** On the two-sided BLS12-381 factor base (`B_f = 2^71.29`, `B_g = 2^71.69`, `lrels = log2(1.09 U') = 63.441701`, `u_f = 548.9 / 71.29`, `u_g = 677.8 / 71.69`), ideals with `lam >= 100.0` form the protected seed. Unprotected ideals satisfy `h = 1 - exp(-sum_(lam < 100) wt * exp(-lam (1 - h)))`, giving **`h = 0.256845`** and **`exc_crit = 0.7625 << exc = 1.09`** (**`+43.0%` margin** above peeling collapse).
- **Zero-Loss Degree-Preserving Excess-Row Trim + Markowitz Guard + LIFO `repurge_stack` (`RANK-BRIDGE.md` §3):** During excess-row trimming from the post-singleton-peeling 2-core down to `D = N + 4`, CADO-NFS degree-2 clique removal drops a candidate row ONLY when every incident column has degree `>= 2` and every cascading degree-1 peel has a unique degree-1 column in its row (pushed onto LIFO `repurge_stack`; rollback if any row would acquire `>= 2` degree-1 columns), and Markowitz merging rejects any pivot-row cancellation that would drop a nonpivot column to degree `0`. Consequently **`post_peel_cols_lost == post_merge_cols_lost == 0` by deterministic invariant** (`0` unrecovered post-peeling or post-merge columns, `0` back-substitution errors, `n_reconstruct_unit == 0` across `12/1`, `12/6`, `13/1`, `13/2`, `14/1`). After BW solves the `N` active core columns, strictly 3-stack triangular back-substitution (`L_rec` forward -> `merge_stack` reverse -> `repurge_stack` reverse, `< 2^64` u in `filtering_and_misc = 2^80`; **ZERO 4th unit-prop scan and ZERO residual-cycle solver**) recovers **100% of all post-peeling 2-core columns**, making `miss_side`'s post-peeling formula exact with zero post-trim/post-merge loss!
- **Exact `<= LX = 2^70` Leaf Discard & Conservative `> 2^70` Children:** Leaf threshold is **`LX = 2^70`** (`[70, 76)` has `n_c = 70, kA = kB = 0`). Any prime `<= 2^70` absent from the post-peeling 2-core discards the candidate cleanly (`exp(-miss)` integrated over `[40, 70)`, charging `exp(-x)(1 + x)` at `x = lam(1 - h)` for `lam < 100` and `exp(-lam)` for `lam >= 100`). Every prime **`> 2^70`** (including `(2^70, B_s]`) is conservatively counted as an unknown Window-A/B child in `ss_eval` and descended strictly downward.

## 3. Root-Only i.i.d. Boot Retries (`B_max = 10`), Zero-Internal-Retry DAG (`k = 1`), Expected Work under `H-EA`, `t = 2^121.4` Meter, & Exploratory `H5` Evidence

Every degree-1 ideal `q > LX = 2^70` yields a 16D lattice `L_q` (`det = q`); constructive 16D LLL extracts `(v_1 = b_1, v_2 = b_k)` generating a free rank-2 submodule `M <= L_q`, and enumerates `alpha v_1 + beta v_2` (`UDUP_NODE = 0.95`; worst-case offsets `(c_f, c_g) = (+131.66, +46.15)`b; `[7e]` shift `(+3.40, +2.35)`b inside `(+10, +4)`b).

1. **Zero-Internal-Retry DAG (`k = 1`) + Explicit Norm Cutoff + Expected Work under `H-EA` + Ideal-Log Cache (`P(repeat) = 0`):** Inside a single root attempt, every reached special-`q` node `v` enforces the deterministic norm cutoff **`|N_f| <= 2^nf_max, |N_g| <= 2^ng_max`** (dropping any lattice point exceeding the node's modeled envelope `nf_max, ng_max`, analogous to boot's `<= 2560`b cutoff), sieves its fixed ball ONCE, and cofactors qualifying survivors (with expected survivor count `cells * (sv + pr)` and expected sieve-hit / ECM / 100-round MR work priced under declared expected-work model inequality **`H-EA`**, wrapped globally inside the hard aborting meter `t = 2^121.4`). Each node takes the **FIRST** qualifying usable relation (`M_usable(v) = 0 =>` abort this attempt; zero internal retries). Strictly decreasing child norms + a per-attempt ideal-log hash table give `P(repeat) = 0` and deterministic **`C_det_1log = 4,554` nodes/log** (`P(node cap) = 0`).
2. **Unconditional Additive Boole Union Bound under `4.6x`-Weaker Reached-Ideal Ensemble `H5` (`NB(1.5)`):** Under declared reached-ideal ensemble production-law model inequality **`H5`** (`Pr(M_usable(q_sigma) = 0 | slot sigma reached) <= p0(mu_model) = (1 + mu_model / 3)^(-1.5)`, part of `NB(r = 1.5, alpha = 2/3)` PGF `(1 + (mu_model / 3)(1 - t))^(-1.5)` for `t in [0, 0.8]`; `p0_max = 0.00185` at `min_mu = 196.31`, **`4.6x` weaker** than `NB(2)`!), Boole's inequality holds unconditionally for arbitrary slot dependence within a single attempt: **`Phi_1 = 0.328283`** (`0.317964` at `384/384`), giving single-attempt failure **`Q_1 = Phi_1 + 10^(-4) + 2^(-254) = 0.328383 < 0.33`**.
3. **Root-Only Strictly i.i.d. Boot Retries (`B_max = 10`) + Exact Wald Work (`H-EA`) + `t = 2^121.4` Meter (Markov):** For each discrete logarithm, run up to **`B_max = 10` strictly independent fresh-uniform-`e_j` root attempts** (clearing the DAG cache between attempts so attempts `1..B_max` are 100% i.i.d. with zero shared state!).
   - **Exact i.i.d. amplification:** Failure across all `B_max` root attempts is `<= Q_1^10` per log -> **`2 * Q_1^10 = 2.92e-05 <= 0.000030`** for BOTH logs!
   - **Exact Wald's identity (`E[W_j] <= boot_own_exp + W_dag_1` independent of `{N_att >= j}`):** Because attempt `j`'s fresh coins `e_j` and expected work `E[W_j] <= boot_own_exp + W_dag_1` (under `H-EA`) are 100% independent of prior attempts `1..j-1`, expected two-log work is `E[D] <= 2 * ((1 - Q_1^10) / (1 - Q_1)) * (boot_own_exp + W_dag_1) = 1.4889x` single-attempt work = **`2^108.841670 <= 2^108.842`** core units (`2^108.817423` at `384/384`).
   - **Hard Combined-Log Meter `t = 2^121.4` (Markov):** Pointwise `min(D, t) <= D` preserves `E[min(D, t)] <= E[D]` with **Markov meter-abort `P(D >= t) <= E[D] / t = 0.00016583 <= 0.000167`** (zero independence needed).
4. **Two Scoped Empirical Compatibility Checks (`q12_43_descent_live` + `prod_sq24_live` in `audit_cm.py`):**
   - **(a) Structurally Matching Asymmetric Degree-`(4, 3)` `Q(zeta_12)` High-Mean Depth-1 + Depth-2 Actual Missing-FB & ECM Usable-Relation Campaign (`224` ideals; `4` 2-hop chains + `4` 1-hop chains + `4` controls = `16` full 2,640-pt balls LIVE!):** Irreducible `f_4(X) = X^4 + X^2 + 2` (deg 16), `g_3(X) = X^3 + w^2` (deg 12), `LX = 2^17`, full 8D LLL sparse basis combinations `B8_w6` (all `2,640` nonzero `{-1,0,+1}^8` vectors of Hamming weight `<= 6` with first nonzero coord `+1`), single large-prime window (`D1`: `[22, 23)`, `nc = 21.5`; `D2` grandchildren `q_0 -> k=1 usable -> q_1 -> k=1 usable -> q_2`: `[20.5, 21.5)`, `nc = 20.0`), **ACTUAL per-leaf missing-FB lookup** (`P_MISS` = product of `437` missing primes in `(2^15, 2^17]`, `gcd(nf, P_MISS) == gcd(ng, P_MISS) == 1`), **ACTUAL per-factor SplitMix64 ECM thinning** (`p_ecm = 0.996` per prime factor `> 2^10`), first-principles Murphy `alpha = (+3.55, +3.04)`, and pre-calibrated **`qbins`/`qbins_sub`/`qbins_lo` usable norm-offset envelope** (NEVER per-cell exact norms!: **`mu_u_d1 = 43.71 g / 106.12 f`**, **`mu_u_d2 = 45.95 g / 96.15 f`** — every slot `>= 43.7`!; `sv_up_d1 = 990.5 / 2417.3`, `sv_up_d2 = 1032.0 / 2353.5`). `audit_cm.py` (`q12_43_descent_live`) reconstructs **4 complete 2-hop chains (`q_0 -> q_1 -> q_2`) + 4 1-hop chains + 4 controls (`16` full balls) LIVE** (`exact_live_match = True`) and audits all 224 ideals against **BOTH `NB(1.5)` AND `NB(1.0)`** PGFs plus actual two-sided survivor upper bounds:

     | Cohort | Ideals | `M_u` range | Empty | `M_u/mu_u` | `var` vs `NB(1.5)` | `sv2/sv_up` (max) | PGF `emp/nb15` (`0/.2/.5/.8`) | PGF `emp/nb10` (`0/.2/.5/.8`) |
     |---|---:|---|---:|---:|---|---:|---:|---:|
     | `ALL_224` | 224 | `15..311` | **`0/224`** | `1.8575` | `0.4122` vs `2.330` | **`0.0832`** (`0.1938`) | `0 / 0 / .000005 / .002655` | `0 / 0 / .000004 / .003104` |
     | `D1_100` | 100 | `23..291` | **`0/100`** | `1.9462` | `0.4103` vs `2.557` | **`0.0875`** (`0.1938`) | `0 / 0 / 0 / .000756` | `0 / 0 / 0 / .000884` |
     | `D2_24` | 24 | `32..311` | **`0/24`** | `1.6434` | `0.3354` vs `1.827` | **`0.0754`** (`0.1334`) | `0 / 0 / 0 / .000398` | `0 / 0 / 0 / .000463` |
     | `CTRL_100` | 100 | `15..255` | **`0/100`** | `1.8201` | `0.4197` vs `2.238` | **`0.0807`** (`0.1747`) | `0 / 0 / .000011 / .005092` | `0 / 0 / .000010 / .005956` |

     Slot-vs-control ratios: **`D1/CTRL = 1.0693`**, **`D2/CTRL = 0.9029`** (with `D2` `M_u/mu_u = 1.6434 >> 1.0`, `M_min = 32`, zero empties!), actual two-sided **`sv2 / sv_upper = 0.0832 < 1.0`** (`max = 0.1938 << 1.0` across all 224 ideals)!
   - **(b) 24-Base Actual Production BLS12-381 `(4, 3)` `Q(zeta_15)` `(f_9, g_3)` Special-`q` 16D Lattice Check (`audit_cm.py`: `prod_sq24_live`):** Across `12 g-side + 12 f-side` production primes spanning `73..258` bits (`384` exact production deg-72/deg-24 norms): `det(v1, v2) % q == 0` and exact `q | N` verified LIVE; all 6 empirical sextiles `p16..max` lie strictly below `qbins` `c_hi + cf, c_hi + cg` (`g-side` headroom `>= 39.68`b `f` / `>= 23.87`b `g`; `f-side` `>= 37.39`b `f` / `>= 27.31`b `g`); single-ball modeled survivor ratio `sv_emp / sv_plan` has mean **`0.0503` (`g`, max `0.1226`)** and **`0.2602` (`f`, max `1.0000`)**.
   - **Explicit Exploratory Scope:** Together, `0/224` high-mean depth-1 & depth-2 actual-missing-FB+ECM usable empties + BOTH `NB(1.5)` & `NB(1.0)` PGF domination + `sv2/sv_upper <= 0.1938 < 1.0` + the 24-base production `(4, 3)` norm-envelope check provide scoped finite-scale compatibility evidence for declared production-law inequalities `H5` and `H-EA` (NOT statistical proof of production tails); production tail and expected-work transfer rest explicitly on declared exploratory model inequalities `H5` and `H-EA`.

## 4. Complete Live-Computed Deterministic-DAG Plan (All `39` Bins / `78` Nodes, Default `MESH = 96/192`, Full de Bruijn ON)

|Bin|S|How|R|n|nA|kA|kB|mu|p0|Phi_dag|Nf/Ng|own|W_det|C_det|
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
|[70,76)|f|sieve|64|70|70.0|0|0|3769.2|2.2e-5|0.00002|823/810|106.59|106.59|1|
|[70,76)|g|sieve|64|70|70.0|0|0|3755.2|2.3e-5|0.00002|893/740|106.59|106.59|1|
|[76,81)|f|sieve|56|76|76.0|1|0|8460.1|6.7e-6|0.00005|826/813|103.59|107.67|3|
|[76,81)|g|sieve|56|76|76.0|1|0|8316.5|6.8e-6|0.00005|902/737|103.60|107.68|3|
|[81,86)|f|sieve|56|76|76.0|1|0|2252.2|4.9e-5|0.00009|843/820|103.53|107.67|3|
|[81,86)|g|sieve|64|76|76.0|1|0|6523.0|9.9e-6|0.00005|938/744|106.60|108.18|3|
|[86,91)|f|sieve|56|76|76.0|2|0|3739.0|2.3e-5|0.00011|861/828|107.35|109.10|5|
|[86,91)|g|sieve|48|76|76.0|2|0|886.9|2.0e-4|0.00029|931/737|104.92|108.70|5|
|[91,96)|f|sieve|56|76|76.0|2|0|1008.9|1.6e-4|0.00025|878/835|106.01|108.81|5|
|[91,96)|g|sieve|56|76|76.0|2|0|895.8|1.9e-4|0.00028|969/744|105.88|108.79|5|
|[96,101)|f|sieve|64|76|76.0|2|0|840.3|2.1e-4|0.00030|910/848|107.43|109.12|5|
|[96,101)|g|sieve|64|76|76.0|2|0|700.7|2.8e-4|0.00037|1006/752|107.26|109.07|5|
|[101,106)|f|sieve|80|76|76.0|2|0|1551.2|8.5e-5|0.00017|950/863|111.83|111.97|5|
|[101,106)|g|sieve|64|81|75.5|1|1|684.5|2.9e-4|0.00044|1028/754|106.97|109.50|9|
|[106,111)|f|sieve|80|81|75.5|1|1|1505.8|8.9e-5|0.00024|968/870|111.79|112.02|9|
|[106,111)|g|sieve|64|81|81.0|2|0|473.3|5.0e-4|0.00071|1051/757|106.81|109.86|13|
|[111,116)|f|sieve|64|81|75.5|2|1|327.8|8.6e-4|0.00106|962/870|112.64|112.81|11|
|[111,116)|g|sieve|64|81|75.5|2|1|235.3|1.4e-3|0.00161|1073/759|111.80|112.09|11|
|[116,121)|f|sieve|80|86|86.0|2|0|1665.9|7.6e-5|0.00037|1003/885|111.76|112.12|13|
|[116,121)|g|sieve|80|96|83.0|1|1|4259.2|1.9e-5|0.00070|1119/769|111.77|112.28|17|
|[121,126)|f|sieve|80|96|83.0|1|1|1750.5|7.1e-5|0.00076|1020/893|111.75|112.27|17|
|[121,126)|g|sieve|80|101|85.5|1|1|3161.3|2.9e-5|0.00085|1141/772|111.76|112.30|17|
|[126,131)|f|sieve|80|96|96.0|2|0|1726.1|7.2e-5|0.00115|1038/900|111.76|112.41|21|
|[126,131)|g|sieve|80|116|100.2|2|0|2629.5|3.8e-5|0.00138|1164/774|111.77|112.47|21|
|[131,136)|f|sieve|64|96|83.0|2|1|513.3|4.4e-4|0.00128|1032/900|110.73|111.84|23|
|[131,136)|g|sieve|64|96|83.0|2|1|280.5|1.1e-3|0.00192|1163/769|109.50|111.40|23|
|[136,141)|f|sieve|64|101|85.5|2|1|598.8|3.5e-4|0.00132|1050/908|110.53|111.80|23|
|[136,141)|g|sieve|64|101|85.5|2|1|301.2|9.8e-4|0.00195|1186/772|109.18|111.38|23|
|[141,146)|f|sieve|64|96|96.0|3|0|402.1|6.4e-4|0.00225|1067/915|111.42|112.48|31|
|[141,146)|g|sieve|80|101|85.5|2|1|504.0|4.6e-4|0.00142|1231/782|112.21|112.73|23|
|[146,151)|f|sieve|64|106|93.6|2|1|446.8|5.4e-4|0.00236|1085/923|110.47|112.99|35|
|[146,151)|g|sieve|80|101|101.0|3|0|1662.1|7.6e-5|0.00209|1254/784|113.91|114.19|31|
|[151,156)|f|sieve|80|106|93.6|2|1|820.1|2.2e-4|0.00203|1125/938|112.86|113.79|35|
|[151,156)|g|sieve|80|106|98.1|2|1|467.6|5.1e-4|0.00259|1276/787|112.50|113.64|35|
|[156,161)|f|sieve|64|116|100.2|2|1|477.7|4.9e-4|0.00450|1120/938|110.79|113.92|45|
|[156,161)|g|sieve|80|106|106.0|3|0|486.2|4.8e-4|0.00270|1299/789|113.16|114.52|43|
|[161,166)|f|sieve|80|116|100.2|2|1|865.9|2.0e-4|0.00421|1160/953|113.04|114.43|45|
|[161,166)|g|sieve|80|116|105.9|2|1|496.2|4.7e-4|0.00461|1321/792|112.67|114.75|53|
|[166,171)|f|sieve|80|126|98.0|2|1|601.0|3.5e-4|0.00435|1178/960|112.56|114.32|55|
|[166,171)|g|sieve|80|131|100.5|2|1|520.5|4.3e-4|0.00453|1344/794|112.11|114.25|63|
|[171,176)|f|sieve|80|126|106.8|2|1|742.4|2.6e-4|0.00494|1195/968|113.31|114.99|79|
|[171,176)|g|sieve|80|131|110.0|2|1|648.9|3.1e-4|0.00508|1366/797|112.61|114.84|87|
|[176,181)|f|sieve|80|141|105.5|2|1|991.0|1.7e-4|0.00491|1213/975|113.29|114.98|75|
|[176,181)|g|sieve|64|131|100.5|3|1|500.5|4.6e-4|0.00522|1366/792|114.61|115.33|73|
|[181,186)|f|sieve|80|146|108.0|2|1|742.0|2.6e-4|0.00647|1230/983|113.40|115.09|99|
|[181,186)|g|sieve|80|131|100.5|3|1|846.4|2.1e-4|0.00497|1411/802|116.51|116.74|73|
|[186,191)|f|sieve|80|151|110.5|2|1|501.7|4.6e-4|0.00692|1248/990|113.44|115.55|111|
|[186,191)|g|sieve|64|141|105.5|3|1|378.3|7.0e-4|0.00618|1411/797|113.84|115.35|89|
|[191,196)|f|sieve|80|161|115.5|2|1|986.9|1.7e-4|0.01269|1265/998|114.80|116.47|137|
|[191,196)|g|sieve|80|141|105.5|3|1|626.9|3.3e-4|0.00581|1456/807|115.68|116.28|89|
|[196,201)|f|sieve|80|166|118.0|2|1|644.4|3.2e-4|0.01474|1283/1005|114.96|116.69|159|
|[196,201)|g|sieve|64|151|110.5|3|1|251.7|1.3e-3|0.00875|1456/802|113.35|115.68|133|
|[201,206)|f|sieve|80|171|120.5|2|1|456.1|5.3e-4|0.01496|1300/1013|115.00|116.70|179|
|[201,206)|g|sieve|64|146|119.9|3|1|380.7|6.9e-4|0.01288|1478/804|115.05|116.36|145|
|[206,211)|f|sieve|80|181|125.5|2|1|771.9|2.4e-4|0.01573|1318/1020|117.09|117.86|235|
|[206,211)|g|sieve|80|161|115.5|3|1|1002.9|1.6e-4|0.01535|1524/814|116.22|117.20|161|
|[211,216)|f|sieve|80|181|125.5|2|1|204.5|1.7e-3|0.01722|1335/1028|115.43|117.12|235|
|[211,216)|g|sieve|80|166|118.0|3|1|791.3|2.3e-4|0.01732|1546/817|115.92|117.18|189|
|[216,221)|f|ecm|20|171|136.3|5|1|128915.8|1.1e-7|0.02543|1209/987|124.24|124.24|349|
|[216,221)|g|sieve|80|171|120.5|3|1|566.3|3.8e-4|0.01747|1569/819|115.70|117.09|209|
|[221,226)|f|ecm|20|176|139.6|5|1|203540.6|5.7e-8|0.02635|1227/995|125.07|125.08|397|
|[221,226)|g|sieve|64|181|125.5|3|1|196.3|1.8e-3|0.02000|1568/814|113.53|116.90|269|
|[226,231)|f|ecm|20|171|148.9|5|1|595310.3|1.1e-8|0.03133|1244/1002|126.59|126.59|449|
|[226,231)|g|sieve|80|181|125.5|3|1|298.8|9.9e-4|0.01914|1614/824|115.35|117.22|269|
|[231,236)|f|ecm|20|176|139.6|6|1|866282.7|6.4e-9|0.02961|1262/1010|125.67|125.68|443|
|[231,236)|g|sieve|64|176|139.6|3|1|216.5|1.6e-3|0.02141|1613/819|115.60|117.20|305|
|[236,241)|f|ecm|20|176|139.6|6|1|350349.0|2.5e-8|0.02961|1279/1017|125.67|125.68|443|
|[236,241)|g|sieve|56|191|149.4|3|1|722.9|2.7e-4|0.02671|1622/817|117.20|118.57|399|
|[241,246)|f|ecm|20|171|148.9|6|1|1372964.7|3.2e-9|0.03578|1297/1025|127.35|127.35|515|
|[241,246)|g|sieve|64|191|149.4|3|1|564.8|3.8e-4|0.02682|1658/824|117.50|118.69|399|
|[246,251)|f|ecm|20|191|149.4|6|1|3991958.1|6.5e-10|0.03978|1314/1032|127.43|127.43|597|
|[246,251)|g|sieve|80|191|149.4|3|1|872.7|2.0e-4|0.02664|1704/834|119.20|119.68|399|
|[251,256)|f|ecm|20|191|149.4|6|1|1698034.8|2.3e-9|0.03978|1332/1040|127.43|127.43|597|
|[251,256)|g|sieve|80|196|152.7|3|1|750.0|2.5e-4|0.03396|1726/837|119.19|119.78|437|
|[256,261)|f|ecm|20|191|149.4|6|1|696640.2|8.9e-9|0.03978|1349/1047|127.43|127.43|597|
|[256,261)|g|sieve|80|201|156.0|3|1|584.1|3.7e-4|0.03870|1749/839|119.18|119.79|503|
|[261,266)|f|ecm|20|186|160.6|6|1|4323408.0|5.8e-10|0.05485|1367/1055|129.57|129.57|715|
|[261,266)|g|sieve|80|206|159.2|3|1|479.9|4.9e-4|0.04990|1771/842|119.41|119.97|589|
