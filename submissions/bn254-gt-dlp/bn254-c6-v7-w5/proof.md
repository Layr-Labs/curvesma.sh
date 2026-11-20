# bn254-gt-dlp: Pure-Integer `C6-m184_w216-rec` `(6,2)` STNFS (Pure-Product `kap=0` Deterministic-Occ Band-1/2 Direct-Apply Sieve, 64-Cell ORR, Radix-4+3 `<2^62` Shoup NTT, 4-Vec `2^24` `Z[omega]` LA, Conditional-Cantelli Zero-Retry DAG + 8-Bit ECM Ladder) — Expected Work `2^86.496544115` Units

Track `bn254-gt-dlp`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `x`. Score unit: `C = 20,000` pure baseline scalar ARMv8.0-A integer instructions (`-mgeneral-regs-only`, `simd = 0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, `p_num = 0`, zero external C/Callgrind measurements). All derivations are self-contained here and in `RANK-BRIDGE.md`, `STORAGE.md`, `KERNELS.md`, and verified live in `audit_bn_cm.py` (no prior verdict/score is a premise). No advice, smaller targets, or amortization. Analytic, exploratory claim.

## Claim

Under declared exploratory STNFS premises `H1`, `H2`, `H3`, `H5`, `H-EA`, `H-SM`, `H-R`, `H-X`, `H-SEP` (§6, `RANK-BRIDGE.md`; zero `H-FFT`, zero Band-3 `H1'`, zero Poisson Band-2 occupancy heuristic), using exact 3-class 2-adic `rho_dB` on pure-product `s = 0.80, kap = 0.0` (`200,000` exact-HNF cells), 100% pure-integer Band-1/2 direct-apply sieve (`W_SB = 128`, 1-bit `occ_bm[]`, deterministic `p_occ = min(1.0, mu_sb)`, 64-cell `ORR <= 0.96`, zero-sort 2-pass lists, Radix-4+3 `<2^62` Shoup NTT, `7`-limb `f`-SRT), 4-vec register-resident `Z[omega]` LA (`2^12..2^24` 7-scale ladder, ZERO TABM), 5-SM with deterministic shifted Beckermann–Labahn $\sigma$-basis extraction (`RANK-BRIDGE.md` §1), and 62-row `LB = 59` zero-retry deterministic-DAG (`k = 1`) priced by distribution-free **Conditional Cantelli** (`kappa_disp = 3.0`) + 8-bit ECM ladder + `B_max = 24` i.i.d. per-trial fresh-`e` boot (`t = 2^112` meter), the attack succeeds w.p. `>= 0.39` (opened-run success `>= 0.999997` via `FIX-U`), and

    E[W] <= 1 + (K / 2^64) * 2^87.854993757 <= 2^86.496544115 units (C = 20,000 Ir),
    K = 7194251771502039637.

Open iff fresh 64-bit coin `u < K` (least integer with `(K/2^64)*0.999997 >= 0.39`; closed run fails at `1` unit). `86.496544115` bounds strictly **expected** work `E[W]` (not every-run). Meters: two-log `t = 2^112` (`E[D] <= 2^85.393743`, Markov `9.789e-9`); split `2^112` (`8.306e-10`). Work rounds up; `success` rounds down; finite prototypes give scoped compatibility evidence, not production proofs.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Exact `sigma(w) == w`, & Live Mod-`r` CRT

BN254: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h` & twist `w`:** Totally imaginary sextic over `Q(sqrt(-2))` (unit rank `2`, order-3 `sigma`): `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -2^9 * 19^2`, `h_K = 1`, `w_K = 2`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`. `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3`, `w(y) = -y^5 + y^4 - y^2 - 2y` (`N(w) = 216`, `|w_j| = sqrt(6)`), `(u*w)^((p^6-1)/2) = -1` in `F_{p^6}`, and **`w(sigma(y)) == w(y)` identically in `Q[y]/(h(y))`** (live).
- **Reciprocal pair & order-6 `tau`:** `f_rec(x) = (-2y^5+2y^4-2y^2-4y+3)x^8 + (3y^5-3y^4+3y^2+6y)x^6 - 12x^4 + (-3y^5+3y^4-3y^2-6y+12)x^2 + (2y^5-2y^4+2y^2+4y-5)`, `g_rec(x) = u*w*x^2 - 1` (`N(lead)=N(const)=729`, `f_rec(x^2->(u*w)^-1) = 0 mod (p,h)`). `tau(y,x) = (sigma(y), -x)` (`Frob^10`, `chi = p^10 mod r`, `chi^2-chi+1 = 0 mod r`); mod-`r` DDF: `h:[1^6]`, `g_rec:[2]`, `f_rec:[2^4]`. Pure-int `G_sym`/`D_fi` Voronoi FD certified for ALL radii via Dedekind `O_lat == O_h`, `Fraction` `LDL^T` (`5,670` pts, `382/382` orbits), and outward-rounded `Decimal(55)` index-1 & convex-hexagon theorem (`KERNELS.md` §0).
- **Exact 3-Class 2-Adic Split & Field Map:** `2 = P_2^2` (`N(P_2)=8, v_{P_2}(w)=1`) partitions primitive `(a,b)` into `(v_2(N_f), v_2(N_g)) = (0,0)` (`7/9`), `(0,3)` (`1/9`), `(6,0)` (`1/9`); `alpha_f_odd = ALPHA_F_RAW + 2/3`, `alpha_g_odd = ALPHA_G_RAW + 1/3`. `phi_iso: F_{p^{12}} -> F_p[y,x]/(h, uwx^2 - 1)` (`wrapper <= 2^53.712288`); uniform `e in [0,r-1]` gives `T = phi_iso(z*g^e) = c_0 + c_1 x`, `R = c_1/c_0 mod (p,h)` (`c_1 = 0` w.p. `1/r` in `FIX-U`), `x = (vlog(T) - e*vlog(g))*vlog(g)^-1 mod r`, checked by `g^x == z` (`verify <= 2^11.000001`).

---

## 2. Pure-Product (`kap = 0.0`) Deterministic-Occ Band-1/2 Direct-Apply Collection (`12.595` Ir/Cell)

- **200k exact-HNF pool (`s = 0.80, kap = 0.0`):** `B_f = 2^61`, `B_g = 2^60` (`lx2 = P_mid = 2^54.5`, `sl = 160, sb = 165`), `U' = 2^53.641526`, target `0.82*U'` (`+0.0666`b). Exact 3-class `rho_dB` + class-matched `Retention` across `136,563` prim + `63,437` nonprim == `200,000` active cells (live `1,500`-sample HNF audit `|z| <= 1.43`): `V = 2^97.8362` (`SE = 0.0042`b), `V_b = 2^52.3194`, `N_a = 2^45.5168`, `N_5 = 2^43.6889`, `L_b = 396.33`, `L_max = 807.40`, `N_f = 402.07`b (**AM-GM `sup log2|N_f| <= 441.66 < 448`b = 7 limbs**, `+6.34`b headroom), `N_g = 472.22`b, `r_g = 0.6254`.
- **Unified Single-Arena Band-1/2 Sieve (`P_min=14, P_cut=45.0, P_mid=54.5, W_SB=128`, ZERO Band 3; `KERNELS.md` §1):**
  - *Band 1 `[2^14, 2^45)` (`rep_lo = 5.024` Ir/cell):* Single-owner `(b_1..b_5)` row arena (`2^{56.18} B << 2^{64} B`), ZERO `M_t[t]`, 8x-unrolled `APPLY8` (`30/8 = 3.75`) + worst-case `C_B1_ROW = 24.0` Ir/(p,b)-row.
  - *Band 2 `[2^45, 2^54.5)` with **Deterministic Occupied-SB Bound `p_occ = min(1.0, mu_sb)`** (`probe=0.417, occ_sb=1.208, hi_ap=3.124` Ir/cell):* `W_SB=128` (`max_sb=8`, `cr2=1.329`), 1-bit `occ_bm[]` 12-op `b2_occ_probe = 12.0`, `b2_sb_enter = 14.0`, 7-op `b2_cand_prod = 7.0`, 6-op `b2_hit_direct <= 7.0`, unified 64-cell `ORR` scan+zero `59/64 = 0.9219 <= 0.96` (`0 / 128,000` errs live). **ZERO Poisson occupancy heuristic:** for ANY prime `p` with `S_p = ceil(p / W_SB)` super-buckets and `N_a` `a`-entries (regardless of residue structure), `#(non-empty SBs) <= min(S_p, N_a)` deterministically (`det_occ_bound_ok = True` across all 24 prototype ideals), so occupied-SB density `<= min(1.0, N_a / S_p) = min(1.0, mu_sb)`; under `H2` equidistribution of `b`-row starting residues mod `p`, entered-SB rate `<= min(1.0, mu_sb)` (mean `0.6067`).
- **Zero-Sort Sequential 2-Pass Lists (`lists <= 2^84.688287`, `KERNELS.md` §1.3):** ZERO `Q(a)` sort, ZERO Pass 3: `C_LIST_B1 = 20.0`, `C_LIST_B2 = 24.0`, `C_LBUCKET = 6.0` (builder RAM `2^52.32 B << 2^64 B`; `zero_sort_exact_ok = True` live).
- **Stage B, Radix-4+3 `< 2^62` Shoup NTT, & `7`-Limb `f`-SRT (`cof_b = 1.862` Ir/cell):** `C_NORMS = 30k`, `417`-chunk `P_14` strip (`C_SB = 24k`), `4` disjoint LP windows (`0` fn), Radix-4+3 `<2^62` 3-prime Shoup NTT (`60*n*k + 220*n` / `60*n_3*(k-2) + 355*n_3`, 43-op R4, live `n=16` & `n=24` convolutions `ntt3_r4_r3_conv_ok = True`), `7`-limb `f`-SRT (`sup_nf = 441.66 < 448`b, `Q_F_SAFETY = 2.0`), nonprim `0.0396` (`6` witnesses). Total: **`12.595` Ir/cell** `=> collection <= 2^87.208049`, global RAM **`2^98.394 B < 2^100 B`**.

---

## 3. `Z[omega]` Filtering (`2^12..2^24` 7-Scale Ladder), 4-Vec `H3` LA, & Live Shifted $\sigma$-Basis Extraction

| Scale `U'` | Engine | `N_sw / U'` | `N_l / U'` | `k5_exp` (`8.00/58`) | `tab_r` | `w_mid` |
|---|---|---:|---:|---:|---:|---:|
| `2^12` | **LIVE Python** (`W <= 161`, `kern=0`, `lost=0`) | `0.014160` | `0.009277` | `20.23` | `0.7760` | `14.07` |
| `2^13` | **LIVE Python** (`W <= 161`, `lost=0`) | `0.006836` | `0.004028` | `19.51` | `0.7806` | `13.71` |
| `2^20` | C `filtsim` (`n = 4` seeds) | `0.001271` | `0.000751` | `13.84` | `0.7482` | `11.54` |
| `2^21` | C `filtsim` (`n = 8` seeds) | `0.001031` | `0.000630` | `13.61` | `0.7436` | `11.25` |
| `2^22` | C `filtsim` (`n = 4` seeds) | `0.000923` | `0.000570` | `13.34` | `0.7411` | `10.96` |
| `2^23` | C `filtsim` (`n = 4`, max `0.000808/0.000510/13.20`) | `0.000761` | `0.000489` | `13.15` | `0.7390` | `10.67` |
| `2^24` | C `filtsim` (`n = 4`, max `0.000695/0.000452/12.92`) | `0.000673` | `0.000441` | `12.90` | `0.7363` | `10.38` |
| **Hard Caps (`H3`)** | **Deterministic Abort (`W <= 161`, ZERO TABM)** | **`0.000850`** | **`0.000550`** | **`13.50`** | **`1.0000`** | **`~0.0`** |

- **Fixed-Instance `H3` & 4-Vec `Z[omega]` LA:** Monotone across all 7 scales; every `2^23..2^24` run lies strictly below the 4 caps (`1.25x–1.26x` headroom above `2^24` mean). `unit4_reg_z6 = 32/4 = 8.00` (20 GPRs `x8..x27`, zero spills), `gen1_reg_z6 = 56 <= 58.0`, `row_tail_3bucket_z6 = 722 <= 740.0` (1-step Barrett, ZERO TABM), `bw = (1024, 35)`, `nslots = 7`: **`la <= 2^83.601174`**, **`gen <= 2^76.777804`**, 1-owner RAM **`2^59.62 B << 2^64 B`** (`lazy = 0`). Degree-preserving trim + LIFO `repurge_stack` loses `0` post-peel cols (live at `2^13` and `185` genuine `C6` rels).
- **5-SM & Live `18x18` Shifted Beckermann–Labahn $\sigma$-Basis (`RANK-BRIDGE.md`, `hx_z6.py`):** $S_{\mathrm{full}} \in \mathbb{F}_r^{R \times 5}$ (`sm <= 2^73.827705`), conditional Schwartz-Zippel `4 N_sw_cap / 2^64 = 2.590e-6`; shifted order-$L$ $\sigma$-basis extraction theorem + 9-limb Kronecker `PM-Basis` bound proved deterministically inside `gen_at_cap` + `la`, verified LIVE on unplanted `F_1187^*` (`18x18` order-19 $\sigma$-basis, `4` SM-nonzero cols, `4/4` `H*f=0`, ZERO dense RREF recovery, `[407]` across all `68/68` purged rels) + `Z[omega]` Schur identity.

---

## 4. 62-Row Zero-Retry DAG (`k = 1`) + Conditional Cantelli (`kappa_disp = 3.0`) + 8-Bit ECM Ladder + `B_max = 24` Boot (`LB = 59`)

- **Zero-Retry Deterministic-DAG (`k = 1`):** Each of the `62` strictly decreasing rows scans its ball once and takes its first usable relation (`k = 1`); `M = 0` fails the root attempt immediately. Unconditional Boole's inequality propagates failure with ZERO internal retries, ZERO after-failure hazard, ZERO intra-DAG independence, deterministic node count (`c_det_one = 411`, `c_det_two_max = 19,728`), and acyclic per-run cache (`repeat_prob = 0.0`).
- **Distribution-Free Conditional Cantelli / Paley-Zygmund Bound (`kappa_disp = 3.0`, ZERO Negative-Binomial Ansatz):**
  For any special-`q` lattice `L_q = bareiss_lll_12d(Q)`, let `M >= 0` be its usable-relation count, `mu(L_q) = E[M | L_q]`, and `D(L_q) = Var(M | L_q) / mu(L_q)` be its within-lattice Fano dispersion. By the exact **distribution-free one-sided Chebyshev / Cantelli inequality** (zero Poisson or Negative-Binomial assumption) and Law of Total Probability:
  $$P(M = 0 \mid L_q) \le \frac{\operatorname{Var}(M \mid L_q)}{\operatorname{Var}(M \mid L_q) + \mu(L_q)^2} = \frac{1}{1 + \mu(L_q) / D(L_q)}, \qquad P(M = 0 \mid \text{slot}) \le E\!\left[\frac{1}{1 + \mu(L_q) / D(L_q)}\;\Big|\;\text{slot}\right].$$
  Moreover, by the Law of Total Variance $\operatorname{Var}(M) = E[\operatorname{Var}(M \mid L_q)] + \operatorname{Var}(E[M \mid L_q]) \ge E[\operatorname{Var}(M \mid L_q)]$, so any empirical cross-lattice Fano ratio $\operatorname{Var}(M)/E[M]$ upper-bounds the average within-lattice Fano dispersion $E[\operatorname{Var}(M \mid L_q)]/E[M]$ (it does not prove the nonlinear expectation bound over rare low-yield reached lattices). The quantitative per-slot inequality `P(M = 0 | slot) <= E[1 / (1 + mu(L_q) / D(L_q)) | slot] <= 1 / (1 + mu_model / 3.0)` (`kappa_disp = 3.0`, `3x` Poisson variance headroom) is therefore an **explicitly declared exploratory STNFS premise (`H5`, §6)**: every node charges **`p_0 = 1 / (1 + mu_model / 3.0)`** (`mu_min = 713.91`, `p0_max = 4.185e-3`).
- **Three-Pillar Live Exploratory Compatibility Support (`audit_bn_cm.py`, `24/24` checks):**
  1. **Above-Cutoff (`B_FB = 2^{10} = 1,024`) Literal Degree-8 `f_rec` (deg 48) + Degree-2 `g_rec` (deg 12) Mini-BN254 `C6` Campaign (`deg8_cantelli_and_prod_yield`, `u_t = 5, p_t = 27631, r_t = 27481` PRIME, 2-window `W_A = (2^{10}, 2^{16}]` `kA <= 6`, `W_B = (2^{16}, 2^{26}]` `kB <= 2`):** ALL `100` special-`q` ideals (`25` `f_ctrl`, `25` `f_d1`, `25` `g_ctrl`, `25` `g_d1`, `q > B_FB = 1,024` strictly above the factor-base cutoff, `k=1` child selection `1,024 < q_1 < q_0`; `bareiss_lll_12d(9/10)`, `is_ideal_primitive_oh`, projective dedup) + **2 live 2-hop depth-2 chains strictly above `B_FB = 1,024`** (`g(6089)->f(5281,M=5)->g(2377,M=3)` and `f(6529)->g(4483,M=5)->f(3187,M=7)`, replayed bit-exact with `2` controls): **`0/100` empty**, cross-lattice Fano `Var(M)/E[M]` (`>= E[Var(M|L_q)]/E[M]`): **`f_ctrl = 0.742`**, **`f_d1 = 0.902`** (`f_all = 0.807`), **`g_ctrl = 0.846`**, **`g_d1 = 1.381`** (`g_all = 1.106`) — **`2.2x–4.0x` below `kappa_disp = 3.0`**, zero adverse `k=1` shift (`f d1/ctrl = 0.961`, `g d1/ctrl = 1.073`).
  2. **Production `60..215`-Bit Between-Lattice Yield Floor (ALL `24/24` Prod `C6` Special-`q` Lattices, LIVE):** Matched to exact `plan59` rows, `mu(L_q) / mu_model` has **`min = 2.43x, med = 6.08x`** under conservative `ALPHA_F = 0.0, ALPHA_G = +3.4`, and **`min = 7.4x, med = 25.1x`** (`mu_true(L_q) >= 5,280`, med `~17,900`) under actual `f_rec`/`g_rec` Murphy `alpha` (`-19.93 / -0.04`).
  3. **End-to-End Strict Depth-2 `224`-Ideal `C6` Campaign (`c6_pgf_campaign`, LIVE 8-ideal + 2-hop replay):** `0/224` empty, `M/mu = 1.2395`, `S_obs/sv_up = 0.1961 < 1`, residual Fano `res_v/m`: **`all = 1.927, d1 = 1.624, d2 = 2.107, ctrl = 2.042 < 3.0`**.
  *(Exploratory scope: finite campaigns and production norm-model audits provide scoped compatibility evidence, not proofs of production slot tails; transfer of the per-slot Conditional Cantelli bound to every reached production slot is an explicit exploratory `H5` premise.)*
- **8-Bit ECM Ladder & Per-Trial Fresh-`e` `i.i.d.` Boot (`boot = [163, 88, 12, 4]`, `ea = [5, 1250]`, `B_max = 24`):**
  Every node and boot trial charges `extraction_ladder(lo, nA, nbits, 0.999, step=8) + (kA + kB)*PRIME`, audited live under `_ecm_best` with upper-endpoint `[lo, hi_end]` enclosure: **`min_dag = 0.6174 > 0.50`** (`208` DAG rungs), **`min_boot = 0.6341 > 0.50`** (`12` Boot rungs to `154.33`b, `tau = 1250`). Boot prime-product tree (`x_EA = 64.4167`, `2*P_x = 2^62.945 B < 2^64 B`) is built once. For each log (`z`, `g`), run `<= B_max = 24` strictly independent root attempts: each trial (`1..T_boot_max`, `eps_boot = 1e-6`) draws a **fresh uniform `e in [0, r-1]`**, runs a **fresh 2-stage `11/40 -> 9/10` Bareiss 12D LLL (`<= 2^40` Ir charged per trial in `trial_b`; joint `2^{39.81} < 2^{40}` live)**, and tests a uniform radius-4 vector (`48*T_boot_max/r = 1.696e-57` in `FIX-U`). `Phi_1 = phi_dag_one + 1e-6 = 0.463929 <= 0.50`: **`2 * Phi_1^{24} = 1.976e-8`**, `E[attempts/log] <= 1.8654`, **`E[D] <= two_logs = 2^85.393743`** (`two_dag = 2^84.478`, `two_boot = 2^84.305`), metered by `t = 2^112` (`Markov 9.789e-9`).

---

## 5. Itemized 15-Row Ledger, `FIX-U`, & Exact Rational Gate

| Row | Comp. (`log2` u) | Decl. (`log2` u) | Row | Comp. (`log2` u) | Decl. (`log2` u) |
|---|---:|---:|---|---:|---:|
| `collection` | `87.208048` | `87.208049` | `misc` | `69.544406` | `69.544407` |
| `chi_setup` | `65.712288` | `65.712288` | `fb` | `65.038073` | `65.038073` |
| `la` | `83.601173` | `83.601174` | `split` | `81.834890` | `81.834891` |
| `gen` | `76.777804` | `76.777804` | `gtree` | `59.185168` | `59.185169` |
| `lists` | `84.688286` | `84.688287` | `cert` | `50.000001` | `50.000001` |
| `two_logs` | `85.393742` | `85.393743` | `verify` | `11.000001` | `11.000001` |
| `canon` | `74.029127` | `74.029127` | `wrapper` | `53.712288` | `53.712288` |
| `sm` | `73.827704` | `73.827705` | **`exp_run` / `E[W]`** | **`87.854993757`** | **`86.496544115`** |

**Strict additive `FIX-U`:** Main Schwartz-Zippel `2.590e-6` + Filtering given `H3` `0.0` + Split `2^112` meter `8.306e-10` + Zero-retry Cantelli(`3.0`)+`B_max=24` boot `2*Phi_1^{24} = 1.976e-8` + Two-log `2^112` Markov `9.789e-9` + Repeat-ideal `0.0` + 64-round MR `1.735e-18` + Group-order/boot-`c_1`/SM `(8+48*T_boot_max)/r = 1.696e-57` + Numerical `0.0` `=>` **`success = 0.99999738 >= 0.999997`**.

---

## 6. Declared Exploratory Heuristics

- **H1:** Full de Bruijn `rho_dB` with odd-part Murphy `alpha`, exact 3-class 2-adic split, and all-`200,000`-cell Euler-product retention/survivor laws (Band-1/2 up to `lx2 = P_mid = 2^54.5`, zero Band 3).
- **H2:** Uniform sampling over certified Voronoi FD (`s = 0.80, kap = 0.0`, `200,000` active exact-HNF cells, `SE = 0.0042`b) and equidistribution of `b`-row starting residues mod `p` (paired with deterministic `#(occ SBs) <= min(S_p, N_a)` -> `p_occ <= min(1.0, mu_sb)`).
- **H3:** Extending the 7-scale monotone `Z[omega]` filtering ladder (`2^12..2^24`) to `U' = 2^53.64` (`exc = 0.82, W <= 161.0`): production matrix satisfies 4 hard caps (`0.00085 / 0.00055 / 13.50 / 1.00`, ZERO TABM); LA/gen/RAM priced deterministically at cap (`lazy = 0`).
- **H-SM / H-R / H-X (`RANK-BRIDGE.md`):** Non-degenerate 5-SM $\chi$-projections mod $r$; `rank M_prod = N_l`, 1D kernel of `[M_prod | S_prod]`; controllability/observability of padded `(B, Y, X)`, yielding conditional Schwartz-Zippel `4 N_sw_cap / 2^64` and deterministic shifted Beckermann–Labahn $\sigma$-basis extraction inside `gen_at_cap` + `la`.
- **H5 / H-EA / H-SEP (Declared Exploratory Model Inequalities — Supported, Not Proved):** `(0) H-EA`: `QBINS`/`BOOT_BINS` lower-bound smoothness and `QBINS_LO`/`BOOT_LO_OFF` upper-bound survivors; `(i) H5`: (a) fresh-`e` i.i.d. boot trial (`bareiss_lll_12d` charged per trial) succeeds w.p. `>= prob_b` (`eps_boot = 1e-6`), and (b) at every DAG slot reached under `k = 1`, the Conditional Cantelli lattice expectation obeys `E[1 / (1 + mu(L_q)/D(L_q)) | slot] <= 1 / (1 + mu_model / 3.0)` (`kappa_disp = 3.0`, compatibly supported by above-cutoff `B_FB=1024` deg-8 `f_rec`/deg-2 `g_rec` `Var/E <= 1.38`, depth-2 `C6` `res_v/m <= 2.11`, and `min mu_true/mu_mod = 7.4x` on `24` prod lattices), propagated via unconditional Boole (`Phi_1 = 0.463929`) and amplified across `B_max = 24` i.i.d. root attempts (`2*Phi_1^{24} = 1.976e-8`); `(ii) H-SEP`: 8-bit ECM ladder achieves per-curve conditional separation `>= 0.50` (`min_dag = 0.6174`, `min_boot = 0.6341`, `63/63` Suyama-12 live; `MONT5 = 367.76`).
