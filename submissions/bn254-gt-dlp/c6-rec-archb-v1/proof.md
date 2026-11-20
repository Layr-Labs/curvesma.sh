# bn254-gt-dlp: Hardened Pure-Scalar Content-Reduced Reciprocal CM-Sextic `(6,2)` STNFS (Full de Bruijn + Direct-`2^21` `Z[omega]` LA) — Expected Work `2^87.372777743` Units

Track `bn254-gt-dlp`. Input: `z = g^x` in `GT` (`g = e(P0, Q0)`), `x` uniform in `0..r-1`; output `x`. The score unit is `C = 20,000` scalar AArch64 instruction equivalents (100% pure baseline scalar AArch64 integer instructions; zero vector/SIMD extensions, zero `FEAT_LSE` atomics). No advice, smaller targets, or amortization across targets are used. Analytic, exploratory claim.

## Claim

Let `W` be the total classical work of one attack on one full-production BN254 target. Under declared heuristics H1-H5, H-EA, H-FFT, H-R, H-X, and H-SEP (Section 7 — explicitly including H3 for production filtering hard-cap conformance and H5 for standard STNFS descent branching and orbit-collision independence), with **full de Bruijn local-density correction `rho_dB` baked into the primary claim**, **deterministic hard-cap linear-algebra pricing at `(N_sw_cap, N_l_cap)`** (zero Jensen/variance premise for expected quadratic LA cost!), and **100% pure scalar AArch64 kernels**, the attack succeeds with probability at least `0.39` over uniform challenges and fresh coins (an opened run succeeds with probability at least `0.969522` by the strict additive union bound `FIX-U`; the run gate opens at most a share `0.402260` of runs), and

    E[W] <= 1 + (K / 2^64) * 2^88.686577254 <= 2^87.372777743 units (C = 20,000 instructions),
    K = 7420388798548898458.

- **Run gate:** Read one fresh 64-bit coin `u` and open the run iff `u < K = 7420388798548898458`, the least integer with `(K / 2^64) * 0.969522 >= 0.39` (`evidence.py` checks in exact rationals that `K` passes and `K - 1` fails). A closed run fails at cost `1` unit.
- **Metered phases:** Both individual logarithms share one combined `2^112`-unit meter; the deterministic split phase has its own `2^112`-unit interruptible budget. Every declared upper bound rounds strictly upward and `success` rounds strictly downward; `evidence.py` recomputes every row, algebraic check, exhaustive fundamental-domain check, retention/survivor law, and the complete 62-row descent recurrence live.

---

## 1. Tower, Content-Reduced Reciprocal Pair, Order-6 Frobenius, & Live Mod-`r` CRT

BN254 parameters: `u = 4965661367192848881`, `p = 36u^4 + 36u^3 + 24u^2 + 6u + 1`, `r = 36u^4 + 36u^3 + 18u^2 + 6u + 1`, `ord_r(p) = 12`.

- **Base field `K_h`:** Rank-2 CM sextic over `K_2 = Q(sqrt(-2))`:
  `h(y) = y^6 - 2y^5 + y^4 + 2y^3 - 4y + 3`, `D_h = -184,832 = -2^9 * 19^2`, signature `(0, 3)`, unit rank `2`, class number `h_K = 1`, `w_K = 2`, `Reg = 2.374167`, `[O_h : Z[y]] = 3`, `1/zeta_K(2) = 0.683542`.
- **Order-3 automorphism & twist:** `sigma(y) = (-2/3)y^5 + y^4 + (1/3)y^3 - (5/3)y^2 - (4/3)y + 3` fixes `K_2` and satisfies `sigma(y) = y^(p^4) = y^(p^10) mod (p, h)`. The `sigma`-invariant twist `w = sqrt(-2)*(1 + sqrt(-2)) = -y^5 + y^4 - y^2 - 2y` has norm `N(w) = 216`, constant modulus `|w_j| = sqrt(6)` across all three complex places (giving exact `sigma`-isometry of `Q_a` and `Q_b`), and `(u*w)^((p^6 - 1)/2) = -1` in `F_{p^6}`.
- **Content-reduced reciprocal pair:** Dividing `w^4 * x^8 * F_BN(1 / (w * x^2))` by `c_0 = 4(1 + sqrt(-2))^2` (`N(c_0) = 2^12 * 3^6`, `21.509775` bits) yields the primitive balanced polynomial:
  `f_rec(x) = (-2y^5 + 2y^4 - 2y^2 - 4y + 3)x^8 + (3y^5 - 3y^4 + 3y^2 + 6y)x^6 - 12x^4 + (-3y^5 + 3y^4 - 3y^2 - 6y + 12)x^2 + (2y^5 - 2y^4 + 2y^2 + 4y - 5)`,
  `g_rec(x) = u * w * x^2 - 1`.
  Both leading and constant coefficients of `f_rec` have norm `729 = 3^6`, `f_rec(x^2 -> (u*w)^-1) = 0 mod (p, h)`, and `N(disc(f_rec)) = 2^96 * 3^96`.
- **Order-6 Frobenius:** `tau(y, x) = (sigma(y), -x)` acts on `F_{p^12}` as `Frob^10` (`(u*w)^((p^10 - 1)/2) = -1`), with eigenvalue `chi = p^10 mod r = 36u^3 + 18u^2 + 6u + 2` (`192` bits, `chi^2 - chi + 1 = 0 mod r`, `chi^3 = -1 mod r`). Every ramified prime of `O_{K_g}` has `tau`-orbit length `< 6` (zero virtual-log corruption; PARI/GP certified `0/1000` valuation mismatches at `{2, 3, 19}`).
- **Live Mod-`r` CRT Certification (`evidence.py [algebra]`):** Distinct-degree factorization modulo `r` is executed live over **all 6 roots** of `h(y) mod r`: `h(y) mod r` has degree pattern `[1, 1, 1, 1, 1, 1]` (squarefree), and at **every single one of the 6 roots**, `g_rec(x) mod r` is irreducible of degree `[2]` and `f_rec(x) mod r` splits as `[2, 2, 2, 2]`. Thus `eps = r^2 - 1` (`508` bits) annihilates every residue-field unit group modulo `r` identically across all 6 components.
- **Exhaustive Small-Radius Fundamental Domain Bijection (`evidence.py [exhaustive_fd]`):** On `T_2(v) <= 50.0`, Fincke-Pohst enumeration finds all `1552` nonzero lattice points in `O_h`, grouping into `107` exact `<+-1, O_h^x, sigma>` orbits. With the 6-facet Voronoi cone + exact integer unit-matrix boundary tie-break (`UMATS[0..5]`), **`107/107`** orbits have exactly 1 canonical representative (`missed_orbits = 0`, `duplicate_orbits = 0`).
- **Small primes & Murphy alpha:** `2 = P_2^2` with `N(P_2) = 8` gives exact 2-adic classes `(v_2(N_f), v_2(N_g))` of `(0,0)` (`7/9`), `(0,3)` (`b in P_2`, `1/9`), and `(6,0)` (`a in P_2`, `1/9`). Exact 100k Murphy alpha plus `+0.08`b uniform safety margin: `alpha(f_rec) = +1.580`b, `alpha(g_rec) = -0.039`b.

---

## 2. Pure-Scalar AArch64 Collection (`Arch B-Hybrid`, Full de Bruijn `rho_dB`)

- **Asymmetric bounds & full de Bruijn yield:** `B_f = 2^60` (batch-tested), `B_g = 2^59` (sieved), `U' = 2^52.666123` order-6 `tau`-orbits, relation target `0.95 * U'` with `+0.0666`b margin. Evaluated under **full de Bruijn local density `rho_dB(u, log2 N) = rho(u) - gamma * rho(u-1) / (log2(N) * ln 2)`** on the bundled `200,000`-sample Voronoi FD pool quadrature (`c6_data.json`, `n_in = 89151`, `n_prim_in = 61463`, `frac_reg = 0.445755` matching the exact 6D ellipsoid integral `0.445421` within `0.001`):
  - Active cells `V = 2^97.7799` (`SE = 0.0077`b), `V_b = 2^54.6717`, canonical `a`-points `N_a = 2^44.2738`, `b`-rows `N_5 = 2^45.6491` (exact 5D Schur complement of `G_b`, average chord `L_b = 520.1`, exact max chord `L_max = 2 * R_b / sqrt(G_b[0][0]) = 1059.52`), relation-weighted norms `N_f = 405.49`b, `N_g = 468.69`b.
- **Clean `(t, a_slot)` active-only score arena & pure-scalar 64-bit SWAR fused scan+zero (`2.50` scalar instr/cell):**
  - Canonical `a`-points are sorted once in increasing order of `Q_a(a)`. Along any `b`-row `b_0 + t * 1` (`0 <= t < L_b`), the active-region predicate `Q_b(b(t)) <= 1 - kap * Q_a(a)` is equivalent to **`a_slot < M(t)`**, a single contiguous prefix `[0, M(t))`! Storing the single-owner score arena in row-major `(t, a_slot)` order (`score[t][0..M(t)-1]`) guarantees that the fused scan+zero touches **strictly the `V_kept = V` active cells** (zero dead cells scanned or zeroed) and fits on ONE single owner (`(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table) = 2^55.33 B << 2^64 B`).
  - Because `|u*w*a^2| / |b^2| > 2^62`, `|log2 N_g,rec(a,b) - log2 N_g,rec(a,0)| < 10^-15` bits (`evidence.py`: `0.0e+00` bits in binary64 across 500 production-scale pairs). Precomputing `bias_table[a_slot] = 0x8000 - T(a_slot)` once for all `N_a` canonical `a`'s costs `< 2^40` units (`width = 0.25`b, `thr = 0.00`).
  - For a fixed `t`-slice, `b(t)` is constant across the inner `a_slot` loop; when `b(t) in P_2` (`v_2 = 3`, probability `1/9`), we add the constant `x7 = 0x0030003000300030` (`+48` quantized bits per lane) to the 4 packed 16-bit scores (`s_i < 2^15`, `bias_i < 2^15`, zero inter-lane carry). Pure baseline scalar AArch64 loop (9 integer instructions per 4 active cells = `2.25` + `0.25` slice setup = **`2.50` scalar instr/cell**, verified `0/10,000` errors with independent varying `T_0..T_3` and `b in P_2` in `evidence.py`):
    `ldr x2,[x0]; ldr x4,[x5],#8; str xzr,[x0],#8; add x2,x2,x7; add x3,x2,x4; tst x3,x6; b.ne .hit; subs x1,x1,#4; b.gt .loop`.
- **End-to-end pure-scalar `Arch B-Hybrid` accounting (`slack = 114`, `slack_b = 119`, live Euler-product `r_g = 0.8587`, `sf = 1.089e-06`, `s0 = 7.149e-08`):**
  - **Band 1 `[2^14, 2^43)` (pointwise `b`-lookup):** Bucketed `a`-table (`0.5` tables/prime, `256` `Q_a`-shells). For each `b`, exact `Q_b(b)` selects the exact shell prefix `[0, s_max(b))` (zero `Q_b(t)` candidate waste; `1.004` partial-shell overhang). Exact 6-instruction pure scalar AArch64 loop (`ldr; ldrh; add; strh; cmp; b.lo`, `C_HIT_LO = 6.0`) + `V_b * (li(2^43) - li(2^14))` lookups at `C_LOOK_LO = 24.0` `=> 7.60` instr/cell.
  - **Band 2 `[2^43, 2^52.5)` (`Q_a`-sorted super-bucket `a`-table, `W_sb = 128`):** Every `b`-row chord (`L_max = 1059.52 < 64 * 128 = 8192`) spans `avg_sb = 5.06` super-buckets and at most 2 adjacent 64-bit bitmap words (exact expectation `bmw = 1.0634` words/probe at `C_PROBE = 18.0` `=> 0.70` instr/cell). Exact live Poisson integration of `1 - exp(-N_a * W_sb / p)` over `[2^43, 2^52.5)` charges `C_OCC_SB = 13.0` scalar instr per occupied super-bucket (`1.55` instr/cell). Exact 6D ellipsoid integration over `r_2 ~ (5/2) r_2^(3/2)` and chord `Q_b(s)` variation gives `E[cands]/E[hits] = 1.8037` at `C_CAND = 13.0` (checking `qa <= qa_max`, `t < L_b`, `a_slot < M[t]`) + `C_HIT_HI = 9.0`/active hit `=> 6.48` instr/cell. Verified `0` mismatches over `1523` hits in `evidence.py`.
  - **Band 3 `[2^52.5, 2^57.5)`:** Pre-walked 7D walk (`E[log2(run_mass)] = 5.52`, `frac(L>=16) = 0.989`, `rec_top = 10.16` instr/cell) + explicit per-`b`-row bucket write (`12.0`), owner-to-owner transfer (`8.0`), and pure scalar apply (`5.0`) = `25.0` instr/hit (`bkt_top = 2.27` instr/cell).
  - **Cofactorization (exact bn15r kernels, zero downward proportional scaling):** `C_NORMS = 27,000`, `C_SB = 80,000`, `gbatch_cost(60)`, `C_TDIV`, `C_BATCH`, `gbatch_cost(59)` `=> cof = 0.96` instr/cell.
  - **Totals:** `32.210` scalar instr/cell (times `1.0033` work allowance) `=> collection <= 2^88.506316` units; 3-pass pure-C `a`-table build (`C_LIST_A = 35.0`) `=> lists <= 2^81.583325` units; peak RAM `2^97.65 B < 2^100 B`.

---

## 3. Deterministic Hard-Cap `Z[omega]` Filtering, Linear Algebra, & Heuristic H3

- **Real C `filtsim_bn254 --ring z6` monotonic scale ladder (`U' = 2^20 -> 2^21`, norms `405/469`, `exc = 0.95`, `W = 200`):**
  - `U' = 2^20` (4 seeds): `N_sw/U' = 0.00211`, `N_l/U' = 0.00118`.
  - `U' = 2^21` (8 independent seeds): `N_sw/U' = 0.00187 +- 0.00003` (`-10.9%` vs `2^20`), `N_l/U' = 0.00107 +- 0.00002` (`-9.8%` vs `2^20`), `unit_share = 97.28%`, explicit-entry unit share `exp_unit = 90.41%`.
- **Deterministic Hard-Cap Abort & Quadratic LA Pricing (Zero Jensen/Variance Premise):**
  - After filtering, the algorithm **deterministically aborts** if either `N_sw > N_sw_cap = 1.75 * m_sw_ub * U' = 2^44.4630` (`nu_sw_cap = 0.003393`, `+81%` above the `2^21` mean) OR `N_l > N_l_cap = 1.75 * m_l_ub * U' = 2^43.6422` (`nu_l_cap = 0.001921`, `+80%` above the `2^21` mean), where `m_sw_ub = 0.001939` and `m_l_ub = 0.001098` are the direct `2^21` mean `+ 2 SE` bounds (no downward `2^21 -> 2^25` mean scaling).
  - With `lazy = 0`, every term of the Block Wiedemann cost (`N_l * N_sw * c_row * (...) + N_l^2 * (...) + N_l * (...) + const`), the generator cost (`gen`), and the 1-owner RAM is **strictly monotonically increasing in both `N_sw` and `N_l`**. Therefore, pricing `la`, `gen`, and 1-owner RAM **at the hard cap `(N_sw_cap, N_l_cap)`** guarantees deterministically on every non-aborted run that `W_LA(N_sw, N_l) <= W_LA(N_sw_cap, N_l_cap)` and `W_gen(N_l) <= W_gen(N_l_cap)` — requiring **no** Jensen inequality or variance premise for expected quadratic LA cost.
  - **Heuristic H3 (Production Filtering Cap Conformance):** At production scale `U' = 2^52.67`, we declare heuristic H3 that `P(N_sw > N_sw_cap or N_l > N_l_cap) <= p_cap_h3 = 0.027819` (calibrated by the 8-seed `2^21` reference-scale `chi^2(7)` Cantelli union bound `1.5974e-02 + 1.1846e-02` and empirically supported by the monotonic `2^20 -> 2^21` shrinkage ladder and `+81%` cap margin; we explicitly declare this production transfer as H3 rather than assuming zero scale drift or normal sampling at `2^52.67`).
- **1-Owner 3-Bucket `Z[omega]` BW (Evaluated Deterministically at Hard Cap `(N_sw_cap, N_l_cap)`):** Conservative explicit-entry `k5 = exp_unit * 7.8 + (1 - exp_unit) * 33.0 = 10.2157` (no implicit-unit discount), conservative `c_tail = 420.0` instr/row, `5` Schirokauer columns, `G4`-tuned dense-band tables (`tab_r = 0.7267`), BW `(m, n) = (768, 35)`: `la <= 2^84.868599`, `gen <= 2^78.661451`. Explicit 1-owner LA RAM at hard cap `N_sw_cap` includes BOTH input block `x` (`160 B/slot`) AND output block `y` (`160 B/slot`) plus matrix rows, dense tables, and 3-bucket buffers: `2^59.70 B < 2^64 B`.
- **Peeling cascade:** Seeded fixed point at `(60, 59, exc=0.95)` gives `e_crit = 0.718` (`> 0.41`b margin) and `h = 0.3471`.

---

## 4. Complete Live 62-Row 12D Descent Plan (`LB = 59`, `ORBIT = 6`, Heuristic H5, `evidence.py [descent]`)

Because `C6-m184_w216-rec` is a `(6, 2)` tower, both boot and special-`q` lattices are 12D. Direct PARI/GP 12D LLL audits (`n=24` boot targets, `n=100` `q~2^100` special-`q` combinations) certify `C6-m184_w216-rec` beats `bn15r` by `-1.71`b on 12D boot and `-8.33`b on special-`q` nodes. Conservatively claiming **zero** of that norm advantage while charging `C6-m184_w216-rec`'s exact `base_primitive = 0.683542` (`counted_frac = 0.644655`), `descent_bn254_c6.py` executes the full 62-row bottom-up recurrence + boot `(163, 90, 16, 4)` live under declared heuristic **H5** (standard STNFS descent branching & orbit-collision model, identical to accepted `bn15r`/`c36`: (i) Poisson/PGF relation-count envelope at each bin-side node, (ii) conditional independence of successive candidate relations at a node giving `psi^k`, and (iii) uniform orbit-collision probability `<= 1/(li(2^59)/6)` per pair of encountered descent ideals):
- `two_logs = 2^83.853597 <= 2^83.853598` units (`phi_max = 0.03268 <= 0.04`, `psi_max = 0.42046 < 0.8`, `boot_psi = 0.46843`, `E[N] = 2^11.8339`, node budget `2^20.977`).
- **Strict `RT` repeated-ideal bound (under H5(iii)):** Uses single-side `tau`-orbit count at the lower bound `2^59`, `N_orbits,min = li(2^59)/6 = 2^51.0977` `=> P(abort) <= 2^-9.1431`, `P(repeat) <= 2^-10.1437`.

<details>
<summary>Complete 62-Row Production Descent Plan Table (`LB = 59`, executed live in `descent_bn254_c6.py`)</summary>

| Bin | Side | `n` | `iA` | `kA` | `kB` | `R` | `thin` |
|---|---|---:|---:|---:|---:|---:|---:|
| `[59,65)` | `f` | `59` | `0` | `0` | `0` | `56` | `1` |
| `[59,65)` | `g` | `59` | `0` | `0` | `0` | `56` | `1` |
| `[65,70)` | `f` | `59` | `0` | `0` | `0` | `64` | `1` |
| `[65,70)` | `g` | `59` | `0` | `0` | `0` | `64` | `1` |
| `[70,75)` | `f` | `59` | `0` | `0` | `0` | `80` | `1` |
| `[70,75)` | `g` | `59` | `0` | `0` | `0` | `80` | `1` |
| `[75,80)` | `f` | `71` | `69` | `2` | `0` | `32` | `1` |
| `[75,80)` | `g` | `67` | `17` | `1` | `1` | `40` | `1` |
| `[80,85)` | `f` | `74` | `35` | `1` | `1` | `32` | `1` |
| `[80,85)` | `g` | `69` | `44` | `2` | `0` | `48` | `0` |
| `[85,90)` | `f` | `72` | `57` | `2` | `0` | `48` | `0` |
| `[85,90)` | `g` | `69` | `49` | `1` | `1` | `48` | `0` |
| `[90,95)` | `f` | `71` | `31` | `1` | `1` | `48` | `0` |
| `[90,95)` | `g` | `69` | `31` | `1` | `1` | `56` | `0` |
| `[95,100)` | `f` | `73` | `26` | `1` | `1` | `56` | `1` |
| `[95,100)` | `g` | `78` | `45` | `2` | `0` | `64` | `0` |
| `[100,105)` | `f` | `73` | `22` | `1` | `1` | `64` | `0` |
| `[100,105)` | `g` | `70` | `30` | `1` | `1` | `80` | `0` |
| `[105,110)` | `f` | `71` | `27` | `1` | `1` | `80` | `0` |
| `[105,110)` | `g` | `76` | `22` | `1` | `1` | `80` | `0` |
| `[110,115)` | `f` | `80` | `48` | `2` | `0` | `96` | `0` |
| `[110,115)` | `g` | `82` | `21` | `1` | `1` | `80` | `0` |
| `[115,120)` | `f` | `84` | `24` | `1` | `1` | `80` | `0` |
| `[115,120)` | `g` | `98` | `21` | `1` | `1` | `64` | `0` |
| `[120,125)` | `f` | `98` | `22` | `1` | `1` | `64` | `0` |
| `[120,125)` | `g` | `99` | `25` | `1` | `1` | `80` | `1` |
| `[125,130)` | `f` | `99` | `25` | `1` | `1` | `80` | `0` |
| `[125,130)` | `g` | `97` | `23` | `1` | `1` | `96` | `0` |
| `[130,135)` | `f` | `108` | `20` | `1` | `1` | `80` | `0` |
| `[130,135)` | `g` | `102` | `31` | `1` | `1` | `96` | `0` |
| `[135,140)` | `f` | `104` | `33` | `1` | `1` | `80` | `0` |
| `[135,140)` | `g` | `96` | `27` | `2` | `1` | `96` | `0` |
| `[140,145)` | `f` | `104` | `43` | `1` | `1` | `96` | `0` |
| `[140,145)` | `g` | `103` | `27` | `2` | `1` | `96` | `0` |
| `[145,150)` | `f` | `111` | `24` | `1` | `1` | `112` | `0` |
| `[145,150)` | `g` | `99` | `46` | `2` | `1` | `96` | `0` |
| `[150,155)` | `f` | `145` | `56` | `2` | `0` | `80` | `0` |
| `[150,155)` | `g` | `104` | `47` | `2` | `1` | `96` | `0` |
| `[155,160)` | `f` | `124` | `56` | `1` | `1` | `80` | `0` |
| `[155,160)` | `g` | `108` | `66` | `2` | `1` | `80` | `0` |
| `[160,165)` | `f` | `129` | `62` | `1` | `1` | `80` | `0` |
| `[160,165)` | `g` | `109` | `64` | `2` | `1` | `96` | `0` |
| `[165,170)` | `f` | `129` | `62` | `1` | `1` | `96` | `0` |
| `[165,170)` | `g` | `134` | `24` | `2` | `1` | `112` | `0` |
| `[170,175)` | `f` | `119` | `72` | `5` | `1` | `14` | `0` |
| `[170,175)` | `g` | `124` | `64` | `2` | `1` | `80` | `0` |
| `[175,180)` | `f` | `118` | `72` | `5` | `1` | `16` | `0` |
| `[175,180)` | `g` | `127` | `64` | `2` | `1` | `96` | `0` |
| `[180,185)` | `f` | `119` | `72` | `5` | `1` | `16` | `0` |
| `[180,185)` | `g` | `142` | `48` | `2` | `1` | `96` | `0` |
| `[185,190)` | `f` | `126` | `72` | `5` | `1` | `14` | `0` |
| `[185,190)` | `g` | `144` | `42` | `2` | `1` | `128` | `0` |
| `[190,195)` | `f` | `125` | `72` | `5` | `1` | `16` | `0` |
| `[190,195)` | `g` | `153` | `48` | `2` | `1` | `96` | `0` |
| `[195,200)` | `f` | `126` | `72` | `6` | `1` | `16` | `0` |
| `[195,200)` | `g` | `160` | `48` | `2` | `1` | `96` | `0` |
| `[200,205)` | `f` | `125` | `72` | `6` | `1` | `16` | `0` |
| `[200,205)` | `g` | `163` | `48` | `2` | `1` | `112` | `0` |
| `[205,210)` | `f` | `132` | `72` | `5` | `1` | `16` | `0` |
| `[205,210)` | `g` | `170` | `48` | `2` | `1` | `112` | `0` |
| `[210,215)` | `f` | `139` | `72` | `5` | `1` | `14` | `0` |
| `[210,215)` | `g` | `151` | `72` | `6` | `1` | `14` | `0` |

</details>

---

## 5. Itemized Ledger, Strict Additive Union Bound (`FIX-U`), & Run Gate

| Ledger Row | Computed (`log2` u) | Declared Upper Bound (`log2` u) |
|---|---:|---:|
| `collection` | `88.506316` | `88.506316` |
| `chi_setup` | `65.712288` | `65.712288` |
| `la` | `84.868599` | `84.868599` |
| `gen` | `78.661450` | `78.661451` |
| `lists` | `81.583325` | `81.583325` |
| `two_logs` | `83.853597` | `83.853598` |
| `canon` | `72.786130` | `72.786130` |
| `sm` | `69.749687` | `69.749688` |
| `misc` | `68.788318` | `68.788318` |
| `fb` | `64.061920` | `64.061920` |
| `bases` | `76.870177` | `76.870177` |
| `split` | `80.275452` | `80.275452` |
| `gtree` | `60.174753` | `60.174754` |
| `cert` | `50.000001` | `50.000001` |
| `verify` | `11.000001` | `11.000001` |
| `wrapper` | `53.712288` | `53.712288` |
| **Expected opened run (`exp_run`)** | **`88.686577253`** | **`88.686577254`** |
| **Gated claim (`E[W]`)** | **`87.372777742`** | **`87.372777743`** |

**Strict additive union bound on opened-run success (`FIX-U`, zero survival-factor multiplication; under H1-H5):**
- Main-phase Schwartz-Zippel (`4 N_sw_cap / 2^64`): `5.258e-06`
- Filtering hard-cap abort (`H3`, `p_cap_h3` for `N_sw > N_sw_cap or N_l > N_l_cap`): `2.782e-02`
- Split `2^112` budget abort: `2.818e-10`
- Log-phase `2^112` meter Markov bound: `3.366e-09`
- Descent node-budget abort (`H5`, `2^-9.1431`): `1.769e-03`
- Descent repeated-ideal collision (`H5(iii)` + `RT`, `2^-10.1437`): `8.840e-04`
- Group-order & numerical (`2/r + 2^-48`): `3.553e-15`
- **Opened-run success:** `1 - sum(failures) = 0.96952272 >= 0.969522`.

---

## 6. Sensitivity Comparison (`evidence.py [sensitivity]`)

| Model Variant (All 100% Pure Scalar AArch64, Hard-Cap LA, Live 62-Row Descent) | Opened (`log2` u) | Gated (`log2` u) | Delta vs #1 (`87.5884`) |
|---|---:|---:|---:|
| **PRIMARY CLAIM (Full de Bruijn `rho_dB` + Direct `2^21` Hard-Cap LA + `c_tail=420`)** | **`88.686577`** | **`87.372778`** | **`-0.2156` b** |
| Full de Bruijn `rho_dB` + Matched-scale `2^21->2^25` Hard-Cap LA + `c_tail=420` | `88.653302` | `87.339501` | `-0.2489` b |
| Plain Dickman `rho` + Direct `2^21` Hard-Cap LA + `c_tail=420` | `88.521620` | `87.207821` | `-0.3806` b |
| Plain Dickman `rho` + Matched-scale Hard-Cap LA + `c_tail=207` | `88.475513` | `87.161712` | `-0.4267` b |

---

## 7. Declared Heuristics

- **H1 / H1':** Standard NFS smoothness under full de Bruijn local density `rho_dB` with 100k Murphy `alpha`, exact 2-adic split, and live Euler-product retention/survivor laws on `N_f`, `N_g`.
- **H2:** Uniform sampling over the certified 6-facet Voronoi fundamental domain (200,000-sample bundled pool quadrature, `SE = 0.0077`b).
- **H3 (Filtering Percolation & Production Cap Heuristic):** At production scale `U' = 2^52.67`, the filtered hexagonal `Z[omega]` matrix dimensions satisfy `N_sw <= N_sw_cap` and `N_l <= N_l_cap` (`CAP_MULT = 1.75` above direct `2^21` `mean + 2 SE`, `+81%` above the `2^21` mean) with failure probability `<= p_cap_h3 = 0.027819` (supported empirically by the monotonic real-C `filtsim_bn254 --ring z6` `2^20 -> 2^21` shrinkage ladder and 8-seed `chi^2(7)` reference-scale check; LA, generator, and 1-owner RAM are priced **deterministically at the hard cap `(N_sw_cap, N_l_cap)`** with `lazy = 0`, requiring no expectation/variance premise for quadratic LA cost).
- **H4 / H-R / H-X:** Full rank of `[M | S]` modulo the virtual-log kernel line and standard Block Wiedemann rank/Krylov regularity on the `Z[omega]` product-form Schur complement.
- **H5 / H-EA / H-FFT / H-SEP (Standard STNFS Descent Branching & Orbit-Collision Heuristic, identical to accepted `bn15r`/`c36`):** At `B_g = 2^59`, (i) usable special-`q` descent relations at each bin-side node follow the standard Poisson/PGF relation-count envelope (`descent_bn254_c6.py`), (ii) successive candidate relations at a node fail conditionally independently (giving the `psi^k` failure bound), and (iii) distinct encountered descent ideals collide in the same `tau`-orbit with probability `<= 1 / (li(2^59)/6)` per pair (`RT`), together with scaled remainder tree FFT pricing (`H-FFT`) and stored-prefix ECM isolation (`H-EA`, `H-SEP`).
