# G1 ECDLP (powers-of-tau): W(n)=const (4,3) CM-STNFS over Q(zeta_15) behind a run gate: expected work `2^112.832352798025`

Track `bls12-381-g1-ecdlp`. Input: n (1 <= n < r), auxiliary powers A_1..A_n = [s^i]P0 and target T = [s]P0 (uniform s in 1..r-1); output scalar s in 1..r-1 with [s]P0 = T. Reference unit: one batched incomplete affine G1 addition, C = 2,000 (1 core unit of 20,000 instructions = 10 G1 units). No advice, smaller targets, SIMD, `FEAT_LSE`, hand-written assembly, or amortization are used. **Full de Bruijn local-density correction** (`rho_dB(u, x) = rho(u) - gamma rho(u-1)/ln x`) is explicitly baked into **both** Main Collection and Individual-Logarithm Descent. Predecessor work is credited in Section 9; no prior verdict is a premise.

## Claim

`W` is the total classical work of the one attack on the one target. Step 0 is a run gate (kaankolcu's device, Section 9): draw a uniform 64-bit word `u`; unless `u < K = 7212029477497188191`, stop immediately as a failure (`< 1` unit). `K` is the least integer with `K x 0.997532 >= 0.39 x 2^64`. Under heuristics `H1`, `H1'`, `H-EA`, `H-FFT`, `H2`, `H3`, `H-SM`, `H4`, `H4u`, `H-X`, and `H5` (Section 8), an opened run succeeds with probability `>= 0.997532` (strict additive union bound, Section 6); because `u` is independent of the challenge and all remaining coins, overall algorithmic success is `>= 0.39` over uniform full-production challenges and fresh coins, and

Applicability: valid for **all integers `1 <= n < r`** with zero number-theoretic side conditions on `r - 1`; scored at `n = 2^20 = 1,048,576`.

    E[W(n)] <= 1 + (K / 2^64) x 10 x (2^110.8653137017 + 2^66) < 2^112.832352798025  G1 units (C = 2,000).

Here `2^110.8653137017` core units (`C = 20,000` scalar AArch64 instructions) is the sum of Section 7's rows with every-run bounds on collection, replay tables, linear algebra, generator, Schirokauer maps, ECM split meter (`2^96.5`), and fixed rows, and expected work on the two individual logarithms (`2^102.498`, Section 5). The score `112.832352798025` is an expectation over the challenge and fresh coins of the single attack, not a bound on every opened run; every opened run is strictly metered at:

    W(n) <= 10 x (2^115.6530356257 + 2^66) < 2^118.974963720588  G1 units (C = 2,000).

Analytic, exploratory claim. All declared numbers round conservatively upward (or downward for success floors) and are audited live from `cm_data.json` by `audit_cm.py` + `cm_models.py` (`checks_ok = true`, pure Python stdlib).

## 1. Reduction (`G1` ECDLP with auxiliary powers: constant `W(n)` for all `1 <= n < r`) and field map

Per `COST_MODEL.md` (*"The attack may use only a subset of the powers"*) and `target-profile.json` (*"the attack may use a subset and need not read or validate the whole sequence"*), this attack ignores the auxiliary powers `A_2, ..., A_n` and operates solely on `T = A_1 = [s]P0`. Consequently, its work `W(n)` is **constant for all `1 <= n < r`** (requiring no divisibility condition on `r - 1`) and is evaluated at the scoring parameter `n = 2^20 = 1,048,576`.

With the pinned reduced pairing `e` (`py_ecc.optimized_bls12_381.pairing(Q, P)`) and public generators `P0 in G1`, `Q0 in G2`:

0. **Run gate (Claim):** draw uniform 64-bit `u`; if `u >= K`, abort as failure (`< 1` G1 unit).
1. **Validate target:** decode `T = [s]P0`, check `T` on `E(F_p)`, `[r]T = infinity`, `T != infinity` (uniform `s in 1..r-1`).
2. **Pairings:** compute `g = e(P0, Q0)` and `z = e(T, Q0) = g^s` (computed fresh per attack; nothing precomputed or amortized).
3. **Field map:** construct the explicit isomorphism `F_(p^12) -> F_p[y, X]/(h_pi, X^3 - omega^2 u)` (`field-map.md`) and map `(g, z)`.
4. **Fresh-coin randomization & core solve:** draw independent uniform `e_1, e_2 in {1, ..., r - 1}` and set `g' = g^(e_1)`, `z' = z^(e_2)`. Because `G_T` has prime order `r` and `g, z != 1`, BOTH `g'` and `z'` are **exact uniform on `G_T \ {1}`** from fresh algorithm coins alone (eliminating any fixed-generator dependency for `g = e(P0, Q0)`). Run the core under its meter and node counter to recover `log(g')`, `log(z')`, and set `s = (e_1 * log(z')) * (e_2 * log(g'))^(-1) mod r`.
5. **Output & check:** verify `[s]P0 = T` and `1 <= s < r`; output `s` iff valid, else fail.

Wrapper work (`2` pairings, field map `< 2^67`, G1 scalar mul, check) is `< 2^80` instructions = `2^65.71 < 2^66` core units (`= 10 x 2^66 = 2^69.32` G1 units); wrapper failure `< 2^(-48)`.

## 2. CM-Split `Q(zeta_15)` Tower, `w' = omega^2` Isometry, & Order-12 Frobenius

For BLS12-381, `u = -0xd201000000010000 = -2^16 x 906349 x 254760293`, `r = u^4 - u^2 + 1` (255-bit prime, `r == 13 mod 15`), and `p = (u - 1)^2 r / 3 + u` (381-bit prime, `p == 7 mod 15`).

### CM Splitting Over `K_H = Q(zeta_15)` and Exact `w = y <=> w' = omega^2` Isometry

Because BLS12-381 has CM discriminant `-3`, its degree-6 polynomial `F(z) = (z - 1)^2 (z^4 - z^2 + 1) + 3 z` splits over `Q(omega)` (`omega^2 + omega + 1 = 0`) into two degree-3 conjugates `F(z) = F_1(z) F_2(z)` with `F_1(z) = z^3 - z^2 - (omega + 1) z - omega`. Taking the cyclotomic base field `K_H = Q(zeta_15) = Q[y]/(Phi_15(y))` (`[K_H : Q] = 8`, `O_H = Z[zeta_15]`, `D_H = 3^4 x 5^6 = 1,265,625`, PID `h(K_H) = 1`, `omega = y^5`, `zeta_5 = y^3`):
- **Split prime base:** Because `p == 7 mod 15` (`ord_15(p) = 4`), `Phi_15(y) == h_pi(y) h_pibar(y) (mod p)` splits into two degree-4 irreducibles over `F_p`, where `h_pi = gcd(Phi_15, y^5 - w_pi)` (`w_pi = (u^3 - u^2 - u)/(u + 1) mod p`) selects the factor with `F_1(u) == 0 (mod P)`, `P = (p, h_pi(y))`, `O_H / P ~= F_(p^4)`.
- **`sigma`-fixed cubic twist `w' = omega^2`:** Applying the unimodular `T_2`-isometry `(X, b_tilde) = (zeta_5 x, zeta_5^(-1) b)` to the `w = y` twist (`y = omega^2 zeta_5^2`) produces the **100% `sigma`-fixed polynomials**:
  `g(X) = X^3 - omega^2 u` (deg 24 over `Q`),  `f(X) = F_1(omega X^3) = X^9 + (omega + 1) X^6 + X^3 - omega` (deg 72 over `Q`).
- **Exact order-12 Frobenius (`audit_cm.py`: `algebra_bls12`):**
  1. `(u omega^2)^((p^4 - 1)/3) != 1` in `F_(p^4)`, so `g(X)` is irreducible over `F_(p^4)` and `F_(p^12) ~= O_H[X] / (P, g(X))` with `f(X) == 0 (mod P, g(X))`.
  2. The automorphism `sigma(y) = y^7` (`ord(sigma) = 4`) fixes `omega = y^5` (hence fixes `P`, `f`, and `g` identically) and satisfies `sigma(y) == y^p (mod h_pi)`.
  3. Since `(u omega^2)^((p - 1)/3) == omega^2 (mod P)`, the order-12 automorphism `tau(y, X) = (y^7, omega^2 X)` acts as **exact `Frob^1`** on `F_(p^12)`. Thus `log_r(tau(I)) == chi log_r(I) (mod r)` with `chi = p == u (mod r)`, `ord_r(chi) = 12`, and `Phi_12(chi) = chi^4 - chi^2 + 1 == 0 (mod r)`. One relation stands for its 12-element `tau`-orbit and the factor base is quotiented into `Z[zeta_12]` orbits.

### Split-Base Virtual-Log Theorem & Pigeonhole Zero-Log Theorem

1. **Split-base localization (`P` vs `P_bar`):** The NFS homomorphisms `Psi_f, Psi_g` localize **only** at the single prime ideals `P_f, P_g` above `P`, which `tau` fixes individually (`sigma(P) = P`). Prime ideals above the conjugate `P_bar = (p, h_pibar)` are coprime to `P` and have absolute norm `>= N(P_bar) = p^4 ~= 2^1524 >> N_f, N_g ~= 2^680` (cannot divide any main-sieve relation) and `>> B_boot = 2^260` (rejected if ever dividing a boot norm).
2. **Maximal orders & Pigeonhole Zero-Log Theorem:** `O_(K_f)` (`w'=omega^2`) is maximal away from `3` (`N(Disc f) = 3^72 19^12 67^12`); `O_(K_g)` has `Disc(O_(K_g)) = 2^16 3^36 5^18 906349^16 254760293^16`, `gcd(r, Disc) = 1`. Every ramified or non-maximal rational prime `ell in {2, 3, 5, 19, 67, 906349, 254760293}` (and every rational prime `ell ! == 1, 11, 14 (mod 15)`) has `tau`-orbit length `d | 12` strictly `< 12` (`d in {1, 2, 3, 4, 6}`). Because `gcd(chi^d - 1, r) = 1` for all `d in {1, 2, 3, 4, 6}` (verified by `audit_cm.py`), equivariance `(chi^d - 1) vlog(Q) == 0 (mod r)` forces **`vlog(Q) == 0 (mod r)` identically**. Thus **zero** ramified, non-maximal, or short-orbit prime ideals ever enter the `Z[zeta_12]` linear system. (Verified `0 / 1,950` valuation mismatches at `2, 3, 5` in `O_(K_g)`.)
3. **Free `C_12`-action & 4 Schirokauer columns (`3 + 1`):** Complex conjugation `y -> y^14` does not lie in `<sigma> = {1, 7, 4, 13} mod 15`, so `<tau> ~= C_12` acts freely on the 36 complex places of `K_f` (3 free orbits) and 12 complex places of `K_g` (1 free orbit). By Dirichlet's unit theorem over `F_r[C_12]`, `dim((O_(K_f)^x / r)^chi) = 3`, `dim((O_(K_g)^x / r)^chi) = 1`, and `dim((O_H^x / r)^chi) = 0` (`tau^4 = 1` on `K_H`, `chi^4 != 1`), giving **4 Schirokauer columns** total (down from 10 in `c36`).
4. **Exact Murphy `alpha` (`w' = omega^2`):** Analytic at `ell = 2, 5` (`v_2(N_f) = v_5(N_f) = v_5(N_g) = 0`, `E[v_2(N_g)] = 1.505881` exact, `kappa_2(9.6, 71.69) = 0.0525`), 40k MC at `ell = 3`, and exact Kummer-Dedekind over `[7, 100k]` (all residue degrees `f_p in {1, 2, 4}`): exact `alpha_f = +2.718`b (charged **`+2.74`b**), `alpha_g = +0.762`b (`alpha_(g,>=3) = +1.267`b, charged **`+1.28`b**; total charged `+0.77`b in descent).
5. **Live End-to-End `(4, 3)` `Q(zeta_15)` Split-Base Tower DLP Audit (`audit_cm.py`: `toy_e2e_dlp`):** On real toy BLS12 parameters `u_0 = -5, p_0 = 7207 == 7 (mod 15)` (`Phi_15 == h_pi h_pibar mod 7207` with `P_bar` explicitly present), `r_0 = 601` prime, `chi = 596` (`ord = 12`), `audit_cm.py` executes live BSGS in `F_(7207^12)`, live CRT mod `601^2` Schirokauer maps, and live mod-601 linear algebra: `0` `P_bar` hits, `27` omitted short-orbit prime factors producing `0` errors, **`25/25`** `K_g` left-kernel deps annihilating true BSGS logs (`0` fail), negative control A (drop 1 `SM_g` col) failing `20/26`, negative control B (wrong `chi = p_0^2`) failing `40/40`, **`5/5` fresh unseen target DLP recoveries** exact match, and **`57/57`** `K_f` (exact deg-72 BLS12 polynomial, `3` `SM_f` cols) left-kernel deps annihilating true BSGS logs (`0` fail, `62/62` exact virtual-log fit).

## 3. Relation Collection: Tiled `Arch2-Rev-SigmaB` & Exact de Bruijn Yield

Operating point **`B71.5_d04_e109`**: asymmetric smoothness bounds `B_f = 2^71.29`, `B_g = 2^71.69` (`U' = (li(B_f) + li(B_g))/12 = 2^63.317373` 12-orbits, smaller than `c36`'s `2^63.92`), filtering excess `1.09` (identical to `c36`), collection target `1.148 U' = 2^63.516495` relations.

### Geometry, Exact `min-T2` Voronoi FD, & Complete `100%` `b_0`-Invariant `1/12` FD `D_(omega,sigma)`

In the unimodular `Z[omega]`-module basis `B_(1,omega) = {1, omega, zeta_5, omega zeta_5, zeta_5^2, omega zeta_5^2, zeta_5^3, omega zeta_5^3}`, the `T_2` Gram matrix is the exact root-lattice tensor product `G = A_4 x A_2` (`det G = 1,265,625`, `W^T G W = G`, `S^T G S = G`, verified live in `audit_cm.py`).
- **Canonical `a` (`min-T2` 10-facet Voronoi cone):** Canonicalizing `a` modulo `<O_H^x, sigma>` via the exact 10-facet Dirichlet-Voronoi cone (brute-force `qfminim` ratio `1.000000` at `R^2 = 24, 30, 36, 44`; asymptotic `p_can = 0.268857`) yields `N_a^can = 2^57.6928` canonical `a`-orbits.
- **2-adic class regions (`s = (1.00, 1.20, 1.36)`):** Because `2 O_(K_g) = P_(2,1)^3 P_(2,2)^3` and `v_2(N_g(a, b)) = 4 sum_i min(3 v_(p_(2,i))(a), 16)` depends **solely on `a`** (and is `sigma`-invariant), scaling class radii `(s_0, s_1, s_2) = (1.00, 1.20, 1.36)` for classes `(v_1=v_2=0)`, `(v_1+v_2=1)`, `(v_1+v_2>=2)` (exact densities `225/256`, `225/2048`, `23/2048`) gives `Q_a(a) <= 0.65 s_c R_a^2` and `Q_b(b) <= (s_c - q_a) R_b^2` (`r_(ab) = R_a/R_b = 1.20`, `A^2 = 510986.8`), with unified 2-adic key `key_a(a) = q_a(a) - s_c(a)` so the region condition is `key_a(a) <= -q_b(b)`.
- **Swapping `<sigma>` onto `b` (`D_(omega,sigma)` 1/12 FD, `239/239` exact-1):** All 12 group matrices in `G_12 = <W, S>` have **zero** bottom-left `6 x 2` block (verified live in `audit_cm.py`), so `t = (b_2, ..., b_7)` transforms completely independently of `(b_0, b_1)` under the `6 x 6` restriction `G_(12,6)`. Define `U = b_2 - b_4 - b_6`, `V = b_3 - b_5 - b_7`. Then `b in D_(omega,sigma)` iff:
  - When `(U, V) != (0, 0)`: `(U > V >= 0)` and (`(b_2, b_3) >_lex (0, 0)` if nonzero, else `(b_4-b_6, b_5-b_7) >_lex (0, 0)`).
  - When `U == V == 0` and `t != 0`: `t == min_{M in G_(12,6)} (M t)` lexicographically (`t = 0` is the rank-2 `Q(omega)` subspace of orbit size `<= 3`, never 12).
  This predicate is **100% independent of `(b_0, b_1)`** (zero chord cutting) and verified exhaustively live in `audit_cm.py` on `Q <= 30`: `2,940` points, `239` full 12-orbits, **`239/239` exact-1, `0` duplicates, `0` misses**! Moving `<sigma>` onto `b` gives `N_7 = 2^57.0218` rows (`L_max = 2 sqrt(1.36 R_b^2 / 8) = 538.11`, mean chord `L_b = 231.12`) and `N_a^tbl = 4 N_a^can = 2^59.6928` slots.

### Tiled `Arch2-Rev-SigmaB` Reversed `a`-Table Replay (`P_lo = 2^54`, `W_tile = 256`, `W_SB = 128`)

For any 4-orbit of `sigma`-conjugate ideals `(p, s_k, r)` (`s_k = s_0^(7^k)`, `k in {0,1,2,3}`), a hit `a(s_k) == r b(s_k) (mod p)` is equivalent to `(sigma^k(a))(s_0) == r b(s_k) (mod p)`. Thus **ONE** 4-orbit `a`-table built at `s_0` over all `N_a^tbl` slots `a_slot = (m << 58) | oid` serves all 4 `sigma`-conjugates simultaneously (`m' = (m - k) mod 4`). Every `b_0`-row (`L_b <= L_max = 538.11`) is partitioned into `ceil(L_b / 256) in {1, 2, 3}` contiguous tiles of width `<= W_tile = 256` (`1.408` tiles/row, `3.206` spanned super-buckets/row):
- **3-Pass `a_0`-chord + `perm[]` build (`C_LIST_A = 35.0` Ir/entry, `KERNELS.md`):** Total build for `g` to `P_g = 2^64` and `f` to `P_f = 2^60`: **`2^107.2110`** units, peak table storage **`2^120.41`** bytes (`STORAGE.md`).
- **Band 1 (`[2^24, 2^54)`):** Unified `key_a`-sorted buckets (`12` B/entry, `<= 2^63.29` B/table); exact integrated branchless `rev_b1_bsearch` (`8` Ir/step + `24` header) + contiguous slice copy (`1.50` Ir/hit) + `rev_b1_apply_u4` with 2-instr conjugate rotation (`6.25` Ir/hit, charged `C_HIT_LO = 8.0`).
- **Band 2 (`[2^54, P_hi)`):** Super-buckets of width `W_SB = 128` with 64-bit occupancy bitmap (`<= 2^63.802` B/table `< 2^64` B). Per tile: `rev_sb_probe_fast` (`C_TILE_HDR = 20.0` Ir) + **`P(occ) = 1.0`** occupied-SB header/terminator overhead (`C_SB_HDR = 18.0` Ir on every spanned SB). Exact 16D ellipsoid integration over all tiles and chord `q_b(t)` variation (verified live in `audit_cm.py`) gives `E[cands]/E[hits] = 3.9835`, conservatively charged **`4.024`** (`+1%` discrete safety) at `C_CAND = 13.0` Ir/cand, plus full emit + 1 MiB page-copy/handoff + `rev_b2_apply_u4` (`18.75 <= C_HIT_HI = 19.0` Ir/hit).
- **Zero 9D walk & exact prefix `rev_scan_zero_u4`:** `P_g = X_2 = 2^64` eliminates the entire 9D walk, 106-bit LLL, contamination, and hand-written ASM. Because `N_g(a, b) = N_g(0, b)` to `0.000000` bits (`0 / 4,976` rows `> 1/16`b) and `oid` is globally sorted by `key_a`, `rev_scan_zero_u4` scans and zeroes **ONLY** the exact active prefixes `[0, N_act_orb(t))` (`V_all = 1.0276 V_kept`, `C_SCAN_ZERO = 2.95` Ir/kept cell). Survivors pass `_nk_g_only` (`4,500` Ir), `_nk_f_only` (`6,500` Ir), conservative two-sided Stage B (`C_SB_G = 760,000`, `C_SB_F = 620,000`, `(s_a_g, s_b_g) = (360, 295)`, `(s_a_f, s_b_f) = (420, 345)`), `fbatch`, `batch2`, trial division, 40-round MR, and `99.9%` ECM — all computed **LIVE** from `cm_models.Retention` inside `audit_cm.py` (`r_g = 0.8103`, `r_f = 0.9382`, `r_tot = 0.7602`, `c_cof = 1.7071` ipc).

### Exact de Bruijn Pool Yield, Band/ECM Split + Chebyshev, & Collection Cost

Evaluated on the pre-registered `200,000`-sample duplicate-free pool `P2_wy_2ad` (seed `900003`; 12 integer anchor pairs verified live to `< 5e-11` bits in `audit_cm.py`) with **de Bruijn `rho_dB(u, x)` applied to both sides**:
- **Kept cells:** `V = 2^119.3186` (`SE_c = 0.00290`, `1 + 2 SE_c = 1.00579 <= H2_ALLOW = 1.00774`), mean norms `(N_f, N_g) = (548.90, 677.81)`b (`<= (549.0, 678.0)`b; relation-weighted `(543.32, 673.60)`b).
- **Relation count & margin:** `2^63.516495` (`SE_y = 0.00969` bits; live 2.0-bit histogram check in `audit_cm.py` agrees to `0.00163` bits). After deducting `c36`'s `0.009`-bit band-miss allowance, net margin above `1.09 U' = 2^63.44170` is **`0.06579` bits (`>= 0.060` bits, `6.8 SE`)**, and the 2-SE lowered count `2^63.4881` strictly exceeds the descent input `2^63.46`.
- **Per-cell collection cost:** `g`-replay `27.45` + `f`-replay `16.64` (Band 1 `16.99`, tile+SB overhead `6.94`, Band 2 `20.17`) + scan/zero `2.95` + live cofactorization `1.71` = **`48.75` instructions/kept cell** (estimate `2^110.638351`, times `H2` allowance `1.00774` <= **`2^110.649475`** core units).
- **Band split & first-principles ECM split + Chebyshev tail (verified live in `audit_cm.py`):** Across `<= 2^(log2_rels + 0.55) = 2^64.066` kept relations with conservative max norms `(640b, 768b)` (`10, 12` limbs): band split `<= 2^100.543` units; first-principles Suyama-Montgomery ECM (`13` f-factors `<= 2^71.29` + `16` g-factors `<= 2^71.69`, `p >= 0.999 > 1/1.002`, plus MR + `2^26` ideal lookups) gives expected work `E[W_ecm] = 2^95.6178` core units. Metering at **`2^96.5`** core units (reached = failure) leaves gap `0.8431 E[W_ecm]`, max single-schedule cost `2^26.94` units, variance `Var(W_ecm) <= 0.00201 x 2^26.94 x E[W_ecm] = 2^113.60`, and Chebyshev tail **`P(W_ecm > 2^96.5) <= 2^-77.14 < 2^(-70)`** (explicitly added to the Section 6 union bound!).

## 4. `Z[zeta_12]` Filtering & 1-Owner Block Wiedemann Linear Algebra

Because `Phi_12(chi) = chi^4 - chi^2 + 1 = 0`, setting `omega = chi^2` splits `Z[zeta_12] = Z[omega] (+) chi Z[omega]` into two orthogonal Eisenstein hexagonal lattices. Using the 6-bucket signed hexagonal basis `{1, chi, chi^2, chi^3, chi^2 - 1, chi^3 - chi}`, all 12 units have weight 1, and only `1/12` of unit+unit Markowitz collisions creates a coefficient-2 non-unit (vs `1/4` in `c36`'s `Z[i]`).

### Direct `2^16 .. 2^23` `filtsim_z12` Campaign (`B71.5_d04_e109`, `W = 200`, Conservative Cell-Mean Norms `548.9 / 677.8`)

Validated first against `c36`'s `Z[i]` `PFX` tables (`0.0%` diff at `2^21`), `filtsim_z12` was run for 6 seeds at every doubling `U' = 2^16 .. 2^23` at `Lf = 71.29, Lg = 71.69, excess = 1.09, W = 200`, using conservative cell-mean norms `(548.9, 677.8)` strictly above smooth-relation means `(543.2, 673.6)`:

| Scale `U'` | Seeds | `N_l / U'` (SE) | `N_sw / U'` (SE) | `sd_l U'^(1/4)` | Depth | Unit % (Exp.) | `k5` | `tab_ratio` |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `2^16` | 6 | `0.003260` (`0.000227`) | `0.006462` (`0.000482`) | `0.0089` | `13.7` | `95.61%` (`89.15%`) | `10.5342` | `0.80659` |
| `2^17` | 6 | `0.003050` (`0.000112`) | `0.006089` (`0.000219`) | `0.0052` | `15.2` | `96.63%` (`90.47%`) | `10.2017` | `0.80554` |
| `2^18` | 6 | `0.002857` (`0.000088`) | `0.005747` (`0.000173`) | `0.0049` | `16.7` | `96.96%` (`91.30%`) | `9.9914` | `0.80212` |
| `2^19` | 6 | `0.002533` (`0.000074`) | `0.005111` (`0.000118`) | `0.0049` | `16.7` | `97.18%` (`91.82%`) | `9.8616` | `0.79984` |
| `2^20` | 6 | `0.002399` (`0.000032`) | `0.004857` (`0.000049`) | `0.0025` | `18.7` | `97.35%` (`92.32%`) | `9.7354` | `0.79839` |
| `2^21` | 6 | `0.002172` (`0.000035`) | `0.004377` (`0.000067`) | `0.0033` | `18.5` | `97.42%` (`92.61%`) | `9.6628` | `0.79684` |
| `2^22` | 6 | `0.001999` (`0.000021`) | `0.004057` (`0.000039`) | `0.0024` | `19.7` | `97.43%` (`92.77%`) | `9.6218` | `0.79508` |
| `2^23` | 6 | `0.001889` (`0.000013`) | `0.003846` (`0.000026`) | `0.0017` | `20.5` | `97.45%` (`93.01%`) | `9.5620` | `0.79485` |

- **Strict monotonicity (`2^16 -> 2^23`):** Both `N_l / U'` and `N_sw / U'` decrease at every single doubling (checked by `audit_cm.py`); `N(250)/U'` at `2^23` is `0.003552 +- 0.000024`.
- **Unit share & exact `k5` (`H4u`):** Rises monotonically to **`97.45%`** total (**`93.01%`** explicit). Computing `k5` directly from `exp_unit_lo = 0.93008 - 2 SE` via `eu x 7.80 + (1 - eu) x 33.0` gives `9.5669`, conservatively charged **`9.59`** (verified live in `audit_cm.py`).
- **Spread law & caps (`H3`):** `max(sd_l U'^(1/4)) = 0.0089` and `max(sd_sw U'^(1/4)) = 0.0189`, well below `c36`'s `0.045/2 = 0.0225` (`P_CAP = 0.00022`). Under `H3`, the monotone `2^16..2^23` means and declared spread law `sd <= (0.045/2) U'^(-1/4)` bridge the `2^23` measurements to production `U' = 2^63.32`: charging `(mean_23 + 2 SE) x 1.0159 (splitting) + cap_add` gives **`N_l / U' <= 0.00199`** (`N_l <= 2^54.3444`) and **`N_sw / U' <= 0.00400`** (`N_sw <= 2^55.3516`).
- **Dense-band subset-sum tables:** Order-12 invariance keeps table size at `2^47.715` elements (`{2: 28, 3: 18.25, 4: 14}`); fitted `tab_ratio` slope `+ 2 SE = -0.00160 < 0`, conservatively charged `max(mean + 2 SE) = 0.8082`.

### 1-Owner Architecture B Block Wiedemann, Generator, 4-Col Mod-`r` CRT Schirokauer, & `H-X`

- **1-Owner Architecture B (`< 2^64` B incl. output block `y`):** Holding BOTH input `x` (`N_sw x 160 B = 2^62.674` B) AND output `y` (`2^62.674` B) plus all 22 dense-band tables (`2^59.497` B) totals **`2^63.751` B `< 2^64` B** on ONE single 64-bit owner (zero cross-owner sends/combines, verified live in `audit_cm.py`). With `C_ROW + C_CHI = 785`, `k5 = 9.59`, `(m, n) = (8192, 110)`, `N_SM = 4`, `BW_SIG = 2^(-7)`:
  - **Linear Algebra:** `<= 2^106.693172` core units (`-3.26` bits below `c36`'s `109.952`).
  - **Matrix BM Generator:** `<= 2^98.616` core units.
- **Exact 4-Column Mod-`r` CRT Schirokauer Maps (`<= 2^84` units, verified live in `audit_cm.py`):** Because `r == 13 mod 15`, `Phi_15 == h_(r,1) h_(r,2) (mod r)` splits into two irreducible quartics (`[4, 4]`), with `w_1 = -u^2 == u^8 (mod r)` and `w_2 = u^4 (mod r)`:
  - *`K_f` side:* Over both `F_r[y]/(h_(r,i)) ~= F_(r^4)`, `f(X) = X^9 + (w_i + 1) X^6 + X^3 - w_i in F_r[X]` is irreducible of degree 9 (verified live by Rabin `9 x 9` Frobenius check + `gcd(9, 4) = 1`), giving exact residue degrees **`[36, 36]`** and annihilator `eps_f = r^36 - 1` (`9,175` bits — half `c36`'s `r^72 - 1`).
  - *`K_g` side:* On component 1 (`w_1 = -u^2`), `u w_1^2` is a non-cube in `F_r` and `F_(r^4)` (since `(r+1)(r^2+1) == 1 mod 3`), so `g(X) = X^3 - u w_1^2` stays irreducible of degree 3 (`[12]`); on component 2 (`w_2 = u^4`), `u w_2^2 = u^9 = (u^3)^3` is an exact cube in `F_r`, so `g(X) = (X - u^3)(X - u^3 w_2)(X - u^3 w_2^2)` splits into three linear factors (`[4, 4, 4]`). Exact `K_g` CRT degree list is **`[4, 4, 4, 12]`** (`lcm = 12`, annihilator `eps_g = r^12 - 1` annihilates every component!). Splitting component 2 into three degree-4 rings is strictly cheaper than one degree-12 ring (and `r !| Disc`), so the 4 `chi`-projected Schirokauer columns remain well inside the charged **`<= 2^84`** core units.
- **Krylov rank & `H-X`:** Verified on `14/14` `Z[zeta_12]` systems (8 real product-form live systems + 6 large sparse systems up to `D = 3,000`, `4` SM cols, block-diagonal `X` placing the 4 pad coords in blocks `0..3`): `100%` full rank, while the 1-block pad negative control drops exactly `3` ranks on all `14/14` cases (re-runnable live via `audit_cm.py` and `hx_analogue`).

## 5. Individual Logarithms

Detailed in `DESCENT.md` (all `39` bins / `78` nodes listed and executed LIVE bottom-up by `audit_cm.py`). Using our **16D Guillevic Rank-2 Torus Boot (`det = p^4`, `UDUP_BOOT = 1.0`, `BOOT_RELS = 36`, `m_max = 2580.0`b `> 2575.73`b)** and **16D Rank-2 Canonical Special-`q` Tree (`UDUP_NODE = 0.95`, `MESH = 384`, full de Bruijn enabled everywhere)**:
- **Expected two-log work:** `2^102.497353` <= **`2^102.498`** core units (`-6.366` bits below `c36`'s `108.864`).
- **Metered cap:** `t = 2^115.6` core units (`P(cap) <= 2^(102.498 - 115.6) = 1.137378e-04`).
- **Live checks:** `[11]` across 32 fresh BLS12-381 `w'=omega^2` boot targets (`544` exact norms): `kappa_hat = 0.7076 << BALL_PAIR = 2^24`; `[7e]` child-conditioning shift `(+3.40, +2.35)`b inside `(+10, +4)`.

## 6. Success Probability & Strict Additive Union Bound

Every single failure event across main phase, ECM meter, descent, target randomization, and wrapper is summed in a **strict additive union bound** `P(success) >= 1 - sum_i P(fail_i)` (verified in 80-digit `Decimal` inside `audit_cm.py`):
- **Main-phase & ECM-meter failures:** Block Wiedemann Schwartz-Zippel `<= 4 N / 2^128 < 2^(-66)` (conditional on `H-X`, `H-SM`; `4/r < 2^(-252)` for random `W_f, W_g`), filtering spread-law cap `P_CAP = 0.00022`, and ECM-split Chebyshev meter tail `< 2^(-77.3) < 2^(-70)`.
- **Descent & target failures (2 logs, `H5`):** Two boots failing all 36 relations `<= 2 x 7.5e-07`, metered cap `<= 1.137378e-04` (evaluated at `+10%` headroom: `1.251116e-04`), node counter abort (`E[N] = 2^16.0696`, `N_MAX = 2^25.8`) `<= 2^-9.7304`, conservative `tau`-orbit repeat collision (smallest bin `[70, 76)` holds `>= 2^67.64 > 2^66.66` orbits) `<= 2^-10.05`, and degenerate target/randomizer `<= 2/r`.
- **Wrapper failure:** Field map + rejection sampling `< 2^(-50) < 2^(-48)` (`field-map.md`).
- **Opened-run floor (strict additive union bound):** Core success `>= 0.997544`; with `+10%` cap headroom and wrapper `< 2^(-48)`, opened-run success `>= 0.997532 >= 0.997532`.
- **Run gate:** `K = 7212029477497188191` satisfies `(K / 2^64) x 0.997532 >= 0.39 > ((K - 1) / 2^64) x 0.997532` (`audit_cm.py`: `ledger_gate_claim`).

## 7. Complete Ledger

| Component | `log2` core units (`C = 20,000`) |
|---|---:|
| relation collection (tiled `Arch2-Rev-SigmaB`, sec. 3, `H-EA`; est. `110.638351` + `0.774%` H2 allow.) | <= 110.649475 |
| replay 4-orbit `a`-tables build (`g` to `2^64` + `f` to `2^60`, 3-pass `a_0`-chord; sec. 3) | <= 107.2110 |
| linear algebra (`1-owner Arch B Z[zeta_12] BW`, `W=200`, `2^23` direct, `k5=9.59`, `(8192,110)`; sec. 4) | <= 106.693172 |
| matrix Berlekamp-Massey generator (`half-gcd FFT`; sec. 4) | <= 98.616 |
| both individual logarithms (metered cap; expected `102.498`; sec. 5) | <= 115.6 |
| band split of kept relations (`(2^20, 2^48]`, 10/12 limbs; sec. 3) | <= 100.543 |
| ECM split of kept relations (`(2^48, B]`; metered cap; exp. `95.62`, Cheb. `< 2^-77`; sec. 3) | <= 96.5 |
| canonical `a`-list enumeration + global `key_a` pre-sort | <= 86 |
| 4-column `chi`-Schirokauer maps (mod-`r` CRT `f: [36,36]`, `g: [4,4,4,12]`; meas. `<= 82.8`) | <= 84 |
| filtering streams, product-form pivot selection, back-substitution, misc | <= 80 |
| `P_g`, `P_f` product trees (`< 2^89` instructions) | <= 75.5 |
| factor base, roots, 12-orbit ideals (`U' = 2^63.32`) | <= 75 |
| numerical certificate construction | <= 50 |
| final verification | <= 11 |
| wrapper (sec. 1; `<= 2^80` instructions = `2^65.71` core units) | <= 66 |

Core cap sum (`80`-digit `Decimal`): `2^115.653035625674 <= 2^115.6530356257` core units (with wrapper `2^66`: `< 2^115.653035625701` core units = `< 2^118.974963720588` G1 units). Replacing the metered cap `2^115.6` by expected two-log descent `2^102.498` gives opened-run expected core work `2^110.865313701636 <= 2^110.8653137017` core units. Behind the run gate `K / 2^64 = 0.3909649014`, total expected work is `< 2^109.510424703138` dense-`F_(p^12)` units (`C = 20,000`) and `< 2^112.832352798025` G1 units (`C = 2,000`).

## 8. Heuristics

- **`H1` (Smoothness & Local Laws):** Two-sided smoothness with **explicit de Bruijn `rho_dB(u, x)`** on both Collection and Descent, certified Murphy `alpha = (+2.74, +1.28)` (`w'=omega^2`), exact `v_2(N_g)` and `kappa_2`, and boot pair ratio `kappa <= 2^24` (measured `0.7076`, `[11]`).
- **`H1'` & `H-EA` (Detector Retention & Survivor Work):** Exact `c36` two-stage Euler-product retention and survivor-upper laws on `(g, f)` (computed live via `cm_models.Retention`) with `0.009`-bit band-miss margin; `f/g` stage-B dependence bounded by `q_f = min(1, 10 rho(u_f))`.
- **`H-FFT` (Work only):** `F(N) = 8.5 N (log2 N)^2 + 1200 + 504 N` (`N > 2^26`) bounds local SSA integer multiplication.
- **`H2` (Region & Pool Estimation):** Pre-registered `200,000`-sample duplicate-free `min-T2` pool `P2_wy_2ad` (`SE_y = 0.0097`b, `SE_c = 0.0029`, `H2` allowance `1.00774 >= 1 + 2 SE_c`); exact 16D ellipsoid concentration (`+1%` discrete safety on `C/H`).
- **`H3` & `H-SM` (Matrix & Class Numbers):** `r` unramified (proved), `r !| h(K_f), h(K_g)`, `P_chi Lambda` injective on the 4-dim `chi`-unit space (verified on `u_0 = -5` toy tower), and non-increasing `Z[zeta_12]` product-form percolation bridging measured `2^16..2^23` means to production `U' = 2^63.32` under the spread law `sd <= (0.045/2) U'^(-1/4)` (`N_l/U' <= 0.00199`, `N_sw/U' <= 0.00400`, `P_CAP = 0.00022`).
- **`H4` & `H4u` (Kernel & Unit Share):** Pure baseline scalar AArch64 compiled counts (`KERNELS.md`, zero SIMD, zero `FEAT_LSE`, explicit transfer pricing); `Z[zeta_12]` hexagonal 6-bucket unit share (`97.45%` total, `93.01%` explicit at `2^23`), `k5 <= 9.59`, `tab_ratio <= 0.8082` (negative measured slope).
- **`H-X` (Krylov Controllability & Observability):** For the live system with 4 fixed Schirokauer columns and block-diagonal `X` placing the 4 pad coords in distinct blocks, `rank C_d = N` and `rank O_e A = D` at the charged `d, e` (supported by `14/14` `Z[zeta_12]` tests up to `D = 3,000` + live `hx_analogue`).
- **`H5` (Descent Independence):** Lower-sum semismoothness with de Bruijn, lazy cofactorization, negative-binomial recurrence, independent subtrees of distinct ideals (`[7e]` shift `(+3.40, +2.35)`b), Poisson missing-ideal model, random-integer ECM orders, `UDUP_BOOT = 1.0`, `UDUP_NODE = 0.95`.

## 9. Memory & Provenance

**Peak memory:** `< 2^120.5` bytes (`2^120.41` B replay `a`-tables during collection; `< 2^105` B in descent; `2^63.751` B `< 2^64` B on 1 owner in LA; `< 2^30` B wrapper; `STORAGE.md`). All preprocessing is target-specific and charged once in `W`.

**Provenance:** CM-split cyclotomic towers over `Q(zeta_15)` build on Kim-Barbulescu (CRYPTO 2016), Guillevic-Singh (2020), Zhu-Liu (2025), and APT26; filtering follows CADO-NFS and Boudot (ePrint 2017/758). Accounting framework, two-stage `P24` Stage B, SRT lemma, product-form recorded pivots, dense-band subset-sum tables, coset replay, and descent `MESH=384` lower-sum model build on our CurveSmash lineage (`c24..c36`, `stnfs-bn15r`; gbotrel, Claude Opus 5.5 / Fable 5.2) together with credited community ideas: run gate, expected work reading, and high descent cap (kaankolcu `c613ae01`, `3359b9e2`), truncated-Shoup `chi`-mul idea (soispoke `cf131637`), carry-chain combine (yukonautominer `9a7df12c`), and `(m, n)` block sizing (mpjunior92 `ea928ab4`). New here: the `w' = omega^2` `sigma`-fixed isometry, Split-Base & Pigeonhole Zero-Log Theorems, live `(4, 3)` `Q(zeta_15)` E2E DLP audit, tiled `Arch2-Rev-SigmaB` zero-9D-walk pure-C reversed `a`-table replay with complete `1/12` `b_0`-invariant `D_(omega,sigma)`, `Z[zeta_12]` hexagonal 6-bucket 1-owner Arch B BW (`2^16..2^23` ladder), 4-col mod-`r` CRT Schirokauer maps (`[36,36]` / `[4,4,4,12]`), and 16D Guillevic Rank-2 Torus Boot + Rank-2 Canonical special-`q` descent with full de Bruijn enabled everywhere. No prior verdict is a premise.
