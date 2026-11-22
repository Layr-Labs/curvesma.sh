# Finite Relation Rank (`H-R`), Schirokauer Compression, Block Wiedemann (`H-X`), & Self-Contained `PM-Basis`

Separates `H-R`, Schirokauer compression (`W_f, W_g`), Block Wiedemann (`H-X`), and a **100% self-contained classic Giorgi-Jeannerod-Villard (2003) / Thomé (2002) `PM-Basis`** generator (`<= 2^93.0` u + live toy execution). No black-box generator citation.

## 1. Schirokauer Compression (`W_f, W_g`), Exact `F_r` Sampling, & `H-R`

- **Data-oblivious Schirokauer (`<= 4/r < 2^-252`):** `P_chi Lambda_f, P_chi Lambda_g` have `F_r`-dims `6, 2` (3D, 1D `chi`-unit images). Pre-collection uniform `3x6` `W_f`, `1x2` `W_g` over `F_r` form `4` columns `S`; under `H-SM`, injectivity loss `<= 4/r < 2^-252` (verified live in `toy_e2e_dlp`).
- **Exact 255b / `K_rej = 64` `F_r` sampling (`< 2^-156`):** Every uniform `F_r` entry of `W_f` (18), `W_g` (2), `Z` (`106*(N+4)`) draws 255b words rejecting `>= r` (`p_rej = 0.0943098873`), capped at `K_rej = 64` (`p_rej^64 = 2^-218.0126`). Over `< 2^61.25` entries (`N <= 2^54.5215`), union bound `< 2^-156.7 < 2^-150` (80-digit `Decimal` in `audit_cm.py`, absorbed in `2^-48` wrapper). `X` draws direct uniform 128b integers (`0` rejection).
- **`H-R` (fixed-instance premise):** Given `(W_f, W_g)`, filtered `A = [M | S]` (`D = N + 4`) is fixed; `H-R` assumes `rank(M) = N`, `rank(A) = D - 1`, so `ker(A)` is the 1D virtual-log line `(w, c)`, `c != 0`.

## 2. `H-X` & Conditional Schwartz-Zippel

Let `B = [M | 0]`, `Y = [S_1..S_4, B Z]` (`n = 112`), `X` partitioning the `4` pad coords into `4` blocks and `N` coords into `m = 8192` balanced blocks with uniform `[0, 2^128)` weights. **`H-X`** asserts at `d = ceil(D/n) + 1 + ceil(D/(128n))`, `e = ceil(D/m) + 1 + ceil(D/(128m))` that an `N x N` minor of `C_d = [B Z, B Y, ..., B^d Y]` and a `(D-1) x (D-1)` minor of `O_e U` are nonzero polynomials. **Conditional Schwartz-Zippel (proved | `H-R, H-X`):** failure `<= 4 N_cap / 2^128 < 2^-66`.

## 3. Self-Contained Classic `PM-Basis` (`<= 2^92.7203 <= 2^93.0` Core Units)

Given `S_k = X B^k Y in F_r^{m x n}` (`0 <= k < L_bm <= 2^47.7706`), compute `[S(z) | -I_m] [q(z); r(z)] == 0 (mod z^L_bm)` (`t = m + n = 8302`, `s = (0^n, 1^m)`) via recursive `PM-Basis(F, sigma, s)` (executed verbatim in `pm_adaptive.py`):
- **Base (`sigma = 1`, `M-Basis-1`):** Column-echelon on `F(0)` (`<= m t^2` `F_r` MACs at `MONT5 = 270` Ir) across `L_bm` leaves: `W_base = L_bm m t^2 * 270` Ir (`= 2^80.5982` u after `/20000`).
- **Recursion (`sigma > 1`, `s1 = floor(sigma/2)`, `s2 = sigma - s1`):** `(P1, s_mid) = PM-Basis(F mod z^s1, s1, s)`; `G = ((F*P1 mod z^sigma) div z^s1) mod z^s2`; `(P2, s_out) = PM-Basis(G, s2, s_mid)`; `P = P1 * P2`.
- **1D 4-Prime 62-Bit NTT Schoolbook `t^3` Matmul (ZERO Strassen/Karatsuba):** `M_top <= 4 L_bm = 2^49.7706 < 2^53`, so EVERY level uses 1D radix-2 NTT over 4 verified primes `p_i = c_i * 2^53 + 1 < 2^62` (`NTT4_PROD = 2^244.70`). Splitting `F_r` into 4 64b limbs gives max convolution `t L_bm 2^128 = 2^188.79 << NTT4_PROD` (`+55.9`b headroom, asserted LIVE). Per `t x t` product (`k = log2 M <= 50`): (1) fwd NTTs (`32 t^2` x `16 M k` = `512 t^2 M k` Ir); (2) pointwise `64 M` scalar `t x t` `F_{p_i}` matmuls (Barrett MAC `17 <= 24` ops -> `1536 M t^3` Ir); (3) inv NTTs + Garner CRT + mod-`r` (`1024 t^2 M k + 6220 t^2 M` Ir). Over `depth = 48` levels (`sum M_j <= 8 L_bm`): `W_rec = 2 * depth * (8 L_bm) * (1536 t^3 + t^2 (1536 k_max + 6220))` Ir (`= 2^92.7199` u after `/20000`). **Total:** `(W_base + W_rec)/20000 = 2^92.7203 <= 2^93.0` u (`< 2^84 B` RAM; Horner in `la_u`).

## 4. Zero-Loss `repurge_stack` & Live Diagnostics

Excess-row clique trim + Markowitz deg-0 guard give `post_peel_cols_lost == post_merge_cols_lost == 0`; 3-stack back-sub (`L_rec` -> `merge_stack` -> `repurge_stack`, `n_reconstruct_unit == 0`) recovers 100% of 2-core cols. Live: `unplanted_two_sided_bw` (unplanted 4-SM NFS `p=467, r=233`, live recursive `pm_basis` + Horner matching `ker(A)`, 38/38 cols, 49/49 rels, exact DLP) and `hx_z12.py` (exact Schur-vs-product `op(x)`, 2,498/2,498 rels, 0 lost cols). Production `H-R, H-X` remain declared exploratory premises.


## Live 100% Unscreened Zero-Oracle Two-Sided `Q(zeta_15)` Filter + UNPLANTED `H-X` + Same-Matrix Check (`w8_verifiers.py`)
- **Tower:** `p_0 = 7, r_0 = 181, chi = 7`, genuinely distinct non-isomorphic two-sided pair `g(X) = X^3 - w^2`, `f(X) = X^3 - (w - 1)w^2` (`N(w-1) = 3`, congruent mod `pi = (7, h_pi)` and fixed by order-12 `tau(y)=y^7, tau(X)=X*w`).
- **ZERO Oracle / ZERO Screening / ZERO RREF:** `dlog_GT` is NEVER computed until post-hoc verification; `bad_cols_screened == 0`; `rref_row_selection == False` (natural first `D = 21` surviving rows).
- **Live Recovery:** Pure 2-core (`34` rels / `29` orbits: `f=15, g=14`, natural `rk_M = 29/29`, `rk_A = 30/31`, `nullity = 1` with data-oblivious compressed `S_f, -S_g`) -> `_dc`-guarded zero-loss merge (`29 -> 19` live, `10` merged, `lost = 0`, `gen_k5_hex6 = 11.6 <= syn_matched = 12.544`) -> UNPLANTED BW `H-X` (`rC = 19`, `rOe_good = 21`, `rOe_bad = 20 < 21`, `rO = 20`, `rH = 20`) + Horner + back-sub -> **POST-HOC ONLY** independent BSGS in `F_(7^12)^*`: `lam = 135 != 0`, **`34/34` (`100%`) exact on BOTH `f` and `g` sides** (`33/34` non-zero).
- **TRUE Same-Underlying-Matrix Full-vs-Proxy Structural Check:** On ONE concrete `Z[zeta_12]` matrix `M_full` (`U = 1024, R = 1689`, zero `0.90` cap, `139` middle/dense columns with `deg >= 21`): filtering 100% of `M_full` (`0` omitted cols) gives `N_l_eff = 236 <= 264` and `N_sw = 14,226 <= 15,942` (`0.8924x`, `lost = 0/0`) vs omitting those exact `139` columns from the EXACT SAME matrix and charging them as static `w_mid`.

## 5. Exact Zero-Cost Factor-Base Blinding Normalization (`gamma_0`, `V_0`, & Exact `lambda`-Cancellation)

Linear algebra recovers a nonzero homogeneous virtual-logarithm functional `V(t) = lambda * log_g0(t) mod r` with unknown nonzero scalar `lambda in (Z/rZ)^*`. To blind target elements `z in {h, T}` uniformly in `G_T` WITHOUT computing `V(g_0)` or performing a third descent:
1. Pick ONE Stage-1 collected relation `(a_0, b_0)` whose `g`-side prime ideals all lie in the factor base and whose known linear-algebra virtual log `V_g(a_0 - b_0 X) != 0 mod r` (probability `1 - 1/r` per relation).
2. Form the **calibrated `G_T` blinding generator** (`1` final exponentiation `<= 2^11` u, ZERO extra descent):
   - `gamma_0 := phi_iso^(-1)(a_0 - b_0 X)^((p^12-1)/r) in G_T \ {{1}}` (exact prime order `r`, verified by `gamma_0 != 1, gamma_0^r = 1`),
   - `V_0 := ((p^12-1)/r) * V_g(a_0 - b_0 X) mod r != 0` (KNOWN EXACTLY from the linear algebra output vector; since `gcd(r, (p^12-1)/r) == 1`, asserted LIVE in `audit_cm.py`, `V(gamma_0) = V_0 = lambda * log_g0(gamma_0) mod r` identically).
3. For each challenge target `z in {h, T}`: each fresh boot trial draws uniform `e in [0, r-1]`, forms `z' := z * gamma_0^e` (EXACT uniform in `G_T`), descends `z'` to recover `V(z') = V(z) + e * V_0 mod r`, and **unblinds using the KNOWN `V_0`**: `V(z) := (V(z') - e * V_0) mod r = lambda * log_g0(z) mod r`.
4. Finally, `V(z) * V(h)^(-1) mod r = (lambda * x) * (lambda)^(-1) = x mod r` (for `bls12-381-g1-ecdlp-plain`, where `h = e(P_0, Q_0)` and `z = e(T, Q_0) = h^x`) cancels the unknown common scale `lambda` IDENTICALLY! Executed and verified **LIVE** with ZERO BSGS oracle inside `w8_verifiers.py` (`gamma0_blind_unblind_lambda_cancel_ok = True`).
