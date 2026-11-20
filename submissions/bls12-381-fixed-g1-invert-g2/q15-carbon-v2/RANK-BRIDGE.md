# Finite Relation Rank (`H-R`), Schirokauer Compression, Block Wiedemann (`H-X`), & Self-Contained `PM-Basis`

Separates `H-R`, Schirokauer compression (`W_f, W_g`), Block Wiedemann (`H-X`), and a **100% self-contained classic Giorgi-Jeannerod-Villard (2003) / Thomé (2002) `PM-Basis`** generator (`<= 2^93.0` u + live toy execution). No black-box generator citation.

## 1. Schirokauer Compression (`W_f, W_g`), Exact `F_r` Sampling, & `H-R`

- **Data-oblivious Schirokauer (`<= 4/r < 2^-252`):** `P_chi Lambda_f, P_chi Lambda_g` have `F_r`-dims `6, 2` (3D, 1D `chi`-unit images). Pre-collection uniform `3x6` `W_f`, `1x2` `W_g` over `F_r` form `4` columns `S`; under `H-SM`, injectivity loss `<= 4/r < 2^-252` (verified live in `toy_e2e_dlp`).
- **Exact 255b / `K_rej = 64` `F_r` sampling (`< 2^-156`):** Every uniform `F_r` entry of `W_f` (18), `W_g` (2), `Z` (`106*(N+4)`) draws 255b words rejecting `>= r` (`p_rej = 0.0943098873`), capped at `K_rej = 64` (`p_rej^64 = 2^-218.0126`). Over `< 2^61.25` entries (`N <= 2^54.5215`), union bound `< 2^-156.7 < 2^-150` (80-digit `Decimal` in `audit_cm.py`, absorbed in `2^-48` wrapper). `X` draws direct uniform 128b integers (`0` rejection).
- **`H-R` (fixed-instance premise):** Given `(W_f, W_g)`, filtered `A = [M | S]` (`D = N + 4`) is fixed; `H-R` assumes `rank(M) = N`, `rank(A) = D - 1`, so `ker(A)` is the 1D virtual-log line `(w, c)`, `c != 0`.

## 2. `H-X` & Conditional Schwartz-Zippel

Let `B = [M | 0]`, `Y = [S_1..S_4, B Z]` (`n = 110`), `X` partitioning the `4` pad coords into `4` blocks and `N` coords into `m = 8192` balanced blocks with uniform `[0, 2^128)` weights. **`H-X`** asserts at `d = ceil(D/n) + 1 + ceil(D/(128n))`, `e = ceil(D/m) + 1 + ceil(D/(128m))` that an `N x N` minor of `C_d = [B Z, B Y, ..., B^d Y]` and a `(D-1) x (D-1)` minor of `O_e U` are nonzero polynomials. **Conditional Schwartz-Zippel (proved | `H-R, H-X`):** failure `<= 4 N_cap / 2^128 < 2^-66`.

## 3. Self-Contained Classic `PM-Basis` (`<= 2^92.7203 <= 2^93.0` Core Units)

Given `S_k = X B^k Y in F_r^{m x n}` (`0 <= k < L_bm <= 2^47.7706`), compute `[S(z) | -I_m] [q(z); r(z)] == 0 (mod z^L_bm)` (`t = m + n = 8302`, `s = (0^n, 1^m)`) via recursive `PM-Basis(F, sigma, s)` (executed verbatim in `pm_adaptive.py`):
- **Base (`sigma = 1`, `M-Basis-1`):** Column-echelon on `F(0)` (`<= m t^2` `F_r` MACs at `MONT5 = 270` Ir) across `L_bm` leaves: `W_base = L_bm m t^2 * 270` Ir (`= 2^80.5982` u after `/20000`).
- **Recursion (`sigma > 1`, `s1 = floor(sigma/2)`, `s2 = sigma - s1`):** `(P1, s_mid) = PM-Basis(F mod z^s1, s1, s)`; `G = ((F*P1 mod z^sigma) div z^s1) mod z^s2`; `(P2, s_out) = PM-Basis(G, s2, s_mid)`; `P = P1 * P2`.
- **1D 4-Prime 62-Bit NTT Schoolbook `t^3` Matmul (ZERO Strassen/Karatsuba):** `M_top <= 4 L_bm = 2^49.7706 < 2^53`, so EVERY level uses 1D radix-2 NTT over 4 verified primes `p_i = c_i * 2^53 + 1 < 2^62` (`NTT4_PROD = 2^244.70`). Splitting `F_r` into 4 64b limbs gives max convolution `t L_bm 2^128 = 2^188.79 << NTT4_PROD` (`+55.9`b headroom, asserted LIVE). Per `t x t` product (`k = log2 M <= 50`): (1) fwd NTTs (`32 t^2` x `16 M k` = `512 t^2 M k` Ir); (2) pointwise `64 M` scalar `t x t` `F_{p_i}` matmuls (Barrett MAC `17 <= 24` ops -> `1536 M t^3` Ir); (3) inv NTTs + Garner CRT + mod-`r` (`1024 t^2 M k + 6220 t^2 M` Ir). Over `depth = 48` levels (`sum M_j <= 8 L_bm`): `W_rec = 2 * depth * (8 L_bm) * (1536 t^3 + t^2 (1536 k_max + 6220))` Ir (`= 2^92.7199` u after `/20000`). **Total:** `(W_base + W_rec)/20000 = 2^92.7203 <= 2^93.0` u (`< 2^84 B` RAM; Horner in `la_u`).

## 4. Zero-Loss `repurge_stack` & Live Diagnostics

Excess-row clique trim + Markowitz deg-0 guard give `post_peel_cols_lost == post_merge_cols_lost == 0`; 3-stack back-sub (`L_rec` -> `merge_stack` -> `repurge_stack`, `n_reconstruct_unit == 0`) recovers 100% of 2-core cols. Live: `unplanted_two_sided_bw` (unplanted 4-SM NFS `p=467, r=233`, live recursive `pm_basis` + Horner matching `ker(A)`, 38/38 cols, 49/49 rels, exact DLP) and `hx_z12.py` (exact Schur-vs-product `op(x)`, 2,498/2,498 rels, 0 lost cols). Production `H-R, H-X` remain declared exploratory premises.
