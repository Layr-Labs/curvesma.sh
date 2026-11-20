# Finite Relation Rank (`H-R`), Schirokauer Compression, Block Wiedemann (`H-X`), & Self-Contained Classic `PM-Basis` Generator

Separates finite-instance relation rank (`H-R`), concrete Schirokauer compression (`W_f, W_g`), Block Wiedemann controllability/observability (`H-X`), and a **100% self-contained classic Giorgi-Jeannerod-Villard (2003) / Thomé (2002) `PM-Basis`** matrix-polynomial generator with explicit finite baseline scalar ARMv8-A instruction bounds (`<= 2^93.0` core units) and live toy `PM-Basis` execution. No external paper black box is relied upon.

## 1. Schirokauer Compression (`W_f, W_g`), Exact Bounded-Retry `F_r` Sampling, & `H-R`

- **Concrete data-oblivious Schirokauer columns (`<= 4/r < 2^(-252)`):** `P_chi Lambda_f` and `P_chi Lambda_g` have ambient `F_r`-dimensions `6` and `2`, containing 3D and 1D `chi`-unit images. Before collection, draw uniform `3x6` `W_f` and `1x2` `W_g` over `F_r`; the `4` Schirokauer columns of `S` are `W_f (P_chi Lambda_f)` and `W_g (P_chi Lambda_g)`. Under `H-SM`, injectivity loss is `<= 4/r < 2^(-252)` (summed in `proof.md §6`, verified live in `toy_e2e_dlp`).
- **Exact 255-bit / `K_rej = 64` uniform `F_r` sampling (`< 2^(-156) << 2^(-150)`):** Every uniform `F_r` entry of `W_f` (18), `W_g` (2), and `Z` (`(n - 4)(N + 4) = 106 (N + 4)`) draws uniform 255-bit words rejecting `>= r` (`p_rej = 1 - r/2^255 = 0.0943098873`), capped at `K_rej = 64` draws/entry: per-entry miss `p_rej^64 = 2^(-218.0126)`. Over `20 + 106 (ceil(N) + 4) < 2^61.25` entries (`N <= 2^54.5215`), the union bound is `< 2^(-156.7) < 2^(-150)` (verified LIVE in 80-digit `Decimal` in `audit_cm.py` and absorbed in `2^(-48)` wrapper). Block-projection weights `X` sample direct uniform 128-bit integers in `[0, 2^128)` with zero rejection.
- **`H-R` (declared fixed-instance premise):** Given `(W_f, W_g)`, the collected, `tau`-quotiented, filtered matrix `A = [M | S]` (`D = N + 4`) is fixed; `H-R` assumes `rank(M) = N` and `rank(A) = D - 1`, so `ker(A)` is the 1D virtual-log line `(w, c)` with `c != 0` (not proved by relation counts or `r !| h`).

## 2. `H-X` (Block Controllability & Observability) & Conditional Schwartz-Zippel

Let `B = [M | 0]` (`D x D`), `V = Im(B) + span(S)` (`dim V = D - 1` under `H-R`), and `Y = [S_1..S_4, B Z]` (`n - 4 = 106` uniform columns `Z`, `n = 110`). Publicly partition the `4` pad coordinates into `4` distinct projection blocks of `X` and the `N` non-pad coordinates into `m = 8192` balanced blocks, sampling one uniform `[0, 2^128)` weight per coordinate. **`H-X`** asserts at block lengths `d = ceil(D/n) + 1 + ceil(D/(128 n))` and `e = ceil(D/m) + 1 + ceil(D/(128 m))` that an `N x N` minor of `C_d = [B Z, B Y, ..., B^d Y]` (deg `<= N` in `Z`) and a `(D-1) x (D-1)` minor of `O_e U` (deg `<= D-1` in `X`) are nonzero polynomials. **Conditional Schwartz-Zippel (proved given `H-R, H-X`):** failure `<= 4 N_cap / 2^128 < 2^(-66)`.

## 3. Self-Contained Classic `PM-Basis` Matrix Generator (`<= 2^92.7203 <= 2^93.0` Core Units)

Given `S_k = X B^k Y in F_r^{m x n}` for `0 <= k < L_bm` where `L_bm = (1 + 2^-7)(N + 4)(1/m + 1/n) + 6 <= 2^47.7706`, we seek `q(z) in F_r[z]^n` (`deg q <= d`) and `r(z) in F_r[z]^m` (`deg r < d`) satisfying `[S(z) | -I_m] [q(z); r(z)] == 0 (mod z^L_bm)` (`t = m + n = 8302` columns, `m = 8192` rows, initial column shift `s = (0^n, 1^m)`). Classic recursive `PM-Basis(F, sigma, s)` (Giorgi-Jeannerod-Villard 2003 / Thomé 2002; executed verbatim in `pm_adaptive.py`):
- **Base case (`sigma = 1`, `M-Basis-1`):** Column-echelon reduction of scalar `m x t` matrix `F(0)` over `F_r` in non-decreasing shift order (`<= m * t^2` `F_r` MACs at `MONT5 = 270` Ir each), multiplying pivot columns by `z` and incrementing their shifts by 1. Across all `L_bm` leaves:
  `W_base = L_bm * m * t^2 * 270 = 2^80.5982` core units.
- **Recursive step (`sigma > 1`, `s1 = floor(sigma/2)`, `s2 = sigma - s1`):**
  1. `(P1, s_mid) = PM-Basis(F mod z^s1, s1, s)` (`t x t`, `deg P1 <= s1`);
  2. `FP1 = F * P1 mod z^sigma` (`m x t` by `t x t`), `G = (FP1 div z^s1) mod z^s2`;
  3. `(P2, s_out) = PM-Basis(G, s2, s_mid)` (`t x t`, `deg P2 <= s2`);
  4. `P = P1 * P2` (`t x t`, `deg P <= sigma`).
- **Explicit 1D 4-Prime 62-Bit NTT Polynomial-Matrix Multiply Cost (Schoolbook `t^3`, ZERO Strassen/Karatsuba):**
  By induction `deg(F * P1) < 2*sigma` and `deg(P1 * P2) <= sigma`, so NTT length `M_sigma = 2^ceil(log2(2*sigma)) <= 4*sigma` suffices at every level. At the top level `M_top <= 4*L_bm = 2^49.7706 < 2^53`, so EVERY level uses **1D radix-2 NTT** across our 4 verified primes `p_i = c_i * 2^53 + 1 < 2^62` (`NTT4_PROD = 2^244.70`). Splitting each `255`-bit `F_r` entry into `4` base-`2^64` limbs gives exact integer convolution ceiling `t * L_bm * (2^64)^2 = 2^188.79 << NTT4_PROD` (`+55.9` bits capacity headroom, asserted LIVE in `audit_cm.py`). For one `t x t` product of NTT length `M` (`k = log2 M <= 50`):
  1. **Forward NTTs** (`2 * t^2` polys x `4` limbs x `4` primes = `32 t^2` length-`M` NTTs, `9` ops/elem/stage in `ntt62_butterfly` charged `16` with loop/twiddle): `512 * t^2 * M * k` Ir.
  2. **Pointwise `t x t` matmuls** (`M` points x `4` primes x `16` limb-pairs = `64 M` scalar `t x t` matmuls over `F_{p_i}`, each MAC `17 <= 24` scalar ARMv8-A Barrett ops verified bit-exactly in `kernels_symbolic`): `1536 * M * t^3` Ir.
  3. **Inverse NTTs + 4-prime Garner CRT + 7-limb mod-`r`** (`64 t^2` inv NTTs at `16 M k` + `16 t^2 M * 220` Garner + `t^2 M * 10 * 270` mod-`r`): `1024 * t^2 * M * k + 6220 * t^2 * M` Ir.
  Summing over `depth = ceil(log2 L_bm) = 48` levels (`2` matmuls/node, `sum M_j <= 8*L_bm` per level):
  `W_rec = 2 * depth * (8 * L_bm) * (1536 * t^3 + t^2 * (1536 * k_max + 6220)) = 2^92.7199` core units.
  **Total:** `(W_base + W_rec) / 20,000 = 2^92.7203 <= 2^93.0` core units (`18.5` bits below total opened-run work!). Top-level 4-prime 4-limb NTT workspace is `4 * 4 * t^2 * (4 L_bm) * 8 = 2^82.81 B < 2^84 B`. Horner evaluation is priced inside `linear_algebra` (`la_u`).

## 4. Zero-Loss `repurge_stack` & Live Diagnostics

Excess-row clique trim + Markowitz deg-0 guard keep `post_peel_cols_lost == post_merge_cols_lost == 0`; 3-stack triangular back-sub (`L_rec` -> `merge_stack` -> `repurge_stack`, `n_reconstruct_unit == 0`) recovers 100% of 2-core columns.
- **`unplanted_two_sided_bw` (live in `audit_cm.py`):** Unplanted 4-SM NFS (`f = x^4 - 5x^2 + x + 1`, `g = x^2 - x - 10`, `p = 467`, `r = 233`), verifying `rank(M)=23, rank(A)=26, ker_dim=1`, 3-SM and 1-block negative controls, **LIVE recursive `pm_basis` (`order=15, cols=8, deg=6`) + Horner matching the exact true kernel `ker(A)`**, and 2-core `reconstructlog` (`38/38` cols, `49/49` rels, exact `F_467^*` DLP).
- **`hx_z12.py` (live inline `12/5` + `12..14` on demand):** Verifies exact Schur-vs-product operator `op(x)`, 4-SM BW + Horner, and zero-loss 3-stack `repurge_stack` (`2,498/2,498` rels, `0` lost cols). Production `H-R` and `H-X` remain declared exploratory premises.
