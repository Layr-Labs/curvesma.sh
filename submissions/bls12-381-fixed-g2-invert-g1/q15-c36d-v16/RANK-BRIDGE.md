# Finite Relation Rank (`H-R`), Schirokauer Compression, and Block Wiedemann (`H-X`)

Separates finite-instance relation rank (`H-R`), concrete Schirokauer compression (`W_f, W_g`), and Block Wiedemann controllability/observability (`H-X`) for the fixed BLS12-381 `Z[zeta_12]` product-form system (`D = N + 4`, `s_sm = 4`).

## 1. Schirokauer Compression (`W_f, W_g`), Bounded-Retry `F_r` Sampling, & `H-R`

- Concrete data-oblivious Schirokauer columns (`<= 4/r < 2^(-252)`): `P_chi Lambda_f` and `P_chi Lambda_g` have ambient `F_r`-dimensions `6` and `2`, containing 3D and 1D `chi`-unit images. Before collection, draw uniform `3x6` `W_f` and `1x2` `W_g` over `F_r`; the `4` Schirokauer columns of `S` are `W_f (P_chi Lambda_f)` and `W_g (P_chi Lambda_g)`. Under `H-SM`, injectivity loss is `<= 4/r < 2^(-252)` (summed in `proof.md §6`, verified live in `toy_e2e_dlp`).
- Bounded-retry uniform `F_r` sampling (`< 2^(-192) << 2^(-100)`): Every uniform `F_r` entry of `W_f, W_g, X, Y` draws uniform 255-bit words rejecting `>= r` (`r > 2^254.85`, reject rate `< 0.10`), capped at `<= 64` draws/entry: failure `< 0.10^64 < 2^(-212)` per entry, `< 2^(-192)` over `< 2^20` entries (absorbed in `2^(-48)` wrapper).
- `H-R` (declared fixed-instance premise): Given `(W_f, W_g)`, the collected, `tau`-quotiented, filtered matrix `A = [M | S]` (`D = N + 4`) is fixed; `H-R` assumes `rank(M) = N` and `rank(A) = D - 1`, so `ker(A)` is the 1D virtual-log line `(w, c)` with `c != 0` (not proved by relation counts or `r !| h`).

## 2. `H-X` (Block Controllability & Observability) & Conditional Schwartz-Zippel

Let `B = [M | 0]` (`D x D`), `V = Im(B) + span(S)` (`dim V = D - 1` under `H-R`), and `Y = [S_1..S_4, B Z]` (`n - 4 = 106` uniform columns `Z`). Publicly partition the `4` pad coordinates into `4` distinct projection blocks of `X` and the `N` non-pad coordinates into `m = 8192` balanced blocks, sampling one uniform `>= 2^128`-support weight per coordinate. `H-X` asserts at block lengths `d = ceil(D/n) + 1 + ceil(D/(128 n))` and `e = ceil(D/m) + 1 + ceil(D/(128 m))` that an `N x N` minor of `C_d = [B Z, B Y, ..., B^d Y]` (deg `<= N` in `Z`) and a `(D-1) x (D-1)` minor of `O_e U` (deg `<= D-1` in `X`) are nonzero polynomials. Conditional Schwartz-Zippel (proved given `H-R, H-X`): failure `<= 4 N_cap / 2^128 < 2^(-66)`.

## 3. Generator Recovery, Zero-Loss `repurge_stack`, & Live Diagnostics

- Generator, Horner, back-substitution: the order basis of proof.md section 4 (PM-basis, <= 2^98.584 units; `bm_toy.py` checks it against the iterative M-basis) yields `f(z)` with `c = f_0[:4] != 0`, and Horner recovers `(w_bar[:N], c) in ker(A)`. Excess-row clique trim and a Markowitz degree-0 guard keep `post_peel_cols_lost == post_merge_cols_lost == 0`; the 3-stack triangular back-substitution (`L_rec`, `merge_stack`, `repurge_stack`) recovers every 2-core column.
- Live diagnostics: (1) `unplanted_two_sided_bw` in `audit_cm.py` runs unplanted 4-SM NFS (`f = x^4-5x^2+x+1, g = x^2-x-10, p = 467, r = 233`), verifying `rank(M)=23, rank(A)=26, ker_dim=1`, 3-SM and 1-block negative controls, Hankel + Horner + 2-core `reconstructlog` (`38/38` cols, `49/49` rels, exact `F_467^*` DLP); (2) `hx_z12.py` (`log2U=12, seed=5` inline + `12..14` on demand) verifies exact Schur-vs-product operator `op(x)`, 4-SM BW + Horner, and zero-loss 3-stack `repurge_stack` (`2,498/2,498` rels, `0` lost cols). Production `H-R` and `H-X` remain declared exploratory premises.
