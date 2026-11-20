# Finite Relation Rank (`H-R`) and Padded Block Wiedemann (`H-X`)

This document separates two distinct finite-instance premises for the fixed BLS12-381 `Z[zeta_12]` product-form system (`D = N + 4`, `s_sm = 4` Schirokauer columns). Neither premise is conflated with class-number coprimality or Schirokauer-map injectivity.

## 1. `H-R` (Finite Relation Spanning — Declared Fixed-Instance Premise)

Collection's region, polynomials, and public ECM seed are challenge- and coin-independent, so the collected, `tau`-quotiented, and filtered relation matrix `A = [M | S]` over `F_r` (`D = N + 4` live rows, `N = N_l` live orbit columns, `4` Schirokauer columns) is **fixed**. **`H-R`** assumes for this actual finite matrix that **`rank(M) = N`** and **`rank(A) = D - 1`**, so `ker(A)` is the one-dimensional line spanned by the virtual-log vector `(w, c)` with `c != 0`.

Counting smooth relations, `r !| h(K_f), h(K_g)`, and Schirokauer-map injectivity on `chi`-units do **not** prove `H-R`, nor do synthetic planted-matrix replays establish saturation of the production relation set. Relation excess (`1.09 U'`), complete generation rules, `filtsim_z12` scaling campaigns, and the standard finite-field NFS framework (Schirokauer §2.2 Step 4; Fried-Gaudry-Heninger-Thome §3.1) provide exploratory support. If `H-R` fails, randomized Block Wiedemann cannot supply missing relations; final verification detects failure but provides no probability bound on `H-R`.

## 2. `H-X` (Finite Block Controllability & Observability) and Conditional Schwartz-Zippel

Let `B = [M | 0]` (`D x D`), so `Im(B) = Im(M)` has dimension `N` and `V = Im(B) + span(S)` has dimension `D - 1` under `H-R`. Form `Y = [S_1, ..., S_4, B Z]` with `n - 4 = 106` independent random columns `Z` over `F_r`. Fix the coordinate partition publicly before sampling weights: the `4` pad coordinates occupy `4` distinct projection blocks of `X`, and the `N` non-pad coordinates are partitioned into `m = 8192` balanced blocks; each coordinate gets one independent uniform `>= 2^128`-support entry on its assigned block (`0` elsewhere).

**`H-X` (declared premise):** At the finite block lengths `d = ceil(D/n) + 1 + ceil(D/(128 n))` and `e = ceil(D/m) + 1 + ceil(D/(128 m))` used by the solver, two specific minors are nonzero polynomials:
1. An `N x N` minor of `C_d = [B Z, B Y, B^2 Y, ..., B^d Y]` (degree `<= N` in `Z`), and
2. A `(D - 1) x (D - 1)` minor of `O_e U`, where `O_e` stacks `X^T, X^T B, ..., X^T B^(e-1)` and `U` is a fixed basis of `V` determined solely by `A` (degree `<= D - 1 = N + 3` in `X`).

**Conditional Schwartz-Zippel theorem (proved given `H-R` and `H-X`):** The product of these two nonzero minors has total degree `<= 2 N + 3 < 4 N`. Sampling independent entries uniformly from supports of size `>= 2^128` yields simultaneous control/observation failure **`<= 4 N_cap / 2^128 < 2^(-66)`**, strictly **conditional on `H-R` and `H-X`** (not a probability that `H-R` or `H-X` hold).

## 3. Generator Selection (Non-Unit First Coord) & Horner Recovery

When both minors are nonzero, `C_d` spans `Im(B)` and `O_e` is injective on `V`, so projected and full relations coincide on `V`. Because `S c in Im(B)` for the true `c != 0`, the fast matrix Berlekamp-Massey generator (`KERNELS.md` §3, `<= 2^98.584` units) produces a basis of the projected relation module containing a column `f(z) = sum_(i=0)^d f_i z^i` of degree `<= d` whose first `4` constant coordinates `c = f_0[:4]` are **not all zero** (we do **not** fix `c_0 = 1`, since the true first Schirokauer coordinate may be `0`). Writing `f_0 = (c, e_0)`, Horner evaluation `w_bar = Z e_0 + sum_(i=1)^d B^(i-1) Y f_i` satisfies `M w_bar[:N] + S c = 0`, so `(w_bar[:N], c)` spans `ker(A)` under `H-R`.

## 4. Live Unscreened 1D-Kernel Product-Form Diagnostic (`hx_z12.py`)

`hx_z12.py` (inline at `log2U = 12` in `audit_cm.py`; on-demand `hx_analogue` for `log2U in 12..14`) takes **ALL** `N` live columns and the first `D = N + 4` live rows **unscreened** from live `filtsim_z12`, plants an exact 1D kernel `(w_true, c_true)` into `S_live` (with `c_true[0] = 0, c_true[1] = 1` on even seeds), evaluates the actual 2-step recorded-pivot product operator `op(x) = (M_live - E_live L_rec^(-1) B_rec) x`, solves the projected Hankel kernel `H f = 0` (`H = O_e [Y, ..., B^d Y]`), and verifies Horner recovery:
- `log2U = 12, seed = 1` (`c_0 = 1`): `N = 38, D = 42, rM = 38, rA = rOeA = rH = 41, rCd = 38, 1-block ctrl = 39 < 41`, `A(w_rec, c_rec) = 0`, exact proportionality `true`.
- `log2U = 12, seed = 2` (`c_0 = 0`): `N = 51, D = 55, rM = 51, rA = rOeA = rH = 54, rCd = 51, 1-block ctrl = 52 < 54`, exact recovery `true`.
- `log2U = 13..14` (`seeds 1, 2`, up to `N = 154, D = 158`): `6/6` pass with exact Horner kernel recovery and negative control failure.

These are algebra/operator diagnostics on constrained synthetic systems; `H-R` and production `H-X` remain declared structural premises under exploratory review.
