# Storage and Addressing

Every object is partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner) with every cross-owner byte transfer priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Prod-ZeroCopy-U16`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 400`, `W_SB = 32`, `P_m4 = 2^56.5`, `P_lo = P_f = 2^59`, `P_g = 2^64`)

- **Canonical `a` packing & exact `2^57` power-of-2 bank stride (ZERO `q_key`, ZERO `x_meta`, ZERO `x_imap`, ZERO global sort):** `N_a^can = V_a = 2^56.7395 < 0.84 * 2^57` canonical `a`-orbits (`16.5%` headroom below `2^57` for finite-lattice boundary points); `N_a^tbl = 4 V_a = 2^58.7395`. Packing `a_slot = (m << 57) | oid` (`59` active bits, exact bank stride `2^57`) matches column-owner accumulators `uint16_t S_g[4 * 2^57]` (`2^60` B) + `S_f[4 * 2^57]` (`2^60` B) = **`2^61.0` B `<< 2^64` B** per column owner (`400` column owners, ZERO stride mismatch!).
- **1-Owner Simultaneous Flat-Array Build (`< 2^63.889 B < 2^64 B` for ALL primes, ZERO cross-owner build transfers, ZERO global sorting invariant):**
  - *Build algorithm (`rev_build_prod_pass12` + `rev_build_prod_prefix`, `KERNELS.md` §1.7):* Because `Arch2-Prod` has NO `q_key` ordering requirement, ONE single 64-bit owner builds each prime's complete replay table in two linear passes in deterministic `oid_m = (m << 57) | oid` counter order (NO `x_meta`, NO `x_imap`, NO multi-owner merge!). Pass 1 writes `x_kbuf` (`8 B * N_a^tbl = 2^61.739 B`) and increments **64-bit `uint64_t` `x_hist`** (`0` overflow at all `p >= 2^24`). Prefix sum writes both immutable `x_off` (`8 B/bkt`) and mutable `x_pos` (`8 B/bkt`). Pass 2 scatters 8-byte records into `x_out` (`8 B * N_a^tbl = 2^61.739 B`).
  - *Simultaneous 1-owner peak RAM (asserted LIVE in `audit_cm.py`):*
    - Band-1a (`p <= 2^56.5`, 4 `m`-partitioned banks): `16 B * N_a^tbl + 64 B * 2^56.5 + 400 MiB =` **`<= 2^63.625` B `< 2^64` B**!
    - Band-1b (`p <= 2^59`, `f`- and `g`-side): `16 B * N_a^tbl + 16 B * 2^59 + 400 MiB =` **`<= 2^63.889` B `< 2^64` B**!
    - Band-2 (`g`-side only, `p <= 2^64`, `W_SB = 32` -> `<= 2^59` SBs): `16 B * N_a^tbl + 16 B * (2^64 / 32) + 400 MiB =` **`<= 2^63.889` B `< 2^64` B**!
  - *Total replay `a`-table storage (`g` to `2^64` + `f` to `2^59`):* **`2^119.074` bytes** (`< 2^120.0` B), built once per target and freed immediately after collection.
- **Full-Chord Rows & `400` Column Owners (`L_max = 393.70 < W_tile = 400`):** Every `b_0`-row (`N_7 = 2^53.8663` rows in `D_(omega,sigma)`, `L_b = 169.10`, `N_a^tbl * N_7 * L_b = V_kept = 2^120.0075` exact!) fits in one full-chord tile (`1` tile/row, `6.784` SBs/row). Each of the `400` column owners holds `S_g` + `S_f` = **`2^61.0` B `<< 2^64` B**.
- **Stage B, SRT, & Gate G `P_70` product trees:** `P24` product (`3.1` MiB), waves `< 2^27` B; `P_g` (`2^62.53` B), `P_f` (`2^57.53` B), and Gate G `P_70` (`2^68.53` B, partitioned across `64` 64-bit owners of `2^62.53 B < 2^64 B` each, 2D NTT transpose transfers covered inside `fft_F`'s `168 M(k1+k2) + 1100 M`) total `< 2^76` B.

### 2. Linear Algebra Storage — Single-Owner `< 2^64 B` 4-Vector Arch B (`n = 112`, `m = 8192`, `te = 0`)

Each Block Wiedemann bundle runs `4` 256-bit vectors (`128 B/slot`) over `N_sw = sw_cap * U' = 0.00505 * 2^63.3174 = 2^55.6879` structured-Gaussian slots (`>= 1.40x` the measured zero-column-loss `2^24` value `0.00360`). The entire double-buffered dense state fits inside **one single 64-bit owner** (`n_owners = 1`, `te = 0`):

- Double-buffered 4-vector state: `2 * N_sw * 128 B = 2^63.8876 B < 2^64 B` (`16.99 EiB < 18.45 EiB`).
- Horner evaluation state (`2 * N_sw * 32 B = 2^61.8876 B < 2^64 B`) also fits comfortably inside the same single owner.
- All `28` independent 4-vector bundles (`n = 112 = 28 * 4`, `m = 8192`, `tab_ratio_charged = 1.00`) stream the read-only sparse matrix with zero cross-owner communication (`c_xown_live = 0.0`).

## 3. Descent, Wrapper, & Aggregate Peak Memory

- **Descent:** Streaming single active CRT special-`q` node `<= 4*|S_16(92)| = 2^103.298 B < 2^104.0 B` (`|S_16(92)| = 2^101.298`, `max_R = 92`, partitioned across `<= 2^42` chunk owners of `<= 2^62 B < 2^64 B` each, asserted live; `0.50 GiB` flat read-only `P_d` `rank_lex` table; `P_70` tree `2^68.53 B` across 64 owners; `192`-ideal blinding pool `B_blind` `< 16 KiB`; `boot < 2^30 B`; memo `<= 3,357` ideals `< 2^22 B`). **Wrapper:** `< 2^30 B`.
- **Peak Memory:** `< 2^120.0 B` (`2^119.074 B` replay `a`-tables during collection; `< 2^104.0 B` during descent).

## Exact Production 2-Adic Class Rule `c(a)` (`w_cls = (225/256, 225/2048, 23/2048)`)

Over `F_2`, `Phi_15(y) = P_1(y) P_2(y) = (y^4+y+1)(y^4+y^3+1)` splits `(2) = P_1 P_2` in `O_H` (`N(P_1)=N(P_2)=16`). Since `v_2(u)=16`, `g(a,b) = a^3 - u w^2 b^3 == a^3 (mod 2^16)`, so `v_2(N_g) = 3 v_2(N(a))` depends solely on `a mod 4` and is `<sigma, omega>`-invariant: `c(a)=0` iff `v_2(N(a))=0` (`P_1, P_2 !| a`, density `225/256`, `s_0=1.0`, `v_2(N_g)=0`); `c(a)=1` iff `v_2(N(a))=4` (`v_{P_1}+v_{P_2}=1`, density `225/2048`, `s_1=1.5`, `v_2(N_g)=12`); `c(a)=2` iff `v_2(N(a))>=8` (density `23/2048`, `s_2=2.0`, `v_2(N_g)>=24`). Pure-product region: `q_a <= s_{c(a)} R2_a`, `q_b <= R2_b`, `a in F_Vor / mu_30`, `b in D_(omega,sigma)`. Verified LIVE in `check_sigmaB_exactness` (`10,878` hits, `0` mismatches) and `pool_2adic_live` (`400` live draws verifying exact `v_2(N_g)`, exact `O_H`-coprimality `383/400`, `z_prim = 0.13`, `prim_yield_ratio = 1.0049`, `nf_mean`, `ng_mean`, class shares, `log2_Va`, `log2_Vb_all`, and `V_kept`).
