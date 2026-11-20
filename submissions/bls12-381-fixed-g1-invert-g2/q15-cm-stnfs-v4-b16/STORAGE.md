# Storage and Addressing

Every object is partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner) with every cross-owner byte transfer priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, `W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

- **Canonical `a` packing & `key_a`-sorted 4-bank layout:** `N_a^can = 2^57.8300` canonical `a`-orbits sorted by unified 2-adic key `key_a(oid) = round(2^16 (q_a(a) - s_c(a)))`. Packing `a_slot = (m << 58) | oid` (`N_a^tbl = 4 N_a^can = 2^59.8300 < 2^60` slots) guarantees that for any column `t`, active cells `key_a(oid) <= -q_b(t)` are four exact contiguous prefixes `[0, N_act_orb(t))`, so `rev_scan_zero_u16` scans and zeroes ONLY active prefixes (`V_all = 1.0276 V_kept`). Global `perm[]` pre-sort arrays are `< 2^63` B each.
- **4-orbit `sigma`-shared `a`-tables (1 owner/table, 64-bit `uint64_t` build histograms, every table `< 2^64` B):**
  - *Build histograms (`rev_build_pass12`, `KERNELS.md` §1):* `x_hist` uses **64-bit `uint64_t` counters** (`ldr/str x2, [x_hist, x1, lsl #3]`, preventing 32-bit overflow near `p = 2^24` where `N_a^tbl / p > 2^32`) and is converted in-place to `uint64_t` `x_pos_sb` prefix offsets (`8` B/bucket, already accounted in both Band-1 and Band-2 peak table budgets below).
  - *Band 1a (`[2^24, 2^48)`, `m`-partitioned) & Band 1b (`[2^48, 2^55)`, `rev_b1_early4_find`):* `12` B/entry SoA (`keys[u32]` + `slots[u64]`) + `u64` heads/hist (`8 * 2^55` B max) + `552` MiB staging = **`<= 2^63.4484` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Band 2 (`[2^55, P_hi)`, `P_g = 2^63`, `P_f = 2^60`):* `W_SB = 64`, `16`-byte `sb_entry_t` + `u64` heads/hist & 1-byte `sb_min_q8` (`9` B/SB) + `552` MiB staging = **`<= 2^63.9399` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Total replay `a`-table storage (`g` to `2^63` + `f` to `2^60`):* **`2^119.655` bytes** (`< 2^120.0` B), built once per target and freed immediately after collection.
- **Full-Chord Rows & `552` Column Owners (`L_max = 550.47 < W_tile = 552`):** Every `b_0`-row (`N_7 = 2^57.2513` rows, `L_b = 236.44`) fits in one full-chord tile (`nt = 1.000`, `4.676` SBs/row). Each of the `552` column owners holds `S_g[t]` + `S_f[t]` = **`2^61.830` B `<< 2^64` B**.
- **Stage B, SRT, & Gate G `P_70` product trees:** `P24` product (`3.1` MiB), waves `< 2^27` B; `P_g` (`2^61.53` B), `P_f` (`2^58.53` B), and Gate G `P_70` (`2^68.53` B, partitioned across `64` 64-bit owners of `2^62.53 B < 2^64 B` each, 2D NTT transpose transfers covered inside `fft_F`'s `168 M(k1+k2) + 1100 M`) total `< 2^76` B.

## 2. Linear Algebra (`1-Owner Architecture B`, `< 2^64` Bytes Total Including Output Block)

At hard caps `N_l <= 2^54.5215`, `N_sw <= 2^55.5215`, `(m, n) = (8192, 110)` (`22` groups of 5 vectors):
- **Input `x` + Output `y` (`160` B/slot each):** `2 x N_sw x 160 = 2^63.843` B.
- **Dense-band subset-sum tables (`{2: 28.0, 3: 18.25, 4: 14.0}`):** `2^47.715` elements (`2^59.497` B across all 22 groups).
- **1-owner fit (verified live in `audit_cm.py`):** `2 x N_sw x 160 + 2^59.497 =` **`2^63.913` B `< 2^64` B** (`0` cross-owner sends, combines, or routing). Stored matrix `< 2^69.5` B, Horner checkpoints `< 2^70.2` B, BM workspace `< 2^84` B.

## 3. Descent, Wrapper, & Aggregate Peak Memory

- **Descent:** Streaming single active special-`q` node `<= 4*|S_16(128)| = 2^110.913 B < 2^111.0 B` (`|S_16(128)| = 2^108.913`, partitioned across `<= 2^49` chunk owners of `<= 2^62 B < 2^64 B` each, asserted live; `0.88 GiB` read-only `P_d` `rank_lex` table; `P_70` tree `2^68.53 B` across 64 owners; `boot < 2^30 B`; memo `<= 642` ideals `< 2^20 B`). **Wrapper:** `< 2^30 B`.
- **Peak Memory:** `< 2^120.0 B` (`2^119.655 B` replay `a`-tables during collection; `< 2^111.0 B` during descent).

## Exact Production 2-Adic Class Rule `c(a)` (`w_cls = (225/256, 225/2048, 23/2048)`)

Over `F_2`, `Phi_15(y) = P_1(y) P_2(y) = (y^4+y+1)(y^4+y^3+1)` splits `(2) = P_1 P_2` in `O_H` (`N(P_1)=N(P_2)=16`). Since `v_2(u)=16`, `g(a,b) = a^3 - u w^2 b^3 == a^3 (mod 2^16)`, so `v_2(N_g) = 3 v_2(N(a))` depends solely on `a mod 4` and is `<sigma, omega>`-invariant: `c(a)=0` iff `v_2(N(a))=0` (`P_1, P_2 !| a`, density `225/256`, `s_0=1.00`, `v_2(N_g)=0`); `c(a)=1` iff `v_2(N(a))=4` (`v_{P_1}+v_{P_2}=1`, density `225/2048`, `s_1=1.20`, `v_2(N_g)=12`); `c(a)=2` iff `v_2(N(a))>=8` (density `23/2048`, `s_2=1.36`, `v_2(N_g)>=24`). Region: `q_a <= 0.65 s_{c(a)}`, `q_a + q_b <= s_{c(a)}`, `a in F_Vor / mu_30`, `b in D_(omega,sigma)`. Verified LIVE in `check_sigmaB_exactness` and `pool_2adic_live` (`400` live draws verifying exact `v_2(N_g)`, exact `O_H`-coprimality `382/400`, `z_prim = 0.42`, `prim_yield_ratio = 1.0066`, `nf_mean`, `ng_mean`, class shares, `log2_Va`, `log2_Vb_all`, and `V_kept`).
