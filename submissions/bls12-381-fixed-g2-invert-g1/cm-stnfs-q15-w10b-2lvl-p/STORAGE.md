# Storage and Addressing

Every object is partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner) with every cross-owner byte transfer priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Prod-ZeroCopy-U32`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 448 = 7 * 64`, `W_SB = 32`, `P_m4s = 2^55.0`, `P_m4 = 2^56.5`, `P_lo = P_f = 2^59`, `P_g = 2^64`)

- **Canonical `a` packing & exact `2^57` power-of-2 bank stride:** `N_a^can = V_a = 2^56.9093 < 0.94 * 2^57` canonical `a`-orbits (`6.1%` headroom below `2^57`); `N_a^tbl = 4 V_a = 2^58.9093`. Packing `a_slot = (m << 57) | oid` matches column-owner accumulators `uint16_t S_g[4 * 2^57]` (`2^60` B) + `S_f[4 * 2^57]` (`2^60` B) = **`2^61.0` B `<< 2^64` B** per column owner (`448` column owners).
- **1-Owner Simultaneous Flat-Array Build (`<= 2^63.956 B < 2^64 B` for ALL primes, ZERO cross-owner build transfers):**
  - Band-1a (`p <= 2^56.5`): `16 B * N_a^tbl + 64 B * 2^56.5 + 448 MiB <= 2^63.694 B < 2^64 B`.
  - Band-1b (`p <= 2^59`) & Band-2 (`g`-side `p <= 2^64`, `W_SB = 32` -> `<= 2^59` SBs): `16 B * N_a^tbl + 16 B * 2^59 + 448 MiB <= 2^63.956 B < 2^64 B`.
  - Total replay `a`-table storage: **`2^119.163` B `< 2^120` B**, built once and freed immediately after collection.
- **Full-Chord Rows (`L_max = 391.70 < W_tile = 448 = 7 * 64`):** `N_7 = 2^53.8149` rows in `D_(omega,sigma)`, `L_b = 168.24`, `N_a^tbl * N_7 * L_b = V_kept = 2^120.1186` exact.
- **Product Trees (`P_f = 2^71.29`, `P_g = 2^71.69`, `P_71 = 2^71`):** `P_g` (`74.41` u `<= 75.5` u) spans `256` 64-bit owners of **`2^61.22 B << 2^64 B`** each (`P_71` spans `128` owners of `2^61.53 B < 2^62 B` each).

## 2. Linear Algebra Storage — 1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (All `< 2^64 B`, Priced Streaming)

Each Block Wiedemann bundle runs `3` 256-bit vectors (`96 B/slot`) over `N_sw = sw_cap * U' = 0.00500 * 2^63.3174 = 2^55.6735` structured-Gaussian slots (`1.463x` above the measured zero-column-loss `2^25` mean `0.003417`):

- **Dense-vector owner (`1` owner):** Double-buffered 3-vector state **`2 * N_sw * 96 B = 2^63.2585 B << 2^64 B`** (`11.03 EiB < 18.45 EiB`); Horner evaluation state **`2 * N_sw * 32 B = 2^61.6735 B << 2^64 B`** (asserted LIVE in `audit_cm.py`).
- **Read-only sparse-matrix shard owners (`n_mat_owners = 16` owners):** The read-only structured sparse matrix (`201` nonzeros/row x `N_sw` rows x `8 B/entry = 2^66.322 B`) is partitioned across **`16` read-only 64-bit shard owners**, each holding **`2^62.322 B << 2^64 B`** (asserted LIVE in `audit_cm.py` as `mat_shard_mem_log2 = 62.322 < 64.0`).
- **Explicit priced cross-owner streaming transfer (`c_stream_mat = 268.0` Ir/row):** During each 3-vector SpMV of all `40` independent bundles (`n = 120 = 40 * 3`, `m = 16,384`, `tab_ratio_charged = 1.00`, ZERO `TABM`), the `16` matrix-shard owners stream each 8-byte matrix entry sequentially into the dense-vector owner (`4.0` Ir per 8-byte streamed entry amortized across `3` vectors = `1.3333` Ir per vector-nonzero -> **`201 * 4 / 3 = 268.0` Ir/row**), explicitly charged inside `c_prod` (`audit_cm.py` live check!).

## 3. Peak Memory & 2-Adic Class Rule

- **Peak Memory:** `< 2^120.0 B` (`2^119.163 B` during collection; `< 2^104.0 B` during descent; known-log `P_f` pool table for `Pure M1` is `1,344 * 24 B = 32 KiB`).
- **2-Adic Class Rule `c(a)` (`w_cls = (225/256, 225/2048, 23/2048)`):** Identical `<sigma, omega>`-invariant rule (`s = (1.0, 1.5, 2.0)` -> `v_2(N_g) = 0, 12, >= 24`), verified live in `sigmaB_exactness` and `pool_2adic_live`.
