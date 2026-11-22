# Storage and Addressing

Every object is partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner) with every cross-owner byte transfer priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Prod-ZeroSlice-1PassB2-U16`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 448 = 7 * 64`, `W_SB = 16`, `P_min = 2^32`, `P_m4s = 2^53.0`, `P_m4 = 2^55.5`, `P_lo = P_f = 2^58.0`, `P_g = 2^63.5`)

- **Canonical `a` packing & exact `2^56` power-of-2 bank stride (`K_a = 1`, `W_SB = 16`, `rec = 62b`):** `N_a^can = V_a = 2^55.8474 < 0.90 * 2^56` canonical `a`-orbits (`10.0%` headroom below `2^56`); `N_a^tbl = 4 V_a = 2^57.8474`. Band-1b packs `a_slot = (m << 56) | oid` (`58`b) and Band-2 packs `rec = (m << 60) | (r << 56) | oid` (`62`b), matching column-owner slab accumulators `uint16_t S_g[4][16][2^56]` (`2^63.0` B, 1-pass `W_SB = 16` direct apply, ZERO emit buffers!) + 1-column `uint16_t S_f[4][1][2^56]` (`2^59.0` B) = **`2^63.087` B `< 2^64` B** per column-slab owner (`448 / 16 = 28` slab owners).
- **1-Owner Simultaneous Flat-Array Build (`<= 2^63.898 B < 2^64 B` for ALL primes, ZERO cross-owner build transfers):**
  - Band-1a (`p <= 2^55.5`): `16 B * N_a^tbl + 64 B * 2^55.5 + 448 MiB <= 2^62.58 B << 2^64 B`.
  - Band-1b (`p <= 2^58.0`): `16 B * N_a^tbl + 16 B * 2^58.0 + 448 MiB <= 2^62.89 B << 2^64 B`.
  - Band-2 (`g`-side `p <= 2^63.5`, `W_SB = 16` -> `<= 2^59.5` SBs): `16 B * N_a^tbl + 16 B * 2^59.5 + 448 MiB <= 2^63.898 B < 2^64 B`.
  - Total replay `a`-table storage: **`2^118.317` B `< 2^120` B**, built once and freed immediately after collection.
- **Full-Chord Rows (`L_max = 439.06 < W_tile = 448 = 7 * 64`):** `N_7 = 2^54.9675` rows in `D_(omega,sigma)`, `L_b = 188.58`, `N_a^tbl * N_7 * L_b = V_kept = 2^120.3739` exact.
- **Product Trees (`P_f = 2^71.29`, `P_g = 2^71.69`, `KH = 0.08`, 3-prime `ntt3_F`):** `P_g` (`<= 75.5` u) spans `256` 64-bit owners of **`<< 2^64 B`** each.

## 2. Linear Algebra Storage — 1 Dense-Vector Owner + 16 Read-Only Matrix-Shard Owners (All `< 2^64 B`, Priced Streaming)

Each Block Wiedemann bundle runs `4` 256-bit vectors (`128 B/slot`) over `N_sw = sw_cap * U' = 0.00498 * 2^63.3174 = 2^55.6677` structured-Gaussian slots (`1.458x` above the measured zero-column-loss `2^25` mean `0.003417`):

- **Dense-vector owner (`1` owner):** Double-buffered 4-vector state **`2 * N_sw * 128 B = 2^63.668 B < 2^64 B`** (`12.71 EiB < 16 EiB`, `20.6%` headroom); Horner evaluation state **`2 * N_sw * 32 B = 2^61.668 B << 2^64 B`** (asserted LIVE in `audit_cm.py`).
- **Read-only sparse-matrix shard owners (`n_mat_owners = 16` owners):** The read-only structured sparse matrix stores an 8-byte pre-scaled column offset for every nonzero plus a 1-byte `uint8_t` magnitude (`1 <= |c_b| <= 255`, loaded via `ldrb w22, [x23], #1`) for each non-unit entry (`< 6.5%` measured; conservatively bounded at `<= 10%` -> `8.10` B/entry x `201` entries/row x `N_sw` rows = `2^66.337 B`), partitioned across **`16` read-only 64-bit shard owners**, each holding **`2^62.337 B << 2^64 B`** (asserted LIVE in `audit_cm.py` as `mat_shard_mem_log2 = 62.337 < 64.0`).
- **Explicit priced cross-owner streaming transfer (`c_stream_mat = 204.0` Ir/row):** During each 4-vector SpMV of all `32` independent bundles (`n = 128 = 32 * 4`, `m = 16,384`, `tab_ratio_charged = 1.00`, ZERO `TABM`), the `16` matrix-shard owners stream every 8-byte offset AND every 1-byte non-unit coefficient sequentially into the dense-vector owner (`<= 201 * 8.10 / 8 = 203.51` 8-byte words/row at `4.0` Ir/8B word amortized across `4` vectors = **`203.51 <= 204.0` Ir/row**), explicitly charged and asserted inside `c_prod` (`audit_cm.py` live check!).

## 3. Peak Memory & 2-Adic Class Rule

- **Peak Memory:** `< 2^120.0 B` (`2^118.317 B` during collection; `< 2^95.0 B` during descent — max node radius `R = 64`, `CELLS_16D[64] = 2^92.91`, `2^40`-entry `uint16_t` score arrays (`2 TiB/side = 4 TiB/owner << 2^64 B`), max score `25,354 < 2^16`, `16D` flat ball-ranking prefix table `16 * (64^2 + 1) * 129 * 16 B = 129.03 MiB/worker < 256 MiB/worker`; wrapper `< 2^30 B`).
- **2-Adic Class Rule `c(a)` (`w_cls = (225/256, 225/2048, 23/2048)`):** Identical `<sigma, omega>`-invariant rule (`s = (1.0, 1.5, 2.0)` -> `v_2(N_g) = 0, 12, >= 24`), verified live in `sigmaB_exactness` and `pool_2adic_live`.
