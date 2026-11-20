# Storage and Addressing

Every object is explicitly bounded and partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner). No single object requires `> 64`-bit indexing, and every cross-owner byte transfer is explicitly priced in `KERNELS.md`.

## 1. Relation Collection (Tiled `Arch2-Rev-SigmaB`, `P_lo = 2^54`, `W_tile = 256`, `W_SB = 128`)

- **Canonical `a` packing & `key_a`-sorted 4-bank layout:** There are `N_a^can = 2^57.6666` canonical `a`-orbits (`oid in [0, N_a^can)`), globally sorted in ascending unified 2-adic key `key_a(oid) = round(2^16 (q_a(a) - s_c(a)))`. Because `sigma` is a `G`-isometry preserving `q_a` and 2-adic class `s_c`, all 4 `sigma`-conjugates `(oid, m)` (`m in {0, 1, 2, 3}`) share the exact same `key_a(oid)`. Packing `a_slot = (m << 58) | oid` (`N_a^tbl = 4 N_a^can = 2^59.6666 < 2^60` slots; top 2 bits = `m`, low 58 bits = `oid < 2^57.69`) gives two exact guarantees:
  1. Rotating conjugate index `m -> (m - k) mod 4` is 2 baseline AArch64 instructions (`sub` + `and #0xfffffffffffffff`).
  2. For **any** `b`-column `t` with quadratic value `q_b(t)`, the active cells `{(oid, m) : key_a(oid) <= -q_b(t)}` are **four exact contiguous prefixes `[0, N_act_orb(t))`** (one per bank `m in 0..3`). Thus `rev_scan_zero_u4` scans and zeroes **ONLY** the active prefix cells (`V_all = 1.0276 V_kept`, covering the `2.69%` `p2skip` cells too with zero stale entries and **zero** bounding-box overhead).
- **Global pre-sort (`perm[]`):** Built once before table construction (`perm[rank]`, `a_slot_sorted[rank]`, `key_sorted[rank]`: `2^62.69` B each, `< 2^64` B per owner).
- **4-orbit `sigma`-shared `a`-tables (1 owner per 4-orbit table, every table `< 2^64` B):**
  - *Band 1 (`[2^24, 2^54)`):* Unified `key_a`-sorted buckets `(u64 a_slot, s32 key_a)` = `12` B/entry (`2^63.27` B) + `u64` heads (`<= 2^57` B) = **`<= 2^63.29` B `< 2^64` B** per table.
  - *Band 2 (`[2^54, P_hi)`, `P_g = 2^64`, `P_f = 2^60`):* Super-buckets of width `W_SB = 128` with `16`-byte `rev_entry_t` (`2^63.69` B) + `u64` heads and 1-bit occupancy bitmap (`<= (2^64 / 128) x 8.125 = 2^60.02` B) + `256` open 1 MiB staging pages (`256` MiB) = **`<= 2^63.778` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Transient build scratch:* `ka_arr[N_a^tbl]` (`2^62.69` B `< 2^64` B, reused across tables).
  - *Total replay `a`-table storage (`g` to `2^64` + `f` to `2^60`, 4-orbit shared):* **`2^120.39` bytes** (`< 2^120.5` B), built once before collection and freed immediately after collection.
- **Tiled `b_0`-rows & `256` Column Owners (exact `L_max = 536.89`):** A `b_0`-row in `D_(omega,sigma)` has exact maximum length `L_max = 2 sqrt(1.36 R_b^2 / 8) = 536.89` (mean `L_b = 230.60`). Every `b_0`-row is partitioned into `ceil(L_b / 256) in {1, 2, 3}` contiguous tiles of width `<= W_tile = 256` (mean `1.405` tiles/row, `3.200` spanned super-buckets/row). During each tile, **256 column owners (`t = 0..255`)** hold the score arrays for the active tile:
  - Column owner `t` holds `S_g[t][0..N_a^tbl-1]` (`2^60.69` B) and `S_f[t][0..N_a^tbl-1]` (`2^60.69` B), total **`2^61.667` B `<< 2^64` B** (verified live in `audit_cm.py`).
  - *Band-1 transfer (priced in `C_HIT_LO = 8.0`):* Table owner binary-searches `k_b = k_0 + t` and copies the contiguous matching slice `[0, hi)` (`12` B/entry, `1.50` Ir/entry) to column owner `t`, which applies locally via `rev_b1_apply_u4` (`6.25` Ir/hit).
  - *Band-2 transfer (priced in `C_HIT_HI = 19.0`):* Table owner scans spanned super-buckets and appends 16-byte `(a_slot, logp)` records (`9` Ir/hit) into column `t`'s open 1 MiB staging page (`256` pages = `256` MiB total per walker); full pages hand off and copy to column owner `t` (`4` Ir/record), which applies locally via `rev_b2_apply_u4` (`5.75` Ir/hit) and returns the page to the free list.
- **Stage B & SRT waves:** `P24` product (`24,197,143` bits = `3.1` MiB), batched waves `< 2^27` B; `P_g, P_f` product trees `< 2^76` B.

## 2. Linear Algebra (`1-Owner Architecture B`, `< 2^64` Bytes Total Including Output Block!)

At our 4-parameter deterministic hard caps (`N_l <= 2^54.5215`, `N_sw <= 2^55.5215`, `k5 <= 10.50`, `tab_r <= 0.85`, `(m, n) = (8192, 110)` in `22` groups of 5 vectors):
- **Input AND Output 5-vector blocks (Architecture B, 160 B/slot):** Each slot stores only the 5-vector slice (`5 x 32 B = 160 B`; no pre-stored `chi^k x_j` needed because the 6-bucket accumulator folds `chi` once per row tail). Holding **BOTH** the input 5-vector block `x` (`N_sw x 160 B = 2^62.843` B) **AND** the output 5-vector block `y` (`2^62.843` B) takes `2^63.843` B.
- **Dense-band subset-sum tables (`{2: 28.0, 3: 18.25, 4: 14.0}`):** Because `12^g / 2 * C(pi(X)/6, g) ~= 4^g / 2 * C(pi(X)/2, g)` (invariant under order-12 quotienting), table size is `2^47.715` elements (`2^55.04` B per 5-vector group, `2^59.497` B across all 22 groups).
- **1-owner fit (verified live in `audit_cm.py`):** `2 x N_sw x 160 + 2^59.497 =` **`2^63.913` B `< 2^64` B** on ONE single 64-bit owner (`0` cross-owner sends, `0` combines, `0` routing!).
- **Sparse matrix, checkpoints, & BM generator:** Stored matrix `< 2^69.5` B, Horner checkpoints (`4 n = 440` vectors) `< 2^70.2` B, Kronecker BM workspace `< 2^84` B.

## 3. Descent Nodes & Aggregate Peak Memory

- **Descent:** Streaming owner layout across the single active special-`q` node (`< 2^105` B peak; `boot` uses `< 2^30` B).
- **Wrapper:** `< 2^30` B.
- **Aggregate Peak Memory:** **`< 2^120.5` bytes** (dominated by the `2^120.39`-byte replay `a`-tables during collection; all subsequent phases use `< 2^105` bytes). All tables are built from scratch for this single target and charged in full in `W`.
