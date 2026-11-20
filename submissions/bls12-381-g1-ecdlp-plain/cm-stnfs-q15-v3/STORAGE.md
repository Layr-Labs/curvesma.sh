# Storage and Addressing

Every object is explicitly bounded and partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner). No single object requires `> 64`-bit indexing, and every cross-owner byte transfer is explicitly priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Rev-SigmaB-v3`, `W_tile = 544`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

- **Canonical `a` packing & `key_a`-sorted 4-bank layout:** There are `N_a^can = 2^57.7874` canonical `a`-orbits (`oid in [0, N_a^can)`), globally sorted in ascending unified 2-adic key `key_a(oid) = round(2^16 (q_a(a) - s_c(a)))`. Because `sigma` is a `G`-isometry preserving `q_a` and 2-adic class `s_c`, all 4 `sigma`-conjugates `(oid, m)` (`m in {0, 1, 2, 3}`) share the exact same `key_a(oid)`. Packing `a_slot = (m << 58) | oid` (`N_a^tbl = 4 N_a^can = 2^59.7874 < 2^60` slots; top 2 bits = `m`, low 58 bits = `oid < 2^57.80`) guarantees that for **any** `b`-column `t` with quadratic value `q_b(t)`, the active cells `{(oid, m) : key_a(oid) <= -q_b(t)}` are **four exact contiguous prefixes `[0, N_act_orb(t))`**, so `rev_scan_zero_u4` scans and zeroes **ONLY** active prefix cells (`V_all = 1.0276 V_kept`, zero bounding-box overhead).
- **Global pre-sort (`perm[]`):** Built once before table construction (`perm[rank]`, `a_slot_sorted[rank]`, `key_sorted[rank]`: `2^62.79` B each, `< 2^64` B per owner).
- **4-orbit `sigma`-shared `a`-tables (1 owner per 4-orbit table, every table `< 2^64` B):**
  - *Band 1a (`[2^24, 2^48)`, `m`-partitioned) & Band 1b (`[2^48, 2^55)`, unified `rev_b1_early4_find`):* `key_a`-sorted buckets `(u64 a_slot, s32 key_a)` = `12` B/entry (`2^63.37` B) + `u64` heads (`<= 2^58` B) + `544` MiB staging = **`<= 2^63.407` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Band 2 (`[2^55, P_hi)`, `P_g = 2^63`, `P_f = 2^60`):* Super-buckets of width `W_SB = 64` with `16`-byte `sb_entry_t` (`2^63.79` B) + `u64` heads and 1-byte `sb_min_q8` (`9` B/SB, `<= (2^63 / 64) x 9 = 2^60.17` B) + `544` open 1 MiB staging pages (`544` MiB) = **`<= 2^63.900` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Transient build scratch:* `ka_arr[N_a^tbl]` (`2^62.79` B `< 2^64` B, reused across tables).
  - *Total replay `a`-table storage (`g` to `2^63` + `f` to `2^60`, 4-orbit shared):* **`2^119.61` bytes** (`< 2^120.0` B), built once before collection and freed immediately after collection.
- **Full-Chord Rows & `544` Column Owners (exact `L_max = 542.53 < W_tile = 544`):** A `b_0`-row in `D_(omega,sigma)` has exact maximum length `L_max = 2 sqrt(1.36 R_b^2 / 8) = 542.53` (mean `L_b = 233.03`), so every `b_0`-row fits in **one single full-chord tile** (`nt = 1.000`, spanning `4.646` SBs/row with per-SB tight `q_min_j`). During each row, **544 column owners (`t = 0..543`)** hold the score arrays:
  - Column owner `t` holds `S_g[t][0..N_a^tbl-1]` (`2^60.79` B) and `S_f[t][0..N_a^tbl-1]` (`2^60.79` B), total **`2^61.787` B `<< 2^64` B** (verified live in `audit_cm.py`).
  - *Band-1a/1b transfers (`C_HIT_LO_M4 = 6.25`, `C_HIT_LO_U = 8.0`):* Table owner copies matching slices (`8` B/hit via `rev_copy_u64_u4` in `m`-partitioned Band 1a; `12` B/hit in Band 1b) to column owner `t`, which applies locally (`rev_apply_8b_u4` / `rev_b1_apply_u4`).
  - *Band-2 transfer (`C_HIT_HI = 14.5`):* Table owner gates each SB via `rev_b2_minkey_gate`, scans passing SBs with per-SB `q_min_j`, and emits 8-byte `oid_m` records grouped by `0.25`-bit `log_p` bucket into column `t`'s open 1 MiB staging page (`544` MiB/walker); full pages copy via `rev_copy_u64_u4` (`2.0` Ir/hit) to column owner `t`, which applies locally via `rev_apply_8b_u4` (`4.25` Ir/hit).
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
- **Aggregate Peak Memory:** **`< 2^120.0` bytes** (dominated by the `2^119.61`-byte replay `a`-tables during collection; all subsequent phases use `< 2^105` bytes). All tables are built from scratch for this single target and charged in full in `W`.
