# Storage and Addressing

Every object is partitioned into 64-bit owner address spaces (`< 2^64` bytes per owner) with every cross-owner byte transfer priced in `KERNELS.md`.

## 1. Relation Collection (`Arch2-Rev-SigmaB-v3`, `W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

- **Canonical `a` packing & `key_a`-sorted 4-bank layout:** `N_a^can = 2^57.8294` canonical `a`-orbits sorted by unified 2-adic key `key_a(oid) = round(2^16 (q_a(a) - s_c(a)))`. Packing `a_slot = (m << 58) | oid` (`N_a^tbl = 4 N_a^can = 2^59.8294 < 2^60` slots) guarantees that for any column `t`, active cells `key_a(oid) <= -q_b(t)` are four exact contiguous prefixes `[0, N_act_orb(t))`, so `rev_scan_zero_u4` scans and zeroes ONLY active prefixes (`V_all = 1.0276 V_kept`). Global `perm[]` pre-sort arrays are `< 2^63` B each.
- **4-orbit `sigma`-shared `a`-tables (1 owner/table, every table `< 2^64` B):**
  - *Band 1a (`[2^24, 2^48)`, `m`-partitioned split-SoA `4B keys[] + 8B oids[]` = `12` B/entry for contiguous 8B `LDP` staging copy) & Band 1b (`[2^48, 2^55)`, interleaved `12` B/entry `rev_b1_early4_find` + 4x `ldur`):* `12` B/entry + `u64` heads + `552` MiB staging = **`<= 2^63.448` B `< 2^64` B** per table.
  - *Band 2 (`[2^55, P_hi)`, `P_g = 2^63`, `P_f = 2^60`):* `W_SB = 64`, `16`-byte `sb_entry_t` + `u64` heads & 1-byte `sb_min_q8` (`9` B/SB) + `552` MiB staging = **`<= 2^63.939` B `< 2^64` B** per table (verified live in `audit_cm.py`).
  - *Total replay `a`-table storage (`g` to `2^63` + `f` to `2^60`):* **`2^119.65` bytes** (`< 2^120.0` B), built once per target and freed immediately after collection.
- **Full-Chord Rows & `552` Column Owners (`L_max = 544.52 < W_tile = 552`):** Every `b_0`-row fits in one full-chord tile (`nt = 1.000`, `4.658` SBs/row). Each of the `552` column owners holds `S_g[t]` + `S_f[t]` = **`2^61.829` B `<< 2^64` B**. Band-1a/1b slice copies (`8`/`12` B) and Band-2 8B `log_p`-grouped staging-page handoffs (`552` MiB/walker) are explicitly priced in `C_HIT_LO_M4`, `C_HIT_LO_U`, and `C_HIT_HI`.
- **Stage B & SRT waves:** `P24` product (`3.1` MiB), waves `< 2^27` B; `P_g, P_f` product trees `< 2^76` B.

## 2. Linear Algebra (`1-Owner Architecture B`, `< 2^64` Bytes Total Including Output Block)

At hard caps `N_l <= 2^54.5215`, `N_sw <= 2^55.5215`, `(m, n) = (8192, 110)` (`22` groups of 5 vectors):
- **Input `x` + Output `y` (`160` B/slot each):** `2 x N_sw x 160 = 2^63.843` B.
- **Dense-band subset-sum tables (`{2: 28.0, 3: 18.25, 4: 14.0}`):** `2^47.715` elements (`2^59.497` B across all 22 groups).
- **1-owner fit (verified live in `audit_cm.py`):** `2 x N_sw x 160 + 2^59.497 =` **`2^63.913` B `< 2^64` B** (`0` cross-owner sends, combines, or routing). Stored matrix `< 2^69.5` B, Horner checkpoints `< 2^70.2` B, BM workspace `< 2^84` B.

## 3. Descent, Wrapper, & Aggregate Peak Memory

- **Descent:** Streaming single active special-`q` node `< 2^105` B (`boot` `< 2^30` B). **Wrapper:** `< 2^30` B.
- **Aggregate Peak Memory:** **`< 2^120.0` bytes** (`2^119.65` B replay `a`-tables during collection; `< 2^105` B thereafter).
