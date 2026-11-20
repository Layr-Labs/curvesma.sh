# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6`)

Every working-set buffer in `BN254 v6` fits comfortably inside a single 64-bit address space (`< 2^64 B` per owner), and total global read-only/bucket storage stays strictly below `2^100 B` (matching accepted `bn15r` `c5-gt`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6c` operating point (`Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 96.9617`, `log2_Na = 43.8648`, `log2_Vb = 54.2626`, `L_max = 1010.62`, `P_cut = 44.0`, `P_mid = 54.0`, `lx2 = 59.0`, `W_SB = 256`, `max_sb = 5 <= 8`), **strict `< 2^64 B` per-owner buffered hit delivery is the PRIMARY model** (`C_HIT_LO = 6.25`, `C_HIT_HI = 14.5`, `C_BKT_B3 = 16.0`):

1. **Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table)` = **`2^54.8501 B`** (`~32.8 PiB`, **`570x` below `2^64 B`**).
   - For reference, even `V_b * 2 B` = **`2^55.2626 B`** (`<< 2^64 B`).
2. **Sharded Band-1 / Band-2 Prime-Table Owners (`< 2^64 B` each, explicit 8B `log_p`-grouped cross-owner hit buffers):**
   - Each prime's Band-2 `a`-table (`16 B * N_a = 2^47.86 B`) + `sb_u7` directory (`9 B * (2^54 / 256) = 2^49.17 B`) fits `< 2^50 B`; independent `< 2^64 B` prime-table owners hold disjoint shards of primes and deliver hit indices to the score-arena owner explicitly via 8-byte `log_p`-grouped buffers (`C_HIT_LO = 6.25`, `C_HIT_HI = 14.5 = 6 emit + 2 copy + 4.25 apply + 2.25 headroom`).
3. **1-Owner `Z[omega]` Block Wiedemann LA (at 4-parameter hard cap `N_sw_cap = 2^44.3538`, `bw = (768, 35)`, `nslots = 7`):**
   - Two `5`-lane `F_r` vector buffers (`2 * 160 B * 7 * N_sw_cap`) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + dense-band TABM tables (`2^55.04 B * 7`) + `65 * 2^53 B` scratch = **`2^59.68 B`** (**`20x` below `2^64 B`**, zero multi-owner LA split needed).

## 2. Global Storage & Conservative `bkt = 25` Sensitivity

- **Global read-only `a`-tables + Band-3 buckets:**
  - Band-1 `[2^14, 2^44)` `a`-tables + `256`-shell directories: `16 * N_a * n_lo + 8 * 256 * sp_lo = 2^87.9 B`.
  - Band-2 `[2^44, 2^54)` `a`-tables + `sb_u7` directories: `16 * N_a * n_hi + 9 * (sp_hi / 256) = 2^96.5 B`.
  - Band-3 `[2^54, 2^59)` 8-byte `log_p`-grouped buckets: `8 B * V * ln(59/54) = 2^96.5 B`.
  - **Total global RAM:** **`2^97.69 B < 2^100 B`** (`audit_bn_cm.py` asserts `ram_log2 < 100.0`).
- **Conservative `C_BKT_B3 = 25.0` Sensitivity (`part="sensitivity"`):**
  Raising Band-3 bucket delivery from `16.0` to `25.0` Ir/hit (keeping strict `< 2^64 B` per-owner buffered `C_HIT_LO = 6.25, C_HIT_HI = 14.5`) moves the gated claim by only `+0.021` bits to **`2^87.1271`** (`-0.461` full bits below `#1`).
