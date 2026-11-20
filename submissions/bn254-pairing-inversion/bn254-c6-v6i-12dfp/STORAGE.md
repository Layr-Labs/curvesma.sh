# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6i`)

Every working-set buffer in `BN254 v6i` fits comfortably inside a single 64-bit address space (`< 2^64 B` per owner), and total global read-only/bucket storage stays strictly below `2^100 B`.

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6i` operating point (`Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 96.9493`, `log2_Na = 43.8586`, `log2_Vb = 54.2564`, `L_max = 1009.90`, `P_cut = 43.0`, `P_mid = 55.0`, `lx2 = 59.0`, `W_SB = 256`, `max_sb = 5 <= 8`), **strict `< 2^64 B` per-owner buffered hit delivery is the PRIMARY model** (`C_HIT_LO = 13.25`, `C_HIT_HI = 14.5`, `C_BKT_B3 = 20.0`):

1. **Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table)` = **`2^54.8429 B`** (`~32.6 PiB`, **`573x` below `2^64 B`**).
   - For reference, even `V_b * 2 B` = **`2^55.2564 B`** (`<< 2^64 B`).
2. **Sharded Band-1 / Band-2 Prime-Table Owners (`< 2^64 B` each, explicit 8B `log_p`-grouped cross-owner hit buffers):**
   - Each prime's Band-1 (`8B keys[] + 8B a_slots[] = 16 B * N_a`) or Band-2 (`16 B * N_a = 2^47.86 B`) `a`-table + `sb_u7` directory (`9 B * (2^55 / 256) = 2^50.17 B`) fits `< 2^51 B` (`b2_1p_ram_log2 = 54.909 << 64.0`); independent `< 2^64 B` prime-table owners hold disjoint shards of primes and deliver direct 54-bit linear arena indices `t * N_a + a_slot` (`ldr x_a, [x_p], #8; add x_h, x_row_base, x_a`) to the score-arena owner explicitly via 8-byte `log_p`-grouped buffers (`C_HIT_LO = 13.25 = 7 emit + 2 copy + 4.25 apply`, `C_HIT_HI = 14.5 >= 13.25`).
3. **1-Owner `Z[omega]` Block Wiedemann LA (at 4-parameter hard cap `N_sw_cap = 2^44.3538`, `bw = (768, 30)`, `nslots = 6`):**
   - Two `5`-lane `F_r` vector buffers (`2 * 160 B * 6 * N_sw_cap`) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + dense-band TABM tables (`2^55.04 B * 6`) + `65 * 2^53 B` scratch = **`2^59.61 B`** (**`21x` below `2^64 B`**, zero multi-owner LA split needed).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.42`):**
   - Product `P_x` has size `2^61.95 B` and even the doubled `2 * P_x` root workspace is **`2^62.95 B < 2^64 B`** (asserted live in `audit_bn_cm.py [ram_ok]`), fitting strictly inside a single 64-bit owner.

## 2. Global Storage & Conservative `bkt = 25` Sensitivity

- **Global read-only `a`-tables + Band-3 buckets:**
  - Band-1 `[2^14, 2^43)` `a`-tables + `256`-shell directories: `16 * N_a * n_lo + 8 * 256 * sp_lo = 2^86.9 B`.
  - Band-2 `[2^43, 2^55)` `a`-tables + `sb_u7` directories: `16 * N_a * n_hi + 9 * (sp_hi / 256) = 2^98.6 B`.
  - Band-3 `[2^55, 2^59)` 16-byte `(a_key, b_idx)` `2^20`-hash-bucketed buffers: `16 B * V * ln(59/55) = 2^97.1 B` (drained with a `256 MB` local open-addressed `a_key -> local slice` table per bucket, ZERO global `2^43.87` lookup table!).
  - **Total global RAM:** **`2^98.92 B < 2^100 B`** (`audit_bn_cm.py` asserts `ram_log2 < 100.0`).
- **Conservative `C_BKT_B3 = 25.0` Sensitivity (`part="sensitivity"`):**
  Raising Band-3 bucket delivery from `20.0` to `25.0` Ir/hit (keeping strict `< 2^64 B` per-owner buffered `C_HIT_LO = 13.25, C_HIT_HI = 14.5`) moves the gated claim by only `+0.007` bits to **`2^87.480`**.
