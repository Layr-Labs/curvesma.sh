# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6z-noband3`)

Every working-set buffer in `BN254 v6z-noband3` fits comfortably inside a single 64-bit address space (`< 2^64 B` per owner), and total global read-only storage stays strictly below `2^100 B` (with **ZERO Band-3 bucket storage**).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6z-noband3` operating point (`Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 97.5260`, `log2_Na = 44.1469`, `log2_Vb = 54.5448`, `L_max = 1044.10`, `P_cut = 43.0`, `P_mid = lx2 = 55.5`, `W_SB = 256`, `max_sb = 6 <= 8`), **strict `< 2^64 B` per-owner buffered hit delivery is the PRIMARY model** (`C_HIT_LO = 13.25`, `C_HIT_HI = 14.5`, zero Band 3):

1. **Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 2 B (bias_table) + N_a * 8 B (Q_a_sorted for M_t[t])` = **`2^55.1846 B`** (`~41.3 PiB`, **`450x` below `2^64 B`**; per-row `M_t[0..L_b-1]` table is `<= 8.4 KiB` in L1 cache!).
   - For reference, even `V_b * 2 B` = **`2^55.5448 B`** (`<< 2^64 B`).
2. **Sharded Band-1 / Band-2 Prime-Table Owners (`< 2^64 B` each, explicit 8B `log_p`-grouped cross-owner hit buffers):**
   - Each prime's Band-1 (`8B keys[] + 8B a_slots[] = 16 B * N_a`) or Band-2 (`16 B * N_a = 2^48.15 B`) `a`-table + `sb_u7` directory (`9 B * (2^55.5 / 256) = 2^50.67 B`) fits `< 2^51 B` (`b2_1p_ram_log2 = 55.252 << 64.0`); independent `< 2^64 B` prime-table owners hold disjoint shards of primes and deliver direct 55-bit linear arena indices `t * N_a + a_slot` (`ldr x_a, [x_p], #8; add x_h, x_row_base, x_a`) to the score-arena owner explicitly via 8-byte `log_p`-grouped buffers (`C_HIT_LO = 13.25 = 7 emit + 2 copy + 4.25 apply`, `C_HIT_HI = 14.5 >= 13.25`).
3. **1-Owner `Z[omega]` Block Wiedemann LA (at 4-parameter hard cap `N_sw_cap = 2^44.3538`, `bw = (1024, 35)`, `nslots = 7`):**
   - Two `5`-lane `F_r` vector buffers (`2 * 160 B * 7 * N_sw_cap`) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + conservative buffer (`2^55.04 B * 7`, ZERO TABM tables required at `TAB_R_CAP = 1.00`) + `65 * 2^53 B` scratch = **`2^59.68 B`** (**`20x` below `2^64 B`**, zero multi-owner LA split needed).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.42`):**
   - Product `P_x` has size `2^61.95 B` and even the doubled `2 * P_x` root workspace is **`2^62.95 B < 2^64 B`** (asserted live in `audit_bn_cm.py [ram_ok]`), fitting strictly inside a single 64-bit owner.

## 2. Global Storage & Conservative `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band-3 buckets!):**
  - Band-1 `[2^14, 2^43)` `a`-tables + `256`-shell directories: `16 * N_a * n_lo + 8 * 256 * sp_lo = 2^87.2 B`.
  - Band-2 `[2^43, 2^55.5)` `a`-tables + `sb_u7` directories: `16 * N_a * n_hi + 9 * (sp_hi / 256) = 2^99.36 B`.
  - **Total global RAM:** **`2^99.36 B < 2^100 B`** (`audit_bn_cm.py` asserts `ram_log2 < 100.0`).
- **Conservative `C_SWAR_PROBE = 22.0` Sensitivity (`part="sensitivity"`):**
  Raising Band-2 SWAR probe cost from `18.0` to `22.0` Ir/row (keeping strict `< 2^64 B` per-owner buffered `C_HIT_LO = 13.25, C_HIT_HI = 14.5`) moves the gated claim by only `+0.031` bits to **`2^87.813`**.
