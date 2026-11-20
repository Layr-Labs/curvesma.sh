# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6z-w3b`)

Every working-set buffer in `BN254 v6z-w3b` fits comfortably inside a single 64-bit address space (`< 2^64 B` per owner), and total global read-only storage stays strictly below `2^100 B` (with **ZERO Band-3 bucket storage**).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6z-w3b` operating point (`skew = 0.80, kap = 0.322`, `Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 97.2722`, `log2_Na = 45.8173`, `log2_Vb = 52.6198`, `L_max = 835.92`, `P_0 = 30.0`, `N_a_blk = 2^35`, `B_a = 1805`, `P_cut = 44.0`, `P_mid = lx2 = 55.5`, `W_SB = 256`, `max_sb = 5 <= 8`), Band-0 runs local-owner direct-apply (`C_HIT_B0 = 4.00`, `2^62.72 B << 2^64 B` per owner), Band 1 uses buffered 16B entries (`C_HIT_LO = 13.25`), and Band 2 uses **8-byte packed `(a_slot << 8) | r_lo`** records (`C_CAND = 12.0`, `C_HIT_HI = 14.5`, zero Band 3):

0. **Band-0 Local-Owner Direct-Apply Block Owners (`B_a = 1805` independent owners of `N_a_blk = 2^35` sorted `a`-slots):**
   - Local `score_b0` arena `(L_max + 2) * 2^35 * 2 B = 2^45.71 B` + local `a_loc[]` tables `8 B * 2^35 * n_0` + residue pointers `8 B * sp_0` = **`2^62.7192 B`** (**`2.43x` below `2^64 B`**; asserted live in `audit_bn_cm.py [ram_ok]`; global sum across `1805` owners `2^73.54 B << 2^99.58 B`).
1. **Band-1/2 Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`, receives explicit 64B streamed `score_b0` chunks priced at `0.4375` Ir/cell inside `C_SCAN_ZERO = 1.95`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 4 B (bias + bias_p2) + N_a * 8 B (Q_a_sorted for M_t[t])` = **`2^56.535 B`** (**`177x` below `2^64 B`**; per-row `M_t_s8[0..L_b-1]` table is `<= 6.7 KiB` in L1 cache!).
   - For reference, even `V_b * 2 B` = **`2^53.6198 B`** (`<< 2^64 B`).
2. **Sharded Band-1 / Band-2 Prime-Table Owners (`< 2^64 B` each, explicit 8B `log_p`-grouped cross-owner hit buffers):**
   - Each prime's Band-1 (`8B keys[] + 8B a_slots[] = 16 B * N_a`) or Band-2 (**8B packed `(a_slot << 8) | r_lo` = `8 B * N_a = 2^48.82 B`**) `a`-table + `sb_u7` directory (`9 B * (2^55.5 / 256) = 2^50.67 B`) fits `< 2^51 B` (`b2_1p_ram_log2 = 56.561 << 64.0`); independent `< 2^64 B` prime-table owners hold disjoint shards of primes and deliver direct 56-bit linear arena indices `t * N_a + a_slot` (`ldr x_a, [x_p], #8; add x_h, x_row_base, x_a`) to the score-arena owner explicitly via 8-byte `log_p`-grouped buffers (`C_HIT_LO = 13.25 = 7 emit + 2 copy + 4.25 apply`, `C_HIT_HI = 14.5 >= 13.25`).
3. **1-Owner `Z[omega]` Block Wiedemann LA (at 4-parameter hard cap `N_sw_cap = 2^44.2118`, `bw = (1024, 35)`, `nslots = 7`):**
   - Two `5`-lane `F_r` vector buffers (`2 * 160 B * 7 * N_sw_cap`) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + conservative buffer (`2^55.04 B * 7`, ZERO TABM tables required at `TAB_R_CAP = 1.00`) + `65 * 2^53 B` scratch = **`2^59.67 B`** (**`20x` below `2^64 B`**, zero multi-owner LA split needed).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.42`):**
   - Product `P_x` has size `2^61.95 B` and even the doubled `2 * P_x` root workspace is **`2^62.95 B < 2^64 B`** (asserted live in `audit_bn_cm.py [ram_ok]`), fitting strictly inside a single 64-bit owner.

## 2. Global Storage & Conservative `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band-3 buckets!):**
  - Band-0 `[2^14, 2^30)` + Band-1 `[2^30, 2^44)` `a`-tables + `256`-shell directories: `<= 2^89.8 B`.
  - Band-2 `[2^44, 2^55.5)` **8B packed** `a`-tables + `sb_u7` directories: `8 * N_a * n_hi + 9 * (sp_hi / 256) = 2^99.58 B`.
  - **Total global RAM:** **`2^99.58 B < 2^100 B`** (`audit_bn_cm.py` asserts `ram_log2 < 100.0`).
- **Conservative `C_SWAR_PROBE = 22.0` Sensitivity (`part="sensitivity"`):**
  Raising Band-2 SWAR probe cost from `18.0` to `22.0` Ir/row (keeping strict `< 2^64 B` per-owner buffered `C_HIT_LO = 13.25, C_HIT_HI = 14.5`) moves the gated claim by only `+0.014` bits to **`2^87.273`**.
