# `STORAGE.md` — Memory & Single-Owner RAM Accounting (`BN254 v6z-w4`)

Every working set in `BN254 v6z-w4d` fits in `< 2^64 B`/owner; global read-only RAM `< 2^100 B` (ZERO Band 3).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At the `v6z-w4` operating point (`skew = 0.80, kap = 0.322`, `Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 97.5887`, `log2_Na = 45.9755`, `log2_Vb = 52.7781`, `L_max = 851.34`, `P_0 = 30.0`, `N_a_blk = 2^35`, `B_a = 2014`, `P_b1_dense = 41.0`, `P_cut = 44.5`, `P_mid = lx2 = 54.75`, `W_SB = 128`, `max_sb = 8 <= 8`):

0. **Band-0 Local-Owner Direct-Apply Block Owners (`B_a = 2014` owners of `N_a_blk = 2^35` sorted `a`-slots):**
   - Local `score_b0` `(L_max + 2) * 2^35 * 2 B` + local `a_loc[]` `8 B * 2^35 * n_0` + residue pointers `8 B * sp_0` = **`2^62.7192 B`** (`2.43x` below `2^64 B`, asserted LIVE in `ram_ok`; global sum `2^73.70 B << 2^99.08 B`).
1. **Band-1/2 Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`):**
   - `(L_max + 2) * N_a * 2 B + N_a * 4 B (bias + bias_p2) + N_a * 8 B (Q_a_sorted)` = **`2^56.714 B`** (`156x` below `2^64 B`; per-row `M_t_s8[0..L_b-1]` table is `<= 6.8 KiB` in L1 cache; `V_b * 2 B = 2^53.78 B << 2^64 B`).
2. **Sharded Band-1 / Band-2 Prime-Table Owners (`< 2^64 B` each, explicit builder & query RAM):**
   - **Query working set (1 Band-2 prime `p <= 2^{54.75}`):** 8B packed `entries[]` (`8 B * N_a`) + `off[]` (`8 B * S_p`) + `sb_u7[]` (`1 B * S_p`) + arena = **`2^56.747 B << 2^64 B`** (`b2_1p_ram_log2`).
   - **3-Pass Flat-Array Builder working set (1 Band-2 prime at a time, temps freed after `p`):** `tmp_sb[]` (`8B*N_a`) + `tmp_rec[]` (`8B*N_a`) + `entries[]` (`8B*N_a`) + `a_lines[]` (`8B*N_a`) + `qa_u7_tbl[]` (`1B*N_a`) + `cnt[]` (`4B*S_p`) + `off[]` (`8B*S_p`) + `cur[]` (`8B*S_p`) + `sb_u7[]` (`1B*S_p`) = `33*N_a + 21*(2^{54.75}/128)` = **`2^52.861 B << 2^64 B`** (`b2_build_1p_ram_log2`, asserted LIVE in `ram_ok`!). 
   - Cross-owner 8B `log_p`-grouped hit buffers deliver direct 56-bit linear indices `t * N_a + a_slot` (`COPY64 = 1.25`, `APPLY8 = 3.75`).
3. **1-Owner `Z[omega]` Block Wiedemann LA (SoA 5-vector layout, at 4-cap `N_sw_cap = 2^43.8774`, `bw = (1024, 35)`, `nslots = 7`):**
   - Two `5`-lane SoA `F_r` vector buffers (`2 * 160 B * 7 * N_sw_cap`, identical total bytes to AoS!) + matrix rows (`W_ROW * 8 B * N_sw_cap`) + conservative buffer (`2^55.04 B * 7`, ZERO TABM tables at `TAB_R_CAP = 1.00`) + `65 * 2^53 B` scratch = **`2^59.646 B`** (`20.4x` below `2^64 B`).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.5208`):**
   - Doubled root workspace `2 * P_x` = **`2^63.050 B < 2^64 B`** (asserted LIVE in `ram_ok`).

## 2. Global Storage & Conservative `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band-3 buckets!):**
  - Band-0 `[2^14, 2^30)` + Band-1 `[2^30, 2^44.5)` `a`-tables + `256`-shell directories: `<= 2^90.4 B`.
  - Band-2 `[2^44.5, 2^54.75)` **8B packed** `entries[]` + `off[]`/`sb_u7[]` directories: `8 * N_a * n_hi + 9 * (sp_hi / 128) = 2^99.08 B`.
  - **Total global RAM:** **`2^99.078 B < 2^100 B`** (`+0.92`b headroom, asserted LIVE in `ram_ok`).
- **Conservative `C_SWAR_PROBE = 22.0` Sensitivity (`part="sensitivity"`):**
  Raising Band-2 SWAR probe cost from `12.0` to `22.0` Ir/row moves the gated claim by only `+0.028` bits (to `2^87.020`).
