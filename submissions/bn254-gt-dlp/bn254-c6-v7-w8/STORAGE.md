# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v7-w8`, Pure Product `kap = 0.0`)

Every single-owner working set fits in `< 2^64 B`; global read-only RAM `< 2^100 B` (ZERO Band 3).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v7-w8` (`s = 0.80, kap = 0.0`, `Lf = 61.0, Lg = 60.0, exc = 0.82, W = 160`, `log2_V = 97.8362`, `log2_Na = 45.5168`, `log2_Vb = 52.3194`, `L_max = 807.40`, `P_min = 14, P_cut = 45.0, P_mid = lx2 = 54.5, W_SB = 128, max_sb = 8`):

1. **Unified Single-Arena `(b_1..b_5)` Row Worker (`uint16_t`, row-major `(t, a_slot)`, ZERO duplicate `score_b0`!):**
   - `(L_max + 2) * N_a * 2 B + N_a * 2 B (bias)` = **`2^56.179 B`** (`226x` below `2^64 B`; `V_b * 2 B = 2^53.32 B`; row-pointer table `L_max * 8 B <= 6.5 KiB` L1).
2. **Band-1 / Band-2 Prime-Table Working Sets (`<< 2^64 B`/owner):**
   - **Band-1 Query (1 prime `p < 2^{45}` + row arena):** `8B * N_a + 8B * p + arena <= 2^{56.19} B << 2^64 B`.
   - **Band-2 Query (1 prime `p <= 2^{54.5}` + row arena):** `8B * N_a` (`entries[]`) + `8.125B * S_p` (`off[]` + 1-bit `occ_bm[]`) + arena = **`2^56.215 B << 2^64 B`** (`b2_1p_ram_log2`, LIVE `ram_ok`).
   - **Zero-Sort Sequential 2-Pass Builder (1 Band-2 prime at a time):** RO `a_lines[]` (`8B*N_a`) + scratch `tmp_sb[]` (`8B*N_a`) + `tmp_rec[]` (`8B*N_a`) + `entries[]` (`8B*N_a`) + `cnt[]` (`4B*S_p`) + `off[]` (`8B*S_p`) + `cur[]` (`8B*S_p`) + 1-bit `occ_bm[]` (`0.125B*S_p`) = `32*N_a + 20.125*(2^{54.5}/128)` = **`2^52.319 B << 2^64 B`** (`b2_build_1p_ram_log2`, LIVE `ram_ok`).
3. **1-Owner `Z[omega]` BW LA (4-vec register-resident, `N_sw_cap = 2^43.9388`, `N_l_cap = 2^43.3173`, `bw = (1024, 35)`, `nslots = 7`):**
   - `2 * 160 B * 7 * N_sw_cap` + matrix `W_ROW * 8 B * N_sw_cap` + `2^55.04 B * 7` (ZERO TABM) + `65 * 2^53 B` = **`2^59.649 B`** (`20.4x < 2^64 B`, LIVE `ram_ok`).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.4167`):** `2 * P_x` = **`2^62.945 B < 2^64 B`** (LIVE `ram_ok`).

## 2. Global Storage & Conservative `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band 3):** Band-1 `[2^14, 2^45)`: `8*N_a*n_lo + 8*sp_lo <= 2^{90.2} B`; Band-2 `[2^45, 2^54.5)`: `8*N_a*n_hi + 8.125*(sp_hi/128) = 2^{98.39} B`; **total `2^98.394 B < 2^100 B`** (`+1.61`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** Raising `C_SWAR_PROBE` from `12.0` to `22.0` moves the gated claim by only `+0.023`b to `2^86.643360`.
