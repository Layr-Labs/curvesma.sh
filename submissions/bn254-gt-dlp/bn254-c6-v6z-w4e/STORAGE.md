# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 v6z-w4e`)

Every single-owner working set fits in `< 2^64 B`; global read-only RAM `< 2^100 B` (ZERO Band 3).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v6z-w4e` (`s=0.80, kap=0.322`, `Lf=61, Lg=60, exc=0.82, W=160`, `log2_V=97.5887`, `log2_Na=45.9755`, `log2_Vb=52.7781`, `L_max=851.34`, `P_0=30, N_a_blk=2^35, B_a=2014`, `P_b1_dense=41, P_cut=44.5, P_mid=lx2=54.75, W_SB=128, max_sb=8`):

0. **Band-0 Block Owners (`B_a = 2014` owners of `2^35` sorted `a`-slots):**
   `score_b0` `(L_max+2)*2^35*2B` + `a_loc[]` `8B*2^35*n_0` + ptrs `8B*sp_0` = **`2^62.7192 B`** (`2.43x < 2^64 B`, LIVE `ram_ok`; global `2^73.70 B`).
1. **Band-1/2 Score-Arena Owner (`uint16_t`, row-major `(t, a_slot)`):**
   `(L_max+2)*N_a*2B + N_a*4B (bias) + N_a*8B (Q_a)` = **`2^56.714 B << 2^64 B`** (`M_t_s8 <= 6.8 KiB` L1; `V_b*2B = 2^53.78 B`).
2. **Sharded Band-1/2 Prime-Table Owners (`< 2^64 B`/owner):**
   - **Query (1 Band-2 prime `p <= 2^{54.75}`):** 8B `entries[]` (`8B*N_a`) + `off[]` (`8B*S_p`) + `sb_u7[]` (`1B*S_p`) + arena = **`2^56.747 B << 2^64 B`** (`b2_1p_ram_log2`).
   - **Zero-Sort 2-Pass Builder (1 prime at a time):** shared RO `slot_of_lex[]` (`8B*N_a`) + `a_line_starts` (`<=1B*N_a`) + `qa_u7_tbl[]` (`1B*N_a`) + scratch `tmp_sb[]` (`8B*N_a`) + `tmp_lo[]` (`1B*N_a`) + `entries[]` (`8B*N_a`) + `cnt[]` (`4B*S_p`) + `off[]` (`8B*S_p`) + `cur[]` (`8B*S_p`) + `sb_u7[]` (`1B*S_p`) = `27*N_a + 21*(2^{54.75}/128)` = **`2^52.603 B << 2^64 B`** (`b2_build_1p_ram_log2`, LIVE `ram_ok`).
   - Cross-owner 8B hit buffers carry 56-bit `t*N_a + a_slot` (`COPY64=1.25, APPLY8=3.75`).
3. **1-Owner `Z[omega]` BW LA (SoA 5-vec, `N_sw_cap = 2^43.8774`, `bw=(1024,35)`, `nslots=7`):**
   `2*160B*7*N_sw_cap` + matrix `W_ROW*8B*N_sw_cap` + `2^55.04B*7` (ZERO TABM) + `65*2^53B` = **`2^59.646 B`** (`20.4x < 2^64 B`).
4. **1-Owner Boot Prime-Product Tree (`x_EA = 64.5208`):** `2*P_x` = **`2^63.050 B < 2^64 B`** (LIVE `ram_ok`).

## 2. Global Storage & `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band 3):** Band-0/1 `<= 2^90.4 B`; Band-2 8B `entries[] + off[] + sb_u7[]` = `2^99.08 B`; **total `2^99.078 B < 2^100 B`** (LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** `C_SWAR_PROBE = 22.0` gives `2^86.960532` (`+0.029`b, still `< 2^87.0`).
