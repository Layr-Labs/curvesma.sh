# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v12`, `kap = 0.0`)

Every single-owner working set `< 2^64 B`; global RO RAM `< 2^100 B` (ZERO Band 3).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v12` (`s=0.80, kap=0.0, Lf=61, Lg=60, exc=0.82, W=160, log2_V=97.5900, log2_Na=45.3937, log2_Vb=52.1962, L_max=796.00, P_min=17, P_cut=41.5, P_mid=lx2=55.0, sl=165, sb=170`, 4-sub-band flat-offset Band 2 `(16d / 64d / 128d / 256s)` with 1-pass `bic`-tagged inline-slot sparse `[52.2, 55.0)`, `max_sb = 5 <= 8`):

1. **Power-of-2 Row-Strided Unified Arena (`(b_1..b_5)` worker, `uint16_t`, `N_a_pad = 2^{46}`, ZERO `x_rtab`, ZERO duplicate `score_b0`):** `(L_max + 2)*2^{47} + 2*N_a = 2^{56.641} B` (`164x < 2^64 B`; `V_b*2 = 2^{53.196} B`; scan touches only `0..N_a-1` per row).
2. **Band-1 / 4-Sub-Band Band-2 Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + arena): `8*N_a + 8*p + arena <= 2^{56.65} B << 2^64 B`.
   - Band-2 (worst-case 1 prime across `(16d/64d/128d/256s)` + arena, charging full `16*N_a` worst-case `ovf_nodes[]`): `16*N_a + 8.125*max(2^{b_hi}/W_k) + arena = 2^{56.665} B << 2^64 B` (`b2_1p_ram_log2`).
   - Lex Zero-Sort / 1-Pass Inline-Slot Builder (worst-case 1 Band-2 prime): `32*N_a + 20.125*max(2^{b_hi}/W_k) = 2^{51.937} B << 2^64 B` (`b2_build_1p_ram_log2`).
3. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^43.9388, N_l_cap=2^43.3173, bw=(1024,35), nslots=7, TAB_R_CAP=1.00`):** `2*160*7*N_sw_cap + W_ROW*8*N_sw_cap + 2^55.04*7 + 65*2^53 = 2^{59.649} B` (`20.4x < 2^64 B`).
4. **1-Owner Boot Tree (`x_EA = 65.3021`, `boot = [180, 76, 12, 4]`, `K_h`-split `0.30`):** `2*P_x = 2^{62.094} B << 2^64 B` (`2^{63.831} B < 2^64 B` even without `K_h`-split).

## 2. Global Storage & Conservative `swar_probe = 22` Sensitivity

- **Global read-only `a`-tables (ZERO Band 3):** Band-1 `[2^17, 2^{41.5})`: `8*N_a*n_lo + 8*sp_lo <= 2^{86.6} B`; 4-sub-band Band-2 `[2^{41.5}, 2^{55.0})`: dense `[41.5, 52.2)` `8*N_a*(0.5*bt_k) + 8*(sp_k/W_k)` + sparse inline-slot `[52.2, 55.0)` `8.125*(sp_k/256) + 32*(0.5*b_coll*N_a)` (charging 2 `ovf_nodes` = 32 B per colliding insertion!) = `2^{98.760} B`; **total `2^{98.760} B < 2^{100} B`** (`+1.24`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** Raising `C_SWAR_PROBE` from `11.0` to `22.0` moves the gated claim by only `+0.047`b to `2^85.905765` (still strictly `< 85.91 << 86.0`).
