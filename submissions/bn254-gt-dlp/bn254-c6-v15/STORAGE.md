# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v15`, `kap = 0.0`)

Every single-owner working set `< 2^64 B`; global RO RAM `< 2^100 B` (ZERO Band 3, ZERO `occ_bm[]`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v15` (`s=0.80, kap=0.0, Lf=61, Lg=60, exc=0.82, W=160, log2_V=97.3839, log2_Na=45.2907, log2_Vb=52.0932, L_max=786.58, P_min=14, P_cut=41.5, P_mid=lx2=55.25, sl=170, sb=140`, 100%-dense 9-sub-band flat-offset Band 2 `(8/16/32/64 x4, 64/128 x2, 256/512/1024 x1)`, `p_occ = 1.0`, ZERO `occ_bm[]`):

1. **Power-of-2 Row-Strided Pre-Initialized Arena (`(b_1..b_5)` worker, `uint16_t`, `N_a_pad = 2^{46}`, ZERO `x_rtab`, ZERO duplicate `score_b0`):** `(L_max + 2)*2^{47} + 4*N_a = 2^{56.625} B` (`166x < 2^64 B`; `4*N_a` holds BOTH `init_a[]` (`shift = 0`, `8/9` of rows) AND `init_a_p2[]` (`shift = 48`, `1/9` of rows); `V_b*2 = 2^{53.093} B`; scan touches only `0..N_a-1` per row).
2. **Band-1 / 100%-Dense Band-2 Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + arena): `8*N_a + 8*p + arena <= 2^{56.64} B << 2^64 B`.
   - Band-2 (worst-case 1 prime across all 9 dense sub-bands `W in {8..1024}` + arena): `8*N_a + 8.0*max(2^{b_hi}/W_k) + arena = 2^{56.635} B << 2^64 B` (`b2_1p_ram_log2`; `max(2^{b_hi}/W_k) = 2^{55.25}/1024 = 2^{45.25}`).
   - Lex Zero-Sort 2-Pass Builder (worst-case 1 Band-2 prime): `32*N_a + 20.0*max(2^{b_hi}/W_k) = 2^{51.075} B << 2^64 B` (`b2_build_1p_ram_log2`).
3. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^43.7461, N_l_cap=2^43.1193, bw=(1024,35), nslots=7, TAB_R_CAP=1.00`):** `2*160*7*N_sw_cap + W_ROW*8*N_sw_cap + 2^55.04*7 + 65*2^53 = 2^{59.638} B` (`20.6x < 2^64 B`).
4. **1-Owner Boot Tree (`x_EA = 65.4062`, `boot = [182, 74, 12, 1]`, `K_h`-split `0.30`):** `2*P_x = 2^{62.198} B << 2^64 B`.

## 2. Global Storage & Sensitivity

- **Global read-only `a`-tables (ZERO Band 3, ZERO `occ_bm[]`):** Band-1 `[2^{14}, 2^{41.5})`: `8*N_a*n_lo + 8*sp_lo <= 2^{86.6} B`; 100%-dense 9-sub-band Band-2 `[2^{41.5}, 2^{55.25})`: `sum_k [8*N_a*(0.5*bt_k) + 8.0*(sp_k/W_k)] = 2^{98.144} B`; **total `2^{98.144} B < 2^{100} B`** (`+1.86`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** Because Band 2 is 100% dense (`probe = 0.0`, ZERO `occ_bm[]`), `conservative_swar22` leaves the gated claim identically `2^85.653930`; `plain_dickman_rho` gives `2^85.519606`; `de_bruijn_tighter_cap` (`cap_scale = 0.80`) gives `2^85.582256`.
