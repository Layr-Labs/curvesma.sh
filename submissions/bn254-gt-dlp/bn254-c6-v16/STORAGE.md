# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v16`, `kap = 0.0`)

Every single-owner working set `< 2^64 B`; global RO RAM `< 2^100 B` (ZERO Band 3, ZERO `occ_bm[]`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v16` (`s=0.80, kap=0.0, Lf=61, Lg=60, exc=0.82, W=160, log2_V=97.2966, log2_Na=45.2470, log2_Vb=52.0496, L_max=782.63, P_min=14, P_cut=41.5, P_mid=lx2=55.50, sl=166, sb=136`, hybrid 2-pass `x4` (`[41.5,49.0)`, `W=8..64`) + 1-pass chained `head[sb]`+`next[ai]` (`[49.0,55.50)`, `W=64..512`) Band 2, `p_occ = 1.0`, ZERO `occ_bm[]`):

1. **Power-of-2 Row-Strided Pre-Initialized Arena (`(b_1..b_5)` worker, `uint16_t`, `N_a_pad = 2^{46}`, ZERO `x_rtab`, ZERO duplicate `score_b0`):** `(L_max + 2)*2^{47} + 4*N_a = 2^{56.618} B` (`167x < 2^64 B`; `4*N_a` holds BOTH `init_a[]` (`shift = 0`, `8/9` of rows) AND `init_a_p2[]` (`shift = 48`, `1/9` of rows); `V_b*2 = 2^{53.050} B`; 1-based `ai in [1..N_a]`, `N_a + 1 < 2^{46}`).
2. **Band-1 / Hybrid Band-2 Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + arena): `8*N_a + 8*p + arena <= 2^{56.63} B << 2^64 B`.
   - Band-2 (worst-case 1 prime across all 8 sub-bands `W in {8..512}` + arena): `8*(N_a+1) + 8.0*max(2^{b_hi}/W_k) + arena = 2^{56.633} B << 2^64 B` (`b2_1p_ram_log2`; `max(2^{b_hi}/W_k) = 2^{55.50}/512 = 2^{46.50}`).
   - Hybrid Zero-Sort Builder (worst-case 1 Band-2 prime): `max(32*N_a + 20.0*max_2p(2^{b_hi}/W_k), 8*(N_a+1) + 8.0*max_ll(2^{b_hi}/W_k)) = 2^{50.425} B << 2^64 B` (`b2_build_1p_ram_log2`; 1-pass chained `[49.0, 55.50)` needs `4x` less build RAM with ZERO `tmp_kta[]`!).
3. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^43.8899, N_l_cap=2^43.2798, bw=(1024,35), nslots=7, TAB_R_CAP=1.00`):** `2*160*7*N_sw_cap + W_ROW*8*N_sw_cap + 2^55.04*7 + 65*2^53 = 2^{59.646} B` (`20.4x < 2^64 B`).
4. **1-Owner Boot Tree (`x_EA = 65.3021`, `boot = [180, 76, 12, 1]`, `K_h`-split `0.30`):** `2*P_x = 2^{62.094} B << 2^64 B`.

## 2. Global Storage & Sensitivity

- **Global read-only `a`-tables (ZERO Band 3, ZERO `occ_bm[]`):** Band-1 `[2^{14}, 2^{41.5})`: `8*N_a*n_lo + 8*sp_lo <= 2^{86.6} B`; Hybrid Band-2 `[2^{41.5}, 2^{55.50})`: `sum_k [8*(N_a+1)*(0.5*bt_k) + 8.0*(sp_k/W_k)] = 2^{98.879} B`; **total `2^{98.879} B < 2^{100} B`** (`+1.12`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** Because Band 2 has ZERO `occ_bm[]` (`probe = 0.0`), `conservative_swar22` leaves the gated claim identically `2^85.577090`; `plain_dickman_rho` gives `2^85.445136`; `de_bruijn_tighter_cap` (`cap_scale = 0.80`) gives `2^85.491164`.
