# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v18`, `kap = 0.0`)

Every single-owner working set `< 2^64 B`; global RO RAM `< 2^100 B` (ZERO Band 3, ZERO `occ_bm[]`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v18` (`s=0.80, kap=0.0, Lf=61.5, Lg=60.5, exc=0.82, W=160, log2_V=97.3839, log2_Na=45.2907, log2_Vb=52.0932, L_max=786.58, P_min=17, P_cut=41.5, P_mid=lx2=55.50, sl=172, sb=136`, hybrid 2-pass `x4` (`[41.5,49.0)`, `W=8..64`) + 1-pass chained `head[sb]`+`next[ai]` (`[49.0,55.50)`, `W=64..512`) Band 2, `p_occ = 1.0`, ZERO `occ_bm[]`):

1. **Power-of-2 Row-Strided Pre-Initialized Arena (`(b_1..b_5)` worker, `uint16_t`, `N_a_pad = 2^{46}`, ZERO `x_rtab`, ZERO duplicate `score_b0`):** `(L_max + 2)*2^{47} + 4*N_a = 2^{56.625} B` (`166x < 2^64 B`; `4*N_a` holds BOTH `init_a[]` (`shift = 0`, `8/9` of rows) AND `init_a_p2[]` (`shift = 48`, `1/9` of rows); `V_b*2 = 2^{53.093} B`; 1-based `ai in [1..N_a]`, `N_a + 1 < 2^{46}`).
2. **Band-1 / Hybrid Band-2 Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + arena): `8*N_a + 8*p + arena <= 2^{56.64} B << 2^64 B`.
   - Band-2 (worst-case 1 prime across all 8 sub-bands `W in {8..512}` + arena): `8*(N_a+1) + 8.0*max(2^{b_hi}/W_k) + arena = 2^{56.640} B << 2^64 B` (`b2_1p_ram_log2`; `max(2^{b_hi}/W_k) = 2^{55.50}/512 = 2^{46.50}`).
   - Hybrid Zero-Sort Builder (worst-case 1 Band-2 prime): `max(32*N_a + 20.0*max_2p(2^{b_hi}/W_k), 8*(N_a+1) + 8.0*max_ll(2^{b_hi}/W_k)) = 2^{50.464} B << 2^64 B` (`b2_build_1p_ram_log2`; 1-pass chained `[49.0, 55.50)` needs `4x` less build RAM with ZERO `tmp_kta[]`!).
3. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^43.5197, N_l_cap=2^43.0433, bw=(1024,40), nslots=8, TAB_R_CAP=1.00`):** `2*160*8*N_sw_cap + W_ROW*8*N_sw_cap + 2^55.04*8 + 65*2^53 = 2^{59.690} B` (`19.8x < 2^64 B`).
4. **1-Owner Split Batch (`SPLIT_Q = 2^{46}, SPLIT_REC = 56 B`):** `56 * 2^{46} = 2^{51.807} B << 2^64 B`.
5. **1-Owner Boot Tree (`x_EA = 65.3021`, `boot = [180, 76, 12, 1]`, `ea = [5, 1250]`, conservative `0.30` cap):** `2*P_{x_EA} = 2^{62.094} B << 2^64 B`.

## 2. Global Storage & Sensitivity

- **Global read-only `a`-tables (ZERO Band 3, ZERO `occ_bm[]`):** Band-1 `[2^{17}, 2^{41.5})`: `8*N_a*n_lo + 8*sp_lo <= 2^{86.6} B`; Hybrid Band-2 `[2^{41.5}, 2^{55.50})`: `sum_k [8*(N_a+1)*(0.5*bt_k) + 8.0*(sp_k/W_k)] = 2^{98.896} B`; **total `2^{98.896} B < 2^{100} B`** (`+1.10`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** Because Band 2 has ZERO `occ_bm[]` (`probe = 0.0`), `conservative_swar22` leaves the gated claim identically `2^85.451403`; `plain_dickman_rho` gives `2^85.316398`; `de_bruijn_tighter_cap` (`cap_scale = 0.80`) gives `2^85.390598`; `flat_shared_arena_ref` gives `2^85.465391`.
