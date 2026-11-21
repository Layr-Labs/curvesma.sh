# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v18-lbs`, `kap = 0.0`)

Every single-owner working set `< 2^64 B`; global RAM `< 2^100 B` (`2^{98.092} B`, ZERO `occ_bm[]`).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v18-lbs` (`s=0.80, kap=0.0, Lf=61, Lg=60, exc=0.82, W=160, log2_V=96.5761, log2_Na=44.8868, log2_Vb=51.6893, L_max=750.72, P_min=14, P_cut=41.5, P_list=54.5, P_mid=lx2=60.0, sl=80, sb=55, B3_SKIP_YIELD_LOSS=0.0080`, hybrid Band 2 `[41.5, 54.5)` + Band-3 $\tau$-orbit 12D lattice bucket-sieve `[54.5, 60.0)`, `p_occ = 1.0`, ZERO `occ_bm[]`):

1. **Power-of-2 Row-Strided Pre-Initialized Single-Chord Arena (`(b_1..b_5)` worker, `uint16_t`, `N_a_pad = 2^{45}`, ZERO `x_rtab`, ZERO duplicate `score_b0`):** `(L_max + 2)*2^{46} + 4*N_a = 2^{55.560} B` (`347x < 2^64 B`; `4*N_a` holds BOTH `init_a[]` (`shift = 0`, `8/9` of rows) AND `init_a_p2[]` (`shift = 48`, `1/9` of rows); `V_b*2 = 2^{52.689} B`; 1-based `ai in [1..N_a]`, `N_a + 1 < 2^{45}`).
2. **Band-1 / Hybrid Band-2 Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + arena): `8*N_a + 8*p + arena <= 2^{55.57} B << 2^64 B`.
   - Band-2 (worst-case 1 prime across all 7 sub-bands `[41.5, 54.5)`, `W in {8..256}` + arena): `8*(N_a+1) + 8.0*max(2^{b_hi}/W_k) + arena = 2^{55.588} B << 2^64 B` (`b2_1p_ram_log2`; `max(2^{b_hi}/W_k) = 2^{54.5}/256 = 2^{46.5}`).
   - Hybrid Zero-Sort Builder (worst-case 1 Band-2 prime): `max(32*N_a + 20.0*max_2p(2^{b_hi}/W_k), 8*(N_a+1) + 8.0*max_ll(2^{b_hi}/W_k)) = 2^{50.112} B << 2^64 B` (`b2_build_1p_ram_log2`; 1-pass chained `[49.0, 54.5)` needs `4x` less build RAM with ZERO `tmp_kta[]`).
3. **Band-3 $\tau$-Orbit 12D Lattice Bucket-Sieve `[54.5, 60.0)` Working Sets (`<< 2^64 B`/owner, `B3_BATCHES = 2^{37}`, `B3_SUPERBUCKETS = 2^{39}`, `B3_CHORDS_PER_SB = 20`):**
   - **Phase-1 Batch Owner (`b3_p1_1own_ram_log2 = 59.600 < 64.0`, `21.1x < 2^64 B`):** Holds read-only padded `canon_a` (`2^{57.269} <= 2^{57.5} B`, `+2` coord padding on all 6 dims), read-only padded 16B `table_b` (`2^{51.080} <= 2^{51.5} B`, `+3` coord padding on all 5 dims), `128.5 KiB` (`65,792 x uint16_t`) `SQRT16` table, per-super-bucket page-chain headers (`16*2^{39} = 2^{43} B`), active `4 KiB` write pages (`4096*2^{39} = 2^{51} B`), and all flushed `4 KiB` pages emitted by one batch (`byt_b3 / 2^{37} = 2^{59.205} B`).
   - **Phase-2 Super-Bucket Owner (`b3_p2_1own_ram_log2 = 60.480 < 64.0`, `11.5x < 2^64 B`):** Holds the flat `20`-chord `uint16_t` score array `20 * 1024 * 2^{46} B = 2^{60.322} B` (`elem_idx = (chord<<55)|(b0_rel<<45)|ai`, `60` bits), `init_a` + `init_a_p2` (`4*N_a`), all `2^{37}` batch headers (`16*2^{37} = 2^{41} B`), all Band-3 bucket pages across all `2^{37}` batches for this super-bucket (`byt_b3 / 2^{39} = 2^{57.205} B`), and one worst-case Band-2 prime table (`8*(N_a+1) + 8*2^{46.5} = 2^{49.9} B`) simultaneously.
4. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^43.8899, N_l_cap=2^43.2798, bw=(1024,35), nslots=7, TAB_R_CAP=1.00`):** `2*160*7*N_sw_cap + W_ROW*8*N_sw_cap + 2^55.04*7 + 65*2^53 = 2^{59.646} B` (`20.4x < 2^64 B`).
5. **1-Owner Boot Tree (`x_EA = 65.3021`, `boot = [180, 76, 12, 1]`, `K_h`-split `0.30`):** `2*P_x = 2^{62.094} B << 2^64 B`.

## 2. Global Storage & Sensitivity

- **Global RAM (`2^{98.092} B < 2^{100} B`, `+1.908`b headroom, LIVE `ram_ok`):** Band-1 `[2^{14}, 2^{41.5})` read-only `a`-tables `2^{83.661} B`; Hybrid Band-2 `[2^{41.5}, 2^{54.5})` `a`-tables `2^{97.642} B`; Band-3 `[2^{54.5}, 2^{60.0})` `4 KiB` (`511`-entry + next ptr) page-linked super-bucket buffers plus `(batch, sb)` header arrays `8*(512/511)*ln(60/54.5)*V + (4096+16)*2^{37}*2^{39} = 2^{96.205} B`; Phase-2 `20`-chord super-bucket score array `2^{60.322} B`; **total `2^{98.092} B < 2^{100} B`** (smaller than `bn254-c6-v17-pi`'s `2^{98.879} B` because `N_a` drops and `[54.5, 55.5)` `a`-lists are eliminated).
- **Sensitivity (`part="sensitivity"`):** `conservative_swar22` leaves the gated claim identically `2^85.449179` (`probe = 0.0`); `band3_c_lbs_85` (`C_LBS = 85.0` vs `81.0`) gives `2^85.473252`; `flat_shared_arena_ref` gives `2^85.458938`; `plain_dickman_rho` gives `2^85.311337`; `de_bruijn_tighter_cap` (`cap_scale = 0.80`) gives `2^85.355022`.
