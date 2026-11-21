# `STORAGE.md` — Memory & Single-Owner RAM (`BN254 C6 v19`, `E12 + JointVor` + `2^{48}` Streaming Super-Tile Lattice Band 2)

Every single-owner working set is **100% deterministically `< 2^64 B`** (ZERO per-chord density cap or overflow assumption); global RAM `< 2^100 B` (`+1.87`b headroom).

## 1. Single-Owner Working Sets (`<< 2^64 B`)

At `v19` (`E12 + JointVor` `W_a=1764, W_b=1000`, `Lf=Lg=61.0, exc=0.76, W=160, log2_V=95.4978, log2_Na=48.1845, B_a=8` equal-count `Q(a)` shells -> `log2_Na_shell=45.1845 < 46, log2_Vb=52.2261, log2_N5=43.6110, L_b=392.08, L_max=798.75, P_min=17, P_cut=41.5, P_lat=51.0, P_top=61.0`, 0-LP Stage B `sl=74.0, sb=None`, `B_1=2^{48}` super-tiles with `<= 6` chords/super-tile deterministically):

1. **Shell Power-of-2 Row-Strided Pre-Initialized Arena (`(b_1..b_5)` shell worker, `uint16_t`, `N_{a,\mathrm{sh},\mathrm{pad}} = 2^{46}`):** `(L_max + 2)*2^{47} + 4*N_{a,\mathrm{sh}} = 2^{56.647} B` (`163x < 2^64 B`; holds BOTH `init_a[]` and `init_a_p2[]`, 1-based `ai in [1..N_{a,\mathrm{sh}}]`).
2. **Band-1 / Lower Band-2 (`[2^{41.5}, 2^{51.0})`, `W in {8, 16, 32, 64, 128}`) Working Sets (`<< 2^64 B`/owner):**
   - Band-1 (1 prime `p < 2^{41.5}` + shell arena): `8*N_{a,\mathrm{sh}} + 8*64*p + arena <= 2^{56.66} B << 2^64 B`.
   - Lower Band-2 (1 prime across `W in {8..128}` + shell arena): `8*(N_{a,\mathrm{sh}}+1) + 8.0*max(2^{b_hi}/W_k) + arena = 2^{56.653} B << 2^64 B` (`b2_1p_ram_log2`; `max(2^{b_hi}/W_k) = 2^{51.0}/128 = 2^{44.0}`).
   - Descending-`ai` 1-Pass Builder (worst-case 1 lower Band-2 prime): `8*(N_{a,\mathrm{sh}}+1) + 8.0*2^{44.0} = 2^{48.711} B << 2^64 B` (`b2_build_1p_ram_log2`).
3. **1-Owner Lattice Band-2 Super-Tile (`B_1 = 2^{48}`, `ubfx(b_packed, 0, 48)`, deterministic `<= 6` chords/super-tile, `64 KiB`-Page Streaming Apply):**
   - Because `max |b_5| <= 376.6` (`753` integer values), fixing `b_5[6:0]` leaves **at most `ceil(753/128) = 6` chords** in ANY super-tile (100% deterministic).
   - **Zero-Overflow `64 KiB`-Page Streaming:** During apply, a super-tile worker NEVER buffers its entire bucket in local memory — it streams its linked chain of `64 KiB` bucket pages **sequentially one `64 KiB` page at a time** (double-buffered `2 * 64 KiB = 2^{17} B`) from the global slab pool into its fixed 6-chord score arena. Consequently, 1-owner RAM is **completely independent of per-chord hit skew** (ZERO per-chord density bound and ZERO overflow procedure needed!). During emit, the emitter's 1-owner working set is only `wptr[2^{48}] = 2^{51.0} B` + current `12x12` LLL basis (`< 4 KiB`).
   - Full-`N_a` 6-chord score arena (`Na_tot_pad = 2^{49}`): `6 * (L_max + 2) * 2^{50} + 4*N_a = 2^{62.23} B`; global `2^{51}`-slot open-addressed `HT` (`16 B`/slot, load `0.14`): `2^{55.00} B`; `wptr[2^{48}]`: `2^{51.00} B`; double-buffered page stream: `2^{17} B`.
   - **Total 1-owner deterministic super-tile RAM:** **`2^{62.241} B` (`3.39x < 2^64 B`, LIVE `b0_1own_ram_log2`)**.
4. **1-Owner `Z[omega]` BW LA (`N_sw_cap=2^{43.1936}, N_l_cap=2^{42.6867}, bw=(1024,35), nslots=8, TAB_R_CAP=1.00`):** `2*160*8*N_sw_cap + W_ROW*8*N_sw_cap + 2^{54.69}*8 + 65*2^{53} = 2^{59.611} B` (`21.0x < 2^64 B`).
5. **1-Owner Split Batch (`SPLIT_Q = 2^{46}, SPLIT_REC = 56 B`):** `56 * 2^{46} = 2^{51.807} B << 2^64 B`.
6. **1-Owner Boot Tree (`x_EA = 65.3021`, `boot = [180, 76, 12, 1]`, `ea = [5, 1250]`, conservative `0.30` cap):** `2*P_{x_EA} = 2^{62.094} B << 2^64 B`.

## 2. Global Storage & Sensitivity

- **Global RAM (ZERO Band 3, ZERO `occ_bm[]`):** Band-1 `[2^{17}, 2^{41.5})`: `8*N_a*n_lo + 8*64*sp_lo <= 2^{89.5} B`; Lower Band-2 `[2^{41.5}, 2^{51.0})`: `2^{97.46} B`; Lattice Band-2 shared `64 KiB` on-demand linked-page bucket slab pool (`lat_bkt_bytes = 16 * V * raw_mult * ln(61/51) = 2^{96.87} B` — shared globally so per-chord skew wastes zero extra pool space, plus `<= N_5 * 64 KiB = 2^{59.61} B` worst-case 1-partial-page tail per non-empty tile + `wptr = 2^{51.0} B`) + `HT` (`2^{55.0} B`) + arena (`2^{56.65} B`); **total `2^{98.128} B < 2^{100} B`** (`+1.87`b headroom, LIVE `ram_ok`).
- **Sensitivity (`part="sensitivity"`):** `plain_dickman_rho` gives `2^{84.983881}` (`exp_run 86.342331`); `conservative_lat_b2_plus10` (`c_lat = 78.5` vs `68.5`) gives `2^{85.234029}` (`exp_run 86.592478`); `de_bruijn_tighter_cap` (`cap_scale = 0.80`) gives `2^{85.090748}` (`exp_run 86.449199`); `flat_shared_arena_ref` gives `2^{85.144723}` (`exp_run 86.503172`).
