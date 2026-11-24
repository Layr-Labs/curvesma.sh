# STORAGE.md: Direct Live uint8_t Score Nodes, 1-Hop 16 KiB Link Pages and the < 2^{100} B Storage Bound

Where the Band 1, Lower Band 2 and Lattice Band 2 data lives, how it moves and what moves cost. Every node is a flat 64-bit scalar machine with `< 2^{64}` B of RAM, addressed by unit-cost `ldr/str/ldp/stp/staddb/ldadd` under COST_MODEL.md. Every cross-node byte goes through an explicit 1-hop 16 KiB link page copied by scalar `ldp/stp` (`1.007 <= 1.05` Ir per 8 B word) or, for Band 1 and Lattice Band 2, written directly into the link page by the producer's counted stores (KERNELS.md §1.1, §1.3). Figures are computed in exact rationals by `stor_bound.py` and checked by `audit_bn_cm.py` (`ram_ok`).

## 1. Node Types and Memory

| Node | Count | Holds | RAM (log2 B) |
|---|---:|---|---:|
| Emitter E_0 | 1 | 16-byte mixed-radix `a_tab = (ai, &B_TAB[k][0])` (`prod_Ma = 79,151,674,016,786,875` entries, `2^{60.136}` B), contiguous 4-plane `B_TAB[0..3]` (`4 * prod_Mb = 1.742e18` 8 B entries, `2^{63.596}` B, full `E_{12}+DELTA` box with zero chord splits and exact `S`-fixed handling on plane `k=3`), 2-op Direct Q12 `isqrt` LUT (`2^{37}` 8 B entries = `1 TiB`, `2^{40}` B), wheel/rot/mod3 tables (`rot_tab` 64 MiB + `mod3_tab` 8 MiB + `STEP1_TAB` 128 KiB + `STEP_K`/`wtab`, `< 2^{27}` B), `2^{35}` direct 16 KiB link pages at 32 KiB stride, `ptrs[0..2^{35}-1]`, `E_a+` build arrays | <= 63.723 (emitter_log2: 63.7219; **2.805 EiB headroom**) |
| Score Nodes C[s] | N_sc = 2^{35} | live `uint8_t` (1 byte/cell) score rows of the sub-bands assigned to `s` (`<= avg + 1 sub-band`), `2^{36}` incoming list link pages, 1 incoming emitter link page, `2^{49}` B dummy sink (`> |E_a+|`), `1 B` `init_a[]` + `16 B` `a(ai)` tables (`17 |E_a+|` B), 1.5 MB Stage-B batch | <= 62.394 (score_node_log2: 62.3926; blocks: 22,024,191 < 2^{26}) |
| List Nodes L[k] | N_list = 2^{36} | Band 1 lists (`uint32_t` delta-packed on `[2^{24}, 2^{29.5})` + 8 B on `[2^{29.5}, 2^{41.5})`, mean `<= 2^{55.37}` B/node) and aligned-`bfi` Lower Band 2 lists (all lists at 8x `<= 2^{56.44}` B/node, node max `<= 2^{58.13} < 2^{60}` B pool), build scratch (`<= 2^{56}` B), `2^{35}` outgoing link pages (`2^{49}` B) + 16 B stream state, sub-band/line descriptors | <= 60.176 (list_node_log2: 60.1754; model 60.1747) |

Score rows (truncated `E_{12}`: `Q(a) <= 17/20, Q(b) <= 17/20, B_a = 2, BLK = 38`, `uint8_t` 1 byte/cell). Each canonical `b` with `0 < Q(b) <= 17/20` owns one row; rows are grouped in `W = 64` sub-bands on `(b_1..b_5)` chords with shared stride `2^{k_s}` bytes, `k_s = max(38, ceil(log2(cap)))`, aligned to `2^{38}`-byte blocks (`W * 2^{k_s}` alignment for `W in {8, 16}` guarantees bits `[k_s .. k_s+3] == 0` for aligned-`bfi` Lower-B2). Greedy least-loaded assignment guarantees each Score Node holds `<= avg + 64 * 2^{48}` B (`<= 2^{62.3926}` B, `blocks = 22,024,191 < 2^{25} < 2^{26}`). `B_TAB[k][pack_b] = (s << 29) | row_block` has 35-bit `s` in `[63:29]`, **bits `[28:26] = 000`**, and 26-bit `row_block` in `[25:0]`; max record address `((blocks - 1) << 38) + |E_a+| + 1 < 2^{64}`. Out-of-domain `b` maps to `0` (C[0] dummy sink blocks `[0, 2^{11})`, `2^{49}` B `> |E_a+|`, disjoint from real rows).

Row capacity: `cap(sb) = #{a in E_a+ : Q(a) <= min(17/20, u_sb) + delta} + D_B1 + 16` (`delta = 0.005`, `eps_or <= 0.0003844 <= delta / 2` proved in exact rationals by `fd_bounds.py`, `D_B1 <= 2^{38.578}`). Every write of every tier lands at `< cap(sb)` of its own row; guard cells `ai > N_a(b)` are dropped by one compare before Stage B.

## 2. Data Movement and Charges

Serial deterministic 1-hop 16 KiB dual-ported link-page model with per-link doorbells (`<= 20` Ir per side handoff, identical to `v31b`):
1. **One-time `uint8_t` score initialisation & read-only `scan8_u8`:** `C[s]` copies 1-byte `init_a[]` (`init_a[0] = 0`) into `[0, L_g]` and scans read-only at 8 cells per 64-bit word (`C_SCAN_G8 = 275/1024`, `scan_zero = 0.4279` Ir/cell, `scan_group_overhead = 1.0179 <= 1.05`, `scan8_u8_errs = 0`).
2. **Lattice Band 2 emission (`E_0 -> link_E_to_C[s] -> C[s]`):** 9-op `LDP+LDADD` direct-to-link-page emit (`ldp x0, x6, [x2]; ldr x3, [x6, x5, lsl #3]; add x8, x17, x3, lsr #26; add x7, x0, x3, lsl #38; ldadd x1, x9, [x8]; str x7, [x9]; tbnz x9, #14`) + 8-op `emit_direct_flush` (`0.014663 <= 0.015` Ir/rec) + x128-unrolled `STADDB` apply (`1.52443 <= 1.55 = C_LAT_APPLY`), all executed instruction by instruction by `v25_verify.py lat_apply_e2e` across all 4 contiguous planes at full production width (`ai` up to `2^{48}-1`, `rb = 22,015,925`). 248 log-p bands on `[2^{45.5}, 2^{61.0}]`; band-end partial pages charged at `4,130` Ir in `misc`.
3. **Band 1 and Lower Band 2 streams (`L[k] -> link_L_to_C[k][s] -> C[s]`):** Band 1 writes 32-bit delta segments on `[2^{24}, 2^{29.5})` (`delta_b1_copy/apply`, `2.82 + 24*P_ov`) and 8 B records `(row_block << 38) + ai` on `[2^{29.5}, 2^{41.5})` (`C_HIT_LO = 3.5699`, `v25_verify.py` check 2) directly into `link_L_to_C[k][s]`. Lower Band 2 streams `(ai << 4) | j_rel` words at 1.05 Ir/8 B inside `C_B2_INT_X4 = 5.30` and inline headers inside `C_SB_ENTER_X4 = 20.05`. Partial pages charged in `misc` (`misc <= 2^{72.504300}`).
4. **Setup broadcasts & precomputation:** `ins_pre` explicitly charges `4 * prod_Ma` (16-B `a_tab`), `2 * (4 * prod_Mb)` (4-plane `B_TAB`), `16 * 2^{37}` (1 TiB Direct Q12 `isqrt` LUT), and `16 * 2^{27}` (`rot_tab`/`mod3_tab`) alongside `init_a[]`/`a(ai)`/`L[k]` broadcasts (`lists <= 2^{78.389935}`).

## 3. Other Single-Worker Working Sets (< 2^{64} B)

1. List builders on `L[k]`: Lower-B2 1-prime build `<= 2^{52.450}` B, Band-1 build `<= 2^{51.905}` B (rounded to `2^{56}` B in `stor_bound.py`; with shell arena `<= 2^{57.566}` B).
2. Stage-B `P_min = 24` `g_rec` SRT: `1.5` MB.
3. 1-owner 16-vec `Z[omega]` BW LA (`bw = (768, 32)`): `<= 2^{63.1574}` B `< 2^{64}` B (`2^{59.3854}` B SpMV + `2^{63.0477}` B PM-Basis workspace).
4. 1-owner deterministic `ntt3_F` Bernstein split (`Q = 2^{46}`): `<= 2^{57.0}` B (`SP` prime-tree level `<= 2^{56.0}` B + batch PT/RT workspace `32 * 29 * 2^{46} <= 2^{55.9}` B; rounded up to `2^{58}` B in `stor_bound.py`).
5. 1-owner boot tree (`boot = [185, 80, 11, 1], ea = [5, 1150], x_EA = 65.5625`, 0.30 cap): `<= 2^{62.3543}` B (wave `<= 2^{53.0}` B).

## 4. Total Machine Storage (< 2^{100} B, pointwise) and Sensitivity

- Exact-rational `stor_bound.py` certificate: `|E_a+| <= 2^{47.6022}`, sub-bands `<= 2^{46.0185}`, `uint8_t` score rows `<= 2^{97.3870}` B; per Score Node `<= 2^{62.3926}` B (`blocks = 22,024,191 < 2^{26}`); `E_0 <= 2^{63.7219}` B (**2.805 EiB headroom**); `L[k] <= 2^{60.1754}` B; **total machine storage `<= 2^{97.9088}` B `< 2^{100}` B** (`ram_log2 = 97.9087`, `storage_cap_abort = 0.0`).
- Sensitivity (`part="sensitivity"`): `plain_dickman_rho` `2^{83.151806}` (`exp_run 84.510257`); `conservative_lat_b2_plus10` (`C_LAT_LEAF + 10` Ir) `2^{83.772528}` (`85.130979`); `conservative_srt_030_cap` `2^{83.317378}` (`84.675829`); `de_bruijn_tighter_cap` (`cap_scale = 0.80`) `2^{83.265168}` (`84.623619`); `flat_shared_arena_ref` `2^{83.266053}` (`84.624504`). `dag_survivors` (2x envelope): `2^{83.355175}`.
