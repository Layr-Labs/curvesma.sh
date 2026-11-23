# Kernels (Raw-Opcode-Verified Baseline Scalar ARMv8.0-A: `Arch2-Prod-ZeroSlice-Buf64K-B2-U16-Streamed`, 3-Prime Radix-4 Lazy Shoup NTT, 1-Dense + 16-Proj + 16-Shard 4-Vec LA, Generous Scalar Cold-Path Bounds)

Every hot-loop kernel cost is derived from explicit baseline scalar ARMv8.0-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero `bl` — all kernels are inline macros or leaf routines returning via `ret`). Cold paths (table build, `.Lhit64` hit-block handler, and 16D descent DFS / `BOUND_CELL`) use generous, manifestly-conservative scalar ARMv8.0-A operation bounds rather than fragile hand-scheduled assembly. `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file (`29` hot-loop blocks, `0` undefined `.L...` labels, per-kernel disjoint physical register maps asserting ZERO symbolic/numbered collisions in every single block), validates opcodes & `LDP`/`STP` addressing, counts every instruction (`cmp`+`b.*`=2 ops, loop control included), simulates Band-1a/1b (`rev_b1_stream_u16` `1.125` Ir/word) + sequential incrementing-pointer `64 KiB`-buffered Band-2 (`W_SB = 16`, `oid < 2^56`, `33/4 = 8.25 <= 10.0 = c_b2_row_k0` Ir/prime/row 4-prime orbit incremental row-step `rev_b2_row_k0_step4` + `19`-op `rev_sb_probe_fast` + `10`-op `rev_b2_direct_sb_hdr` + `15`-op prime-shard `rev_b2_emit_1` writing `[x7], #16` with `tst x7, #0xffff` full-buffer drain + `6`-op slab-owner `rev_b2_recv_apply_1` reading `[x_rx], #16` to exact `x_rx_end` + `<= 16/4096 = 0.004` Ir/hit flush wrapper = `21.004 <= 22.0` Ir/hit, plus `28` row-end partial drains/row) + `cbz`-guarded 64-cell `ORR` scan (`47/64` incl. entry guard and explicit `.Lscan_done: ret` before `.Lhit64:` + blunt `1,000` Ir/hit-block cold bound `<= 0.84`) + 4-vec/1-vec `Z[zeta_12]` LA + AAPCS64 `MONT5 = 280` (`279` opcodes) + Barrett `row_red5_z12` + 3-prime Radix-4/3 Lazy Shoup NTT + 4-prime Garner CRT + exact Python verification of counting-sort build, 4D Fincke-Pohst DFS vs brute force with ambient sign canonicalization `v -> -v`, 55-bit `SQ_LUT`, 16D half-ball bijective ranking, `+6` outward-rounded `bmul_6`, and input-uniform Bareiss LLL bounds, asserting `opcodes/items <= charged`.

## 1. Collection Kernels (`Arch2-Prod-ZeroSlice-Buf64K-B2-U16-Streamed`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 448 = 7 * 64`, `W_SB = 16`, `P_min = 2^32`, `P_m4s = 2^52.0`, `P_m4 = 2^57.5`, `P_lo = P_f = 2^58.0`, `P_g = 2^63.5`)

Because the sieving region is a pure product of two 8D balls (`qa(a) <= s_c R2_a` with 2-adic radii `s_c = (1.0, 1.5, 2.0)` and `qb(b) <= R2_b`, `A^2 = 321,500`, `r_ab = 0.83`), every canonical `a`-table entry (`N_a^can = V_a = 2^55.8744 < 0.92 * 2^56`, `N_a^tbl = 4 V_a = 2^57.8744`) is unconditionally in-region for every `b` (`V_b = 2^64.5535`, `L_max = 440.08 < W_tile = 448`, `L_b = 189.02`, `N_7 = 2^54.9910`, `V_kept = 2^120.4279`). There is **ZERO `q_key`, ZERO `q_tab`, ZERO `sb_min_q8`, ZERO binary search, ZERO global sorting invariant, and ZERO `f`-side Band 2** (`P_f = P_lo = 2^58.0`).
- **Bias-initialized `S_g` & zero-wrap `uint16_t` accumulators (`S_g, S_f`) with exact `2^56` bank stride (`K_a = 1`, `W_SB = 16`):** Quarter-bit weights on `<= 705`b norms reach `s <= 3,000`. Initializing `S_g` lanes to `bias16 = (32768 - T) & 0x7FFF` (`1 <= T <= 3,000`) keeps `(32768 - T) + s in [29768, 35767] subset [0, 65535]` with **ZERO 16-bit wraparound/carry between lanes** and `s >= T <==> bit 15 (0x8000) == 1`. Each column-slab owner owns `16` consecutive columns `[16*S, 16*S+15]` (`28` slab owners cover `W_tile = 448 = 28 * 16`) and allocates `uint16_t S_g[16][4][2^56]` (`2^63.0` B, indexed by `62b` `(u << 58) | (m' << 56) | oid`) + 1-column `uint16_t S_f[4][1][2^56]` (`2^59.0` B) = **`2^63.087` B `< 2^64` B** (exact power-of-2 bank stride `2^56`, `oid < 2^56`).
- **Zero-copy 8B records with explicit cross-owner streaming across Band-1a/1b and sequential incrementing-pointer `64 KiB`-buffered `(x4, w_lp)` transfer & drain in Band-2:**
  - Band-1a-long (`[2^32, 2^52.0)`): `m`-partitioned `oid < 2^56` -> explicit 16-rec burst stream `rev_b1_stream_u16` (`18/16 = 1.125` Ir/hit) + 16x-unrolled `rev_b1a_apply_u16` (`58/16 = 3.625` Ir/hit) = **`4.750 <= 4.78` Ir/hit**, with 4-bank wrap/tail/setup overhead (`4` banks x `<= 86.75 + 9.25 = 96.0` = **`384.0` ops/prime/col** across ALL `k >= 0`) charged inside `b1a_long_lk = 384.0`.
  - Band-1a-short (`[2^52.0, 2^57.5)`): guarded 1-record stream `rev_b1_short_entry_guard` + `rev_b1_stream_tail1` (`4.0`) + `rev_b1a_apply_1` (`6.0`) = **`10.0` Ir/hit** (`b1a_short_lk = 40.0` = `4 * 10.0` ops/prime/col across 4 banks incl. `cbz` `k = 0` entry guards!).
  - Band-1b (`[2^57.5, 2^58.0)`): guarded 1-record stream `rev_b1_short_entry_guard` + `rev_b1_stream_tail1` (`4.0`) + 8B `(m << 56) | oid` `rev_b1b_apply_1` (`8.0`) = **`12.0` Ir/hit** (`b1b_lk = 12.0` ops/bkt incl. `cbz` `k = 0` entry guards!).
  - Band-2 (`g`-side `[2^58.0, 2^63.5)` only — **Sequential Incrementing-Pointer `64 KiB`-Buffered `(x4, w_lp)` Transfer & Exact-Bound Drain**): Each Band-2 prime orbit is assigned its own dedicated build/replay owner (`~2^56.05` owners; each builds its single prime locally at peak `<= 14.981 EiB < 16.00 EiB` with ZERO cross-owner build transfers, frees `kbuf`/`pos`, and holds `<= 7.491 EiB` max / `4.696 EiB` mean resident + `416 B` incremental row-step state). For each `D_(omega, sigma)` row, `rev_b2_row_k0_step4` (`33/4 = 8.25` Ir/prime/row + `<= 1.354` outer-chord reset <= `10.0 = c_b2_row_k0` Ir/prime/row, priced explicitly in `audit_cm.py`) produces the exact `0 <= x_k0 < x_p` entry inputs for all `4` orbit primes before `rev_sb_probe_fast`. Under sequential single-threaded row execution across the prime owners (one active emitter at a time -> ZERO multi-producer race or tearing), each of the `28` slab owners maintains a `64 KiB`-aligned `64 KiB` (`4,096`-record, `16 B/record`) receive buffer with live write cursor `x_wr[S]` in the row controller's `224 B` table `x_tx_tbl[0..27]`. The active prime owner runs `rev_b2_emit_1` (`15` ops/hit: post-index `stp x4, x_wlp, [x7], #16`, saves `x7`, checks `tst x7, #0xffff; b.eq .Lb2_flush_slab`); when full (or once per slab at row end — only `28` partial drains/row, ZERO `w_lp` class-boundary flushes and ZERO `(0,0)` EOF sentinels), the slab owner executes `rev_b2_recv_apply_1` (`6` ops/hit: post-index `ldp x4, x5, [x_rx], #16` up to exact `x_rx_end`). Total = `15 + 6 + 16/4096 = 21.004 <= 22.0 = C_HIT_HI` Ir/hit (`+0.996` Ir/hit slack), plus `28 * 32 / (N_a_tbl * L_b)` row-end flush overhead priced explicitly in `audit_cm.py`.

### 1.1a `rev_b2_row_k0_step4` (Band-2 4-prime orbit incremental `D_(omega, sigma)` row-step `k_0` init: `33/4 = 8.25` Ir/prime/row + `<= 1.354` outer-chord reset <= `C_B2_ROW_K0 = 10.0` Ir/prime/row)

Because `sigma(omega) = omega` (`7*5 === 4*5 === 13*5 === 5 mod 15`), all `4` primes `s_k = s_0^(7^k)` (`k = 0..3`) in a Band-2 prime orbit share the EXACT same modulus `p` (`2^58 <= p < 2^63.5`), the EXACT same `b_0`-coefficient `pwk[0] = 1`, and the EXACT same `b_1`-coefficient `pwk[1] = w = s_0^5 mod p`. Moreover, the `D_(omega, sigma)` fundamental-domain predicate `in_D12(b)` depends ONLY on `(b_2..b_7)`, so for every fixed `(b_2..b_7)`, `b_1` sweeps an uninterrupted contiguous integer chord of mean length `L_b1 = (96/105) * sqrt(R2_b / 6) = 232.30` rows (`> L_b = 189.02`). Maintaining `k_tail[0..3]` (`32 B` in `[x_st]`) across consecutive `b_1 -> b_1 + 1` rows and converting signed `x_b0` (`|x_b0| <= L_max = 440.08 << 2^58 <= p`, precondition `-p < x_b0 < p`) ONCE per 4-prime orbit computes overflow-proof `0 <= x10..x13 < x_p` (`k_0[0..3]`) in `33/4 = 8.25` Ir/prime/row (disjoint GPR map `x_st=x19, x_p=x20, x_d1=x21, x_pms=x22, x_b0=x23`, zero labels, verified bit-exact at `p = 2^64 - 59` and on `Q(zeta_15)` 4-prime orbits in `cm_models.py`). Outer `(b_2..b_7)` chord resets (`<= 25` ops/prime per chord: shared `(b1_lo - 1)*w mod p` + `4` modular adds from `kt_lvl[2..7][0..3]`, `416 B` total per owner) are conservatively bounded at `256 / 189.0 = 1.354` Ir/prime/row -> `8.25 + 1.354 = 9.604 <= 10.0 = C_B2_ROW_K0` Ir/prime/row:

```asm
; kernel: rev_b2_row_k0_step4 items=4 charged=10.0 key=c_b2_row_k0
add x0, x_b0, x_p
cmp x_b0, #0
csel x0, x0, x_b0, lt
sub x1, x_p, x0
ldp x2, x3, [x_st]
ldp x4, x5, [x_st, #16]
subs x6, x2, x_pms
add x7, x2, x_d1
csel x2, x6, x7, hs
subs x6, x3, x_pms
add x7, x3, x_d1
csel x3, x6, x7, hs
subs x6, x4, x_pms
add x7, x4, x_d1
csel x4, x6, x7, hs
subs x6, x5, x_pms
add x7, x5, x_d1
csel x5, x6, x7, hs
stp x2, x3, [x_st]
stp x4, x5, [x_st, #16]
subs x6, x2, x1
add x7, x2, x0
csel x10, x6, x7, hs
subs x6, x3, x1
add x7, x3, x0
csel x11, x6, x7, hs
subs x6, x4, x1
add x7, x4, x0
csel x12, x6, x7, hs
subs x6, x5, x1
add x7, x5, x0
csel x13, x6, x7, hs
ret
```

### 1.1 `rev_sb_probe_fast` (Band-2 row start & wraparound: `19 <= C_TILE_HDR = 20.0` Ir/row — internal custom-ABI leaf-frame helper whose row/tile caller saves callee-saved `x19..x25` once per row/tile, zero inner-loop overhead; precondition `0 <= x_k0 < x_p` supplied by `rev_b2_row_k0_step4`; disjoint GPR map `x_k0=x7, x_p=x8, x_lb=x9, x_col0=x10, x_krot56=x11, x_off=x12, x_rot=x23`; exact `L_b == 0` guard, `sbi_start`, `x_rot`, signed starting slab `x22 = ((sbi_start<<4) + x_col0 - x_k0) >> 4`, `n_sb`, and `x19 = &off[sbi_start]` for BOTH normal `x_col0 = 0` and wraparound `x_col0 = p - k_0`)

```asm
; kernel: rev_sb_probe_fast items=1 charged=20.0 key=c_tile_hdr
mov x21, #0
cbz x_lb, .Lempty_chord
lsr x4, x_k0, #4
sub x1, x_col0, x_k0
and x0, x1, #15
orr x_rot, x_krot56, x0, lsl #60
lsl x22, x4, #4
add x22, x22, x1
asr x22, x22, #4
sub x5, x_p, x_k0
cmp x5, x_lb
csel x5, x_lb, x5, hi
add x6, x_k0, x5
sub x6, x6, #1
lsr x6, x6, #4
sub x21, x6, x4
add x21, x21, #1
add x19, x_off, x4, lsl #3
.Lempty_chord:
ret
```

### 1.2 `rev_b2_direct_sb_hdr` (Band-2 super-bucket header with `cbz x21` empty-chord guard: `10 <= C_SB_HDR = 10.0` Ir/SB)

```asm
; kernel: rev_b2_direct_sb_hdr items=1 charged=10.0 key=c_sb_hdr
cbz x21, .Lb2_sb_done
.Lsb:
ldp x5, x6, [x19]
add x19, x19, #8
cmp x5, x6
b.eq .Lnext_sb
add x2, x20, x5
add x3, x20, x6
; --- inline rev_b2_emit_1 (.Lb2_emit) runs on prime-shard owner for [x2, x3), preserving x19..x22 ---
.Lnext_sb:
add x22, x22, #1
subs x21, x21, #1
b.ne .Lsb
.Lb2_sb_done:
```

### 1.3 Band-1a-Long 16x-Unrolled + Band-1a-Short 1-Record Apply (`rev_b1a_apply_u16`: `58/16 = 3.625` Ir/hit + `rev_b1_stream_u16` `1.125` = `4.750 <= C_HIT_LO_M4 = 4.78` Ir/hit; `rev_b1a_apply_1`: `6.0` + `4.0` stream = `10.0` Ir/hit)

Flag-safe pre/post wrap (`4` ops once per bucket) + 16-hit unrolled body (`58` ops) + guarded 1-record loop (`6` ops, with explicit `cbz` `k = 0` entry guard `rev_b1_short_entry_guard`), paired with explicit cross-owner bucket streaming `rev_b1_stream_wrap16` + `rev_b1_stream_u16` + `rev_b1_stream_tail1` (§1.3b). Combined worst-case overhead above `4.750` Ir/hit for ANY bucket length `k >= 0` is `<= 39.625 + 47.125 = 86.75 <= 96.0` ops/bkt (`>= 9.25` ops/bank setup margin even at `k = 15`; 4 banks x `96.0` = `384.0` ops/prime/col charged in `b1a_long_lk = 384.0`, simulated flag-exact across ALL `k = 0..63` in `audit_cm.py` with `0` overreads!):

```asm
; kernel: rev_b1a_apply_wrap16 items=1 charged=4.0 key=c_hit_lo_m4_wrap
subs x_cnt, x_cnt, #16
b.lo .Ltail_fix16
; --- inline rev_b1a_apply_u16 (.Lb1a16) runs here while C == 1 ---
.Ltail_fix16:
adds x_cnt, x_cnt, #16
b.eq .Ldone
; --- inline rev_b1a_apply_1 (.Lb1a1) runs here while x_cnt != 0 ---
.Ldone:
```

```asm
; kernel: rev_b1_short_entry_guard items=1 charged=2.0 key=c_b1_short_guard
cbz x_cnt, .Lb1_empty
; --- inline rev_b1a_apply_1 / rev_b1_stream_tail1 / rev_b1b_apply_1 runs here only when x_cnt > 0 ---
.Lb1_empty:
```

```asm
; kernel: rev_b1a_apply_u16 items=16 charged=3.65 key=c_hit_lo_m4
.Lb1a16:
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
subs x_cnt, x_cnt, #16
b.hs .Lb1a16
```

```asm
; kernel: rev_b1a_apply_1 items=1 charged=6.0 key=c_hit_lo_m4_s
.Lb1a1:
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1a1
```

### 1.3b Explicit Cross-Owner 16-Record Burst Stream (`rev_b1_stream_u16`: `18/16 = 1.125 <= 1.15` Ir/word; `rev_b1_stream_tail1`: `4.0` Ir/word)

Because each Band-1a/1b/2 prime's read-only replay table (`<= 2^63.90 B < 2^64 B`) resides in its prime-shard owner while the `28` column-slab accumulators (`S_g[16][4][2^56] + S_f[4][1][2^56] = 2^63.087 B < 2^64 B` each) reside in the `28` slab owners (`STORAGE.md` §1), every active Band-1a-long bucket slice `[off[k_t], off[k_t+1])` (and every LA sparse-matrix shard row in §4, aggregate mean `<= 213.644` 8B words/row incl. 64B 24-bucket header & packing slack!) is streamed in `128B` (`16`-word) bursts via `rev_b1_stream_wrap16` (`4` ops) + `rev_b1_stream_u16` (`18/16 = 1.125 <= 1.15` Ir/word) + `rev_b1_stream_tail1` (`4.0` Ir/rem word; total overhead `<= 4 + 15*(4 - 1.125) = 47.125 <= 48` ops/bkt, 4 banks x `<= 96.0` charged inside `b1a_long_lk = 384.0`, simulated flag-exact across `k = 0..63` in `verify_kernels_symbolic` with `0` overreads!). On short Band-1a-short and Band-1b buckets, `rev_b1_short_entry_guard` (`cbz` `k = 0` guard) + direct 1-record `rev_b1_stream_tail1` (`4.0` Ir/hit) streams directly with zero unroll overhead:

```asm
; kernel: rev_b1_stream_wrap16 items=1 charged=4.0 key=c_b1_stream_wrap
subs x_cnt, x_cnt, #16
b.lo .Lb1s_fix16
; --- inline rev_b1_stream_u16 (.Lb1s16) runs here while C == 1 ---
.Lb1s_fix16:
adds x_cnt, x_cnt, #16
b.eq .Lb1s_done
; --- inline rev_b1_stream_tail1 (.Lb1s1) runs here while x_cnt != 0 ---
.Lb1s_done:
```

```asm
; kernel: rev_b1_stream_u16 items=16 charged=1.15 key=c_b1_stream
.Lb1s16:
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
ldp x2, x3, [x_src], #16
stp x2, x3, [x_dst], #16
ldp x4, x5, [x_src], #16
stp x4, x5, [x_dst], #16
ldp x6, x7, [x_src], #16
stp x6, x7, [x_dst], #16
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
ldp x2, x3, [x_src], #16
stp x2, x3, [x_dst], #16
ldp x4, x5, [x_src], #16
stp x4, x5, [x_dst], #16
ldp x6, x7, [x_src], #16
stp x6, x7, [x_dst], #16
subs x_cnt, x_cnt, #16
b.hs .Lb1s16
```

```asm
; kernel: rev_b1_stream_tail1 items=1 charged=4.0 key=c_b1_stream_tail
.Lb1s1:
ldr x0, [x_src], #8
str x0, [x_dst], #8
subs x_cnt, x_cnt, #1
b.ne .Lb1s1
```

### 1.4 Band-1b Direct 1-Record 8B `{m:2, oid:56}` Apply (`rev_b1b_apply_1`: `8.0` + `rev_b1_stream_tail1` `4.0` = `12.0 <= C_HIT_LO_U = 12.0` Ir/hit, `b1b_lk = 12.0` ops/bkt)

Every Band-1b hit packs `(m << 56) | oid` into one 64-bit word, streamed via `rev_b1_stream_tail1` (`4.0` Ir/hit) and rotating bank `m -> (m - k) & 3` via `sub x4, x0, x_k56; and x4, x4, x_mask58` directly into the `2^56`-stride accumulator (`8.0 + 4.0 = 12.0` Ir/hit):

```asm
; kernel: rev_b1b_apply_1 items=1 charged=8.0 key=c_hit_lo_u
.Lb1b1:
ldr x0, [x_p], #8
sub x4, x0, x_k56
and x4, x4, x_mask58
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1b1
```

### 1.5 Sequential Incrementing-Pointer `64 KiB`-Buffered Band-2 `(x4, w_lp)` Transfer & Exact-Bound Drain (`15.0 + 6.0 + 0.004 = 21.004 <= 22.0 = C_HIT_HI` Ir/hit)

Under sequential single-threaded row execution across the `~2^56.05` 1-prime-orbit Band-2 build/replay owners (`<= 14.981 EiB` peak build, `<= 7.491 EiB` max / `4.696 EiB` mean resident, ZERO cross-owner build transfers, processed one owner at a time so only one emitter is active at any instant — ZERO multi-producer race or `stp`/`ldp` tearing), each of the `28` slab owners maintains a `64 KiB`-aligned `64 KiB` (`4,096`-record, `16 B/record`) receive buffer `[buf_base[S], buf_base[S] + 65536)` with live write cursor `x_wr[S]` stored in the row controller's `28 * 8 = 224 B` cursor table `x_tx_tbl[0..27]`. The active prime-shard executes `rev_b2_emit_1` (`15 <= 15.0` Ir/hit; disjoint GPR map `x_rot=x23, x_mask58=x24, x_tx_tbl=x25, x_wlp=x26`): computes `x4` and `x6 = x22 + carry`, skips `x6 >= 28`, loads `x7 = x_tx_tbl[x6]`, writes `stp x4, x_wlp, [x7], #16` (advancing `x7`), stores updated `x7` back, and tests `tst x7, #0xffff; b.eq .Lb2_flush_slab`. When a `64 KiB` buffer fills (`4,096` records -> `rev_b2_flush_wrap` using disjoint GPR map `x_rx=x10, x_rx_end=x11, x_tx_tbl=x25` matching `cm_models.py:532`, completely disjoint from `rev_b2_direct_sb_hdr`'s live `x19..x22`; `<= 16 / 4096 = 0.0039` Ir/hit incl. `x_acc` handoff load inside `charged = 0.01` = `40.96` ops/flush) OR once per slab at row end (`28` partial drains/row, skipping empty buffers, priced explicitly as `28 * 32 / (N_a_tbl * L_b)` in `audit_cm.py`; ZERO `w_lp` class-boundary flushes because every record carries its own `w_lp`, and ZERO `(0, 0)` EOF sentinel needed because `x_rx_end` is exact), the slab owner executes `rev_b2_recv_apply_1` (`6 <= 6.0` Ir/hit; separate slab-owner GPR map `x_rx=x20, x_rx_end=x21, x_acc=x22`). Total = **`15.0 + 6.0 + 0.004 = 21.004 <= 22.0 = C_HIT_HI`** Ir/hit (`+0.996` Ir/hit headroom). Out-of-chord padded columns `[L_b, 448)` accumulate harmlessly inside `S_g[16][4][2^56]` (priced via `cand_bd = 1 + 16/L_b`), while `rev_scan_zero_u64` scans active columns `t < L_b`:

```asm
; kernel: rev_b2_emit_1 items=1 charged=15.0 key=c_b2_emit
.Lb2_emit:
ldr x0, [x2], #8
adds x4, x0, x_rot
adc x6, x22, xzr
lsr x9, x4, #60
and x4, x4, x_mask58
orr x4, x4, x9, lsl #58
cmp x6, #28
b.hs .Lb2_skip
ldr x7, [x_tx_tbl, x6, lsl #3]
stp x4, x_wlp, [x7], #16
str x7, [x_tx_tbl, x6, lsl #3]
tst x7, #0xffff
b.eq .Lb2_flush_slab
.Lb2_skip:
cmp x2, x3
b.lo .Lb2_emit
```

```asm
; kernel: rev_b2_flush_wrap items=4096 charged=0.01 key=c_b2_flush_wrap
; GPR map (disjoint from rev_b2_direct_sb_hdr x19..x22): x_rx=x10, x_rx_end=x11, x_tx_tbl=x25
.Lb2_flush_slab:
sub x_rx, x7, #16, lsl #12
mov x_rx_end, x7
str x_rx, [x_tx_tbl, x6, lsl #3]
; --- hand [x_rx, x_rx_end) + slab owner x_acc base load to slab owner x6 running rev_b2_recv_apply_1 (.Lb2_recv .. .Lb2_rx_done); <= 16 ops/flush = 0.004 Ir/hit (charged = 0.01 = 40.96 ops/flush) ---
b .Lb2_skip
```

```asm
; kernel: rev_b2_recv_apply_1 items=1 charged=6.0 key=c_b2_recv_apply
.Lb2_recv:
ldp x4, x5, [x_rx], #16
ldrh w9, [x_acc, x4, lsl #1]
add w9, w9, w5
strh w9, [x_acc, x4, lsl #1]
cmp x_rx, x_rx_end
b.lo .Lb2_recv
.Lb2_rx_done:
```

### 1.6 Bias-Initialized 64-Cell `uint16_t` `ORR` Scan + Reset (`rev_scan_zero_u64`: `47/64 = 0.734375 <= 0.74` incl. `cbz` Entry Guard & Explicit `.Lscan_done: ret` + Blunt `1,000`-Ir/Hit-Block Cold-Path Bound `<= C_SCAN_ZERO = 0.84`)

Uses `x0..x15` for the 16 loaded `S_g` words (`64` `uint16_t` cells), `x16, x17, x19, x20` for the 15-op `orr` tree (preserving all 64 `S_g` scores on entry to `.Lhit64`), and `x_g, x_f, x_b16, x_8000, x_cnt` in `x21..x25` (`<= 28` GPRs, ZERO `x18/x29/x30`, signed-offset `ldp`/`stp` only; internal custom-ABI leaf-frame helper whose column/tile caller saves callee-saved `x19..x25` once per column/tile, zero inner-loop overhead). Canonical `a`-table entries are partitioned strictly by their static 2-adic class **`c(a) in {0, 1, 2}`** (`3` contiguous partitions, each padded to an exact multiple of 64 entries so `x_cnt % 64 == 0` is invariant) with ONE constant `g`-side threshold **`T_g[c] = (299, 287, 275)` bits** (`bias16 = (32768 - 4*T_g[c]) & 0x7FFF`, ZERO `Q_CM(a)` radial shells, ZERO `b`-norm predictor!). `cbz x_cnt, .Lscan_done` guards `x_cnt == 0` on entry, and explicit `.Lscan_done: ret` terminates the loop cleanly before the out-of-line `.Lhit64:` cold handler (zero fall-through). When `tst x16, x_8000` fires (probability `<= 64 * sga_mean ~ 1.7e-3` per 64-cell block), `.Lhit64` runs a plain scalar C cold-path handler that checks each flagged lane in `x0..x15` against constant `S_f >= T_f_const = 118` bits (`472` quarter-bits), emits global `oid`, and returns to **`.Lscan_zero_sf`** (so `S_f` zeroing across the 64 `uint16_t` cells / 128 B, `add x_f, x_f, #128`, `subs x_cnt, x_cnt, #64`, and exit on the final 64-cell block ALWAYS execute!). Conservatively charging all `47` opcodes (even though `cbz` and `ret` execute once per column) plus a blunt **`1,000` Ir per hit-block** cold bound gives `scan_tot = 0.734375 * (1 + 16/189.0) + sga_mean * 1000.0 = 0.8231 <= 0.84 = C_SCAN_ZERO`:

```asm
; kernel: rev_scan_zero_u64 items=64 charged=0.74 key=c_scan_zero_raw
cbz x_cnt, .Lscan_done
.Lscan64:
ldp x0, x1, [x_g]
ldp x2, x3, [x_g, #16]
ldp x4, x5, [x_g, #32]
ldp x6, x7, [x_g, #48]
ldp x8, x9, [x_g, #64]
ldp x10, x11, [x_g, #80]
ldp x12, x13, [x_g, #96]
ldp x14, x15, [x_g, #112]
stp x_b16, x_b16, [x_g]
stp x_b16, x_b16, [x_g, #16]
stp x_b16, x_b16, [x_g, #32]
stp x_b16, x_b16, [x_g, #48]
stp x_b16, x_b16, [x_g, #64]
stp x_b16, x_b16, [x_g, #80]
stp x_b16, x_b16, [x_g, #96]
stp x_b16, x_b16, [x_g, #112]
add x_g, x_g, #128
orr x16, x0, x1
orr x17, x2, x3
orr x16, x16, x17
orr x17, x4, x5
orr x19, x6, x7
orr x17, x17, x19
orr x16, x16, x17
orr x19, x8, x9
orr x20, x10, x11
orr x19, x19, x20
orr x20, x12, x13
orr x17, x14, x15
orr x20, x20, x17
orr x19, x19, x20
orr x16, x16, x19
tst x16, x_8000
b.ne .Lhit64
.Lscan_zero_sf:
stp xzr, xzr, [x_f]
stp xzr, xzr, [x_f, #16]
stp xzr, xzr, [x_f, #32]
stp xzr, xzr, [x_f, #48]
stp xzr, xzr, [x_f, #64]
stp xzr, xzr, [x_f, #80]
stp xzr, xzr, [x_f, #96]
stp xzr, xzr, [x_f, #112]
add x_f, x_f, #128
subs x_cnt, x_cnt, #64
b.ne .Lscan64
.Lscan_done:
ret
.Lhit64:
; --- out-of-line cold-path scalar C handler bounded conservatively by 1,000 Ir/hit-block (processes 1..64 hits + emits global oid + branches back to .Lscan_zero_sf so S_f zeroing, x_f+=128, x_cnt-=64, and final-block exit ALWAYS execute) ---
```

### 1.7 1-Owner Flat-Array Counting-Sort / Prefix-Sum / Scatter Build (Generous Scalar Bounds: `C_LIST_A = 32.0`, `C_LBKT = 10.0`)

Replay-table construction is `3.15%` of total work (`~2^105.3` u), so rather than hand-scheduling build assembly we charge generous scalar ARMv8.0-A bounds: **`C_LIST_A = 32.0` Ir/entry** (2-pass overflow-proof modular step `p = 2^64 - 59` + histogram increment + scatter of 8B `oid` or `(r << 60) | (m << 56) | oid`, vs `<= 22.5` raw ops -> `+9.5` Ir/entry slack!) and **`C_LBKT = 10.0` Ir/bucket** (`+1.0` histogram zero-init + all 4 per-bank terminal sentinels charged in `bld`, vs `6.0` raw ops -> `+4.0` Ir/bucket slack!). Exact counting-sort orbit `0..V_a-1` mapping and overflow-proof 64-bit modular stepping at `p = 2^64 - 59` are verified LIVE in `cm_models.py`, and 1-owner build RAM fits strictly `< 2^64 B` (`<= 2^63.905 B` for ALL primes, asserted LIVE in `audit_cm.py`).

## 2. Exact 5-Limb Baseline Scalar ARMv8-A CIOS Montgomery Multiplication (`MONT5 = 280.0`, `279` Opcodes, AAPCS64-Compliant)

One complete carry-exact, shift-included two-chain CIOS outer round executes **`48`** baseline scalar ARMv8-A opcodes (`5 x 48 = 240` across 5 unrolled rounds) + **`39`** AAPCS64 callee-saved `x19..x27` prologue/epilogue (`sub sp, sp, #80` / `add sp, sp, #80` explicit frame allocation; exact fixed argument register assignment **`x_a = x5, x_b = x6, x_n0 = x7, x_out = x8, x_p = x9`**, so `x_a = x5` and `x_p = x9` are consumed before the 5 rounds reuse `x5` and `x9` as scratch while `x_b = x6, x_n0 = x7, x_out = x8` remain untouched!) + load/zero-init/final-sub/store opcodes = **`279 <= 280.0 = MONT5`** Ir (inlined macro without `ret` uses `278` ops, or `279` ops as an internal leaf routine; bit-exactness of this exact 279-opcode register machine is verified LIVE on 32 random 320-bit `(a, b, p)` triples inside `audit_cm.py`):

```asm
; kernel: mont5_cios_round items=1 charged=48.0 key=mont5_round
ldr x5, [x_b], #8
mul x10, x0, x5
mul x11, x1, x5
mul x12, x2, x5
mul x13, x3, x5
mul x14, x4, x5
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x27, x27, x13
adcs x19, x19, x14
adcs x20, x20, xzr
adc x21, xzr, xzr
umulh x10, x0, x5
umulh x11, x1, x5
umulh x12, x2, x5
umulh x13, x3, x5
umulh x14, x4, x5
adds x16, x16, x10
adcs x17, x17, x11
adcs x27, x27, x12
adcs x19, x19, x13
adcs x20, x20, x14
adc x21, x21, xzr
mul x9, x15, x_n0
mul x10, x9, x22
mul x11, x9, x23
mul x12, x9, x24
mul x13, x9, x25
mul x14, x9, x26
adds xzr, x15, x10
adcs x15, x16, x11
adcs x16, x17, x12
adcs x17, x27, x13
adcs x27, x19, x14
adcs x19, x20, xzr
adc x20, x21, xzr
umulh x10, x9, x22
umulh x11, x9, x23
umulh x12, x9, x24
umulh x13, x9, x25
umulh x14, x9, x26
adds x15, x15, x10
adcs x16, x16, x11
adcs x17, x17, x12
adcs x27, x27, x13
adcs x19, x19, x14
adc x20, x20, xzr
```

```asm
; kernel: mont5_cios_epilogue items=1 charged=39.0 key=mont5_epi
sub sp, sp, #80
stp x19, x20, [sp, #0]
stp x21, x22, [sp, #16]
stp x23, x24, [sp, #32]
stp x25, x26, [sp, #48]
str x27, [sp, #64]
ldp x0, x1, [x_a, #0]
ldp x2, x3, [x_a, #16]
ldr x4, [x_a, #32]
ldp x22, x23, [x_p, #0]
ldp x24, x25, [x_p, #16]
ldr x26, [x_p, #32]
mov x15, #0
mov x16, #0
mov x17, #0
mov x27, #0
mov x19, #0
mov x20, #0
; --- 5x unrolled mont5_cios_round (5 x 48 = 240 ops) inserted here ---
subs x10, x15, x22
sbcs x11, x16, x23
sbcs x12, x17, x24
sbcs x13, x27, x25
sbcs x14, x19, x26
sbcs xzr, x20, xzr
csel x15, x10, x15, hs
csel x16, x11, x16, hs
csel x17, x12, x17, hs
csel x27, x13, x27, hs
csel x19, x14, x19, hs
stp x15, x16, [x_out, #0]
stp x17, x27, [x_out, #16]
str x19, [x_out, #32]
ldp x19, x20, [sp, #0]
ldp x21, x22, [sp, #16]
ldp x23, x24, [sp, #32]
ldp x25, x26, [sp, #48]
ldr x27, [sp, #64]
add sp, sp, #80
ret
```

Total exact AAPCS64-compliant `MONT5 = 5 * 48 + 39 = 279 <= 280.0` Ir.

## 2b. Exact 2-Limb Baseline Scalar ARMv8.0-A CIOS Montgomery Multiplication (`mont2_cios`, `47 <= 48` Opcodes <= `58.84` Ir Charged for `L = 2`)

For small composite ECM branches (`<= 128` bits, `L_min = 2` in `q15_ecm16.py`, which removes the `5`-limb floor inflation), `mulmod_instr(nbits, l_min=2)` charges `MONT5_desc * (2/5)^2 = 367.76 * 4 / 25 = 58.8416` Ir (even under `MONT5 = 280.0`, `L = 2` CIOS executes in **`47 <= 48.0 < 58.84` Ir**!). Inputs `(a0, a1, b0, b1, n0, n1, n0inv)` reside in AAPCS64 caller-saved registers `x0..x6` (zero stack spills, zero callee-saved registers, zero banned ops, disjoint GPRs `x0..x15`), returning `a * b * 2^{-128} mod N` in `x0:x1` (`0 <= t < 2N < 2^129` after both rounds for every 128-bit odd modulus `N <= 2^128 - 1`, verified bit-exactly LIVE in `cm_models.verify_kernels_symbolic` across edge moduli `2^128 - 1, 2^127 + 1, 3` and `64` random 128-bit odd moduli):

```asm
; kernel: mont2_cios items=1 charged=48.0 key=mont2_cios
mul x7, x0, x2
umulh x8, x0, x2
mul x9, x1, x2
umulh x10, x1, x2
adds x8, x8, x9
adc x9, x10, xzr
mul x11, x7, x6
mul x12, x11, x4
umulh x13, x11, x4
mul x14, x11, x5
umulh x15, x11, x5
adds xzr, x7, x12
adcs x7, x8, x14
adcs x8, x9, xzr
adc x9, xzr, xzr
adds x7, x7, x13
adcs x8, x8, x15
adc x9, x9, xzr
mul x10, x0, x3
umulh x11, x0, x3
mul x12, x1, x3
umulh x13, x1, x3
adds x7, x7, x10
adcs x8, x8, x12
adcs x9, x9, xzr
adc x10, xzr, xzr
adds x8, x8, x11
adcs x9, x9, x13
adc x10, x10, xzr
mul x11, x7, x6
mul x12, x11, x4
umulh x13, x11, x4
mul x14, x11, x5
umulh x15, x11, x5
adds xzr, x7, x12
adcs x7, x8, x14
adcs x8, x9, xzr
adc x9, x10, xzr
adds x7, x7, x13
adcs x8, x8, x15
adc x9, x9, xzr
subs x10, x7, x4
sbcs x11, x8, x5
sbcs xzr, x9, xzr
csel x0, x10, x7, hs
csel x1, x11, x8, hs
ret
```

## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Exact 2nd-Moment ECM Chebyshev

- **Explicit Schoolbook Multi-Limb Scalar ARMv8-A Tower Norm Derivation (`Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`, ZERO Karatsuba shortcut):**
  Explicit schoolbook scalar ARMv8-A `L_1 x L_2` limb multiplication uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)` (`mul/umulh/adds/adcs` carry chains). On the production sieving region `Q_CM(a) <= s_c * 266,845 <= 533,690` (all 2-adic classes `s_c in {1, 1.5, 2}`), `Q_CM(b) <= 387,349.4`, exact Cauchy-Schwarz via `max diag(G_8^{-1}) = 2/5` proves `|a_j| <= 462`, `|b_j| <= 394`, `|a^3_j| <= 2^{26.90} < 2^{27}` and `|b^3_j| <= 2^{26.20} < 2^{27}` (1 signed limb; `4` `8x8` 1-limb poly muls + `Phi_15` reduction = `400` Ir), `|g(a, b)_j| = |a^3 - u w^2 b^3|_j <= 2^{89.9} < 2^{90}` (2 signed limbs; `+96` Ir), and `|f(a, b)_j| <= 2^{82.2} < 2^{83}` (fits in 2 signed limbs, conservatively priced as 3 limbs: `a^6, b^6` = `2 x 928 = 1,856` Ir, `a^9, b^9, a^6 b^3, a^3 b^6` = `4 x 1,152 = 4,608` Ir, combine = `128` Ir). Using pure schoolbook `mul_LL` across the 3-step tower norm `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:
  1. **`_nk_g_only` (`2 -> 4 -> 8 -> 16` limbs):** `32` `2x2` (`1,536`) + `6` `4x4` (`960`) + `3` `8x8` (`1,728`) + adds/subs/setup/cold-scan (`1,544`) = **`5,768 <= 6,500` Ir charged**!
  2. **`_nk_f_only` (`3 -> 6 -> 12 -> 24` limbs):** reuses `a^3, b^3` -> `a^6, b^6, a^9, b^9, a^6 b^3, a^3 b^6` + combine (`6,592`) + 3-limb schoolbook tower norm (`4,608 + 4,032 + 3,744 = 12,384`) = **`18,976 <= 20,000` Ir charged**!
- **Parsed 62-Bit Shoup Radix-2 NTT Butterfly (`ntt62_butterfly`: `18` opcodes/2-elems/prime -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 fwd + 1 inv across 4 primes) & 4-Prime Garner CRT (`fft_F`):**
  `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`log2(NTT4_PROD) = 244.70`, all 4 verified prime LIVE). One exact baseline scalar ARMv8-A radix-2 Shoup butterfly (`w_s = floor(w * 2^64 / p)` precomputed, `t = w*v - umulh(v, w_s)*p` reduced to `[0, p)`) takes **`18` opcodes per 2 elements (`9.0` Ir/elem/prime -> `36.0` Ir/elem/transform across 4 primes -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 forward + 1 inverse NTTs in `fft_F`, `1.56x` headroom!)**; scalar `F_{p_i}` Barrett MAC in `PM-Basis` uses `mul/umulh` + 1-step Barrett `umulh/msub` + 2 conditional subs + add/sub/csel + 2 loads + loop = **`17 <= 24` opcodes** (verified bit-exactly LIVE on all 4 NTT primes in `kernels_symbolic`); 4-prime Garner CRT uses 3 Shoup mod-muls (`15` ops) + 3 multi-limb mul-adds (`<= 75` ops) = `<= 90 << 220` linear term in `fft_F` (`168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`; 2D `168 M_tot log2(M_tot) + 1100 M_tot + 1200`). Both `ntt62_butterfly` (64 random triples) and 4-prime Garner CRT (32 random 240-bit integers) are verified bit-exactly LIVE in `audit_cm.py`, alongside `srt_toy` (128 planted leaves, `0` mismatches):

```asm
; kernel: ntt62_butterfly items=2 charged=10.0 key=ntt62_elem_prime
.Lntt:
ldr x0, [x_u]
ldr x1, [x_v]
ldp x2, x3, [x_w], #16
mul x4, x1, x2
umulh x5, x1, x3
msub x4, x5, x_p, x4
subs x6, x4, x_p
csel x4, x6, x4, hs
add x7, x0, x4
subs x8, x7, x_p
csel x7, x8, x7, hs
subs x9, x0, x4
add x10, x9, x_p
csel x9, x9, x10, hs
str x7, [x_u], #8
str x9, [x_v], #8
subs x_cnt, x_cnt, #1
b.ne .Lntt
```

- **Parsed 3-Prime `< 2^62` Radix-4 Lazy Shoup NTT (`ntt62_r4_shoup`: `49` opcodes / `8` elem-stages = `6.125 <= 6.2` Ir/elem/stage/prime -> `3 primes x 3 transforms = 55.125 <= 55.25` Ir/elem/stage in `ntt3_F`) & 5-Transform Middle-Product `RT` (`gbatch_ntt3`):**
  Three `< 2^62` primes divisible by `3 * 2^55`: `p_1 = 69 * 2^55 + 1`, `p_2 = 57 * 2^55 + 1`, `p_3 = 54 * 2^55 + 1` (`prod = 2^182.696 > 2^54 * (2^64 - 1)^2`, all 3 verified prime + 3-prime Garner CRT + bit-exact `n = 16` & `n = 24` radix-4/radix-3 Shoup convolutions using exact `subs`/`lo` ARMv8-A semantics LIVE in `cm_models.verify_kernels_symbolic`). `ntt3_F(M) = min(55.25*n2*k + (26*n2 if k odd) + 220*n2, 55.25*n3*(k-2) + (26*n3 if k-2 odd) + 364*n3)` (`+ 512*n` if `n > 2^26`; the explicit `+26*n` term prices the final radix-2 pass whenever the power-of-2 exponent is odd, and the linear `+220*n` absorbs inter-group pointer resets). Inputs live in `[0, 2p)` (`x_p2 = 2p < 2^63`, guaranteeing Shoup outputs `< 2p` and exact unsigned borrow via `subs` + `csel ..., lo`): `4` loads + `3` twiddle `ldp`s + `4` Lazy Shoup muls (`mul/umulh/msub`) + `8` lazy mod-`2p` add/subs/csel (`3` ops each) + `4` stores + `2` loop (`subs`/`b.ne`) = **`49` opcodes / `8` elem-stages = `6.125 <= 6.2` Ir/elem/stage/prime** (`9 * 6.125 = 55.125 <= 55.25`):

```asm
; kernel: ntt62_r4_shoup items=8 charged=6.2 key=ntt62_r4_elem_stage
.Lntt_r4:
ldr x0, [x_a0]
ldr x1, [x_a1]
ldr x2, [x_a2]
ldr x3, [x_a3]
ldp x4, x5, [x_tw], #16
mul x6, x1, x4
umulh x7, x1, x5
msub x1, x7, x_p, x6
ldp x4, x5, [x_tw], #16
mul x6, x2, x4
umulh x7, x2, x5
msub x2, x7, x_p, x6
ldp x4, x5, [x_tw], #16
mul x6, x3, x4
umulh x7, x3, x5
msub x3, x7, x_p, x6
add x8, x0, x2
subs x9, x8, x_p2
csel x8, x9, x8, hs
subs x10, x0, x2
add x9, x10, x_p2
csel x10, x9, x10, lo
add x11, x1, x3
subs x9, x11, x_p2
csel x11, x9, x11, hs
subs x12, x1, x3
add x9, x12, x_p2
csel x12, x9, x12, lo
mul x6, x12, x_j
umulh x7, x12, x_jp
msub x12, x7, x_p, x6
add x0, x8, x11
subs x9, x0, x_p2
csel x0, x9, x0, hs
subs x2, x8, x11
add x9, x2, x_p2
csel x2, x9, x2, lo
add x1, x10, x12
subs x9, x1, x_p2
csel x1, x9, x1, hs
subs x3, x10, x12
add x9, x3, x_p2
csel x3, x9, x3, lo
str x0, [x_a0], #8
str x1, [x_a1], #8
str x2, [x_a2], #8
str x3, [x_a3], #8
subs x_cnt, x_cnt, #1
b.ne .Lntt_r4
```
- **Stage B (`P_min = 32`, `KH = 0.08`, Static 3-Class `T_g[c] = (299, 287, 275)` + Constant `T_f_const = 118`b Pre-Gates, `step = 0.25`):** `gbatch_ntt3(pbits_split(32, 0.08), leaf=L, c_leaf=10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16))` with AAPCS64 `MONT5 = 280.0`. `S_g` is gated solely by static 2-adic class constant `T_g[c] = (299, 287, 275)` (ZERO `Q_CM(a)` radial shells, ZERO per-`bgo` predictor!); inside `.Lhit64`, after `S_g` fires, the cold-path handler checks the **single constant `f`-side threshold `S_f >= T_f_const = 118` bits (`472` quarter-bits)** (ZERO pre-`N_f` knowledge of `bf`!) before computing `N_g` (`6,500` Ir) + `Stage-B(g)` (`sbg = 315`) -> `N_f` (`20,000` Ir) + `Stage-B(f)` (`sbf = 375`) -> `fbatch`/`gbatch`, giving exact **`c_cof_2d_upper = 1.0557` Ir/cell** (evaluating `cst_g`/`cst_f.survivor_upper` at `max(bin_lo, bin_hi)` across each 1-bit bin; analytic Chebotarev split density for `g_3` over `Q(zeta_15)` is `5/72 = 0.06944 < 0.08`, 50k-prime census `0.06995 / 0.05951`, live `2^14` check `0.07197 / 0.05649 <= 0.08`).
- **Cauchy-Schwarz Correlated Random-Sum Collection ECM Chebyshev (`cm_models.ecm_mean_var`) & Branch-Tree-Capped `L_min = 2` Geometric-Ladder Descent ECM (`q15_ecm16.py`):** Before Suyama ECM, an exact prime-exponent Newton integer perfect-power precheck (`_is_perfect_power(N)`, explicitly charged at `NEWTON_PRECHECK_INSTR = 2^26` Ir > 32.06M schoolbook ops for 2,590-bit boot inputs) eliminates all `N = p^k` singular redraws. Collection and descent retry regimes are cleanly decoupled: (i) for **post-`2^48` collection cofactors** (`cm_models.ecm_mean_var`), conservative exploratory full-schedule retry `p_r = 1e-3` with `E[T] = E_raw / (1 - p_r)` and the **exact Cauchy-Schwarz upper bound on correlated random-sum variance** combined with deterministic truncation at `n_rels_upper = 2^63.530` relations (`Pr[N_rels > n_rels_upper] == 0`) yield `ew_l2 = 94.5038 << 95.3025` and Chebyshev **`2^-70.62 < 2^-70`**; (ii) for **descent cofactors** (`q15_ecm16.py`), `tree_branch_caps(nbits, b_min=32, l_min=2)` computes the **exact dynamic-programming worst-case active-branch limb-square cap `sum L_i^2`** and **`sum gcd_instr(b_i)`** over ALL partitions `b_1 + ... + b_t <= nbits` into composite branches `b_i >= b_min = 32` (every prime `<= 2^16` stripped by `TDIV`), using the `47`-op `mont2_cios` kernel (`L_min = 2`, charged `58.84` Ir) to close Codex `V8-02` completely; every rung runs a **`12`-bit geometric ECM ladder (`ladder_v9_instr` from `16`b up to `fbits`)** with conservative separation penalty **`SEPARATION = 0.5`** (`curves = ceil(log(1-s)/log(1-0.5*p1))`) and **explicit failed-replay retry pricing** `cost = curves * (per_curve + nf * p1 * replay_1curve + (p_r / (1 - p_r)) * (per_curve + replay_1curve)) + 2^26` at `p_r = 0.0043` (every rejected simultaneous-hit draw pays BOTH `per_curve` AND full `replay_1curve`, closing Codex `V8-03`).

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Dense + 16-Shard `< 2^64 B` 4-Vec 4-Col Arch B Sign-Split Scalar ARMv8.0-A + `rev_b1_stream_u16` Matrix Streaming + Exact 1-Step Barrett `row_red5_z12`)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`, `tab_ratio_charged = 1.00` plain sparse with ZERO `TABM`, `32` sequential 4-vector bundles of `128 B/slot`, `n = 128 = 32 * 4`, `m = 32,768`, peak single-owner dense-vector RAM `2 * N_sw * 128 B = 2^63.368 B < 2^64 B`, `16` read-only sparse-matrix shard owners `2^62.107 B << 2^64 B` each (byte-balanced across 16 shards; each row stores a `64 B` = `8`-word header holding `(d_i, ceil(nu_i/8))` plus all **`24` `uint16_t` sub-list lengths** for the 12 signed unit + 12 signed non-unit `+-chi^b` buckets, `d_i` 8B pre-scaled column offsets, and `ceil(nu_i/8)` packed `uint8_t` magnitudes -> `<= 1,709.15 B/row`), explicit priced cross-owner burst streaming via `rev_b1_stream_wrap16` + `rev_b1_stream_u16` + `rev_b1_stream_tail1` (§1.3b): `w_i = 8 + d_i + ceil(nu_i/8) <= 8 + d_i + nu_i/8 + 7/8`, aggregate mean `E[w_i] <= 213.64375` words/row, `E[ops(w_i)] <= 1.125 * 213.64375 + 47.125 = 287.474 <= 291.0 = c_stream_mat_1vec` (1-vec Horner) and `<= 73.0 * 4` (`c_stream_mat = 73.0` Ir/vec-row in 4-vec Krylov)!):
- **`unit_run_4vec4_z12` & `unit_run_4vec4_neg_z12` (`122` opcodes / `16` `(entry, vec)` pairs = `7.625 <= 7.65` Ir/vec-entry, `28 <= 29` GPRs `x0..x17, x19..x28`, ZERO `x18/x29/x30`, `0` spills):** Each column holds 4 contiguous 256-bit vectors (`128 B/col`), addressed via buffer-independent pre-scaled `uint64_t` byte offset `col * 128 < 2^63` (`N_sw * 128 = 2^62.368 B < 2^64 B`) added explicitly to `x6 = x_base` (`ldr x7, [x4], #8; add x7, x6, x7` — works natively on both ping-pong buffers with ZERO vector copying!). Each loop iteration processes 4 columns across all 4 vectors (`[x7, #0/#16]`, `[x7, #32/#48]`, `[x7, #64/#80]`, `[x7, #96/#112]`) into 4 signed 5-limb accumulators (`v0`: `x8..x12`; `v1`: `x13..x17`; `v2`: `x19..x23`; `v3`: `x24..x28`), `4 * (2 + 4 * 7) + 2 = 122 / 16 = 7.625 <= 7.65` Ir/vec-entry. The `0..3` leftover tail entries per sub-list are handled by a 1-entry 4-vec scalar tail (worst-case `12` signed unit sub-lists x `<= 3` entries x `8` ops/vec-entry = `<= 288` Ir/vec, explicitly priced inside `row_tail_6bucket_z12 = 504 <= c_row_chi = 960` Ir):

```asm
; kernel: unit_run_4vec4_wrap items=1 charged=4.0 key=unit_run_4vec4_wrap
subs x5, x5, #4
b.lo .Lrun4_tail
; --- inline unit_run_4vec4_z12 / _neg (.Lrun4v4 / .Lrun4v4n) runs here while C == 1 ---
.Lrun4_tail:
adds x5, x5, #4
b.eq .Lrun4_done
; --- 1..3 leftover entries run 1-entry 4-vec tail (<= 288 Ir/vec across 12 sub-lists, priced inside row_tail_6bucket_z12 = 504 <= c_row_chi = 960) ---
.Lrun4_done:
```

```asm
; kernel: unit_run_4vec4_z12 items=16 charged=7.65 key=unit_run_4vec4
.Lrun4v4:
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
adds x13, x13, x0
adcs x14, x14, x1
adcs x15, x15, x2
adcs x16, x16, x3
adc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
adds x19, x19, x0
adcs x20, x20, x1
adcs x21, x21, x2
adcs x22, x22, x3
adc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
adds x24, x24, x0
adcs x25, x25, x1
adcs x26, x26, x2
adcs x27, x27, x3
adc x28, x28, xzr
subs x5, x5, #4
b.hs .Lrun4v4
```

```asm
; kernel: unit_run_4vec4_neg_z12 items=16 charged=7.65 key=unit_run_4vec4_neg
.Lrun4v4n:
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
ldr x7, [x4], #8
add x7, x6, x7
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldp x0, x1, [x7, #32]
ldp x2, x3, [x7, #48]
subs x13, x13, x0
sbcs x14, x14, x1
sbcs x15, x15, x2
sbcs x16, x16, x3
sbc x17, x17, xzr
ldp x0, x1, [x7, #64]
ldp x2, x3, [x7, #80]
subs x19, x19, x0
sbcs x20, x20, x1
sbcs x21, x21, x2
sbcs x22, x22, x3
sbc x23, x23, xzr
ldp x0, x1, [x7, #96]
ldp x2, x3, [x7, #112]
subs x24, x24, x0
sbcs x25, x25, x1
sbcs x26, x26, x2
sbcs x27, x27, x3
sbc x28, x28, xzr
subs x5, x5, #4
b.hs .Lrun4v4n
```

- **`row_red5_z12` (`67 <= 68.0` Ir/bucket-vec — Exact 1-Step Signed 5-Limb Barrett Reduction for ALL `|a_4| < 2^60`):**
  Handles ANY signed 5-limb bucket accumulator `A = (x10..x14)` with `|a_4| < 2^60` (`|A| < 2^316`, covering EVERY `|c| <= 255` row up to `+-2^59 * (R_BLS - 1) > 10^17 * (R_BLS - 1)` with ZERO carry overflow!). Computes branchless `U = |A|`, `top = (u_4 << 2) | (u_3 >> 62)`, `q_hat = umulh(top, mu_r)` where `mu_r = floor(2^318 / R_BLS) = 0x8d54253b7fb78ddf < 2^64` (`0 <= floor(U / R_BLS) - q_hat <= 2` unconditionally), `rem = U - q_hat * R_BLS in [0, 3 * R_BLS)`, two branchless conditional `- R_BLS` steps, and branchless sign fix:

```asm
; kernel: row_red5_z12 items=1 charged=68.0 key=row_red5_barrett
asr x15, x14, #63
eor x0, x10, x15
eor x1, x11, x15
eor x2, x12, x15
eor x3, x13, x15
eor x4, x14, x15
subs x0, x0, x15
sbcs x1, x1, x15
sbcs x2, x2, x15
sbcs x3, x3, x15
sbc x4, x4, x15
lsl x5, x4, #2
orr x5, x5, x3, lsr #62
umulh x5, x5, x_mu
mul x6, x5, x22
mul x7, x5, x23
mul x8, x5, x24
mul x9, x5, x25
umulh x16, x5, x22
umulh x17, x5, x23
umulh x20, x5, x24
umulh x19, x5, x25
adds x7, x7, x16
adcs x8, x8, x17
adcs x9, x9, x20
adc x19, x19, xzr
subs x0, x0, x6
sbcs x1, x1, x7
sbcs x2, x2, x8
sbcs x3, x3, x9
sbc x4, x4, x19
subs x6, x0, x22
sbcs x7, x1, x23
sbcs x8, x2, x24
sbcs x9, x3, x25
sbcs x16, x4, xzr
csel x0, x6, x0, hs
csel x1, x7, x1, hs
csel x2, x8, x2, hs
csel x3, x9, x3, hs
csel x4, x16, x4, hs
subs x6, x0, x22
sbcs x7, x1, x23
sbcs x8, x2, x24
sbcs x9, x3, x25
sbcs xzr, x4, xzr
csel x0, x6, x0, hs
csel x1, x7, x1, hs
csel x2, x8, x2, hs
csel x3, x9, x3, hs
orr x6, x0, x1
orr x7, x2, x3
orr x6, x6, x7
cmp x6, #0
csel x15, xzr, x15, eq
subs x6, x22, x0
sbcs x7, x23, x1
sbcs x8, x24, x2
sbc x9, x25, x3
tst x15, x15
csel x0, x6, x0, ne
csel x1, x7, x1, ne
csel x2, x8, x2, ne
csel x3, x9, x3, ne
stp x0, x1, [x_out, #0]
stp x2, x3, [x_out, #16]
ret
```

- **Sign-Split 2-Vec `gen_run1_2v_z12` (`c_b > 0`) & `gen_run1_2v_neg_z12` (`c_b < 0`) (`43/2 = 21.50 <= 22.0` Ir/vec-entry, `24 <= 29` GPRs avoiding `x18/x29/x30`, `0` spills; caller-side non-empty `cbz x24` precondition guarantees `x24 >= 1` on entry):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`, `< 6.5%` of entries, stored as 1-byte `uint8_t` magnitudes loaded and zero-extended to `x22` via `ldrb w22, [x23], #1`) into positive and negative sub-lists, executed as two spill-free 2-vector passes (`(v0, v1)` at offsets `#0..#48`, then `(v2, v3)` at offsets `#64..#112`, explicit `add x8, x9, x8` from `x_base = x9`): `3` ptr/coeff loads + `2 * 19` mul/umulh/carry ops + `2` loop = **`43 / 2 = 21.50 <= 22.0` Ir/vec-entry** (`k5_live = eu_lo * 7.65 + (1 - eu_lo) * 22.0 = 8.5581 <= k5_charged = 10.10` Ir, `+1.542` Ir headroom!):

```asm
; kernel: gen_run1_2v_z12 items=2 charged=22.0 key=gen_run1_2v_coord
.Lgen2v:
ldr x8, [x21], #8
add x8, x9, x8
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adc x14, x14, x7
ldp x0, x1, [x8, #32]
ldp x2, x3, [x8, #48]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x15, x15, x4
adcs x16, x16, x5
adcs x17, x17, x6
adcs x19, x19, x7
adc x20, x20, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x16, x16, x4
adcs x17, x17, x5
adcs x19, x19, x6
adc x20, x20, x7
subs x24, x24, #1
b.gt .Lgen2v
```

```asm
; kernel: gen_run1_2v_neg_z12 items=2 charged=22.0 key=gen_run1_2v_neg_coord
.Lgen2vn:
ldr x8, [x21], #8
add x8, x9, x8
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x10, x10, x4
sbcs x11, x11, x5
sbcs x12, x12, x6
sbcs x13, x13, x7
sbc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x11, x11, x4
sbcs x12, x12, x5
sbcs x13, x13, x6
sbc x14, x14, x7
ldp x0, x1, [x8, #32]
ldp x2, x3, [x8, #48]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x15, x15, x4
sbcs x16, x16, x5
sbcs x17, x17, x6
sbcs x19, x19, x7
sbc x20, x20, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x16, x16, x4
sbcs x17, x17, x5
sbcs x19, x19, x6
sbc x20, x20, x7
subs x24, x24, #1
b.gt .Lgen2vn
```

- **Explicit 1-Vector Horner Kernels (`unit_run_1vec4_z12` / `_neg`: `38/4 = 9.50 <= 9.60` Ir/entry; `gen_run1_1v_z12` / `_neg`: `24.0 <= 24.0` Ir/entry with caller-side non-empty `cbz x24` precondition `x24 >= 1` -> `k5_horner_charged = 11.80`):** Reuses the EXACT SAME read-only sparse-matrix shards (`128*j` stored byte offsets) on the `32 B/slot` 1-vector Horner dense buffer via **`add x7, x6, x7, lsr #2`** (`(128*j) >> 2 = 32*j` exact byte offset, ZERO matrix re-encoding, ZERO offset-shift bug!), simulated bit-exact in `cm_models.py`:

```asm
; kernel: unit_run_1vec4_z12 items=4 charged=9.60 key=unit_run_1vec4
.Lrun1v4:
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
adds x8, x8, x0
adcs x9, x9, x1
adcs x10, x10, x2
adcs x11, x11, x3
adc x12, x12, xzr
subs x5, x5, #4
b.hs .Lrun1v4
```

```asm
; kernel: unit_run_1vec4_neg_z12 items=4 charged=9.60 key=unit_run_1vec4_neg
.Lrun1v4n:
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
ldr x7, [x4], #8
add x7, x6, x7, lsr #2
ldp x0, x1, [x7, #0]
ldp x2, x3, [x7, #16]
subs x8, x8, x0
sbcs x9, x9, x1
sbcs x10, x10, x2
sbcs x11, x11, x3
sbc x12, x12, xzr
subs x5, x5, #4
b.hs .Lrun1v4n
```

```asm
; kernel: gen_run1_1v_z12 items=1 charged=24.0 key=gen_run1_1v_coord
.Lgen1v:
ldr x8, [x21], #8
add x8, x9, x8, lsr #2
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adc x14, x14, x7
subs x24, x24, #1
b.gt .Lgen1v
```

```asm
; kernel: gen_run1_1v_neg_z12 items=1 charged=24.0 key=gen_run1_1v_neg_coord
.Lgen1vn:
ldr x8, [x21], #8
add x8, x9, x8, lsr #2
ldrb w22, [x23], #1
ldp x0, x1, [x8, #0]
ldp x2, x3, [x8, #16]
mul x4, x0, x22
mul x5, x1, x22
mul x6, x2, x22
mul x7, x3, x22
subs x10, x10, x4
sbcs x11, x11, x5
sbcs x12, x12, x6
sbcs x13, x13, x7
sbc x14, x14, xzr
umulh x4, x0, x22
umulh x5, x1, x22
umulh x6, x2, x22
umulh x7, x3, x22
subs x11, x11, x4
sbcs x12, x12, x5
sbcs x13, x13, x6
sbc x14, x14, x7
subs x24, x24, #1
b.gt .Lgen1vn
```

- **`chi_mul_z12_fast` (`50` Ir), `row_tail_6bucket_z12` (`504` Ir = `12` signed unit sub-lists x max `4` tail entries x `8` ops + `2 * 50` `chi_mul` + `20` dispatch), 1-Dense + 16-Projection + 16-Shard 4-Vec Arch B Row Overhead (`c_row_chi = 960.0` Ir, `c_stream_mat = 73.0` Ir/vec-row Krylov, `c_stream_mat_1vec = 291.0` Ir/row 1-vec Horner, `k5_horner_charged = 11.80`, `c_xproj = 70.0` Ir/vec-row incl. `4.5` Ir/vec-row cross-owner vector stream to the `16` projection owners + explicit `4.5*n = 576` Ir/row Horner `Y = [S, BZ]` read stream, `n_owners = 1`, `n_proj_owners = 16`, `n_mat_owners = 16`, `te = 0`):** `6 * 67 + 504 + 16 = 922 <= 6 * 68 + 504 + 16 = 928 <= c_row_chi = 960.0` Ir worst-case; aggregate matrix-shard `u16`-burst streaming (`E[w_i] <= 8 + 201 + 0.15*201/8 + 7/8 = 213.644` words/row incl. 64B 24-bucket header with `W_max = 4,096` row cap -> `1.125*213.644 + 47.125 = 287.47 <= 291.0 = c_stream_mat_1vec`, `<= 73.0 * 4` Krylov) charged inside `c_prod` / `c_prod_horner` (`tab_ratio_charged = 1.00`, ZERO `TABM`).

## 5. Descent Kernels (Generous Scalar Bounds for Uncapped `dsieve16` & Exact-Integer `BOUND_CELL`, Plus Input-Uniform Bareiss LLL)

Descent is `5.85%` of opened-run work (`2^106.2192` u). Rather than hand-scheduling fragile cold-path assembly for 16D Fincke-Pohst DFS or interval ball arithmetic, we charge generous, manifestly-conservative scalar ARMv8.0-A operation bounds and verify every mathematical property end-to-end in exact integer Python:
- **Uncapped `dsieve16` on `S_16(R) = B_16(R)/{+-1}` (`SIEVE_HIT = 2,600` Ir/hit, ZERO per-prime skip caps):** Every relative-degree-1 prime ideal with `p in [2^8, 2^70]` is reduced via Bareiss LLL (`PAIR_SETUP = 2^35.2`) and enumerated completely over `B_16(R)/{+-1}` using plain compiled scalar C (`-march=armv8-a -mgeneral-regs-only`) 16D Fincke-Pohst Q32 DFS + leaf vector `v = V_stk[0] + z_0 * b_0` + `v = 0` exclusion + ambient sign canonicalization (`v -> -v` if first nonzero ambient coordinate `< 0`) + 128-bit bijective half-ball rank `[0, C_R - 1]` via `129.03 MiB` flat prefix table + owner routing + `uint16_t` score add. Explicit 16D Fincke-Pohst measurements on 99/100-LLL-reduced prime lattices confirm `edge/hit = 2.24` (`1.24` internal nodes/hit) at `R = 12, p ~ 2^28` (`8.0e7` hits) and `edge/hit = 2.13` (`1.13` internal nodes/hit) at `R = 16, p ~ 2^38` (`7.9e6` hits), rising to `4.64` only in the sparse small-hit edge regime (`R = 20, p ~ 2^52`, `1.7e4` hits; all `76` nodes use `R in [12..64]`, where `54/67` sieve nodes sit in the dense `R in [40..64]` regime and the `13/67` smaller-radius sieve nodes at `R in {14, 16, 20, 28, 32}` have empty/sparse DFS edge traversals `~2,000` Ir/prime absorbed by `> 10^6x` inside `PAIR_SETUP = 2^35.2` Ir/prime, which carries `3.65 billion` Ir/prime slack over the `2^35.06` Bareiss LLL bound!):
  - **Per-hit breakdown (`<= 914–1,244 << 2,600` Ir/hit):** DFS traversal (`<= 1.25` internal nodes/hit x `120` Ir + leaf check `80` Ir = `<= 230` Ir; `<= 560` Ir even at sparse edge `4.64`) + leaf vector materialization (`64` Ir) + ambient sign canonicalization (`40` Ir) + 16D half-ball bijective ranking (`16` table lookups + 128-bit adds = `320` Ir) + owner routing / `uint16_t` update / compiler spill headroom (`260` Ir) = **`<= 914–1,244 << 2,600` Ir/hit** (`> 2.0x` headroom). Exact 4D Fincke-Pohst half-ball DFS (100% match vs brute force on 18,174 points), ambient sign canonicalization (explicitly verifying witness `z = (-1, 1, 0, 0) -> v = (-3, 3, 1, 0) -> (3, -3, -1, 0)`), 55-bit `SQ_LUT`, and bijective 16D half-ball ranking (`552/552` canonical + sign-flipped `R = 2` hits + `v = 0` exclusion) are verified LIVE in `cm_models.py`.
- **Certified Outward-Rounded Integer `BOUND_CELL (+6)` (`BOUND_CELL = 5,000` Ir/cell, `q15_bound_cell.py`):** Evaluates strict signed-64b outward-rounded complex ball arithmetic (`bmul_6` with `+6 >= ceil(sqrt(5) + 3)` rounding margin, `13` `bmul_6` calls/place across `4` complex embeddings, `<= +0.5343b < +0.625b` `SQ_LUT` envelope). **Per-cell breakdown:** `3,792` Ir exact integer ball-mul core (`16 + 4*(104 + 13*48 + 128 + 88)`) + `<= 1,208` Ir headroom for vector materialization, `SQ_LUT` integer square root, `log2`, `(nf_max, ng_max)` cap compare, and compiler spills = **`5,000` Ir/cell**. Deterministically rejects any cell with `ub_f_rem > nf_max` or `ub_g_rem > ng_max` (guaranteeing ZERO over-norm cells ever enter ECM). Exact `bmul_6` (`+6`) outward rounding (including 50-digit `Decimal` worst-case witness) and 12-cell real-lattice exact integer `BOUND_CELL` certification (`bc_gmin >= 0.0`) are verified LIVE in `cm_models.py` and `verify_q15_descent.py`.
- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every special-`q` node lattice `L_Q` (`det = q`) and `G_T` boot lattice runs exact 2-phase `(11/40 Euclidean -> 99/100 Q_CM = diag(G_8, G_8))` Bareiss LLL (`lll_16d_int(..., use_cm=True)` in `cm_models.py`, verified live on both `Bsq` and `I_16`; `dsieve16` `L_p` uses Euclidean `lll_16d_int(..., use_cm=False)`), where `G_8` has `det(G_8) = 3^4 * 5^6` and `lambda_max(G_8) = 15` (proved live via Sylvester minors of `15*I_8 - G_8`, so `d_k^(CM) <= 15^k * d_k^(Euc)` unconditionally!), yielding explicit derived uniform bounds **`PAIR_SETUP = 2^35.2`** (70-bit `L_p` Euclidean `11/40 -> 99/100`, `2^35.06 < 2^35.2`), **`LLL_NODE = 2^40.0`** (special-`q` `q <= 2^262` `11/40 Euc -> 99/100 Q_CM`, `2^39.41 < 2^40.0`), and **`LLL_BOOT = 2^47.0`** (381-bit `G_T` boot WITH `S_boot = 2^21` scaling on the 8 `b`-coords BEFORE Phase 1: `det = p^4 * 2^168 <= 2^1692`, `L_d <= 58 <= 59` limbs -> **`2^46.75 < 2^47.0`**), all derived and asserted LIVE in `verify_kernels_symbolic` alongside `BOUND_CELL = 5,000`, `SIEVE_HIT = 2,600`, `SCAN = 10`, `ZERO = 2`, `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`.
