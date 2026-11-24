# Kernels (Raw-Opcode-Verified Baseline Scalar ARMv8.0-A: `Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8`, 3-Prime Unrolled Radix-4 Lazy Shoup NTT + Lehmer32 GCD, 1-Dense + 16-Proj + 16-Shard 4-Vec LA, Explicit Build & Single-Sub Barrett Kernels)

Every hot-loop kernel cost is derived from explicit baseline scalar ARMv8.0-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero `bl` — all kernels are inline macros or leaf routines returning via `ret`). Cold paths (`.Lhit64` hit-block handler, `.Lb2_chunk_slow`, and 16D descent DFS / `BOUND_CELL`) use generous, manifestly-conservative scalar ARMv8.0-A operation bounds. `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file (`50` hot-loop blocks, `0` undefined `.L...` labels, per-kernel disjoint physical register maps asserting ZERO symbolic/numbered collisions in every single block), validates opcodes & `LDP`/`STP` addressing, counts every instruction (`cmp`+`b.*`=2 ops, loop control included), and simulates: (i) explicit 2-pass counting-sort build (`rev_b2_build_hist_u4` `8.0` + `rev_b2_build_scatter_u4` `7.0` + chord resets `<= 16.0 = C_LIST_A`, `rev_build_prefix_u4` `3.0 <= 4.0 = C_LBKT`), (ii) `32x -> 4x -> bin21` Band-1a (`4.625 <= 4.63` / `5.50 <= 5.52` Ir/hit, `b1a_long_lk = 203.0`, `b1a_short_lk = 70.0`) & Galois-orbit-specialized Band-1b `u4 -> bin21` (`6.75` Ir/hit `u4` / `11.25` Ir/hit `u1`, `b1bl_lk = 20.0`, `b1bs_lk = 8.0`), (iii) deferred-flush `128 KiB`-buffered `8 B` Band-2 (`W_SB = 16`, `rev_b2_row_k0_step4` `8.25 + 0.339 <= 8.60`, `25 + 5 = 30.0` tile probe/guard, `9 + <=1 = 10.0` `tbnz`-guarded SB hdr, flush-free **`8`-op** `bfi` `rev_b2_emit_c0/c1` + `3.5625` `rev_b1a_apply_u32` recv + `<= 0.01` tail/chunk = **`11.58 = C_HIT_HI`**), (iv) `2`-op wrap + `43`-op post-indexed 64-cell `ORR` scan (`<= 0.76`), (v) `86/16 = 5.375 <= 5.38` 2-butterfly `ntt62_r4_u2` (`48.50` Ir/elem/stage), schoolbook small-`M` crossover, `100*n` 3-prime Garner CRT, and `29/2 = 14.5 <= 16.0` `lehmer32_step` 64-bit Lehmer GCD (`c_cof_2d_upper = 0.7445`), (vi) 4-vec/1-vec `Z[zeta_12]` LA with single-sub **`57`-op** Barrett `row_red5_z12` (`|A| <= 2^{290}`, `c_row_chi = 600.0`, `c_stream_mat = 63.0 / 250.0`, `k5 = 9.45 / 11.20`, `c_xproj = 65.0`), AAPCS64 `MONT5 = 280` (`279` ops), `47`-op `mont2_cios`, and `BOUND_CELL = 4,200`, asserting `opcodes/items <= charged`.

## 1. Collection Kernels (`Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 448 = 7 * 64`, `W_SB = 16`, `P_min = 2^32`, `P_m4s = 2^50.5`, `P_m4 = 2^52.5`, `P_1bs = 2^56.5`, `P_lo = P_f = 2^58.0`, `P_g = 2^63.5`)

Because the sieving region is a pure product of two 8D balls (`qa(a) <= s_c R2_a` with 2-adic radii `s_c = (1.0, 1.5, 2.0)` and `qb(b) <= R2_b`, `A^2 = 320,500`, `r_ab = 0.83`), every canonical `a`-table entry (`V_a = 2^55.8564 < 0.91 * 2^56`, `N_a^tbl = 4 V_a = 2^57.8564`) is in-region for every `b` (`V_b = 2^64.5355`, `L_max = 439.40 < W_tile = 448`, `L_b = 188.72`, `N_7 = 2^54.9752`, `V_kept = 2^120.3919`). There is **ZERO `q_key`, ZERO `q_tab`, ZERO `sb_min_q8`, ZERO binary search, ZERO global sorting invariant, and ZERO `f`-side Band 2** (`P_f = P_lo = 2^58.0`, `28` slab owners unchanged at `8.50 EiB < 16.00 EiB` each).

### 1.1a `rev_b2_row_k0_step4` (Band-2 4-prime orbit incremental `D_(omega, sigma)` row-step `k_0` init: `33/4 = 8.25` Ir/prime/row + `<= 0.339` outer-chord reset <= `C_B2_ROW_K0 = 8.60` Ir/prime/row)

Because `sigma(omega) = omega` (`7*5 === 4*5 === 13*5 === 5 mod 15`), all `4` primes `s_k = s_0^(7^k)` (`k = 0..3`) in a Band-2 prime orbit share the EXACT same modulus `p` (`2^58 <= p < 2^63.5`), the EXACT same `b_0`-coefficient `pwk[0] = 1`, and the EXACT same `b_1`-coefficient `pwk[1] = w = s_0^5 mod p`. Moreover, the `D_(omega, sigma)` fundamental-domain predicate `in_D12(b)` depends ONLY on `(b_2..b_7)`, so for every fixed `(b_2..b_7)`, `b_1` sweeps an uninterrupted contiguous integer chord of mean length `L_b1 = (96/105) * sqrt(R2_b / 6) = 232.30` rows (`> L_b = 188.72`). Maintaining `k_tail[0..3]` (`32 B` in `[x_st]`) across consecutive `b_1 -> b_1 + 1` rows and converting signed `x_b0` (`|x_b0| <= L_max = 439.40 << 2^58 <= p`, precondition `-p < x_b0 < p`) ONCE per 4-prime orbit computes overflow-proof `0 <= x10..x13 < x_p` (`k_0[0..3]`) in `33/4 = 8.25` Ir/prime/row (disjoint GPR map `x_st=x19, x_p=x20, x_d1=x21, x_pms=x22, x_b0=x23`, zero labels, verified bit-exact at `p = 2^64 - 59` and on `Q(zeta_15)` 4-prime orbits in `cm_models.py`). Outer `(b_2..b_7)` chord resets (`<= 25` ops/prime per chord: shared `(b1_lo - 1)*w mod p` + `4` modular adds from `kt_lvl[2..7][0..3]`, `416 B` total per owner, exact `<= 25 / 232.30 = 0.108`) are conservatively bounded at `64 / 188.72 = 0.339` Ir/prime/row (`> 3.1x` headroom!) -> `8.25 + 0.339 = 8.589 <= 8.60 = C_B2_ROW_K0` Ir/prime/row:

```asm
; kernel: rev_b2_row_k0_step4 items=4 charged=8.60 key=c_b2_row_k0
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

### 1.1 Band-1 Shared `k_t` Step & Pre-Based Bucket Header (`b1a_long_lk = 203.0`, `b1a_short_lk = 70.0`, `b1bl_lk = 20.0`, `b1bs_lk = 8.0`)

Pre-adding `flat_buf_base` into `off_ptr[k] = flat_buf_base + 8 * prefix_sum[k]` at table-build time (ZERO extra RAM, ZERO extra build ops) and stepping `k_t` ONCE per prime across all 4 banks eliminates per-bank `add` base adjustments and `lsr #3` shifts (all byte counts `x_cnt = x_end - x_p = 8 * k` remain exact multiples of 8).

```asm
; kernel: rev_b1_kt_step_shared items=1 charged=3.0 key=c_b1_kt_step
subs x1, x0, x_pms
add x2, x0, x_st
csel x0, x1, x2, hs
```

```asm
; kernel: rev_b1_bkt_hdr_u1 items=1 charged=4.0 key=c_b1_bkt_hdr_u1
add x1, x_off, x0, lsl #3
ldp x_p, x_end, [x1]
cmp x_p, x_end
b.eq .Lb1_empty
```

```asm
; kernel: rev_b1_bkt_hdr_cnt items=1 charged=4.0 key=c_b1_bkt_hdr_cnt
add x1, x_off, x0, lsl #3
ldp x_p, x_end, [x1]
subs x_cnt, x_end, x_p
b.eq .Lb1_empty
```

---

### 1.2 Band-1a-Long (`[2^32, 2^50.5)`) & Band-1a-Short (`[2^50.5, 2^52.0)`): `32x -> 4x -> 1x` Stream & Apply (`c_hit_lo_m4 = 4.63`, `c_hit_lo_m4_s = 5.52`)

`rev_b1_stream_u32` (`34/32 = 1.0625` Ir/hit) + `rev_b1a_apply_u32` (`114/32 = 3.5625` Ir/hit) = **`4.625 <= 4.63` Ir/hit** (`rev_b1a_apply_u32` is ALSO reused verbatim as the Band-2 slab-owner receive kernel on `64 KiB = 256 * 256 B` full buffers with ZERO tail!).
`rev_b1_stream_u4` (`6/4 = 1.50` Ir/hit) + `rev_b1a_apply_u4` (`16/4 = 4.00` Ir/hit) = **`5.50 <= 5.52` Ir/hit** (reused as BOTH the `Band-1a-short` main loop AND the `Band-1a-long` mid-tail!).

```asm
; kernel: rev_b1_stream_wrap32 items=1 charged=4.0 key=c_b1_stream_wrap32
subs x_cnt, x_cnt, #256
b.lo .Lb1s_fix32
; --- inline rev_b1_stream_u32 (.Lb1s32) runs here while C == 1 ---
.Lb1s_fix32:
adds x_cnt, x_cnt, #256
b.eq .Lb1s_done
; --- falls into rev_b1_stream_wrap4 ---
```

```asm
; kernel: rev_b1_stream_wrap4 items=1 charged=4.0 key=c_b1_stream_wrap4
subs x_cnt, x_cnt, #32
b.lo .Lb1s_fix4
; --- inline rev_b1_stream_u4 (.Lb1s4) runs here while C == 1 ---
.Lb1s_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb1s_done
; --- falls into rev_b1_stream_tail1 (.Lb1s1) while x_cnt != 0 ---
.Lb1s_done:
```

```asm
; kernel: rev_b1_stream_u32 items=32 charged=1.065 key=c_b1_stream
.Lb1s32:
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
subs x_cnt, x_cnt, #256
b.hs .Lb1s32
```

```asm
; kernel: rev_b1_stream_u4 items=4 charged=1.50 key=c_b1_stream_u4
.Lb1s4:
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
ldp x2, x3, [x_src], #16
stp x2, x3, [x_dst], #16
subs x_cnt, x_cnt, #32
b.hs .Lb1s4
```

```asm
; kernel: rev_b1_stream_tail1 items=1 charged=4.0 key=c_b1_stream_tail
.Lb1s1:
ldr x0, [x_src], #8
str x0, [x_dst], #8
subs x_cnt, x_cnt, #8
b.ne .Lb1s1
```

```asm
; kernel: rev_b1_stream_bin21 items=1 charged=6.0 key=c_b1_stream_bin21
tbz x_cnt, #4, .Lb1s_b1
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
.Lb1s_b1:
tbz x_cnt, #3, .Lb1s_b0done
ldr x0, [x_src], #8
str x0, [x_dst], #8
.Lb1s_b0done:
```

```asm
; kernel: rev_b1a_apply_wrap32 items=1 charged=4.0 key=c_hit_lo_m4_wrap32
subs x_cnt, x_cnt, #256
b.lo .Ltail_fix32
; --- inline rev_b1a_apply_u32 (.Lb1a32) runs here while C == 1 ---
.Ltail_fix32:
adds x_cnt, x_cnt, #256
b.eq .Lb1_empty
; --- falls into rev_b1a_apply_wrap4 ---
```

```asm
; kernel: rev_b1a_apply_wrap4 items=1 charged=4.0 key=c_hit_lo_m4_wrap4
subs x_cnt, x_cnt, #32
b.lo .Ltail_fix4
; --- inline rev_b1a_apply_u4 / rev_b1b_eor_u4 / rev_b1b_sub_u4 runs here while C == 1 ---
.Ltail_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb1_empty
; --- falls into 1-record tail while x_cnt != 0 ---
.Lb1_empty:
```

```asm
; kernel: rev_b1a_apply_u32 items=32 charged=3.565 key=c_hit_lo_m4
.Lb1a32:
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
subs x_cnt, x_cnt, #256
b.hs .Lb1a32
```

```asm
; kernel: rev_b1a_apply_u4 items=4 charged=4.0 key=c_hit_lo_m4_s
.Lb1a4:
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
subs x_cnt, x_cnt, #32
b.hs .Lb1a4
```

```asm
; kernel: rev_b1a_apply_1 items=1 charged=6.0 key=c_hit_lo_m4_t1
.Lb1a1:
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
subs x_cnt, x_cnt, #8
b.ne .Lb1a1
```

```asm
; kernel: rev_b1a_apply_bin21 items=1 charged=13.0 key=c_hit_lo_m4_bin21
tbz x_cnt, #4, .Lb1a_b1
ldp x0, x1, [x_p], #16
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
.Lb1a_b1:
tbz x_cnt, #3, .Lb1a_b0done
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
.Lb1a_b0done:
```

---

### 1.3 Galois-Orbit-Specialized Band-1b (`[2^52.5, 2^56.5)` `u4 -> bin21` at `6.75` Ir/hit, `b1bl_lk = 20.0`; `[2^56.5, 2^58.0)` `u1` at `11.25` Ir/hit, `b1bs_lk = 8.0`)

Every rational prime `p === 1 (mod 15)` splits into 4 equal-density Galois-conjugate ideals `k_c in {0, 1, 2, 3}` sharing one replay table:
- `k_c = 0`: `(m - 0) & 3 = m` -> reuses `rev_b1a_apply_u4` (`4.0`) / `rev_b1a_apply_1` (`6.0`).
- `k_c = 2`: `(m - 2) & 3 = m ^ 2` -> 1-op `eor x0, x0, x_bit57` (`x_bit57 = 1 << 57`): `rev_b1b_eor_u4` (`5.0`) / `rev_b1b_eor_1` (`7.0`).
- `k_c = 1, 3`: `(m - k_c) & 3` -> 2-op `sub + and`: `rev_b1b_sub_u4` (`6.0`) / `rev_b1b_apply_1` (`8.0`).
- Exact orbit average: `u4`: `(4 + 6 + 5 + 6)/4 + 1.50 = 6.75` Ir/hit; `u1`: `(6 + 8 + 7 + 8)/4 + 4.00 = 11.25` Ir/hit.

```asm
; kernel: rev_b1b_eor_u4 items=4 charged=5.0 key=c_hit_1bl_k2
.Lb1b_eor4:
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
eor x0, x0, x_bit57
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
eor x1, x1, x_bit57
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
eor x2, x2, x_bit57
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
eor x3, x3, x_bit57
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
subs x_cnt, x_cnt, #32
b.hs .Lb1b_eor4
```

```asm
; kernel: rev_b1b_sub_u4 items=4 charged=6.0 key=c_hit_1bl_k13
.Lb1b_sub4:
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
sub x0, x0, x_k56
and x0, x0, x_mask58
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
sub x1, x1, x_k56
and x1, x1, x_mask58
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
sub x2, x2, x_k56
and x2, x2, x_mask58
ldrh w6, [x_acc, x2, lsl #1]
add w6, w6, w_lp
strh w6, [x_acc, x2, lsl #1]
sub x3, x3, x_k56
and x3, x3, x_mask58
ldrh w7, [x_acc, x3, lsl #1]
add w7, w7, w_lp
strh w7, [x_acc, x3, lsl #1]
subs x_cnt, x_cnt, #32
b.hs .Lb1b_sub4
```

```asm
; kernel: rev_b1b_eor_1 items=1 charged=7.0 key=c_hit_1bs_k2
.Lb1b_eor1:
ldr x0, [x_p], #8
eor x4, x0, x_bit57
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #8
b.ne .Lb1b_eor1
```

```asm
; kernel: rev_b1b_apply_1 items=1 charged=8.0 key=c_hit_1bs_k13
.Lb1b1:
ldr x0, [x_p], #8
sub x4, x0, x_k56
and x4, x4, x_mask58
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #8
b.ne .Lb1b1
```

```asm
; kernel: rev_b1b_eor_bin21 items=1 charged=16.0 key=c_hit_1bl_k2_bin21
tbz x_cnt, #4, .Lb1be_b1
ldp x0, x1, [x_p], #16
eor x0, x0, x_bit57
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
eor x1, x1, x_bit57
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
.Lb1be_b1:
tbz x_cnt, #3, .Lb1be_done
ldr x0, [x_p], #8
eor x0, x0, x_bit57
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
.Lb1be_done:
```

```asm
; kernel: rev_b1b_sub_bin21 items=1 charged=19.0 key=c_hit_1bl_k13_bin21
tbz x_cnt, #4, .Lb1bs_b1
ldp x0, x1, [x_p], #16
sub x0, x0, x_k56
and x0, x0, x_mask58
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
sub x1, x1, x_k56
and x1, x1, x_mask58
ldrh w5, [x_acc, x1, lsl #1]
add w5, w5, w_lp
strh w5, [x_acc, x1, lsl #1]
.Lb1bs_b1:
tbz x_cnt, #3, .Lb1bs_done
ldr x0, [x_p], #8
sub x0, x0, x_k56
and x0, x0, x_mask58
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
.Lb1bs_done:
```

---

### 1.4 Band-2 (`g`-Side `[2^58.0, 2^63.5)` Only): Deferred-Flush `128 KiB`-Aligned `8 B` Class Drains + Deterministic Tile Guard + Flush-Free `8`-Op `bfi` Emit (`c_hit_hi = 11.58`, `c_sb_hdr = 10.0`, `c_tile_hdr = 30.0`, `c_b2_row_k0 = 8.60`)

- `23` distinct `w_lp` values (`232..254`) -> `23 * 28 = 644` partial class drains/row (`~2^-51` Ir/cell); each slab buffer `x_tx_tbl[-1..28]` (`30` buffers = `3.84 MiB` total, incl. sentinel discard sinks `[-1, 28]` + read-only guard slot `[29]`) is `128 KiB`-aligned (`[base_s, base_s + 64 KiB)` primary + `[base_s + 64 KiB, base_s + 128 KiB)` overrun).
- **100% Deterministic Overflow-Free Flush Deferral (ZERO Heuristics):** `rev_sb_probe_fast` (`25` ops) computes exact `x7 = off_ptr[x6+1] - off_ptr[x4]` (total bytes across all `x21` SBs in the tile); `rev_b2_sb_wrap` (`5` ops, `c_tile_hdr = 25 + 5 = 30.0`) checks `x7 <= 65,536 B` (any giant tile with `> 8,192` hits falls back to `.Lb2_chunk_slow`, splitting each SB into `<= 32 KiB` chunks using the exact same `8`-op emit loop with `<= 0.004` Ir/hit chunk overhead). Inside `.Lb2_emit`, `tst`/`b.eq` is eliminated completely -> **`8`-op** `rev_b2_emit_c0` (`C==0`) / `rev_b2_emit_c1` (`5 + 3 = 8` ops on `C==1`, `+1` `b .Lnext_sb` at SB end counted in `c_sb_hdr = 10.0`). At `.Lnext_sb`, single-instruction `tbnz x10, #16, .Lb2_flush_sb` drains all `[base_s, x10)` (`>= 8,192` records, `0`-copy!) via `rev_b1a_apply_u32 -> u4 -> bin21` (`3.5625` Ir/hit + `<= 0.005` tail/flush + `<= 0.004` chunk = **`11.572 <= 11.62 = C_HIT_HI`**) and resets `x10 = base_s` (`bit 16 == 0`).

```asm
; kernel: rev_sb_probe_fast items=1 charged=25.0 key=c_tile_hdr_probe
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
ldr x2, [x19], #8
add x7, x6, #1
ldr x7, [x_off, x7, lsl #3]
sub x7, x7, x2
add x14, x_tx_tbl, x22, lsl #3
ldp x10, x11, [x14]
.Lempty_chord:
ret
```

```asm
; kernel: rev_b2_sb_wrap items=1 charged=5.0 key=c_b2_sb_wrap
cbz x21, .Lb2_sb_done
cmp x7, #16, lsl #12
b.hi .Lb2_chunk_slow
; --- inline rev_b2_direct_sb_hdr (.Lsb) runs here ---
tbnz x10, #16, .Lb2_flush_last
.Lb2_after_last:
str x10, [x14]
.Lb2_sb_done:
.Lb2_chunk_slow:
; --- out-of-line <=32 KiB chunked 8-op emit fallback (<=0.004 Ir/hit) for giant tiles >64 KiB ---
```

```asm
; kernel: rev_b2_direct_sb_hdr items=1 charged=9.0 key=c_sb_hdr
.Lsb:
ldr x3, [x19], #8
cmp x2, x3
b.eq .Lnext_sb
; --- inline 8-op rev_b2_emit_c0 / rev_b2_emit_c1 runs on [x2, x3), finishing with x2 == x3 ---
.Lnext_sb:
tbnz x10, #16, .Lb2_flush_sb
.Lb2_after_flush:
str x10, [x14], #8
mov x10, x11
ldr x11, [x14, #8]
subs x21, x21, #1
b.ne .Lsb
```

```asm
; kernel: rev_b2_emit_c0 items=1 charged=8.0 key=c_b2_emit
.Lb2_emit:
ldr x0, [x2], #8
adds x4, x0, x_rot
lsr x9, x4, #60
bfi x4, x9, #58, #6
b.hs .Lb2_c1
str x4, [x10], #8
cmp x2, x3
b.lo .Lb2_emit
```

```asm
; kernel: rev_b2_emit_c1 items=1 charged=4.0 key=c_b2_emit_c1
.Lb2_c1:
str x4, [x11], #8
cmp x2, x3
b.lo .Lb2_emit
b .Lnext_sb
```

```asm
; kernel: rev_b2_flush_sb items=8192 charged=0.01 key=c_b2_flush_sb
.Lb2_flush_sb:
and x_rx, x10, #0xfffffffffffe0000
sub x_cnt, x10, x_rx
mov x10, x_rx
sub x6, x14, x_tx_tbl
cmp x6, #0
b.lt .Lb2_fret
cmp x6, #224
b.ge .Lb2_fret
; --- hand [x_rx, x_rx + x_cnt) (>= 8,192 recs) to slab owner (x6 >> 3) running rev_b1a_apply_u32 -> u4 -> bin21 ---
.Lb2_fret:
b .Lb2_after_flush
```

```asm
; kernel: rev_b2_flush_last items=8192 charged=0.01 key=c_b2_flush_last
.Lb2_flush_last:
and x_rx, x10, #0xfffffffffffe0000
sub x_cnt, x10, x_rx
mov x10, x_rx
sub x6, x14, x_tx_tbl
cmp x6, #0
b.lt .Lb2_flret
cmp x6, #224
b.ge .Lb2_flret
; --- hand [x_rx, x_rx + x_cnt) (>= 8,192 recs) to slab owner (x6 >> 3) running rev_b1a_apply_u32 -> u4 -> bin21 ---
.Lb2_flret:
b .Lb2_after_last
```

---

### 1.5 Post-Indexed `43`-Op 64-Cell `ORR` Scan + `2`-Op Row Wrap (`rev_scan_zero_u64`: `43/64 = 0.6719 <= 0.675`, `c_scan_zero: 0.79 -> 0.76`)

Moving `cbz x_cnt, .Lscan_done` and `ret` into the once-per-row `2`-op wrapper `rev_scan_zero_wrap` (`2 / 188.72 = 0.0106` Ir/cell) leaves `.Lscan64` at `43` opcodes (`35` `S_g` ops x `(1 + 16/188.72)` + `8` `S_f` ops + `2/188.72` + `1,000 * sga_mean` = `0.7554 <= 0.76` Ir/cell):

```asm
; kernel: rev_scan_zero_wrap items=1 charged=2.0 key=c_scan_zero_wrap
cbz x_cnt, .Lscan_done
; --- inline 43-op rev_scan_zero_u64 (.Lscan64) runs here ---
.Lscan_done:
ret
```

```asm
; kernel: rev_scan_zero_u64 items=64 charged=0.675 key=c_scan_zero_raw
.Lscan64:
ldp x0, x1, [x_g]
ldp x2, x3, [x_g, #16]
ldp x4, x5, [x_g, #32]
ldp x6, x7, [x_g, #48]
ldp x8, x9, [x_g, #64]
ldp x10, x11, [x_g, #80]
ldp x12, x13, [x_g, #96]
ldp x14, x15, [x_g, #112]
stp x_b16, x_b16, [x_g, #16]
stp x_b16, x_b16, [x_g, #32]
stp x_b16, x_b16, [x_g, #48]
stp x_b16, x_b16, [x_g, #64]
stp x_b16, x_b16, [x_g, #80]
stp x_b16, x_b16, [x_g, #96]
stp x_b16, x_b16, [x_g, #112]
stp x_b16, x_b16, [x_g], #128
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
stp xzr, xzr, [x_f, #16]
stp xzr, xzr, [x_f, #32]
stp xzr, xzr, [x_f, #48]
stp xzr, xzr, [x_f, #64]
stp xzr, xzr, [x_f, #80]
stp xzr, xzr, [x_f, #96]
stp xzr, xzr, [x_f, #112]
stp xzr, xzr, [x_f], #128
subs x_cnt, x_cnt, #64
b.ne .Lscan64
.Lhit64:
; --- out-of-line cold-path scalar C handler bounded conservatively by 1,000 Ir/hit-block; returns to .Lscan_zero_sf ---
```

### 1.7 Explicit Pure-Scalar 1-Owner Counting-Sort / Prefix-Sum / Scatter Build (`C_LIST_A = 16.0`, `C_LBKT = 4.0`, `build_u = 104.2560 <= 104.2563`)

Pass 1 (`rev_b2_build_hist_u4`: `32/4 = 8.0` Ir/entry) executes `4x` overflow-proof modular steps (`p = 2^64 - 59`), stores `kbuf`, and accumulates byte counts (`+8`) into `hist[]`; `rev_build_prefix_u4` (`12/4 = 3.0 <= 4.0 = C_LBKT` Ir/bucket, with `+1.0` extra in `audit_cm.py`) computes `pos[s] = flat_base + sum_{t<s} hist[t]` while zeroing `hist[]` in-place via `2x stp xzr, xzr` (prepending `pos[-1] = flat_base` turns `pos[-1..n_sb-1]` directly into final `off_ptr[0..n_sb]` after Pass 2 with ZERO copy!); Pass 2 (`rev_b2_build_scatter_u4`: `28/4 = 7.0` Ir/entry) packs `(ka & 15) << 60 | (m_b << 56) | oid` in 1 `orr ..., lsl #60` and scatters `8 B` records. Total `8.0 + 7.0 + <=0.5` chord resets & tail = `15.5 <= 16.0 = C_LIST_A` (Band-1 omits `lsr #4` & `orr`, strictly faster!). Peak 1-owner build RAM `<= 14.93 EiB < 16.00 EiB` for ALL primes:

```asm
; kernel: rev_b2_build_hist_u4 items=4 charged=8.0 key=c_bld_hist
.Lbld_h4:
subs x4, x_k, x_pms
add x5, x_k, x_d
csel x0, x4, x5, hs
subs x4, x0, x_pms
add x5, x0, x_d
csel x1, x4, x5, hs
subs x4, x1, x_pms
add x5, x1, x_d
csel x2, x4, x5, hs
subs x4, x2, x_pms
add x5, x2, x_d
csel x_k, x4, x5, hs
stp x0, x1, [x_kbuf], #16
stp x2, x_k, [x_kbuf], #16
lsr x4, x0, #4
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
lsr x4, x1, #4
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
lsr x4, x2, #4
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
lsr x4, x_k, #4
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
subs x_cnt, x_cnt, #4
b.ne .Lbld_h4
```

```asm
; kernel: rev_build_prefix_u4 items=4 charged=4.0 key=c_lbkt
.Lbld_pfx4:
ldp x0, x1, [x_hist]
ldp x2, x3, [x_hist, #16]
stp xzr, xzr, [x_hist], #16
stp xzr, xzr, [x_hist], #16
add x4, x_cur, x0
add x5, x4, x1
add x6, x5, x2
stp x_cur, x4, [x_pos], #16
stp x5, x6, [x_pos], #16
add x_cur, x6, x3
subs x_cnt, x_cnt, #4
b.ne .Lbld_pfx4
```

```asm
; kernel: rev_b2_build_scatter_u4 items=4 charged=7.0 key=c_bld_scat
.Lbld_s4:
ldp x0, x1, [x_kbuf], #16
ldp x2, x3, [x_kbuf], #16
lsr x4, x0, #4
orr x5, x_om, x0, lsl #60
add x_om, x_om, #1
ldr x6, [x_pos, x4, lsl #3]
str x5, [x6], #8
str x6, [x_pos, x4, lsl #3]
lsr x4, x1, #4
orr x5, x_om, x1, lsl #60
add x_om, x_om, #1
ldr x6, [x_pos, x4, lsl #3]
str x5, [x6], #8
str x6, [x_pos, x4, lsl #3]
lsr x4, x2, #4
orr x5, x_om, x2, lsl #60
add x_om, x_om, #1
ldr x6, [x_pos, x4, lsl #3]
str x5, [x6], #8
str x6, [x_pos, x4, lsl #3]
lsr x4, x3, #4
orr x5, x_om, x3, lsl #60
add x_om, x_om, #1
ldr x6, [x_pos, x4, lsl #3]
str x5, [x6], #8
str x6, [x_pos, x4, lsl #3]
subs x_cnt, x_cnt, #4
b.ne .Lbld_s4
```

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

For small composite ECM branches (`<= 128` bits, `L_min = 2` in `q15_ecm16.py`, eliminating the `5`-limb floor inflation), `mulmod_instr(nbits, l_min=2)` charges `MONT5_desc * (2/5)^2 = 367.76 * 4 / 25 = 58.8416` Ir (even under `MONT5 = 280.0`, `L = 2` CIOS executes in **`47 <= 48.0 < 58.84` Ir**!). Inputs `(a0, a1, b0, b1, n0, n1, n0inv)` reside in AAPCS64 caller-saved registers `x0..x6` (zero stack spills, zero callee-saved registers, zero banned ops, disjoint GPRs `x0..x15`), returning `a * b * 2^{-128} mod N` in `x0:x1` (`0 <= t < 2N < 2^129` after both rounds for every 128-bit odd modulus `N <= 2^128 - 1`, verified bit-exactly LIVE in `cm_models.verify_kernels_symbolic` across edge moduli `2^128 - 1, 2^127 + 1, 3` and `64` random 128-bit odd moduli):

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

## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Deterministic Stepped-GCD `H-ECM`

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

- **Parsed 3-Prime `< 2^62` Unrolled Radix-4 Lazy Shoup NTT (`ntt62_r4_u2`: `86` opcodes / `16` elem-stages = `5.375 <= 5.38` Ir/elem/stage/prime -> `9 * 5.375 = 48.375 <= 48.50` in `ntt3_F`), Schoolbook Small-`M` Crossover (`<= 256` limbs), `100*n` Garner CRT, & 64-Bit Lehmer GCD (`lehmer32_step`: `29/2 = 14.5 <= 16.0`):**
  Three `< 2^62` primes divisible by `3 * 2^55`: `p_1 = 69 * 2^55 + 1`, `p_2 = 57 * 2^55 + 1`, `p_3 = 54 * 2^55 + 1` (`4*p_1 < 2^{63.11} < 2^64`, `prod = 2^182.696 > 2^54 * (2^64 - 1)^2`). Pairing 2 butterflies via `4x ldp`/`4x stp` and computing `x12 = (x1 + 2p) - x3 in [0, 4p) < 2^{63.11}` in 2 ops (Lazy Shoup output `< (1 + 4p/2^64)p < 1.54p < 2p`, verified LIVE across all 3 primes and `n = 64` convolution!) cuts `ntt3_F(M)` to `min(2*m^2 + 4*m if m<=256 else inf, 48.50*n2*k + (26*n2 if k odd) + 100*n2, 48.50*n3*(k-2) + (26*n3 if k-2 odd) + 244*n3)` (`+ 512*n` if `n > 2^26`; schoolbook `mul_LL(m/2, m/2) = 8*(m/2)^2 + 8*(m/2) = 2*m^2 + 4*m`). Stage-B `c_leaf` uses 64-bit Lehmer GCD (`lehmer32_step` `29/2 = 14.5 <= 16.0` ops/step + `L`-limb linear combination `<= 64*L` -> `(nbits/24)*(64*L + 544)`, `>= 1.48x` over worst-case Fibonacci & `2.2x` over random 680b pairs, verified LIVE!), dropping **`c_cof_2d_upper` from `1.0485 -> 0.7445` Ir/cell**:

```asm
; kernel: ntt62_r4_u2 items=16 charged=5.38 key=ntt62_r4_elem_stage
.Lntt_r4_u2:
ldp x0, x14, [x_a0]
ldp x1, x15, [x_a1]
ldp x2, x16, [x_a2]
ldp x3, x17, [x_a3]
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
ldp x4, x5, [x_tw], #16
mul x6, x15, x4
umulh x7, x15, x5
msub x15, x7, x_p, x6
ldp x4, x5, [x_tw], #16
mul x6, x16, x4
umulh x7, x16, x5
msub x16, x7, x_p, x6
ldp x4, x5, [x_tw], #16
mul x6, x17, x4
umulh x7, x17, x5
msub x17, x7, x_p, x6
add x8, x0, x2
subs x9, x8, x_p2
csel x8, x9, x8, hs
subs x10, x0, x2
add x9, x10, x_p2
csel x10, x9, x10, lo
add x11, x1, x3
subs x9, x11, x_p2
csel x11, x9, x11, hs
add x12, x1, x_p2
sub x12, x12, x3
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
add x8, x14, x16
subs x9, x8, x_p2
csel x8, x9, x8, hs
subs x10, x14, x16
add x9, x10, x_p2
csel x10, x9, x10, lo
add x11, x15, x17
subs x9, x11, x_p2
csel x11, x9, x11, hs
add x12, x15, x_p2
sub x12, x12, x17
mul x6, x12, x_j
umulh x7, x12, x_jp
msub x12, x7, x_p, x6
add x14, x8, x11
subs x9, x14, x_p2
csel x14, x9, x14, hs
subs x16, x8, x11
add x9, x16, x_p2
csel x16, x9, x16, lo
add x15, x10, x12
subs x9, x15, x_p2
csel x15, x9, x15, hs
subs x17, x10, x12
add x9, x17, x_p2
csel x17, x9, x17, lo
stp x0, x14, [x_a0], #16
stp x1, x15, [x_a1], #16
stp x2, x16, [x_a2], #16
stp x3, x17, [x_a3], #16
subs x_cnt, x_cnt, #2
b.ne .Lntt_r4_u2
```

```asm
; kernel: lehmer32_step items=2 charged=16.0 key=lehmer32_step_op
.Lleh_loop:
cbz x1, .Lleh_done
udiv x9, x0, x1
msub x10, x9, x1, x0
madd x11, x9, x4, x2
madd x12, x9, x5, x3
cmp x10, x12
b.lo .Lleh_done
sub x8, x1, x10
add x13, x4, x11
cmp x8, x13
b.lo .Lleh_done
cbz x10, .Lleh_done
udiv x9, x1, x10
msub x14, x9, x10, x1
madd x15, x9, x11, x4
madd x16, x9, x12, x5
cmp x14, x15
b.lo .Lleh_done
sub x8, x10, x14
add x13, x12, x16
cmp x8, x13
b.lo .Lleh_done
mov x0, x10
mov x1, x14
mov x2, x11
mov x3, x12
mov x4, x15
mov x5, x16
b .Lleh_loop
.Lleh_done:
```
- **Stage B (`P_min = 32`, `KH = 0.08`, Static 3-Class `T_g[c] = (299, 287, 275)` + Constant `T_f_const = 118`b Pre-Gates, `step = 0.25`):** `gbatch_ntt3(pbits_split(32, 0.08), leaf=L, c_leaf=10 * MONT5 * (L/5)^2 + (nb/24) * (64 L + 544))` with AAPCS64 `MONT5 = 280.0`, `ntt62_r4_u2`, schoolbook `M <= 256` crossover, and `lehmer32_step` 64-bit Lehmer GCD -> exact **`c_cof_2d_upper = 0.7445` Ir/cell** (analytic Chebotarev split density `5/72 = 0.06944 < 0.08`, live `2^14` check `0.07197 / 0.05649 <= 0.08`).
- **Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree `H-ECM` Ladder (`q15_ecm16.py` & `cm_models.ecm_mean_var`):** Before Suyama ECM, an exact prime-exponent Newton integer perfect-power precheck (`_is_perfect_power(N)`, explicitly charged at `NEWTON_PRECHECK_INSTR = 2^26` Ir > 32.06M schoolbook ops for 2,590-bit boot inputs) eliminates all `N = p^k` inputs. In Stage 1, every curve steps individually by each prime factor `q` across all prime powers `q^k <= B1` (`<= pi_pp_upper(B1) = pi_upper(B1) + 2*sqrt(B1)` steps, checking `gcd(Z, N)` after every step so distinct valuations `e_1 < e_2` of `q` split cleanly); in Stage 2, every giant step `G` builds a `J`-leaf binary product tree (`h_tree = ceil(log2(max(2, J)))` levels) replacing any `t_j == 0 mod N` leaf with `1` BEFORE multiplying into the tree so every leaf satisfies `gcd(t_j, N) < N`, giving deterministic per-curve cost `(ops + 1.1*J*G) * mul + (pi_pp_upper(B1) + (1 + 2*h_tree)*G) * gcd` (ZERO curve replay, `var_t = 0.0`): whenever `gcd(root, N) == N` in a giant step, descending the binary product tree (`<= 2*h_tree` gcds, `0` extra muls) is guaranteed to isolate a proper factor `1 < d < N` (even when another leaf `t_j == 0 mod N` exists). Under the standard random-integer ECM group-order heuristic `H-ECM` (verified live on 160 composite proper-split trials `N = p * q_comp` counting only `1 < g < N`, plus a live Stage-2 `W = 210, J = 24, h_tree = 5` binary-product-tree giant-step bisection on `N = 100003 * 100043` and `[0, p, q, 1]` zero-leaf-skip test in `q15_ecm16.verify_suyama_stepped`), each rung `(r, lo_r)` of the `12`-bit geometric ladder targets `s_rung = 1 - (1 - success) / n_rungs` with `curves = ceil(log1p(-s_rung) / log1p(-p1))`. Exact DP `tree_branch_caps(nbits, b_min, l_min=2)` bounds `sum L_i^2` and `sum gcd_instr(b_i)` across ALL active composite-branch partitions (`mont2_cios` `47 <= 48` ops, charged `58.84` Ir). (i) For **post-`2^48` collection cofactors** (`cm_models.ecm_mean_var`, `lo_r = 48`, `b_min = 96`, `success = 0.999999` per cofactor), deterministic truncation at `n_rels_upper = 2^63.565` yields `ew_ecm_live = 96.9769 <= 96.9930` u, **`var_t = 0.0`**, with per-relation yield retention `1 - 29*10^-6 = 0.999971` folded directly into collection yield; (ii) for **descent cofactors** (`q15_ecm16.py`, `lo_bits = 16`, `b_min = 32`, `l_min = 2`, `success = 0.999`), total ladder cost is 100% deterministic.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Dense + 16-Shard `< 2^64 B` 4-Vec 4-Col Arch B Sign-Split Scalar ARMv8.0-A + `u32->u4->bin21` Matrix Streaming + Exact Single-Sub `57`-Op Barrett `row_red5_z12`)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`, `tab_ratio_charged = 1.00` plain sparse with ZERO `TABM`, `32` sequential 4-vector bundles of `128 B/slot`, `n = 128 = 32 * 4`, `m = 32,768`, peak single-owner dense-vector RAM `2 * N_sw * 128 B = 2^63.368 B = 10.323 EiB < 2^64 B`, `16` read-only sparse-matrix shard owners `4.287 EiB << 2^64 B` each, explicit priced cross-owner burst streaming via `rev_b1_stream_u32 -> u4 -> bin21`: `E[w_i] <= 212.63875` words/row, `E[ops(w_i)] <= 1.0625 * 212.63875 + 24.0 = 249.93 <= 250.0 = c_stream_mat_1vec` and `<= 63.0 * 4` (`c_stream_mat = 63.0` Ir/vec-row in 4-vec Krylov)!):
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

- **`row_red5_z12` (`57 <= 58.0` Ir/bucket-vec — Exact Single-Sub Signed 5-Limb Barrett Reduction for ALL `|A| <= 2^{290}`):**
  Handles ANY signed 5-limb bucket accumulator `A = (x10..x14)` with `|A| <= 2^{290}` (covering EVERY `W_max = 4,096, |c| <= 255` row up to `4,096 * 255 * (R_BLS - 1) < 2^{275}` by `> 10^4x`!). Since `mu_r = floor(2^318 / R_BLS) = 0x8d54253b7fb78ddf` and `top = floor(|A| / 2^{254}) <= 2^{36}`, exact error `|A|/R_BLS - q_hat < 2^{254}/R_BLS + top/2^{64} + 1 < 1.5522 < 2` unconditionally -> `rem = U - q_hat * R_BLS in [0, 2 * R_BLS)` ALWAYS, so ONE single conditional `- R_BLS` + branchless sign fix reduces to exact `[0, R_BLS)` in **`57 <= 58.0`** opcodes (`6 * 57 = 342` ops/row):

```asm
; kernel: row_red5_z12 items=1 charged=58.0 key=row_red5_barrett
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

- **Sign-Split 2-Vec `gen_run1_2v_z12` (`c_b > 0`) & `gen_run1_2v_neg_z12` (`c_b < 0`) (`43/2 = 21.50 <= 22.0` Ir/vec-entry, `24 <= 29` GPRs avoiding `x18/x29/x30`, `0` spills; caller-side non-empty `cbz x24` precondition guarantees `x24 >= 1` on entry):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`, provisioned up to `nonunit_ub = 0.11` covering every `2^18..2^25` ladder point `0.10915 -> 0.08179` and `1 - (eu_25 - 2*se) = 0.08212`, stored as 1-byte `uint8_t` magnitudes loaded and zero-extended to `x22` via `ldrb w22, [x23], #1`) into positive and negative sub-lists, executed as two spill-free 2-vector passes (`(v0, v1)` at offsets `#0..#48`, then `(v2, v3)` at offsets `#64..#112`, explicit `add x8, x9, x8` from `x_base = x9`): `3` ptr/coeff loads + `2 * 19` mul/umulh/carry ops + `2` loop = **`43 / 2 = 21.50 <= 22.0` Ir/vec-entry** (`k5_live = eu_lo * 7.65 + (1 - eu_lo) * 22.0 = 8.8284 <= k5_charged = 9.45` Ir, `+0.62` Ir headroom; even at `nonunit_ub = 0.11`, `0.89 * 7.65 + 0.11 * 22.0 = 9.2285 <= 9.45`!):

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

- **Explicit 1-Vector Horner Kernels (`unit_run_1vec4_z12` / `_neg`: `38/4 = 9.50 <= 9.60` Ir/entry; `gen_run1_1v_z12` / `_neg`: `24.0 <= 24.0` Ir/entry with caller-side non-empty `cbz x24` precondition `x24 >= 1` -> `k5_horner_charged = 14.00`):** Reuses the EXACT SAME read-only sparse-matrix shards (`128*j` stored byte offsets) on the `32 B/slot` 1-vector Horner dense buffer via **`add x7, x6, x7, lsr #2`** (`(128*j) >> 2 = 32*j` exact byte offset, ZERO matrix re-encoding, ZERO offset-shift bug!), simulated bit-exact in `cm_models.py`:

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

- **`chi_mul_z12_fast` (`50` Ir), Non-Double-Counted `row_tail_6bucket_z12` (`240` Ir = `12` unit sub-lists x (`4/4` wrap + `<=3*0.35` tail delta above `7.65`) + `12*0.5` non-unit header + `2*50` `chi_mul` + `56` header/zero/store + `53` slack), 1-Dense + 16-Projection + 16-Shard 4-Vec Arch B Row Overhead (`c_row_chi = 600.0` Ir, `c_stream_mat = 63.0` Ir/vec-row Krylov, `c_stream_mat_1vec = 250.0` Ir/row 1-vec Horner via `u32->u4->bin21` `1.0625*212.639 + 24.0 = 249.93 <= 250.0`, `nonunit_ub = 0.11`, `k5_charged = 9.45`, `k5_horner_charged = 11.20`, `c_xproj = 65.0` Ir/vec-row `>= 64.5`, `n_owners = 1`, `n_proj_owners = 16`, `n_mat_owners = 16`, `te = 0`):** `6 * 57 + 240 + 16 = 598 <= 600.0 = c_row_chi` Ir worst-case!

## 5. Descent Kernels (Generous Scalar Bounds for Uncapped `dsieve16` & Exact-Integer `BOUND_CELL`, Plus Input-Uniform Bareiss LLL)

Descent is charged at `2^107.5362` u (`live = 2^107.536063` u, `MODEA_V03_A3_Q25_C4`, `76` nodes). Rather than hand-scheduling fragile cold-path assembly for 16D Fincke-Pohst DFS or interval ball arithmetic, we charge generous, manifestly-conservative scalar ARMv8.0-A operation bounds and verify every mathematical property end-to-end in exact integer Python:
- **Uncapped `dsieve16` on `S_16(R) = B_16(R)/{+-1}` (`SIEVE_HIT = 1,400` Ir/hit, ZERO per-prime skip caps):** Every relative-degree-1 prime ideal with `p in [2^8, 2^68.5]` is reduced via Bareiss LLL (`PAIR_SETUP = 2^35.2`) and enumerated completely over `B_16(R)/{+-1}` using plain compiled scalar C (`-march=armv8-a -mgeneral-regs-only`) 16D Fincke-Pohst Q32 DFS + leaf vector `v = V_stk[0] + z_0 * b_0` + `v = 0` exclusion + ambient sign canonicalization (`v -> -v` if first nonzero ambient coordinate `< 0`) + 128-bit bijective half-ball rank `[0, C_R - 1]` via `129.03 MiB` flat prefix table + owner routing + `uint16_t` score add. Explicit 16D Fincke-Pohst measurements on 99/100-LLL-reduced prime lattices confirm `edge/hit = 2.24` (`1.24` internal nodes/hit) at `R = 12, p ~ 2^28` (`8.0e7` hits) and `edge/hit = 2.13` (`1.13` internal nodes/hit) at `R = 16, p ~ 2^38` (`7.9e6` hits), rising to `4.64` only in the sparse small-hit edge regime (`R = 20, p ~ 2^52`, `1.7e4` hits; all `76` nodes use `R <= 104`, where sparse DFS edge traversals `~2,000` Ir/prime at smaller radii are absorbed by `> 10^6x` inside `PAIR_SETUP = 2^35.2` Ir/prime, which carries `3.65 billion` Ir/prime slack over the `2^35.06` Bareiss LLL bound!):
  - **Per-hit breakdown (`<= 914–1,244 <= 1,400` Ir/hit):** DFS traversal (`<= 1.25` internal nodes/hit x `120` Ir + leaf check `80` Ir = `<= 230` Ir; `<= 560` Ir even at sparse edge `4.64`) + leaf vector materialization (`64` Ir) + ambient sign canonicalization (`40` Ir) + 16D half-ball bijective ranking (`16` table lookups + 128-bit adds = `320` Ir) + owner routing / `uint16_t` update / compiler spill headroom (`260` Ir) = **`<= 914–1,244 <= 1,400` Ir/hit** (`>= 156` Ir headroom even at worst sparse edge `4.64`, `+486` Ir in dense regime). Exact 4D Fincke-Pohst half-ball DFS (100% match vs brute force on 18,174 points), ambient sign canonicalization (explicitly verifying witness `z = (-1, 1, 0, 0) -> v = (-3, 3, 1, 0) -> (3, -3, -1, 0)`), 55-bit `SQ_LUT`, and bijective 16D half-ball ranking (`552/552` canonical + sign-flipped `R = 2` hits + `v = 0` exclusion) are verified LIVE in `cm_models.py`.
- **Certified Outward-Rounded Integer `BOUND_CELL (+6)` (`BOUND_CELL = 4,200` Ir/cell, `q15_bound_cell.py`):** Evaluates strict signed-64b outward-rounded complex ball arithmetic (`bmul_6` with `+6 >= ceil(sqrt(5) + 3)` rounding margin, `13` `bmul_6` calls/place across `4` complex embeddings, `<= +0.5343b < +0.625b` `SQ_LUT` envelope). **Per-cell breakdown:** `3,792` Ir exact integer ball-mul core (`16 + 4*(104 + 13*48 + 128 + 88)`) + `128` Ir vector materialization + `80` Ir `SQ_LUT` integer square root + `64` Ir `log2` & `(nf_max, ng_max)` cap compare + `136` Ir compiler spill headroom = **`4,200` Ir/cell**. Deterministically rejects any cell with `ub_f_rem > nf_max` or `ub_g_rem > ng_max` (guaranteeing ZERO over-norm cells ever enter ECM). Exact `bmul_6` (`+6`) outward rounding (including 50-digit `Decimal` worst-case witness) and 12-cell real-lattice exact integer `BOUND_CELL` certification (`bc_gmin >= 0.0`) are verified LIVE in `cm_models.py` and `verify_q15_descent.py`.
- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every single special-`q` node lattice `L_Q` (`det = q <= 2^260 < 2^276`) and `G_T` boot lattice runs exact 2-phase `(11/40 Euclidean -> 99/100 Q_CM = diag(G_8, G_8))` Bareiss LLL (`lll_16d_int(..., use_cm=True)` in `cm_models.py`, verified live on both `Bsq` and `I_16`; `dsieve16` `L_p` uses Euclidean `lll_16d_int(..., use_cm=False)`), where `G_8` has `det(G_8) = 3^4 * 5^6` and `lambda_max(G_8) = 15` (proved live via Sylvester minors of `15*I_8 - G_8`, so `d_k^(CM) <= 15^k * d_k^(Euc)` unconditionally!), yielding explicit derived uniform bounds **`PAIR_SETUP = 2^35.2`** (70-bit `L_p` Euclidean `11/40 -> 99/100`, `2^35.06 < 2^35.2`), **`LLL_NODE = 2^40.0`** (conservative 276-bit envelope covering every single special-`q` `q <= 2^260` under `11/40 Euc -> 99/100 Q_CM`, `L <= 10` limbs `<= 640`b, `2^39.46 < 2^40.0`), and **`LLL_BOOT = 2^47.0`** (381-bit `G_T` boot WITH `S_boot = 2^21` scaling on the 8 `b`-coords BEFORE Phase 1: `det = p^4 * 2^168 <= 2^1692`, `L_d <= 58 <= 59` limbs -> **`2^46.75 < 2^47.0`**), all derived and asserted LIVE in `verify_kernels_symbolic` alongside `BOUND_CELL = 4,200`, `SIEVE_HIT = 1,400`, `SCAN = 10`, `ZERO = 2`, `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`.
