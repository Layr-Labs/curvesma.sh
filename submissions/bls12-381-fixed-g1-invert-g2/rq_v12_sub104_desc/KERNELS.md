# Kernels (Raw-Opcode-Verified Scalar ARMv8.1-A with Disclosed `FEAT_LSE STADDH`: `Arch2-Prod-ZeroSlice-BP2-W28-2Pass-U4-Buf128K-B2-U32-Bin21-Defer8-LSE`, 3-Prime Unrolled Radix-4 Lazy Shoup NTT + Lehmer32 GCD, 1-Dense + 16-Proj + 16-Shard 4-Vec LA, Explicit Build & Single-Sub Barrett Kernels)
All hot-loop kernels use scalar general-purpose-register ARMv8.1-A integer instructions (`-march=armv8.1-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `DC ZVA`, zero `fp64`, zero `x18`, zero `bl`/`blr` — all kernels are inline macros or leaf routines returning via `ret`), plus **one disclosed scalar `FEAT_LSE` atomic integer instruction** from ARMv8.1-A, priced at `1` scalar instruction = `1 / 20,000` units under `COST_MODEL.md`:
- **`STADDH Ws, [Xn]`** (32-bit alias of `LDADDH Ws, WZR, [Xn]`; scalar registers; base-register `[Xn]` addressing only), used exclusively for `uint16_t` score updates on Symmetric Two-Pass 2-bank-pair `{0,2}/{1,3}` slab owners (`15.375 EiB < 16.000 EiB`): shared fixed-base Band-1 accumulator `S_b1` at `0x0000_0000_0000_0000` (`0.500 EiB`), `1`-bit `g`-survivor bitmap `M_g[0..27]` at `0x0800_0000_0000_0000` (`0.875 EiB`), and shared Band-2 superbucket slab `S_b2[0..27]` at `0x2000_0000_0000_0000..0xFFFF_FFFF_FFFF_FFFF` (`14.000 EiB`, slot `col_in_sb in 0..27` at `(col_in_sb + 4) << 59`); bit-identical to `ldrh; add; strh` on every `uint16_t` cell with zero wrap-around (`<= 35,767 <= 65,535`, verified live in `audit_cm.py`). All other `FEAT_LSE` opcodes (`LDADD*`, `SWP*`, `CAS*`) remain forbidden.

Cold paths (`.Lhit64_g`, `.Lhit128_f`, `.Lb2_chunk_slow`, 16D descent DFS / `BOUND_CELL`) use conservative scalar bounds. `audit_cm.py` (`kernels_symbolic`) parses all `73` fenced `asm` blocks (`0` undef labels, disjoint GPR maps), validates opcodes & `LDP`/`STP`/`STADDH`, counts instructions (`cmp`+`b.*`=2 ops), and simulates: (i) counting-sort build (`rev_b2_build_hist_u4` `12.0` + `rev_b2_build_scatter_u4` `7.0` (`bfi`) + resets `<= 19.5 = C_LIST_A_B2`, Band-1 `<= 16.0`, `rev_build_prefix_u4` `3.0 <= 4.0`), (ii) Two-Pass `{0,2}/{1,3}` fixed-base `S_b1` `STADDH` Band-1a (`2.625 <= 2.63`, `3.50 <= 3.52` on both `f`/`g`; `b1a_long_lk = 206.0`, `b1a_short_lk = 73.0`) & 1-bit-XOR `STADDH` Band-1b (`4.00` on both `f`/`g`, `b1bl_lk = 40.0`), (iii) `W_SB = 28` `+4<<59` `128 KiB`-buffered Band-2 on both `f` (`[55.5, 60.5)`) & `g` (`[55.5, 63.5)`) (`rev_b2_row_k0_step4` `8.25 + 0.442 <= 8.70`, exact 64-bit div-by-28 `M7 = 0x924924924924924a` `32 + 5 = 37.0`, `u4` `<= 7.08` / `c_sb_hdr_u4 = 16.0`, `u2` `<= 7.58` / `14.5`, `u1` `<= 9.08` / `10.0`), (iv) Two-Pass `M_g` scan+zero with closed periodic `S_b1` handoff (`1.0625 + 0.2891 + 0.15625` pad `<= 1.50`), (v) `ntt62_r4_u2` (`5.375 <= 5.38`, `48.50`/stage), schoolbook `<= 256`, `100*n` CRT, `lehmer32_step` (`14.5 <= 16.0`), (vi) 4-vec/1-vec LA with `57`-op `row_red5_z12` (`|A| <= 2^{290}`, `c_row_chi = 600.0`, `c_stream_mat = 63.0 / 251.0`, `k5 = 9.45 / 11.25`, `c_xproj = 65.0`), AAPCS64 `MONT5 = 280` (`279` ops), `47`-op `mont2_cios`, and `BOUND_CELL = 4,200`.

## 1. Collection Kernels (`Arch2-Prod-ZeroSlice-BP2-W28-2Pass-U4-Buf128K-B2-U32-Bin21-Defer8-LSE`, Product `L2(a) x L1(b)`, `W_tile = 448 = 7 * 64`, `W_SB = 28`, `P_min = 2^32`, `P_m4s = 2^51.5`, `P_m4 = 2^52.5`, `P_lo = 2^55.5`, `P_4u = 2^60.5`, `P_2u = 2^61.5`, `P_f = 2^60.5`, `P_g = 2^63.5`)
Product `L2(a) x L1(b)` (`A^2=299,400`, `r_ab=1.065`, `R_1=670.79`, `V_a=2^56.9021=0.9344*2^57`, `Vol_8=2^64.289`, `V_{b,D12}=Vol_8/12=2^60.704`, `L_max=335.39`, `+27 pad<14*28=392`, `L_b=144.89`, `V_kept=4V_a*V_{b,D12}=2^119.6059`): ZERO `q_key`/`q_tab`/bsearch. Both `g` (`[55.5,63.5)`) and `f` (`[55.5,60.5)`) reuse `14` slab owners/pair (`x_tx_tbl[0..13]`, `15.375 EiB < 16 EiB`): `S_b1` (`0.500 EiB` at `0x0`), `M_g[28]` (`0.875 EiB` at `0x0800...`), `S_b2[28]` (`14.000 EiB` at `(col_in_sb+4)<<59`). Pass 1 runs `g` (`rev_scan_g_u64` writes `M_g`, re-arms `S_b1` to `bias_g` on `r<r_last`, `0` on `r==r_last`); Pass 2 runs `f` (`rev_scan_f_u128` skips `>99.6%` zero `M_g`, restores `bias_g` on `r==r_last`).

### 1.1a `rev_b2_row_k0_step4` (Band-2 4-prime orbit incremental `D_(omega, sigma)` row-step `k_0` init: `33/4 = 8.25` Ir/prime/row + `<= 0.442` outer-chord reset <= `C_B2_ROW_K0 = 8.70` Ir/prime/row)
Because `sigma(omega) = omega` (`7*5 === 4*5 === 13*5 === 5 mod 15`), all `4` primes `s_k = s_0^(7^k)` (`k = 0..3`) in a Band-2 prime orbit share the EXACT same modulus `p` (`2^55.5 <= p < 2^63.5`), `pwk[0] = 1`, and `pwk[1] = w = s_0^5 mod p`. Moreover, `in_D12(b)` depends ONLY on `(b_2..b_7)`, and `L_1(b)` is convex along `b_1` with `|sigma_k(b_1 y)| = |b_1|`, so `b_1` sweeps an uninterrupted contiguous integer chord of the same harmonic-mean length `L_b1 = L_b = 144.89` rows. Maintaining `k_tail[0..3]` (`32 B` in `[x_st]`) across consecutive `b_1 -> b_1 + 1` rows and converting signed `x_b0` (`|x_b0| <= L_max = 335.39 << 2^55.5 <= p`, precondition `-p < x_b0 < p`) ONCE per 4-prime orbit computes overflow-proof `0 <= x10..x13 < x_p` (`k_0[0..3]`) in `33/4 = 8.25` Ir/prime/row (disjoint GPR map `x_st=x19, x_p=x20, x_d1=x21, x_pms=x22, x_b0=x23`, zero labels, verified bit-exact at `p = 2^64 - 59` and on `Q(zeta_15)` 4-prime orbits in `cm_models.py`). Outer `(b_2..b_7)` chord resets (`<= 25` ops/prime per chord: shared `(b1_lo - 1)*w mod p` + `4` modular adds from `kt_lvl[2..7][0..3]`, `416 B` total per owner, exact `<= 25 / 144.89 = 0.173`) are conservatively bounded at `64 / 144.89 = 0.442` Ir/prime/row (`> 2.5x` headroom!) -> `8.25 + 0.442 = 8.692 <= 8.70 = C_B2_ROW_K0` Ir/prime/row:
```asm
; kernel: rev_b2_row_k0_step4 items=4 charged=8.70 key=c_b2_row_k0
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
### 1.1 Band-1 Shared `k_t` Step & Pre-Based Bucket Header (`b1a_long_lk = 206.0`, `b1a_short_lk = 73.0`, `b1bl_lk = 40.0` (summed across both `{0,2}`/`{1,3}` owners))
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
### 1.2 Shared Fixed-Base `S_b1` Band-1a (`[2^32, 2^51.5)` & `[2^51.5, 2^52.5)`): Disclosed `STADDH` `32x -> 4x -> bin21` Stream & Apply (`2.63 / 3.52` Ir/hit on BOTH `f` and `g`!)
Because Symmetric Two-Pass runs `g`-side Band-1 (`Pass 1`) and `f`-side Band-1 (`Pass 2`) sequentially into the **EXACT SAME** fixed-base accumulator `S_b1` at `base = 0x0000_0000_0000_0000` (`2` banks x `2^57` x `2 B = 0.500 EiB`), **BOTH** `f`-side and `g`-side Band-1 `a`-table entries store the exact 64-bit `uint16_t*` pointer `(m_bit << 58) | (2 * oid)` (`m_bit in {0, 1}`, `oid < 2^57`, `0` extra build ops, ZERO per-hit `add x_acc` on `g`!). `rev_b1_stream_u32` (`1.0625`) + `rev_b1a_apply_f_u32` (`1.5625` Ir/hit, **ALSO reused verbatim as the Band-2 `S_b2` slab-owner receive kernel!**) = **`2.625 <= 2.63` Ir/hit** on both `f` and `g`; `u4`: `1.50 + 8/4 = 3.50 <= 3.52` Ir/hit on both `f` and `g`.
```asm
; kernel: rev_b1_stream_wrap32 items=1 charged=4.0 key=c_b1_stream_wrap32
subs x_cnt, x_cnt, #256
b.lo .Lb1s_fix32
.Lb1s_fix32:
adds x_cnt, x_cnt, #256
b.eq .Lb1s_done
```
```asm
; kernel: rev_b1_stream_wrap4 items=1 charged=4.0 key=c_b1_stream_wrap4
subs x_cnt, x_cnt, #32
b.lo .Lb1s_fix4
.Lb1s_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb1s_done
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
.Ltail_fix32:
adds x_cnt, x_cnt, #256
b.eq .Lb1_empty
```
```asm
; kernel: rev_b1a_apply_wrap4 items=1 charged=4.0 key=c_hit_lo_m4_wrap4
subs x_cnt, x_cnt, #32
b.lo .Ltail_fix4
.Ltail_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb1_empty
.Lb1_empty:
```
```asm
; kernel: rev_b1a_apply_f_u32 items=32 charged=1.565 key=c_hit_lo_m4_f
.Lb1af32:
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
staddh w_lp, [x2]
staddh w_lp, [x3]
ldp x4, x5, [x_p], #16
staddh w_lp, [x4]
staddh w_lp, [x5]
ldp x6, x7, [x_p], #16
staddh w_lp, [x6]
staddh w_lp, [x7]
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
staddh w_lp, [x2]
staddh w_lp, [x3]
ldp x4, x5, [x_p], #16
staddh w_lp, [x4]
staddh w_lp, [x5]
ldp x6, x7, [x_p], #16
staddh w_lp, [x6]
staddh w_lp, [x7]
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
staddh w_lp, [x2]
staddh w_lp, [x3]
ldp x4, x5, [x_p], #16
staddh w_lp, [x4]
staddh w_lp, [x5]
ldp x6, x7, [x_p], #16
staddh w_lp, [x6]
staddh w_lp, [x7]
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
staddh w_lp, [x2]
staddh w_lp, [x3]
ldp x4, x5, [x_p], #16
staddh w_lp, [x4]
staddh w_lp, [x5]
ldp x6, x7, [x_p], #16
staddh w_lp, [x6]
staddh w_lp, [x7]
subs x_cnt, x_cnt, #256
b.hs .Lb1af32
```
```asm
; kernel: rev_b1a_apply_f_u4 items=4 charged=2.0 key=c_hit_lo_m4_s_f
.Lb1af4:
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
staddh w_lp, [x2]
staddh w_lp, [x3]
subs x_cnt, x_cnt, #32
b.hs .Lb1af4
```
```asm
; kernel: rev_b1a_apply_f_bin21 items=1 charged=7.0 key=c_hit_lo_m4_bin21_f
tbz x_cnt, #4, .Lb1af_b1
ldp x0, x1, [x_p], #16
staddh w_lp, [x0]
staddh w_lp, [x1]
.Lb1af_b1:
tbz x_cnt, #3, .Lb1af_b0done
ldr x0, [x_p], #8
staddh w_lp, [x0]
.Lb1af_b0done:
```
### 1.3 Shared Fixed-Base 2-Bank-Pair `{0,2}/{1,3}` 1-Bit-XOR `STADDH` Band-1b (`[2^52.5, 2^55.5)` `u4 = 4.00` on BOTH `f` and `g`, `b1bl_lk = 40.0` across 2 owners)
Partitioning the `4` Galois banks into two disjoint bank-pairs `Owner E = {0, 2}` (`m_bit = 0` for bank `0`, `1` for bank `2`) and `Owner O = {1, 3}` (`m_bit = 0` for bank `1`, `1` for bank `3`) guarantees that under EVERY Galois conjugate `k_c in {0, 1, 2, 3}`, `(m - k_c) mod 4` maps `{0, 2}` and `{1, 3}` bijectively to one single destination bank-pair owner with ZERO cross-pair mixing:
- **Identity (`50%` of conjugates: `k_c = 0` on E/O, `k_c = 1` on O->E, `k_c = 3` on E->O):** reuses `rev_b1a_apply_f_u4` (`2.0`) on both `f` and `g`.
- **1-Bit XOR (`50%` of conjugates: `k_c = 2` on E/O, `k_c = 1` on E->O, `k_c = 3` on O->E):** `m_bit ^ 1` via 1-op `eor x0, x0, x_bit58` (`x_bit58 = 1 << 58`, ZERO carry into bit `59`!): `3.0` on both `f` and `g`.
- Exact orbit average: `u4`: `(2+3)/2 + 1.50 = 4.00` Ir/hit on **BOTH** `f` and `g` (Band-1 stops at `2^55.5` on both sides, eliminating `1bs` completely!).
```asm
; kernel: rev_b1b_eor_f_u4 items=4 charged=3.0 key=c_hit_1bl_k2_f
.Lb1bf_eor4:
ldp x0, x1, [x_p], #16
eor x0, x0, x_bit58
staddh w_lp, [x0]
eor x1, x1, x_bit58
staddh w_lp, [x1]
ldp x2, x3, [x_p], #16
eor x2, x2, x_bit58
staddh w_lp, [x2]
eor x3, x3, x_bit58
staddh w_lp, [x3]
subs x_cnt, x_cnt, #32
b.hs .Lb1bf_eor4
```
```asm
; kernel: rev_b1b_eor_f_bin21 items=1 charged=10.0 key=c_hit_1bl_k2_bin21_f
tbz x_cnt, #4, .Lb1bfe_b1
ldp x0, x1, [x_p], #16
eor x0, x0, x_bit58
staddh w_lp, [x0]
eor x1, x1, x_bit58
staddh w_lp, [x1]
.Lb1bfe_b1:
tbz x_cnt, #3, .Lb1bfe_done
ldr x0, [x_p], #8
eor x0, x0, x_bit58
staddh w_lp, [x0]
.Lb1bfe_done:
```
### 1.4 Symmetric Two-Pass Band-2 (`f`: `[2^55.5, 2^60.5)`, `g`: `[2^55.5, 2^63.5)`, `W_SB = 28`, `15.375 EiB` Slab Owners, `+4<<59` Zero-`orr/bfi` Pointer Emit, Exact 64-Bit Div-by-28 `M7 = ceil(2^66/7) = 0x924924924924924a`, 100%-Worst-Case `c1` `u4` `[55.5, 60.5)` + `u2` `[60.5, 61.5)` + `u1` `[61.5, 63.5)`)
- **Slab Layout (`14` slab owners/pair `x_tx_tbl[0..13]`, `L_max = 335.39`, active+pad `<= 362 < 14*28 = 392 < 448`):** `S_b1[2][2^57]` at `0x0000_0000_0000_0000` (`0.500 EiB`); `M_g[0..27][2][2^57]/8` at `0x0800_0000_0000_0000` (`0.875 EiB`); `S_b2[0..27][2][2^57]` at `0x2000_0000_0000_0000..0xFFFF_FFFF_FFFF_FFFF` (`14.000 EiB`, slot `col_in_sb in 0..27` at `(col_in_sb + 4) << 59`, total `15.375 EiB < 16.000 EiB`).
- **Exact `+4<<59` Zero-`orr/bfi` Emit (`3,136/3,136` verified live):** Records pack `(ka << 59) | (m_bit << 58) | (2*oid)` with `ka in 0..27`. With `x_t0 = (t0_c + 4) << 59`, `adds x4, x0, x_t0` leaves bits `0..58` intact: carry `C == 1` occurs iff `ka + t0_c >= 28` (since `28 + 4 = 32 = 2^5`!). On `c0` (`C == 0`), bits `59..63` equal `(ka + t0_c + 4) in [4, 31]` (exact `uint16_t*` into `S_b2[ka + t0_c]`, zero `orr`/`bfi`); on `c1` (`C == 1`), `add x4, x4, x_off4` (`4 << 59`) gives next slab's `S_b2[ka + t0_c - 28]`. XOR conjugates apply `eor x0, x0, x_bit58`.
- **Exact 64-Bit Unsigned Div-by-28 (`M7 = ceil(2^66 / 7) = 0x924924924924924a`, `7 * M7 = 2^66 + 6`):** With `y = x >> 2 = 7q + r < 2^62`, `0 <= 6q + 6*M7 < 0.911 * 2^66 < 2^66`, so `lsr #2; umulh M7; lsr #2` gives exact `floor(x / 28)` and `msub` gives exact `x mod 28 in [0, 27]` for all `x in [0, 2^64)` (`32 + 5 = 37.0 <= C_TILE_HDR`).
- **100%-Worst-Case `c1` Pricing (ZERO `P(c1)` assumption) & Complete Executable Binary Remainder Tails:** `[55.5, 60.5)` `u4->bin3` runs `rev_b2_u4_wrap_bin3` (`id`, `3` ops -> `.Lb2_u4_tail`) / `rev_b2_u4_wrap_eor` (`eor`, `3` ops -> `.Lb2_u4_tail_eor`) + 4-record `rev_b2_emit_id/eor_u4` (`id` `20/4 = 5.00`, `eor` `24/4 = 6.00` -> `<= 7.08 = C_HIT_HI_U4`, `+1` exit `b` if `k >= 4`) + executable `rev_b2_u4_tail_id` / `rev_b2_u4_tail_eor` (`.Lt2_*1_c1` falls through directly to `.Lb2_u4_t1*`: `bit 4` `<= 10` `id` / `<= 12` `eor` incl. `tbz`, `bit 3` `<= 7` `id` / `<= 8` `eor` incl. `tbz` & `b .Lnext_sb`); for EVERY `k >= 0` (incl. `k = 7` all-carry: `3 + 21 + 10 + 7 = 41 <= 42`), total wrap+emit+exit+tail is `<= 7.0 + 5.00*k` (`id`) / `<= 7.0 + 6.00*k` (`eor`), plus `9.0` header = `16.0 <= 16.0 = c_sb_hdr_u4` (verified for `k = 0..15` in `cm_models.py` and `k = 0..64` in `w8_verifiers.py`). `[60.5, 61.5)` `u2->bin1` (`rev_b2_u2_wrap_bin1`/`_eor` `3` ops + `rev_b2_emit_id/eor_u2` `11/13` ops + `+1` exit `b` + executable `rev_b2_u2_tail_id/eor` ending with `b .Lnext_sb`: `<= 5.5 + 5.5*k` / `<= 5.5 + 6.5*k`, `+ 9.0 <= 14.5 = c_sb_hdr_u2`, `<= 7.58 = C_HIT_HI_U2`); `[61.5, 63.5)` `u1` (ending with `b .Lnext_sb`, `+ 9.0 = 10.0 = c_sb_hdr`): `<= 9.08`, `10.0`.
- `.Lnext_sb` drains `>=8,192`-ptr bursts (`2x128 KiB` ping-pong rx bufs/slab, `<4 MiB<<0.625 EiB`) via `rev_b1a_apply_f_u32->u4->bin21` (`1.5625` `STADDH`, `14` slabs, `x_tx_tbl[-1..15]`). `rev_b2_flush_sb`/`_last` (`9`+`<=16` sync desc/ack `<=25/8192=0.0031<=0.01`) flush `>=64 KiB` bursts; at each `w_lp`-bucket end, an unconditional 14-slab epilogue drains any `<64 KiB` partial buf (`<=14*25=350` ops/bucket, `<2^-11` Ir/ptr) before weight change/scan. Padding writes in `[L_b, ceil(L_b/28)*28)` are excluded by active scan and zeroed by `rev_pad_zero_u64` (`0.15625`), charged via `cand_bd = 1 + 28/L_b` and `c_scan_zero = 1.50`. Wrapped chords (`k0+L_b>p`, `<2^-46.8`) run a 2nd probe (`<2^-41` Ir/row, verified in `cm_models.py` & `sigmaB_exactness`).
```asm
; kernel: rev_sb_probe_fast items=1 charged=32.0 key=c_tile_hdr_probe
mov x21, #0
cbz x_lb, .Lempty_chord
lsr x4, x_k0, #2
umulh x4, x4, x_m7
lsr x4, x4, #2
msub x0, x4, x_28, x_k0
sub x1, x_col0, x0
add x1, x1, #28
lsr x22, x1, #2
umulh x22, x22, x_m7
lsr x22, x22, #2
msub x0, x22, x_28, x1
sub x22, x22, #1
add x0, x0, #4
lsl x_t0, x0, #59
sub x5, x_p, x_k0
cmp x5, x_lb
csel x5, x_lb, x5, hi
add x6, x_k0, x5
sub x6, x6, #1
lsr x6, x6, #2
umulh x6, x6, x_m7
lsr x6, x6, #2
add x21, x6, #1
ldr x7, [x_off, x21, lsl #3]
sub x21, x21, x4
add x19, x_off, x4, lsl #3
ldr x2, [x19], #8
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
tbnz x10, #16, .Lb2_flush_last
.Lb2_after_last:
str x10, [x14]
.Lb2_sb_done:
.Lb2_chunk_slow:
```
```asm
; kernel: rev_b2_direct_sb_hdr items=1 charged=9.0 key=c_sb_hdr
.Lsb:
ldr x3, [x19], #8
cmp x2, x3
b.eq .Lnext_sb
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
; kernel: rev_b2_u4_wrap_bin3 items=1 charged=3.0 key=c_b2_u4_wrap
sub x_cnt, x3, x2
subs x_cnt, x_cnt, #32
b.lo .Lb2_u4_tail
```
```asm
; kernel: rev_b2_emit_id_u4 items=4 charged=5.00 key=c_b2_emit_id_u4
.Lb2_idu4:
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lb2_id4_0_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_1_c1
str x5, [x10], #8
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lb2_id4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_idu4
```
```asm
; kernel: rev_b2_u4_tail_id items=1 charged=26.0 key=c_b2_u4_tail_id
.Lb2_u4_tail:
tbz x_cnt, #4, .Lb2_u4_t1
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lt2_id0_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lt2_id1_c1
str x5, [x10], #8
b .Lb2_u4_t1
.Lt2_id0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lt2_id1_c1
str x5, [x10], #8
b .Lb2_u4_t1
.Lt2_id1_c1:
add x5, x5, x_off4
str x5, [x11], #8
.Lb2_u4_t1:
tbz x_cnt, #3, .Lnext_sb
ldr x0, [x2], #8
adds x4, x0, x_t0
b.hs .Lt1_id_c1
str x4, [x10], #8
b .Lnext_sb
.Lt1_id_c1:
add x4, x4, x_off4
str x4, [x11], #8
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_id_u4_c1a items=1 charged=15.0 key=c_b2_emit_id_u4_c1a
.Lb2_id4_0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_1_c1
str x5, [x10], #8
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lb2_id4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_idu4
b .Lb2_u4_tail
```
```asm
; kernel: rev_b2_emit_id_u4_c1b items=1 charged=12.0 key=c_b2_emit_id_u4_c1b
.Lb2_id4_1_c1:
add x5, x5, x_off4
str x5, [x11], #8
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lb2_id4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_idu4
b .Lb2_u4_tail
```
```asm
; kernel: rev_b2_emit_id_u4_c1c items=1 charged=8.0 key=c_b2_emit_id_u4_c1c
.Lb2_id4_2_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_id4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_idu4
b .Lb2_u4_tail
```
```asm
; kernel: rev_b2_emit_id_u4_c1d items=1 charged=5.0 key=c_b2_emit_id_u4_c1d
.Lb2_id4_3_c1:
add x5, x5, x_off4
str x5, [x11], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_idu4
b .Lb2_u4_tail
```
```asm
; kernel: rev_b2_u4_wrap_eor items=1 charged=3.0 key=c_b2_u4_wrap_eor
sub x_cnt, x3, x2
subs x_cnt, x_cnt, #32
b.lo .Lb2_u4_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u4 items=4 charged=6.00 key=c_b2_emit_eor_u4
.Lb2_eoru4:
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor4_0_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_1_c1
str x5, [x10], #8
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_eoru4
```
```asm
; kernel: rev_b2_u4_tail_eor items=1 charged=29.0 key=c_b2_u4_tail_eor
.Lb2_u4_tail_eor:
tbz x_cnt, #4, .Lb2_u4_t1_eor
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lt2_eor0_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lt2_eor1_c1
str x5, [x10], #8
b .Lb2_u4_t1_eor
.Lt2_eor0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lt2_eor1_c1
str x5, [x10], #8
b .Lb2_u4_t1_eor
.Lt2_eor1_c1:
add x5, x5, x_off4
str x5, [x11], #8
.Lb2_u4_t1_eor:
tbz x_cnt, #3, .Lnext_sb
ldr x0, [x2], #8
eor x0, x0, x_bit58
adds x4, x0, x_t0
b.hs .Lt1_eor_c1
str x4, [x10], #8
b .Lnext_sb
.Lt1_eor_c1:
add x4, x4, x_off4
str x4, [x11], #8
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_eor_u4_c1a items=1 charged=17.0 key=c_b2_emit_eor_u4_c1a
.Lb2_eor4_0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_1_c1
str x5, [x10], #8
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_eoru4
b .Lb2_u4_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u4_c1b items=1 charged=14.0 key=c_b2_emit_eor_u4_c1b
.Lb2_eor4_1_c1:
add x5, x5, x_off4
str x5, [x11], #8
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor4_2_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_eoru4
b .Lb2_u4_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u4_c1c items=1 charged=8.0 key=c_b2_emit_eor_u4_c1c
.Lb2_eor4_2_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_eor4_3_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_eoru4
b .Lb2_u4_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u4_c1d items=1 charged=5.0 key=c_b2_emit_eor_u4_c1d
.Lb2_eor4_3_c1:
add x5, x5, x_off4
str x5, [x11], #8
subs x_cnt, x_cnt, #32
b.hs .Lb2_eoru4
b .Lb2_u4_tail_eor
```
```asm
; kernel: rev_b2_u2_wrap_bin1 items=1 charged=3.0 key=c_b2_u2_wrap
sub x_cnt, x3, x2
subs x_cnt, x_cnt, #16
b.lo .Lb2_u2_tail
```
```asm
; kernel: rev_b2_emit_id_u2 items=2 charged=5.50 key=c_b2_emit_id_u2
.Lb2_idu2:
ldp x0, x1, [x2], #16
adds x4, x0, x_t0
b.hs .Lb2_id0_c1
str x4, [x10], #8
.Lb2_id_h1:
adds x5, x1, x_t0
b.hs .Lb2_id1_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_idu2
```
```asm
; kernel: rev_b2_u2_tail_id items=1 charged=9.0 key=c_b2_u2_tail_id
.Lb2_u2_tail:
tbz x_cnt, #3, .Lnext_sb
ldr x0, [x2], #8
adds x4, x0, x_t0
b.hs .Lt1_u2_id_c1
str x4, [x10], #8
b .Lnext_sb
.Lt1_u2_id_c1:
add x4, x4, x_off4
str x4, [x11], #8
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_id_u2_c1 items=1 charged=8.0 key=c_b2_emit_id_u2_c1
.Lb2_id0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_id1_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_idu2
b .Lb2_u2_tail
```
```asm
; kernel: rev_b2_emit_id_u2_c1b items=1 charged=5.0 key=c_b2_emit_id_u2_c1b
.Lb2_id1_c1:
add x5, x5, x_off4
str x5, [x11], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_idu2
b .Lb2_u2_tail
```
```asm
; kernel: rev_b2_u2_wrap_eor items=1 charged=3.0 key=c_b2_u2_wrap_eor
sub x_cnt, x3, x2
subs x_cnt, x_cnt, #16
b.lo .Lb2_u2_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u2 items=2 charged=6.50 key=c_b2_emit_eor_u2
.Lb2_eoru2:
ldp x0, x1, [x2], #16
eor x0, x0, x_bit58
eor x1, x1, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor0_c1
str x4, [x10], #8
adds x5, x1, x_t0
b.hs .Lb2_eor1_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_eoru2
```
```asm
; kernel: rev_b2_u2_tail_eor items=1 charged=10.0 key=c_b2_u2_tail_eor
.Lb2_u2_tail_eor:
tbz x_cnt, #3, .Lnext_sb
ldr x0, [x2], #8
eor x0, x0, x_bit58
adds x4, x0, x_t0
b.hs .Lt1_u2_eor_c1
str x4, [x10], #8
b .Lnext_sb
.Lt1_u2_eor_c1:
add x4, x4, x_off4
str x4, [x11], #8
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_eor_u2_c1 items=1 charged=8.0 key=c_b2_emit_eor_u2_c1
.Lb2_eor0_c1:
add x4, x4, x_off4
str x4, [x11], #8
adds x5, x1, x_t0
b.hs .Lb2_eor1_c1
str x5, [x10], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_eoru2
b .Lb2_u2_tail_eor
```
```asm
; kernel: rev_b2_emit_eor_u2_c1b items=1 charged=5.0 key=c_b2_emit_eor_u2_c1b
.Lb2_eor1_c1:
add x5, x5, x_off4
str x5, [x11], #8
subs x_cnt, x_cnt, #16
b.hs .Lb2_eoru2
b .Lb2_u2_tail_eor
```
```asm
; kernel: rev_b2_emit_id_c0 items=1 charged=7.0 key=c_b2_emit_id_u1
.Lb2_id1:
ldr x0, [x2], #8
adds x4, x0, x_t0
b.hs .Lb2_id1_u1c1
str x4, [x10], #8
cmp x2, x3
b.lo .Lb2_id1
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_id_c1 items=1 charged=5.0 key=c_b2_emit_id_u1c1
.Lb2_id1_u1c1:
add x4, x4, x_off4
str x4, [x11], #8
cmp x2, x3
b.lo .Lb2_id1
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_eor_c0 items=1 charged=8.0 key=c_b2_emit_eor_u1
.Lb2_eor1:
ldr x0, [x2], #8
eor x0, x0, x_bit58
adds x4, x0, x_t0
b.hs .Lb2_eor1_u1c1
str x4, [x10], #8
cmp x2, x3
b.lo .Lb2_eor1
b .Lnext_sb
```
```asm
; kernel: rev_b2_emit_eor_c1 items=1 charged=5.0 key=c_b2_emit_eor_u1c1
.Lb2_eor1_u1c1:
add x4, x4, x_off4
str x4, [x11], #8
cmp x2, x3
b.lo .Lb2_eor1
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
cmp x6, #112
b.ge .Lb2_fret
; xfer (x_rx,x_cnt) to slab owner x6>>3 (runs recv u32->u4->bin21 in own GPRs; sender x2 preserved)
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
cmp x6, #112
b.ge .Lb2_flret
; xfer (x_rx,x_cnt) to slab owner x6>>3 (runs recv u32->u4->bin21 in own GPRs; sender x2 preserved)
.Lb2_flret:
b .Lb2_after_last
```
### 1.5 Two-Pass `M_g` Scan+Zero (`rev_scan_g_u64` `1.0625` + `rev_scan_f_u128` `0.2891` + `rev_pad_zero_u64` `0.15625` + `2`-op Wrap <= `1.50`)
- **Closed Periodic `S_b1` Handoff (ZERO Extra Ops):** At startup (`<2^-55` Ir/cell), `S_b1 = bias_g[c] = 0x8000 - T_g[c]` (`T_g=(302,290,278)`) and `S_b2[0..27]=0`. For active rows `r=0..r_last`:
  - **Pass 1 (`rev_scan_g_u64`, `68` ops / `64` cells = `1.0625` Ir/cell):** Enters every active `g` row `r` with `S_b1 == bias_g[c]`. Loads `16` `uint64_t` words from `S_b1` and `S_b2[r]`, sums lane-wise into `x0..x15` (zero 16-bit carry since `bias_g + S_g <= 35,767 < 65,536`), writes `x_b16` to `S_b1` (`8x stp`) and `0` to `S_b2[r]` (`8x stp`), runs the 15-op `orr` tree (preserving `x0..x15` for `.Lhit64_g`), tests `x_8000`, and writes `M_g[r]` (`0` on `> 99.8%` fast path). For `r < r_last`, the outer row loop passes `x_b16 = bias_g[c]` (re-arming `S_b1` for `g` row `r + 1`); **on the last active `g` row `r == r_last`, the outer row loop passes `x_b16 = 0` (`xzr`)**, leaving `S_b1 == 0` ready for Pass 2 (`f`) row `r = 0` with ZERO extra instructions! Any padding rows `r = r_last + 1..27` run `rev_pad_zero_u64` on `S_b2[r]` only (`S_b1` stays `0`).
  - **Pass 2 (`rev_scan_f_u128`, `37` ops / `128` cells = `0.2891` Ir/cell):** Enters every active `f` row `r` with `S_b1 == 0` (unbiased `f` accumulator). Loads `128` bits of `M_g[r]` (`ldp x0, x1, [x_mg], #16; orr x0, x0, x1; cbnz x0, .Lhit128_f`). On the `> 99.7%` all-zero fast path (`128 * sga_mean < 0.003`), skips ALL `S_b1`/`S_b2` loads and stores `x_init_f` across `256 B` of `S_b1` (`16x stp x_init_f, x_init_f`) and `0` across `256 B` of `S_b2[r]` (`16x stp xzr, xzr`). On `.Lhit128_f`, only set `M_g` lanes test `S_b1 + S_b2[r] >= T_f` (`122`), emit, then likewise store `x_init_f` into `S_b1` and `0` into `S_b2[r]`. For `r < r_last`, the outer row loop passes `x_init_f = 0` (leaving `S_b1 == 0` for `f` row `r + 1`); **on the last active `f` row `r == r_last`, the outer row loop passes `x_init_f = bias_g[c]`**, re-arming `S_b1 == bias_g[c]` for the NEXT superbucket's Pass 1 (`g`) row `r = 0` with ZERO extra instructions! Any padding rows `r = r_last + 1..27` run `rev_pad_zero_u64` on `S_b2[r]` only (`S_b1` stays `bias_g[c]`). Verified end-to-end across 3 consecutive full & partial superbuckets on ONE single shared `S_b1` array in `cm_models.py`.
  - **Padding-Only Zeroing (`rev_pad_zero_u64`, `10` ops / `64` cells = `0.15625` Ir/cell):** Zeroes ONLY `S_b2[r]` on padding rows (`2 * (28 / L_b) * (10 / 64) = 0.0604` Ir/cell). Total `1.0625 + 0.2891 + 0.0604 + 4/144.89 + (1000 + 1200) * 2.311e-05 = 1.4904 <= 1.50 = c_scan_zero`!
```asm
; kernel: rev_scan_zero_wrap items=1 charged=2.0 key=c_scan_zero_wrap
cbz x_cnt, .Lscan_done
.Lscan_done:
ret
```
```asm
; kernel: rev_scan_g_u64 items=64 charged=1.065 key=c_scan_g_raw
.Lscan64_g:
ldp x0, x1, [x_b1]
ldp x16, x17, [x_b2]
add x0, x0, x16
add x1, x1, x17
ldp x2, x3, [x_b1, #16]
ldp x16, x17, [x_b2, #16]
add x2, x2, x16
add x3, x3, x17
ldp x4, x5, [x_b1, #32]
ldp x16, x17, [x_b2, #32]
add x4, x4, x16
add x5, x5, x17
ldp x6, x7, [x_b1, #48]
ldp x16, x17, [x_b2, #48]
add x6, x6, x16
add x7, x7, x17
ldp x8, x9, [x_b1, #64]
ldp x16, x17, [x_b2, #64]
add x8, x8, x16
add x9, x9, x17
ldp x10, x11, [x_b1, #80]
ldp x16, x17, [x_b2, #80]
add x10, x10, x16
add x11, x11, x17
ldp x12, x13, [x_b1, #96]
ldp x16, x17, [x_b2, #96]
add x12, x12, x16
add x13, x13, x17
ldp x14, x15, [x_b1, #112]
ldp x16, x17, [x_b2, #112]
add x14, x14, x16
add x15, x15, x17
stp x_b16, x_b16, [x_b1, #16]
stp x_b16, x_b16, [x_b1, #32]
stp x_b16, x_b16, [x_b1, #48]
stp x_b16, x_b16, [x_b1, #64]
stp x_b16, x_b16, [x_b1, #80]
stp x_b16, x_b16, [x_b1, #96]
stp x_b16, x_b16, [x_b1, #112]
stp x_b16, x_b16, [x_b1], #128
stp xzr, xzr, [x_b2, #16]
stp xzr, xzr, [x_b2, #32]
stp xzr, xzr, [x_b2, #48]
stp xzr, xzr, [x_b2, #64]
stp xzr, xzr, [x_b2, #80]
stp xzr, xzr, [x_b2, #96]
stp xzr, xzr, [x_b2, #112]
stp xzr, xzr, [x_b2], #128
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
b.ne .Lhit64_g
str xzr, [x_mg], #8
subs x_cnt, x_cnt, #64
b.ne .Lscan64_g
.Lhit64_g:
```
```asm
; kernel: rev_scan_f_u128 items=128 charged=0.290 key=c_scan_f_raw
.Lscan128_f:
ldp x0, x1, [x_mg], #16
orr x0, x0, x1
cbnz x0, .Lhit128_f
stp x_init_f, x_init_f, [x_b1, #16]
stp x_init_f, x_init_f, [x_b1, #32]
stp x_init_f, x_init_f, [x_b1, #48]
stp x_init_f, x_init_f, [x_b1, #64]
stp x_init_f, x_init_f, [x_b1, #80]
stp x_init_f, x_init_f, [x_b1, #96]
stp x_init_f, x_init_f, [x_b1, #112]
stp x_init_f, x_init_f, [x_b1, #128]
stp x_init_f, x_init_f, [x_b1, #144]
stp x_init_f, x_init_f, [x_b1, #160]
stp x_init_f, x_init_f, [x_b1, #176]
stp x_init_f, x_init_f, [x_b1, #192]
stp x_init_f, x_init_f, [x_b1, #208]
stp x_init_f, x_init_f, [x_b1, #224]
stp x_init_f, x_init_f, [x_b1, #240]
stp x_init_f, x_init_f, [x_b1], #256
stp xzr, xzr, [x_b2, #16]
stp xzr, xzr, [x_b2, #32]
stp xzr, xzr, [x_b2, #48]
stp xzr, xzr, [x_b2, #64]
stp xzr, xzr, [x_b2, #80]
stp xzr, xzr, [x_b2, #96]
stp xzr, xzr, [x_b2, #112]
stp xzr, xzr, [x_b2, #128]
stp xzr, xzr, [x_b2, #144]
stp xzr, xzr, [x_b2, #160]
stp xzr, xzr, [x_b2, #176]
stp xzr, xzr, [x_b2, #192]
stp xzr, xzr, [x_b2, #208]
stp xzr, xzr, [x_b2, #224]
stp xzr, xzr, [x_b2, #240]
stp xzr, xzr, [x_b2], #256
subs x_cnt, x_cnt, #128
b.ne .Lscan128_f
.Lhit128_f:
```
```asm
; kernel: rev_pad_zero_u64 items=64 charged=0.160 key=c_pad_zero_raw
.Lpad64:
stp xzr, xzr, [x_b2, #16]
stp xzr, xzr, [x_b2, #32]
stp xzr, xzr, [x_b2, #48]
stp xzr, xzr, [x_b2, #64]
stp xzr, xzr, [x_b2, #80]
stp xzr, xzr, [x_b2, #96]
stp xzr, xzr, [x_b2, #112]
stp xzr, xzr, [x_b2], #128
subs x_cnt, x_cnt, #64
b.ne .Lpad64
```
### 1.7 Explicit Pure-Scalar 2-Bank-Pair Counting-Sort / Prefix-Sum / Scatter Build (`C_LIST_A = 16.0` Band-1, `C_LIST_A_B2 = 19.5` Band-2, `C_LBKT = 4.0`, `build_u <= 105.4240`)
Pass 1 for Band-2 (`rev_b2_build_hist_u4`: `48/4 = 12.0` Ir/entry) executes `4x` overflow-proof modular steps (`p = 2^64 - 59`), computes exact 64-bit `s = floor(kv / 28)` (`lsr #2; umulh M7; lsr #2` with `M7 = ceil(2^66/7) = 0x924924924924924a`, exact for all `kv in [0, 2^64)`) and `ka = kv - 28*s` (`msub`), packs `(ka << 59) | s` directly into `kbuf`, and accumulates byte counts (`+8`) into `hist[s]`; `rev_build_prefix_u4` (`12/4 = 3.0 <= 4.0 = C_LBKT` Ir/bucket, `+1.0` extra in `audit_cm.py`) computes `pos[s] = flat_base + sum_{t<s} hist[t]` while zeroing `hist[]` in-place via `2x stp xzr, xzr`; Pass 2 (`rev_b2_build_scatter_u4`: `28/4 = 7.0` Ir/entry, ZERO `umulh`/`msub`) extracts `s = w & MASK59` (`x_m59 = (1<<59)-1`), inserts `x_om = (m_bit << 58) | (2*oid)` into bits `0..58` via 1-op pure ARMv8.0-A `bfi x0, x_om, #0, #59` (preserving `(ka << 59)` in bits `59..63`), increments `x_om += 2`, and scatters `8 B` records. Total Band-2 `12.0 + 7.0 + <=0.5` chord resets & tail = `19.5 <= C_LIST_A_B2 = 19.5` (Band-1 omits div-by-28 & packing, keeping `<= 16.0 = C_LIST_A`!). Peak per-owner build RAM `<= 10.204 EiB < 16.00 EiB` across ALL primes:
```asm
; kernel: rev_b2_build_hist_u4 items=4 charged=12.0 key=c_bld_hist
.Lbld_h4:
subs x4, x_k, x_pms
add x5, x_k, x_d
csel x0, x4, x5, hs
lsr x4, x0, #2
umulh x4, x4, x_m7
lsr x4, x4, #2
msub x5, x4, x_28, x0
orr x6, x4, x5, lsl #59
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
subs x4, x0, x_pms
add x5, x0, x_d
csel x1, x4, x5, hs
lsr x4, x1, #2
umulh x4, x4, x_m7
lsr x4, x4, #2
msub x5, x4, x_28, x1
orr x7, x4, x5, lsl #59
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
stp x6, x7, [x_kbuf], #16
subs x4, x1, x_pms
add x5, x1, x_d
csel x2, x4, x5, hs
lsr x4, x2, #2
umulh x4, x4, x_m7
lsr x4, x4, #2
msub x5, x4, x_28, x2
orr x6, x4, x5, lsl #59
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
subs x4, x2, x_pms
add x5, x2, x_d
csel x_k, x4, x5, hs
lsr x4, x_k, #2
umulh x4, x4, x_m7
lsr x4, x4, #2
msub x5, x4, x_28, x_k
orr x7, x4, x5, lsl #59
ldr x5, [x_hist, x4, lsl #3]
add x5, x5, #8
str x5, [x_hist, x4, lsl #3]
stp x6, x7, [x_kbuf], #16
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
and x4, x0, x_m59
bfi x0, x_om, #0, #59
add x_om, x_om, #2
ldr x6, [x_pos, x4, lsl #3]
str x0, [x6], #8
str x6, [x_pos, x4, lsl #3]
and x4, x1, x_m59
bfi x1, x_om, #0, #59
add x_om, x_om, #2
ldr x6, [x_pos, x4, lsl #3]
str x1, [x6], #8
str x6, [x_pos, x4, lsl #3]
and x4, x2, x_m59
bfi x2, x_om, #0, #59
add x_om, x_om, #2
ldr x6, [x_pos, x4, lsl #3]
str x2, [x6], #8
str x6, [x_pos, x4, lsl #3]
and x4, x3, x_m59
bfi x3, x_om, #0, #59
add x_om, x_om, #2
ldr x6, [x_pos, x4, lsl #3]
str x3, [x6], #8
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
  Schoolbook `L_1 x L_2` limb mul uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)`. On `L_2(a) x L_1(b)` (`Q_CM(a) <= s_c * R2_a <= 2 * 318,861 = 637,722`, `sum_{k=1}^4 |sigma_k(b)| <= R_1 = 670.79`, `Q_CM(b) <= 2 R_1^2 = 899,911`), Cauchy-Schwarz (`max diag(G_8^{-1}) = 2/5`) and Vandermonde inverse (`C_lin = max_{j,p} 2|V^{-1}_{j,p}| = 0.6529`) prove `|a_j| <= 506`, `|b_j| <= C_lin*R_1 <= 438`, `|a^3_j| <= 2^{26.81}`, `|b^3_j| <= 2^{27.55}` (1 signed limb `<= 2^{63}`; `4` `8x8` 1-limb poly muls + `Phi_15` mod = `400` Ir), `|g(a,b)_j| <= 2^{91.27}` (2 signed limbs `<= 2^{127}`; `+96` Ir), and `|f(a,b)_j| <= 2^{85.01}` (2 limbs, priced as 3: `|a^6_j|,|b^6_j| <= 2^{55.72}` 1 limb priced as 2 = `1,856` Ir; `|a^9_j|,|b^9_j| <= 2^{83.89}` 2 limbs priced as 3 `a^9,b^9,a^6b^3,a^3b^6` = `4,608` Ir; combine = `128` Ir). Intermediate tower norms fit comfortably (`|N_1(g)|<=2^{178}<<2^{255}`, `|N_2(g)|<=2^{356}<<2^{511}`, `|N_Q(g)|<=2^{688}<<2^{1023}`; `|N_1(f)|<=2^{172}<<2^{383}`, `|N_2(f)|<=2^{343}<<2^{767}`, `|N_Q(f)|<=2^{685}<<2^{1535}`):
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
  Three `< 2^62` primes divisible by `3 * 2^55`: `p_1 = 69 * 2^55 + 1`, `p_2 = 57 * 2^55 + 1`, `p_3 = 54 * 2^55 + 1` (`4*p_1 < 2^{63.11} < 2^64`, `prod = 2^182.696 > 2^54 * (2^64 - 1)^2`). Pairing 2 butterflies via `4x ldp`/`4x stp` and computing `x12 = (x1 + 2p) - x3 in [0, 4p) < 2^{63.11}` in 2 ops (Lazy Shoup output `< (1 + 4p/2^64)p < 1.54p < 2p`, verified LIVE across all 3 primes and `n = 64` convolution!) cuts `ntt3_F(M)` to `min(2*m^2 + 4*m if m<=256 else inf, 48.50*n2*k + (26*n2 if k odd) + 100*n2, 48.50*n3*(k-2) + (26*n3 if k-2 odd) + 244*n3)` (`+ 512*n` if `n > 2^26`; schoolbook `mul_LL(m/2, m/2) = 8*(m/2)^2 + 8*(m/2) = 2*m^2 + 4*m`). Stage-B `c_leaf` uses 64-bit Lehmer GCD (`lehmer32_step` `29/2 = 14.5 <= 16.0` ops/step + `L`-limb linear combination `<= 64*L` -> `(nbits/24)*(64*L + 544)`, `>= 1.48x` over worst-case Fibonacci & `2.2x` over random 680b pairs, verified LIVE!), dropping **`c_cof_2d_upper` to `0.7760` Ir/cell**:
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
- **Stage B (`P_min = 32`, `KH = 0.08`, Static 3-Class `T_g[c] = (302, 290, 278)` + Constant `T_f_const = 122`b Pre-Gates, `step = 0.25`):** `gbatch_ntt3(pbits_split(32, 0.08), leaf=L, c_leaf=10*MONT5*(L/5)^2 + (nb/24)*(64*L+544))` (`MONT5 = 280`, `ntt62_r4_u2`, schoolbook `M <= 256`, `lehmer32_step` 64b Lehmer GCD) -> exact **`c_cof_2d_upper = 0.7760` Ir/cell** (Chebotarev `5/72 = 0.0694 < 0.08`, live `2^14` `0.07197/0.05649 <= 0.08`).
- **Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree `H-ECM` Ladder (`q15_ecm16.py`, `cm_models.ecm_mean_var`):** Exact Newton perfect-power precheck (`2^26` Ir > 32.06M ops for 2,590b boot) removes all `N = p^k`. Stage 1 steps per prime factor `q` across all `q^k <= B1` (`<= pi_pp_upper(B1)` gcds); Stage 2 builds a `J`-leaf binary product tree (`h_tree = ceil(log2(max(2, J)))`) replacing any `t_j == 0 mod N` leaf with `1` before tree build -> deterministic `(ops + 1.1*J*G)*mul + (pi_pp_upper(B1) + (1 + 2*h_tree)*G)*gcd` (ZERO curve replay, `var_t = 0.0`, verified live on 160 composite proper splits + Stage-2 tree-bisect & zero-leaf-skip in `q15_ecm16.verify_suyama_stepped`). `tree_branch_caps(nbits, b_min, l_min=2)` uses exact `s_i = b_i - 1` DP (`mont2_cios` `47 <= 48` ops). Collection (`lo_r = 48, b_min = 96, success = 0.999999`, `n_rels_upper = 2^(log2_rels+0.55)`) has `var_t = 0.0` with `1 - 29e-6` yield factor in `det_cred_pen`; descent (`lo_bits = 16, b_min = 32, success = 0.999`) is 100% deterministic.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Dense + 16-Shard `< 2^64 B` 4-Vec 4-Col Arch B + `u32->u4->bin21` Streaming + Single-Sub `57`-Op Barrett `row_red5_z12`)
Over `Z/rZ` (`r` 255b -> 4 limbs, `phi(12) = 4`, `tab_ratio = 1.00` ZERO-`TABM`, `32` 4-vec bundles `128 B/slot`, `n = 128, m = 32,768`, 1-owner dense RAM `< 16 EiB`, 16 sparse shards `<< 16 EiB`, explicit cross-owner streaming `c_stream_mat = 63.0 / 251.0`):
- **`unit_run_4vec4_z12` & `_neg` (`122/16 = 7.625 <= 7.65` Ir/vec-entry, `28` GPRs, ZERO `x18/x29/x30`, `0` spills):** Pre-scaled byte offset `col * 128` added to `x6 = x_base` (ping-pong native, 0 copies); 4 columns x 4 vectors into 4 signed 5-limb accs (`x8..x28`), tails priced inside `c_row_chi = 600.0`:
```asm
; kernel: unit_run_4vec4_wrap items=1 charged=4.0 key=unit_run_4vec4_wrap
subs x5, x5, #4
b.lo .Lrun4_tail
.Lrun4_tail:
adds x5, x5, #4
b.eq .Lrun4_done
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
- **Sign-Split 2-Vec `gen_run1_2v_z12` (`c_b > 0`) & `gen_run1_2v_neg_z12` (`c_b < 0`) (`43/2 = 21.50 <= 22.0` Ir/vec-entry, `24 <= 29` GPRs avoiding `x18/x29/x30`, `0` spills; caller-side non-empty `cbz x24` precondition guarantees `x24 >= 1` on entry):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`, provisioned up to `nonunit_ub = 0.113` covering every `2^18..2^25` ladder point and `1 - (eu_25 - 2*se) = 0.08217`, stored as 1-byte `uint8_t` magnitudes loaded and zero-extended to `x22` via `ldrb w22, [x23], #1`) into positive and negative sub-lists, executed as two spill-free 2-vector passes (`(v0, v1)` at offsets `#0..#48`, then `(v2, v3)` at offsets `#64..#112`, explicit `add x8, x9, x8` from `x_base = x9`): `3` ptr/coeff loads + `2 * 19` mul/umulh/carry ops + `2` loop = **`43 / 2 = 21.50 <= 22.0` Ir/vec-entry** (`k5_live = eu_lo * 7.65 + (1 - eu_lo) * 22.0 = 8.8291 <= k5_charged = 9.45` Ir, `+0.62` Ir headroom; even at `nonunit_ub = 0.113`, `0.887 * 7.65 + 0.113 * 22.0 = 9.2716 <= 9.45`!):
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
- **Explicit 1-Vector Horner Kernels (`unit_run_1vec4_z12` / `_neg`: `38/4 = 9.50 <= 9.60` Ir/entry; `gen_run1_1v_z12` / `_neg`: `24.0 <= 24.0` Ir/entry with caller-side non-empty `cbz x24` precondition `x24 >= 1` -> `k5_horner_charged = 11.25`):** Reuses the EXACT SAME read-only sparse-matrix shards (`128*j` stored byte offsets) on the `32 B/slot` 1-vector Horner dense buffer via **`add x7, x6, x7, lsr #2`** (`(128*j) >> 2 = 32*j` exact byte offset, ZERO matrix re-encoding, ZERO offset-shift bug!), simulated bit-exact in `cm_models.py`:
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
- **`chi_mul_z12_fast` (`50` Ir), Non-Double-Counted `row_tail_6bucket_z12` (`240` Ir = `12` unit sub-lists x (`4/4` wrap + `<=3*0.35` tail delta above `7.65`) + `12*0.5` non-unit header + `2*50` `chi_mul` + `56` header/zero/store + `53` slack), 1-Dense + 16-Projection + 16-Shard 4-Vec Arch B Row Overhead (`c_row_chi = 600.0` Ir, `c_stream_mat = 63.0` Ir/vec-row Krylov, `c_stream_mat_1vec = 251.0` Ir/row 1-vec Horner via `u32->u4->bin21` `1.0625*212.714 + 24.0 = 250.01 <= 251.0`, `nonunit_ub = 0.113`, `k5_charged = 9.45`, `k5_horner_charged = 11.25`, `c_xproj = 65.0` Ir/vec-row `>= 64.5`, `n_owners = 1`, `n_proj_owners = 16`, `n_mat_owners = 16`, `te = 0`):** `6 * 57 + 240 + 16 = 598 <= 600.0 = c_row_chi` Ir worst-case!

## 5. Descent Kernels (Generous Scalar Bounds for Uncapped `dsieve16`, Exact-Integer `BOUND_CELL`, & Input-Uniform Bareiss LLL)
- **Uncapped `dsieve16` on `S_16(R) = B_16(R)/{+-1}` (`SIEVE_HIT = 1,400` Ir/hit, ZERO skip caps):** Every degree-1 ideal `p in [2^8, 2^69.0]` runs Bareiss LLL (`PAIR_SETUP = 2^35.2`) and full 16D Fincke-Pohst Q32 DFS + leaf materialization (`64`) + ambient sign canonicalization (`40`) + bijective 16D half-ball ranking via `129.03 MiB` prefix table (`320`) + routing/update/spill (`260`) = `<= 914–1,244 <= 1,400` Ir/hit (4D DFS `18,174` pts, sign-flip witness, 55b `SQ_LUT`, and 16D `R = 2` `552/552` bijection verified live in `cm_models.py`).
- **Certified Outward-Rounded `BOUND_CELL (+6)` (`4,200` Ir/cell, `q15_bound_cell.py`):** Signed-64b `bmul_6` (`+6 >= ceil(sqrt(5)+3)`, `<= +0.5343b < +0.625b` `SQ_LUT` envelope): `3,792` exact ball-mul core + `128 + 80 + 64 + 136 = 4,200` Ir/cell (50-digit `Decimal` witness & 12-cell real-lattice `bc_gmin >= 0` verified live).
- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every `L_Q` (`det = q <= 2^260`) and `G_T` boot lattice (`S_boot = 2^21` scaling) runs `11/40` Euclidean -> `99/100` `Q_CM = diag(G_8, G_8)` (`det(G_8) = 3^4 * 5^6`, `lambda_max(G_8) = 15` proved live via Sylvester minors): **`PAIR_SETUP = 2^35.2`** (`2^35.06`), **`LLL_NODE = 2^40.0`** (`2^39.46`), **`LLL_BOOT = 2^47.0`** (`2^46.75`), alongside `SCAN = 10`, `ZERO = 2`, `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`.
