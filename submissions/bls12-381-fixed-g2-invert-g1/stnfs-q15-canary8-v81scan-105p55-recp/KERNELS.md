# Kernels (Raw-Opcode-Verified Pure-Scalar ARMv8.1-A with Disclosed `FEAT_LSE STADDB`/`STADD`: `Arch4-K2-Lq075-ColdNP2D-NoDir1-Rot555-2D56bNTT3-M15-W1W8-Bit58-U8Fused-Shoup61-Hex6`)

All hot-loop kernels use pure-scalar general-purpose-register ARMv8.1-A integer instructions (`-march=armv8.1-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `DC ZVA`, zero `fp64`, zero `x18`, zero `x29`, zero `bl`/`blr` — all kernels are inline macros or leaf routines returning via `ret`), plus **two disclosed scalar `FEAT_LSE` atomic integer instructions** from ARMv8.1-A, each priced at `1` scalar instruction = `1 / 20,000` units under `COST_MODEL.md`:
- **`STADDB Ws, [Xn]`** (32-bit byte alias of `LDADDB Ws, WZR, [Xn]`; scalar registers; base-register `[Xn]` addressing only), used for `uint8_t` score updates on fused `S_b1` (`0.25 EiB`) and `S_b2[0..31]` (`8.00 EiB` per slab owner across `M = 15` collision-free chord slab owners) with zero wraparound (deterministically proved `<= 248 <= 255` on `g` at `4.00`b `round` with per-column odd-norm `B_{g,odd,max}(b, c)`-adaptive `T_g_lo(b, c) in [30, 115]` — realized `<= 245` at `sl_g_tgt = 428.0` — and `<= 190 <= 240 <= 255` on `f` at `4.00`b `round` from the Vandermonde limb ceiling `log2 |N_f| <= 720`, verified live in `audit_cm.py`).
- **`STADD Xs, [Xn]`** (64-bit alias of `LDADD Xs, XZR, [Xn]`; `[Xn]` addressing only), used exclusively for in-place `uint64_t pos[1..B]` byte-count accumulation during dense table build across ALL active primes `2^{33} <= p < 2^{61.3}` (`8 * N_{a,k} < 2^{61} << 2^{64}` -> overflow mathematically impossible!). All other `FEAT_LSE` opcodes (`LDADD*`, `STADDH`, `SWP*`, `CAS*`) remain forbidden.

## 1. Collection Kernels (`K = 2` Equal-`V_a` Radial Shells `bp = [1, 1]` + Complex-Place `L_{3/4}(b)` Ball + `M = 15` Full-Chord Slab Owners `480` Cols + Split-`p_rot` / 2-Tier Band-2 `W_B=1`/`W_B=8` + Assembled Cold-Path 2D NP 2x/4x SWAR Scan)

Each shell `k in {0, 1}` owns `V_a,k = 0.7431 * 2^56 < 0.95 * 2^56` (`N_a,k = 4 * V_a,k = 2^57.572`, `stride = 56`, `bp = 1`). Under `K_{3/4}(R_{q,k})` (`R_{q,0} = 982.55`, `L_max = 389.92`, `|g_j| < 2^92.92`, `|f_j| < 2^89.10`), `M = 15` slab owners (`32` cols each, `480 >= 390 + 7`) hold an entire chord with zero cross-column collisions (`9.250 EiB` main arrays + `<= 0.0374 EiB` model max-bin `G_val` `<= 9.287 EiB`, or `<= 13.250 EiB < 16.000 EiB` with fixed `4.000 EiB` `G_val`; Direct-1 Tier 2c inactive at `p_g = [61.3, 61.3] == p_dir`). Band-2 `g`-primes emit via `tx_tbl[-1..15]` before each slab applies Band-1a `g`-primes, runs `rev_scan_g_u8_fused_256` + `.Lhit128_g`, and clears `<= 7` overhang cols (`rev_pad_zero_u8_128`); Pass 2 (`f`, zero overhang) runs `rev_scan_f_u8_512` + `.Lhit128_f` testing `S_f + S_g >= T_g_lo(b, c) + 55` before refilling `S_b1` and clearing `S_b2`. Sequential shells keep peak tables `2^{116.662} B < 2^{120} B` and peak 1-owner build `11.028 EiB < 16.000 EiB`.

### 1.1a `rev_b2_row_k0_step4` (Band-2 4-prime orbit incremental `D_(omega, sigma)` row-step `k_0` init: `33/4 = 8.25` Ir/prime/row + `<= 0.691` outer-chord reset <= `C_B2_ROW_K0 = 8.95` Ir/prime/row)
Because `sigma(omega) = omega`, all `4` primes in a Band-2 orbit share `p`, `pwk[0] = 1`, and `pwk[1] = w = s_0^5 mod p`. Every `(b_1..b_7)` fiber sieves ONE contiguous bounding span `[b_{0,min}, b_{0,max}]` of length `<= 390 < 449` (zero multi-component resets; interior non-convex gaps `< 0.1%` cleared by `rev_pad_zero_u8_128`, covered by `H2_ALLOW = 1.005`). Maintaining `k_tail[0..3]` (`32 B`) computes overflow-proof `0 <= x10..x13 < x_p` in `33/4 = 8.25` Ir/prime/row (disjoint GPRs `x19..x23`, verified at `p = 2^64 - 59`) plus `<= 64/92.65 = 0.691` outer-chord reset -> `8.941 <= 8.95 = C_B2_ROW_K0`:
```asm
; kernel: rev_b2_row_k0_step4 items=4 charged=8.95 key=c_b2_row_k0
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
### 1.1 Band-1 Shared `k_t` Step & Pre-Based Bucket Header (`b1a_long_lk = 31.5`)
Each Band-1a prime `[2^33, 2^52.5)` owns its own dedicated pre-based flat table with `(1 << 58) | (((m - k_c) & 3) << 56) | oid` pre-baked across all 4 Galois conjugates `k_c in {0..3}` (ZERO runtime rotation ops in `rev_b1a_apply_f_u64`, ONE contiguous slice per residue `k_t`, exact on all 4/4 conjugates!). Pre-adding `flat_buf_base` into `off_ptr[k] = flat_buf_base + 8 * prefix_sum[k]` at table-build time (ZERO extra RAM, ZERO extra build ops) eliminates per-bank base adjustments and `lsr #3` shifts (all byte counts `x_cnt = x_end - x_p = 8 * k` remain exact multiples of 8).
```asm
; kernel: rev_b1_kt_step_shared items=1 charged=3.0 key=c_b1_kt_step
subs x1, x0, x_pms
add x2, x0, x_st
csel x0, x1, x2, hs
```
```asm
; kernel: rev_b1_bkt_hdr_cnt items=1 charged=4.0 key=c_b1_bkt_hdr_cnt
add x1, x_off, x0, lsl #3
ldp x_p, x_end, [x1]
subs x_cnt, x_end, x_p
b.eq .Lb1_empty
```
### 1.2 Shared `uint8_t S_b1` Band-1a (`[2^33, 2^52.5)`): `64x` + Straight-Line `32-16-8-4-2-1` Binary Tail `STADDB` (`2.5625 <= 2.57` Ir/hit, Exact Worst-Case `b1a_long_lk = 31.44 <= 31.5`)

Each Band-1a prime owns its own flat table storing the exact pre-rotated 64-bit `uint8_t*` pointer `(1 << 58) | (((m - k_c) & 3) << 56) | oid` into `S_b1` (`1` slice per residue `k_t`). `rev_b1_stream_u64` (`66/64 = 1.03125`) + `rev_b1a_apply_f_u64` (`98/64 = 1.53125` `STADDB`; Band-2 TX-buffer flushes use the register-disjoint `rev_b2_rx_drain_*` blocks in §1.4 with the exact same `98/64 = 1.53125` rate) = **`2.5625 <= 2.57` Ir/hit**. After the `64x` bulk loop, straight-line binary tail dispatch (`rev_b1_stream_bin32` + `rev_b1a_apply_f_bin32` testing bits `#8..#3` of `x_cnt = 8 * r`, `0 <= r < 64`) executes every remainder `r in 0..63` with exact deterministic worst-case excess `mx1al = 21.4375` (at `kl = 1`, verified exhaustively for all `kl in 0..255` in `w8_verifiers.py`): `6.0 + 1.0 * (4.0 + 21.4375) = 31.4375 <= 31.5 = b1a_long_lk`.

```asm
; kernel: rev_b1_stream_pre64 items=1 charged=2.0 key=c_b1_stream_pre64
subs x_cnt, x_cnt, #512
b.lo .Lb1s_fix64
```
```asm
; kernel: rev_b1_stream_u64 items=64 charged=1.035 key=c_b1_stream
.Lb1s64:
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
subs x_cnt, x_cnt, #512
b.hs .Lb1s64
```
```asm
; kernel: rev_b1_stream_bin32 items=1 charged=72.0 key=c_b1_stream_bin32
.Lb1s_fix64:
adds x_cnt, x_cnt, #512
b.eq .Lb1s_done
tbz x_cnt, #8, .Lb1sb16
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
.Lb1sb16:
tbz x_cnt, #7, .Lb1sb8
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
.Lb1sb8:
tbz x_cnt, #6, .Lb1sb4
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
ldp x2, x3, [x_src], #16
stp x2, x3, [x_dst], #16
ldp x4, x5, [x_src], #16
stp x4, x5, [x_dst], #16
ldp x6, x7, [x_src], #16
stp x6, x7, [x_dst], #16
.Lb1sb4:
tbz x_cnt, #5, .Lb1sb2
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
ldp x2, x3, [x_src], #16
stp x2, x3, [x_dst], #16
.Lb1sb2:
tbz x_cnt, #4, .Lb1sb1
ldp x0, x1, [x_src], #16
stp x0, x1, [x_dst], #16
.Lb1sb1:
tbz x_cnt, #3, .Lb1s_done
ldr x0, [x_src], #8
str x0, [x_dst], #8
.Lb1_empty:
.Lb1s_done:
```
```asm
; kernel: rev_b1a_apply_pre64 items=1 charged=2.0 key=c_hit_lo_m4_pre64
subs x_cnt, x_cnt, #512
b.lo .Lb1a_fix64
```
```asm
; kernel: rev_b1a_apply_f_u64 items=64 charged=1.535 key=c_hit_lo_m4_f
.Lb1af64:
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
subs x_cnt, x_cnt, #512
b.hs .Lb1af64
```
```asm
; kernel: rev_b1a_apply_f_bin32 items=1 charged=103.0 key=c_hit_lo_m4_bin32
.Lb1a_fix64:
adds x_cnt, x_cnt, #512
b.eq .Lb1a_done
tbz x_cnt, #8, .Lb1af_b16
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
.Lb1af_b16:
tbz x_cnt, #7, .Lb1af_b8
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
.Lb1af_b8:
tbz x_cnt, #6, .Lb1af_b4
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
ldp x4, x5, [x_p], #16
staddb w_lp, [x4]
staddb w_lp, [x5]
ldp x6, x7, [x_p], #16
staddb w_lp, [x6]
staddb w_lp, [x7]
.Lb1af_b4:
tbz x_cnt, #5, .Lb1af_b2
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
ldp x2, x3, [x_p], #16
staddb w_lp, [x2]
staddb w_lp, [x3]
.Lb1af_b2:
tbz x_cnt, #4, .Lb1af_b1
ldp x0, x1, [x_p], #16
staddb w_lp, [x0]
staddb w_lp, [x1]
.Lb1af_b1:
tbz x_cnt, #3, .Lb1a_done
ldr x0, [x_p], #8
staddb w_lp, [x0]
.Lb1a_done:
```

### 1.4 `M = 15` Full-Chord Slab-Routed Band-2 (Pre-Rotated `W_B = 1` on `[2^52.5, 2^55.5)`, Shared `W_B = 1` on `[2^55.5, 2^60.0)`, `W_B = 8` on `[2^60.0, 2^61.3)`; Direct 1-Slot Tier 2c Inactive at `p_g = [61.3, 61.3] == p_dir`)

All Band-2 records store `rec = ((k & 31) << 59) | (0 << 58) | (m << 56) | oid` (`ubfiz x3, x0, #59, #5`) across `M = 15` slab owners (`480 > 396.92` cols, zero ring collisions). Active Band-2 primes span `3` constant-weight bands per side (`w_lp in {13, 14, 15}` at `q = 4.0`b `round`) sharing `15` `512 KiB`-aligned TX buffers `tx_tbl[0..14]` (`tx_tbl[-1], tx_tbl[15]` = sinks). At `128 KiB` (`tbnz #17`), register-disjoint `rev_b2_rx_drain_{pre64,u64,bin32}` drains at `1.5379 <= 1.54` (`1.53125 + 0.015`) using dead scratch `{x0, x4..x7, x9, x16, x17}`; end-of-band partial drains cost `<= 4,905 <= 20,700` ops/side/chord.

1. **Tier 2a (`W_B = 1`: Pre-Rotated `[2^52.5, 2^55.5)` + Shared `[2^55.5, 2^60.0)` on `g`, `f`):** Pre-rotated 4-table split (`<= 4.442 EiB`) calls `rev_b2_emit_id_u16` on all 4/4 conjugates (`3.73375 <= 3.735 = c_hit_hi_rot`); shared `[2^55.5, 2^60.0)` gives `4.48375 <= 4.485 = c_hit_hi_u4`. `rev_w1_probe_fast` (`12`) + `rev_b2_w1_wrap` (`5`) + `tail_slab` (`2`) = `19 <= 20.0 = c_tile_hdr_w1`; `slab_hdr` (`11`) + `bin8421` (`<= 8.875`) = `19.875 <= 19.9 = c_sb_hdr_u4` (zero `f` overhang).
2. **Tier 2b (`W_B = 8`, `[2^60.0, 2^61.3)` `g`-only):** `rev_w8_probe_fast` (`13`) + `rev_b2_w8_wrap` (`6`) = `19 <= 20.0 = c_tile_hdr_w8`; `slab_hdr` (`16`) loads `ldp x11, x10, [x_tx]`, splits boundary bucket (`1/4`) via `emit_split_*` (`tbnz #17` every `32 B`), emits `3` interior buckets (`3/4`) via `emit_*` into `x10` (`5.311875 <= 5.313 = c_hit_hi_u2`; `30.875 <= 30.9 = c_sb_hdr_u2`; `cand_bd = 1 + 8/Lb_k`; `occ(tx[s]) <= 131,216 B << 256 KiB` unconditionally, verified in `w8_verifiers.py`).
3. **Tier 2c (`Direct 1-Slot Sparse u2`, inactive since `p_g = [61.3, 61.3] == p_dir`):** Preserved ASM blocks below (`rev_sb_probe_fast`, `rev_b2_sb_wrap`, `rev_b2_sieve_direct1_{u2,hit2,tail1}`, `c_tile_hdr = 18.0`, `3.5 + 13.0*p_ne + 11.0*ovf_m + 1.54*lam`) remain unused at runtime.

```asm
; kernel: rev_w1_probe_fast items=1 charged=12.0 key=c_tile_hdr_w1_probe
cbz x_lb, .Lw1_done
sub x5, x_p, x_k0
cmp x5, x_lb
csel x5, x5, x_lb, lo
add x6, x_k0, x5
sub x0, x_col0, x_k0
lsl x0, x0, #59
orr x_t0, x0, x_crot
lsr x1, x_col0, #5
add x_tx, x_tx_base, x1, lsl #3
add x_oc, x_off, x_k0, lsl #3
add x_oe, x_off, x6, lsl #3
```
```asm
; kernel: rev_b2_w1_wrap items=1 charged=11.0 key=c_b2_w1_wrap
and x2, x_col0, #31
cbz x2, .Lw1_aligned
mov x3, #32
sub x3, x3, x2
add x_on, x_oc, x3, lsl #3
cmp x_on, x_oe
csel x_on, x_oe, x_on, hi
b .Lw1_prefix_emit
.Lw1_aligned:
sub x7, x_oe, x_oc
lsr x21, x7, #8
cbz x21, .Lw1_tail_slab
```
```asm
; kernel: rev_b2_w1_slab_hdr items=1 charged=11.0 key=c_sb_hdr_w1
.Lw1_next_slab:
ldr x2, [x_oc]
ldr x3, [x_oc, #256]
add x_oc, x_oc, #256
subs x_cnt, x3, x2
b.eq .Lw1_step
ldr x10, [x_tx]
tbnz x10, #17, .Lw1_flush_x10
.Lw1_do_emit:
str x10, [x_tx]
.Lw1_step:
add x_tx, x_tx, #8
subs x21, x21, #1
b.ne .Lw1_next_slab
```
```asm
; kernel: rev_b2_w1_tail_slab items=1 charged=9.0 key=c_b2_w1_tail
.Lw1_tail_slab:
cmp x_oc, x_oe
b.hs .Lw1_done
ldr x2, [x_oc]
ldr x3, [x_oe]
subs x_cnt, x3, x2
b.eq .Lw1_done
ldr x10, [x_tx]
tbnz x10, #17, .Lw1_flush_t10
.Lw1_tail_emit:
str x10, [x_tx]
.Lw1_done:
```
```asm
; kernel: rev_b2_w1_prefix_emit items=1 charged=10.0 key=c_b2_w1_prefix
.Lw1_prefix_emit:
ldr x2, [x_oc]
ldr x3, [x_on]
mov x_oc, x_on
subs x_cnt, x3, x2
b.eq .Lw1_pre_step
ldr x10, [x_tx]
tbnz x10, #17, .Lw1_flush_p10
.Lw1_pre_emit:
str x10, [x_tx]
.Lw1_pre_step:
add x_tx, x_tx, #8
b .Lw1_aligned
```
```asm
; kernel: rev_b2_flush_x10 items=16384 charged=0.01 key=c_b2_flush_x10
.Lw1_flush_x10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw1_do_emit
.Lw1_flush_t10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw1_tail_emit
.Lw1_flush_p10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw1_pre_emit
.Lb2i_fl16:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lb2i_fl16_ret
.Lb2c_fl16:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lb2c_fl16_ret
```
```asm
; kernel: rev_b2_flush_w8_x10 items=16384 charged=0.01 key=c_b2_flush_w8_x10
.Lw8_flush_x10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw8_chk_x11
.Lw8_flush_t10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw8_tchk_x11
.Lw8_flush_p10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw8_pre_emit
.Lw8_flush_pa10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw8_pa_emit
.Lw8_flush_pat10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lw8_pat_emit
.Lspi_fl10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lspi_fl10_r
.Lspc_fl10:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lspc_fl10_r
```
```asm
; kernel: rev_b2_flush_x11 items=16384 charged=0.01 key=c_b2_flush_x11
.Lw8_flush_x11:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lw8_do_emit
.Lw8_flush_t11:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lw8_tdo_emit
.Lspi_fl11:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lspi_fl11_r
.Lspc_fl11:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lspc_fl11_r
```
```asm
; kernel: rev_b2_rx_drain_pre64 items=1 charged=2.0 key=c_b2_rx_drain_pre64
subs x_fcnt, x_fcnt, #512
b.lo .Lb2d_fix64
```
```asm
; kernel: rev_b2_rx_drain_u64 items=64 charged=1.535 key=c_b2_rx_drain_u64
.Lb2d64:
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
subs x_fcnt, x_fcnt, #512
b.hs .Lb2d64
```
```asm
; kernel: rev_b2_rx_drain_bin32 items=1 charged=103.0 key=c_b2_rx_drain_bin32
.Lb2d_fix64:
adds x_fcnt, x_fcnt, #512
b.eq .Lb2d_done
tbz x_fcnt, #8, .Lb2d_b16
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
.Lb2d_b16:
tbz x_fcnt, #7, .Lb2d_b8
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
.Lb2d_b8:
tbz x_fcnt, #6, .Lb2d_b4
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
ldp x7, x9, [x_rx], #16
staddb w_lp, [x7]
staddb w_lp, [x9]
ldp x16, x17, [x_rx], #16
staddb w_lp, [x16]
staddb w_lp, [x17]
.Lb2d_b4:
tbz x_fcnt, #5, .Lb2d_b2
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
ldp x5, x6, [x_rx], #16
staddb w_lp, [x5]
staddb w_lp, [x6]
.Lb2d_b2:
tbz x_fcnt, #4, .Lb2d_b1
ldp x0, x4, [x_rx], #16
staddb w_lp, [x0]
staddb w_lp, [x4]
.Lb2d_b1:
tbz x_fcnt, #3, .Lb2d_done
ldr x0, [x_rx], #8
staddb w_lp, [x0]
.Lb2d_done:
```
```asm
; kernel: rev_b2_emit_id_pre16 items=1 charged=2.0 key=c_b2_emit_id_pre16
subs x_cnt, x_cnt, #128
b.lo .Lb2i_fix16
```
```asm
; kernel: rev_b2_emit_id_u16 items=16 charged=2.20 key=c_b2_emit_id_u16
.Lb2_id16:
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
stp x8, x9, [x10], #16
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
stp x8, x9, [x10], #16
tbnz x10, #17, .Lb2i_fl16
.Lb2i_fl16_ret:
subs x_cnt, x_cnt, #128
b.hs .Lb2_id16
```
```asm
; kernel: rev_b2_emit_id_bin8421 items=1 charged=37.0 key=c_b2_emit_id_bin8421
.Lb2i_fix16:
adds x_cnt, x_cnt, #128
b.eq .Lb2i_emit_done
tbz x_cnt, #6, .Lb2i_b4
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
stp x8, x9, [x10], #16
.Lb2i_b4:
tbz x_cnt, #5, .Lb2i_b2
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
stp x4, x5, [x10], #16
.Lb2i_b2:
tbz x_cnt, #4, .Lb2i_b1
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
stp x0, x1, [x10], #16
.Lb2i_b1:
tbz x_cnt, #3, .Lb2i_emit_done
ldr x0, [x2], #8
add x0, x0, x_t0
str x0, [x10], #8
.Lb2i_emit_done:
```
```asm
; kernel: rev_b2_emit_c123_pre16 items=1 charged=2.0 key=c_b2_emit_c123_pre16
subs x_cnt, x_cnt, #128
b.lo .Lb2c_fix16
```
```asm
; kernel: rev_b2_emit_c123_u16 items=16 charged=3.20 key=c_b2_emit_c123_u16
.Lb2_c123_16:
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
and x0, x0, #0xfbffffffffffffff
and x1, x1, #0xfbffffffffffffff
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
and x4, x4, #0xfbffffffffffffff
and x5, x5, #0xfbffffffffffffff
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
and x6, x6, #0xfbffffffffffffff
and x7, x7, #0xfbffffffffffffff
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
and x8, x8, #0xfbffffffffffffff
and x9, x9, #0xfbffffffffffffff
stp x8, x9, [x10], #16
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
and x0, x0, #0xfbffffffffffffff
and x1, x1, #0xfbffffffffffffff
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
and x4, x4, #0xfbffffffffffffff
and x5, x5, #0xfbffffffffffffff
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
and x6, x6, #0xfbffffffffffffff
and x7, x7, #0xfbffffffffffffff
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
and x8, x8, #0xfbffffffffffffff
and x9, x9, #0xfbffffffffffffff
stp x8, x9, [x10], #16
tbnz x10, #17, .Lb2c_fl16
.Lb2c_fl16_ret:
subs x_cnt, x_cnt, #128
b.hs .Lb2_c123_16
```
```asm
; kernel: rev_b2_emit_c123_bin8421 items=1 charged=52.0 key=c_b2_emit_c123_bin8421
.Lb2c_fix16:
adds x_cnt, x_cnt, #128
b.eq .Lb2c_emit_done
tbz x_cnt, #6, .Lb2c_b4
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
and x0, x0, #0xfbffffffffffffff
and x1, x1, #0xfbffffffffffffff
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
and x4, x4, #0xfbffffffffffffff
and x5, x5, #0xfbffffffffffffff
stp x4, x5, [x10], #16
ldp x6, x7, [x2], #16
add x6, x6, x_t0
add x7, x7, x_t0
and x6, x6, #0xfbffffffffffffff
and x7, x7, #0xfbffffffffffffff
stp x6, x7, [x10], #16
ldp x8, x9, [x2], #16
add x8, x8, x_t0
add x9, x9, x_t0
and x8, x8, #0xfbffffffffffffff
and x9, x9, #0xfbffffffffffffff
stp x8, x9, [x10], #16
.Lb2c_b4:
tbz x_cnt, #5, .Lb2c_b2
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
and x0, x0, #0xfbffffffffffffff
and x1, x1, #0xfbffffffffffffff
stp x0, x1, [x10], #16
ldp x4, x5, [x2], #16
add x4, x4, x_t0
add x5, x5, x_t0
and x4, x4, #0xfbffffffffffffff
and x5, x5, #0xfbffffffffffffff
stp x4, x5, [x10], #16
.Lb2c_b2:
tbz x_cnt, #4, .Lb2c_b1
ldp x0, x1, [x2], #16
add x0, x0, x_t0
add x1, x1, x_t0
and x0, x0, #0xfbffffffffffffff
and x1, x1, #0xfbffffffffffffff
stp x0, x1, [x10], #16
.Lb2c_b1:
tbz x_cnt, #3, .Lb2c_emit_done
ldr x0, [x2], #8
add x0, x0, x_t0
and x0, x0, #0xfbffffffffffffff
str x0, [x10], #8
.Lb2c_emit_done:
```
```asm
; kernel: rev_w8_probe_fast items=1 charged=13.0 key=c_tile_hdr_w8_probe
cbz x_lb, .Lw8_done
sub x5, x_p, x_k0
cmp x5, x_lb
csel x5, x5, x_lb, lo
add x6, x_k0, x5
add x6, x6, #7
lsr x6, x6, #3
add x_oe, x_off_base, x6, lsl #3
lsr x21, x_k0, #3
add x_off, x_off_base, x21, lsl #3
sub x0, x_col0, x_k0
lsl x1, x0, #59
orr x_t0, x1, x_crot
```
```asm
; kernel: rev_b2_w8_wrap items=1 charged=6.0 key=c_b2_w8_wrap
cbnz x_col0, .Lw8_r2_pre
cbz x1, .Lw8_pure_align
sub x_tx, x_tx_base, #8
.Lw8_aligned:
sub x7, x_oe, x_off
lsr x21, x7, #5
cbz x21, .Lw8_tail_slab
```
```asm
; kernel: rev_b2_w8_slab_hdr items=1 charged=16.0 key=c_sb_hdr_w8
.Lw8_next_slab:
ldp x2, x3, [x_off]
ldr x8, [x_off, #32]
add x_off, x_off, #32
ldp x11, x10, [x_tx]
tbnz x10, #17, .Lw8_flush_x10
.Lw8_chk_x11:
tbnz x11, #17, .Lw8_flush_x11
.Lw8_do_emit:
subs x_cnt, x3, x2
b.eq .Lw8_int
.Lw8_int:
subs x_cnt, x8, x3
b.eq .Lw8_step
.Lw8_step:
stp x11, x10, [x_tx]
add x_tx, x_tx, #8
subs x21, x21, #1
b.ne .Lw8_next_slab
```
```asm
; kernel: rev_b2_w8_tail_slab items=1 charged=15.0 key=c_b2_w8_tail
.Lw8_tail_slab:
cmp x_off, x_oe
b.hs .Lw8_done
ldp x2, x3, [x_off]
ldr x8, [x_oe]
ldp x11, x10, [x_tx]
tbnz x10, #17, .Lw8_flush_t10
.Lw8_tchk_x11:
tbnz x11, #17, .Lw8_flush_t11
.Lw8_tdo_emit:
subs x_cnt, x3, x2
b.eq .Lw8_tint
.Lw8_tint:
subs x_cnt, x8, x3
b.eq .Lw8_tstep
.Lw8_tstep:
stp x11, x10, [x_tx]
.Lw8_done:
```
```asm
; kernel: rev_b2_w8_prefix_pure items=1 charged=63.0 key=c_b2_w8_prefix_pure
.Lw8_r2_pre:
sub x2, xzr, x0
lsr x2, x2, #3
sub x2, x2, x21
and x2, x2, #3
cbz x2, .Lw8_r2_settx
add x_on, x_off, x2, lsl #3
cmp x_on, x_oe
csel x_on, x_oe, x_on, hi
add x3, x_col0, x21, lsl #3
sub x3, x3, x_k0
asr x3, x3, #5
add x_tx, x_tx_base, x3, lsl #3
add x21, x21, x2
ldr x2, [x_off]
ldr x3, [x_on]
mov x_off, x_on
subs x_cnt, x3, x2
b.eq .Lw8_r2_settx
ldr x10, [x_tx]
tbnz x10, #17, .Lw8_flush_p10
.Lw8_pre_emit:
str x10, [x_tx]
.Lw8_r2_settx:
cmp x_off, x_oe
b.hs .Lw8_done
lsr x1, x_t0, #59
cbz x1, .Lw8_pure_align
add x2, x_col0, x21, lsl #3
sub x2, x2, x_k0
sub x2, x2, #1
asr x2, x2, #5
add x_tx, x_tx_base, x2, lsl #3
b .Lw8_aligned
.Lw8_pure_align:
add x2, x_col0, x21, lsl #3
sub x2, x2, x_k0
asr x2, x2, #5
add x_tx, x_tx_base, x2, lsl #3
sub x7, x_oe, x_off
lsr x21, x7, #5
cbz x21, .Lw8_pa_tail
.Lw8_pa_lp:
ldr x2, [x_off]
ldr x3, [x_off, #32]
add x_off, x_off, #32
subs x_cnt, x3, x2
b.eq .Lw8_pa_step
ldr x10, [x_tx]
tbnz x10, #17, .Lw8_flush_pa10
.Lw8_pa_emit:
str x10, [x_tx]
.Lw8_pa_step:
add x_tx, x_tx, #8
subs x21, x21, #1
b.ne .Lw8_pa_lp
.Lw8_pa_tail:
cmp x_off, x_oe
b.hs .Lw8_done
ldr x2, [x_off]
ldr x3, [x_oe]
subs x_cnt, x3, x2
b.eq .Lw8_done
ldr x10, [x_tx]
tbnz x10, #17, .Lw8_flush_pat10
.Lw8_pat_emit:
str x10, [x_tx]
b .Lw8_done
```
```asm
; kernel: rev_b2_emit_split_id_pre4 items=1 charged=2.0 key=c_b2_split_id_pre4
subs x_cnt, x_cnt, #32
b.lo .Lb2spi_fix4
```
```asm
; kernel: rev_b2_emit_split_id_u4 items=4 charged=6.5 key=c_b2_split_id_u4
.Lb2sp_id4:
ldp x0, x1, [x2], #16
adds x0, x0, x_t0
b.hs .Lspi0
str x0, [x11], #8
b .Lspi0_n
.Lspi0:
str x0, [x10], #8
.Lspi0_n:
adds x1, x1, x_t0
b.hs .Lspi1
str x1, [x11], #8
b .Lspi1_n
.Lspi1:
str x1, [x10], #8
.Lspi1_n:
ldp x4, x5, [x2], #16
adds x4, x4, x_t0
b.hs .Lspi2
str x4, [x11], #8
b .Lspi2_n
.Lspi2:
str x4, [x10], #8
.Lspi2_n:
adds x5, x5, x_t0
b.hs .Lspi3
str x5, [x11], #8
b .Lspi3_n
.Lspi3:
str x5, [x10], #8
.Lspi3_n:
tbnz x10, #17, .Lspi_fl10
.Lspi_fl10_r:
tbnz x11, #17, .Lspi_fl11
.Lspi_fl11_r:
subs x_cnt, x_cnt, #32
b.hs .Lb2sp_id4
```
```asm
; kernel: rev_b2_emit_split_id_bin21 items=1 charged=21.0 key=c_b2_split_id_bin21
.Lb2spi_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb2spi_done
tbz x_cnt, #4, .Lspi_b1
ldp x0, x1, [x2], #16
adds x0, x0, x_t0
b.hs .Lspi_t0
str x0, [x11], #8
b .Lspi_t0n
.Lspi_t0:
str x0, [x10], #8
.Lspi_t0n:
adds x1, x1, x_t0
b.hs .Lspi_t1
str x1, [x11], #8
b .Lspi_t1n
.Lspi_t1:
str x1, [x10], #8
.Lspi_t1n:
.Lspi_b1:
tbz x_cnt, #3, .Lb2spi_done
ldr x0, [x2], #8
adds x0, x0, x_t0
b.hs .Lspi_t2
str x0, [x11], #8
b .Lb2spi_done
.Lspi_t2:
str x0, [x10], #8
.Lb2spi_done:
```
```asm
; kernel: rev_b2_emit_split_c123_pre4 items=1 charged=2.0 key=c_b2_split_c123_pre4
subs x_cnt, x_cnt, #32
b.lo .Lb2spc_fix4
```
```asm
; kernel: rev_b2_emit_split_c123_u4 items=4 charged=7.5 key=c_b2_split_c123_u4
.Lb2sp_c4:
ldp x0, x1, [x2], #16
adds x0, x0, x_t0
and x0, x0, #0xfbffffffffffffff
b.hs .Lspc0
str x0, [x11], #8
b .Lspc0_n
.Lspc0:
str x0, [x10], #8
.Lspc0_n:
adds x1, x1, x_t0
and x1, x1, #0xfbffffffffffffff
b.hs .Lspc1
str x1, [x11], #8
b .Lspc1_n
.Lspc1:
str x1, [x10], #8
.Lspc1_n:
ldp x4, x5, [x2], #16
adds x4, x4, x_t0
and x4, x4, #0xfbffffffffffffff
b.hs .Lspc2
str x4, [x11], #8
b .Lspc2_n
.Lspc2:
str x4, [x10], #8
.Lspc2_n:
adds x5, x5, x_t0
and x5, x5, #0xfbffffffffffffff
b.hs .Lspc3
str x5, [x11], #8
b .Lspc3_n
.Lspc3:
str x5, [x10], #8
.Lspc3_n:
tbnz x10, #17, .Lspc_fl10
.Lspc_fl10_r:
tbnz x11, #17, .Lspc_fl11
.Lspc_fl11_r:
subs x_cnt, x_cnt, #32
b.hs .Lb2sp_c4
```
```asm
; kernel: rev_b2_emit_split_c123_bin21 items=1 charged=24.0 key=c_b2_split_c123_bin21
.Lb2spc_fix4:
adds x_cnt, x_cnt, #32
b.eq .Lb2spc_done
tbz x_cnt, #4, .Lspc_b1
ldp x0, x1, [x2], #16
adds x0, x0, x_t0
and x0, x0, #0xfbffffffffffffff
b.hs .Lspc_t0
str x0, [x11], #8
b .Lspc_t0n
.Lspc_t0:
str x0, [x10], #8
.Lspc_t0n:
adds x1, x1, x_t0
and x1, x1, #0xfbffffffffffffff
b.hs .Lspc_t1
str x1, [x11], #8
b .Lspc_t1n
.Lspc_t1:
str x1, [x10], #8
.Lspc_t1n:
.Lspc_b1:
tbz x_cnt, #3, .Lb2spc_done
ldr x0, [x2], #8
adds x0, x0, x_t0
and x0, x0, #0xfbffffffffffffff
b.hs .Lspc_t2
str x0, [x11], #8
b .Lb2spc_done
.Lspc_t2:
str x0, [x10], #8
.Lb2spc_done:
```
```asm
; kernel: rev_sb_probe_fast items=1 charged=13.0 key=c_tile_hdr_probe
cbz x_lb, .Lsd1_done
lsr x4, x_k0, #5
and x0, x_k0, #31
sub x0, x_col0, x0
asr x1, x0, #5
add x_tx, x_tx_base, x1, lsl #3
lsl x0, x0, #59
orr x_t0, x0, x_crot
sub x5, x_p, x_k0
cmp x5, x_lb
csel x5, x5, x_lb, lo
add x6, x_k0, x5
add x_slot, x_slot_base, x4, lsl #3
```
```asm
; kernel: rev_b2_sb_wrap items=1 charged=5.0 key=c_b2_sb_wrap
sub x6, x6, #1
lsr x6, x6, #5
sub x21, x6, x4
subs x21, x21, #1
b.lo .Lsd1_tail1
```
**Direct 1-Slot Sparse Band-2 (inactive at `p_g = [61.3, 61.3] == p_dir`; preserved ASM blocks below):** Overflow heads store `(1 << 58) | node_idx` decoded via `add x4, x_ovf_adj, x1, lsl #4`; `rev_b2_sieve_direct1_{u2,hit2,tail1}` (`3.5 + 13.0*p_ne + 11.0*ovf_m + 1.54*lam`) remain unused at runtime.

```asm
; kernel: rev_b2_sieve_direct1_u2 items=2 charged=3.5 key=c_b2_dir_empty2
.Lsd1_u2:
ldp x1, x2, [x_slot], #16
and x3, x1, x2
cmn x3, #1
b.ne .Lsd1_hit2
.Lsd1_next2:
add x_tx, x_tx, #16
subs x21, x21, #2
b.hs .Lsd1_u2
```
```asm
; kernel: rev_b2_sieve_direct1_tail1 items=1 charged=57.0 key=c_b2_dir_tail1
.Lsd1_tail1:
adds x21, x21, #2
cbz x21, .Lsd1_done
ldr x1, [x_slot]
cmn x1, #1
b.eq .Lsd1_done
tbnz x1, #58, .Lsd1_ovf_t
adds x1, x1, x_t0
cset x6, hs
and x1, x1, #0xfbffffffffffffff
ldr x7, [x_tx, x6, lsl #3]
str x1, [x7], #8
str x7, [x_tx, x6, lsl #3]
tbnz x7, #17, .Lsd1_flt_s
b .Lsd1_done
.Lsd1_flt_s:
and x_rx, x7, #0xfffffffffff80000
sub x_fcnt, x7, x_rx
str x_rx, [x_tx, x6, lsl #3]
b .Lsd1_done
.Lsd1_ovf_t:
ldp x10, x11, [x_tx]
.Lsd1_ovf_tlp:
add x4, x_ovf_adj, x1, lsl #4
ldp x5, x1, [x4]
adds x5, x5, x_t0
and x5, x5, #0xfbffffffffffffff
b.hs .Lsd1_ovf_thi
str x5, [x10], #8
tbnz x10, #17, .Lsd1_flt_lo
b .Lsd1_ovf_tc
.Lsd1_ovf_thi:
str x5, [x11], #8
tbnz x11, #17, .Lsd1_flt_hi
.Lsd1_ovf_tc:
tbnz x1, #58, .Lsd1_ovf_tlp
adds x1, x1, x_t0
and x1, x1, #0xfbffffffffffffff
b.hs .Lsd1_ovf_t2hi
str x1, [x10], #8
tbnz x10, #17, .Lsd1_flt_tlo
b .Lsd1_ovf_tend
.Lsd1_ovf_t2hi:
str x1, [x11], #8
tbnz x11, #17, .Lsd1_flt_thi
.Lsd1_ovf_tend:
stp x10, x11, [x_tx]
b .Lsd1_done
.Lsd1_flt_lo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf_tc
.Lsd1_flt_hi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf_tc
.Lsd1_flt_tlo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf_tend
.Lsd1_flt_thi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf_tend
.Lsd1_done:
ret
```
```asm
; kernel: rev_b2_sieve_direct1_hit2 items=1 charged=130.0 key=c_b2_dir_hit2
.Lsd1_hit2:
cmn x1, #1
b.eq .Lsd1_b1
tbnz x1, #58, .Lsd1_ovf0
adds x1, x1, x_t0
cset x6, hs
and x1, x1, #0xfbffffffffffffff
ldr x7, [x_tx, x6, lsl #3]
str x1, [x7], #8
str x7, [x_tx, x6, lsl #3]
tbnz x7, #17, .Lsd1_fl0
.Lsd1_b1:
cmn x2, #1
b.eq .Lsd1_next2
tbnz x2, #58, .Lsd1_ovf1
adds x2, x2, x_t0
csinc x6, x_one, x_one, lo
and x2, x2, #0xfbffffffffffffff
ldr x7, [x_tx, x6, lsl #3]
str x2, [x7], #8
str x7, [x_tx, x6, lsl #3]
tbnz x7, #17, .Lsd1_fl1
b .Lsd1_next2
.Lsd1_fl0:
and x_rx, x7, #0xfffffffffff80000
sub x_fcnt, x7, x_rx
str x_rx, [x_tx, x6, lsl #3]
b .Lsd1_b1
.Lsd1_fl1:
and x_rx, x7, #0xfffffffffff80000
sub x_fcnt, x7, x_rx
str x_rx, [x_tx, x6, lsl #3]
b .Lsd1_next2
.Lsd1_ovf0:
ldp x10, x11, [x_tx]
.Lsd1_ovf0_lp:
add x4, x_ovf_adj, x1, lsl #4
ldp x5, x1, [x4]
adds x5, x5, x_t0
and x5, x5, #0xfbffffffffffffff
b.hs .Lsd1_ovf0_hi
str x5, [x10], #8
tbnz x10, #17, .Lsd1_fl0_lo
b .Lsd1_ovf0_c
.Lsd1_ovf0_hi:
str x5, [x11], #8
tbnz x11, #17, .Lsd1_fl0_hi
.Lsd1_ovf0_c:
tbnz x1, #58, .Lsd1_ovf0_lp
adds x1, x1, x_t0
and x1, x1, #0xfbffffffffffffff
b.hs .Lsd1_ovf0_thi
str x1, [x10], #8
tbnz x10, #17, .Lsd1_fl0_tlo
b .Lsd1_ovf0_end
.Lsd1_ovf0_thi:
str x1, [x11], #8
tbnz x11, #17, .Lsd1_fl0_thi
.Lsd1_ovf0_end:
stp x10, x11, [x_tx]
b .Lsd1_b1
.Lsd1_fl0_lo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf0_c
.Lsd1_fl0_hi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf0_c
.Lsd1_fl0_tlo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf0_end
.Lsd1_fl0_thi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf0_end
.Lsd1_ovf1:
ldp x10, x11, [x_tx, #8]
.Lsd1_ovf1_lp:
add x4, x_ovf_adj, x2, lsl #4
ldp x5, x2, [x4]
adds x5, x5, x_t0
and x5, x5, #0xfbffffffffffffff
b.hs .Lsd1_ovf1_hi
str x5, [x10], #8
tbnz x10, #17, .Lsd1_fl1_lo
b .Lsd1_ovf1_c
.Lsd1_ovf1_hi:
str x5, [x11], #8
tbnz x11, #17, .Lsd1_fl1_hi
.Lsd1_ovf1_c:
tbnz x2, #58, .Lsd1_ovf1_lp
adds x2, x2, x_t0
and x2, x2, #0xfbffffffffffffff
b.hs .Lsd1_ovf1_thi
str x2, [x10], #8
tbnz x10, #17, .Lsd1_fl1_tlo
b .Lsd1_ovf1_end
.Lsd1_ovf1_thi:
str x2, [x11], #8
tbnz x11, #17, .Lsd1_fl1_thi
.Lsd1_ovf1_end:
stp x10, x11, [x_tx, #8]
b .Lsd1_next2
.Lsd1_fl1_lo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf1_c
.Lsd1_fl1_hi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf1_c
.Lsd1_fl1_tlo:
and x_rx, x10, #0xfffffffffff80000
sub x_fcnt, x10, x_rx
mov x10, x_rx
b .Lsd1_ovf1_end
.Lsd1_fl1_thi:
and x_rx, x11, #0xfffffffffff80000
sub x_fcnt, x11, x_rx
mov x11, x_rx
b .Lsd1_ovf1_end
```

### 1.5 Fused `uint8_t S_b1 + uint8_t S_b2` 2x/4x SWAR Scan & Assembled Pure-ARMv8.1-A Cold-Path 2D NP Gate (`134/256 + 78/512 + pad/cold = 0.72190 <= 0.7222 = c_scan_zero`)

Each column's active bank slice is padded to a multiple of `512` cells (`<= 511` zero pad cells out of `2^{55.57}`, `< 2^{-46}` overhead; `M_g` base `x23` is 64B-aligned). Shared `4.00`b `round` weights with column-constant `T_g_lo(b, c) in [30, 115]` (`sl_g_tgt = 428.0`, `bias_g = 128 - T_g_lo`, `delta_jw = 55`) guarantee `S_b1 + S_b2 <= 248 <= 255` on `g` (realized `<= 245`) and `<= 190 <= 240 <= 255` on `f` (zero carry in 64-bit SWAR `add`). **Complete Assembled Pure-ARMv8.1-A Scan Suite (4 blocks below; zero unassembled cold-path padding, zero new ISA, live-emulated instruction-by-instruction in `canary1_verifiers.py::verify_scan_asm_live`):**
- **Pass 1 (`g`) `rev_scan_g_u8_fused_256` (`134/256 = 0.5234375` Ir/cell) + `rev_scan_g_cold_128` (`63` static ops):** Hot loop (`x21=S_b1, x22=S_b2, x23=M_g, x24=b8_nxt, x25=0x8080..., x26=cnt, x27=ctx`; 96B `[x27]` layout: `#0=G_val_ptr, #8=bias_g, #16..55=saved x21,x22,x24,x26,x23, #64..95=32B quarter spill`) preserves the 4-quarter ORs in `x16 (x0..x3)`, `x17 (x4..x7)`, `x19 (x8..x11)`, `x20 (x12..x15)` at zero extra hot-path ops. On any bit-7 hit, `.Lhit128_g` tests `x17, x16, x20, x19`, spills only triggered 4-word quarters to `[x27, #64..95]` (`x4..x7` are dead after Q1), runs one shared strict MSB-first byte-extraction loop pushing `S_g = byte - bias_g` to `G_val` (`[x27]`) and setting `M_g`, stores `M_g`, and dispatches back to `.Lsg_mid`/`.Lsg_tail` via `tbnz x23, #4, .Lsg_mid`. Exact per-path op count for any non-empty 128-cell bitmask (`k >= 1` hits, `q_trig <= min(4, k)` triggered quarters): fixed save/dispatch/restore minus skipped hot `stp` `<= 22`, per-triggered-quarter spill+word-scan delta `<= 34 < 35` (`Q1:+31, Q0:+32, Q2:+33, Q3:+34`), and per-hit extraction `= 13`, giving **`<= 22 + 35*q_trig + 13*k <= 22 + 48*k`** ops/block (checked across `808` structured + random cases in `verify_scan_asm_live`).
- **Pass 2 (`f`) `rev_scan_f_u8_512` (`78/512 = 0.15234375` Ir/cell) + `rev_scan_f_cold_128` (`37` static ops):** Hot loop (`x21=S_b1, x22=S_b2, x23=M_g, x24=init_f, x26=cnt`) runs `ldp x0, x1, [x23], #16; orr x0, x0, x1; cbnz x0, .Lhit128_f`. `.Lhit128_f` reloads `(x0, x1)` via **`ldp x0, x1, [x23, #-16]`**, pops `G_val` (`x19`) in exact 1-to-1 MSB-first FIFO order (`x14=63, x15=127, x25=b1_base, w27=T_jw, x28=1<<63`), emits 2D-NP survivors (`S_f + S_g >= T_jw`) to `x20`, and dispatches to `.Lsf_ref0..3` via `tbnz x23, #5/#4` BEFORE refilling `S_b1`/`S_b2`. Exact per-path op count for any non-empty 128-cell bitmask (`k >= 1` hits): reload + two half-word checks + dispatch `<= 7` ops, plus `<= 14` ops per hit (`14` on survivor, `12` on reject), giving **`<= 7 + 14*k`** ops/block.
- **Combined Cold + Pad + Wrap:** By linearity of expectation, union bound, and concavity of `min(1, 128*s)` (Jensen's inequality on `sga_mean = 2.0476e-4`), `p_hit128 = min(1, 128 * sga_mean) = 0.026209` -> `scan_cold = p_hit128 * (29/128) + 62 * sga_mean = 0.018633` Ir/cell; `pad/wrap = (7/92.654)*(10/128) + 2/92.654 = 0.027488` Ir/cell (where the conservative `2/lb_live` per-cell term also absorbs per-column pointer/constant setup for `x14..x28`); total **`0.721902 <= 0.7222 = c_scan_zero`**!

```asm
; kernel: rev_scan_zero_wrap items=1 charged=2.0 key=c_scan_zero_wrap
cbz x_cnt, .Lscan_done
.Lscan_done:
ret
```
```asm
; kernel: rev_scan_g_u8_fused_256 items=256 charged=0.525 key=c_scan_g_raw
.Lsg_loop256:
ldp x0, x1, [x21]
ldp x16, x17, [x22]
add x0, x0, x16
add x1, x1, x17
ldp x2, x3, [x21, #16]
ldp x16, x17, [x22, #16]
add x2, x2, x16
add x3, x3, x17
ldp x4, x5, [x21, #32]
ldp x16, x17, [x22, #32]
add x4, x4, x16
add x5, x5, x17
ldp x6, x7, [x21, #48]
ldp x16, x17, [x22, #48]
add x6, x6, x16
add x7, x7, x17
ldp x8, x9, [x21, #64]
ldp x16, x17, [x22, #64]
add x8, x8, x16
add x9, x9, x17
ldp x10, x11, [x21, #80]
ldp x16, x17, [x22, #80]
add x10, x10, x16
add x11, x11, x17
ldp x12, x13, [x21, #96]
ldp x16, x17, [x22, #96]
add x12, x12, x16
add x13, x13, x17
ldp x14, x15, [x21, #112]
ldp x16, x17, [x22, #112]
add x14, x14, x16
add x15, x15, x17
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
orr x16, x0, x1
orr x28, x2, x3
orr x16, x16, x28
orr x17, x4, x5
orr x28, x6, x7
orr x17, x17, x28
orr x19, x8, x9
orr x28, x10, x11
orr x19, x19, x28
orr x20, x12, x13
orr x28, x14, x15
orr x20, x20, x28
orr x28, x16, x17
orr x28, x28, x19
orr x28, x28, x20
tst x28, x25
b.ne .Lhit128_g
stp xzr, xzr, [x23], #16
.Lsg_mid:
ldp x0, x1, [x21]
ldp x16, x17, [x22]
add x0, x0, x16
add x1, x1, x17
ldp x2, x3, [x21, #16]
ldp x16, x17, [x22, #16]
add x2, x2, x16
add x3, x3, x17
ldp x4, x5, [x21, #32]
ldp x16, x17, [x22, #32]
add x4, x4, x16
add x5, x5, x17
ldp x6, x7, [x21, #48]
ldp x16, x17, [x22, #48]
add x6, x6, x16
add x7, x7, x17
ldp x8, x9, [x21, #64]
ldp x16, x17, [x22, #64]
add x8, x8, x16
add x9, x9, x17
ldp x10, x11, [x21, #80]
ldp x16, x17, [x22, #80]
add x10, x10, x16
add x11, x11, x17
ldp x12, x13, [x21, #96]
ldp x16, x17, [x22, #96]
add x12, x12, x16
add x13, x13, x17
ldp x14, x15, [x21, #112]
ldp x16, x17, [x22, #112]
add x14, x14, x16
add x15, x15, x17
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
orr x16, x0, x1
orr x28, x2, x3
orr x16, x16, x28
orr x17, x4, x5
orr x28, x6, x7
orr x17, x17, x28
orr x19, x8, x9
orr x28, x10, x11
orr x19, x19, x28
orr x20, x12, x13
orr x28, x14, x15
orr x20, x20, x28
orr x28, x16, x17
orr x28, x28, x19
orr x28, x28, x20
tst x28, x25
b.ne .Lhit128_g
stp xzr, xzr, [x23], #16
.Lsg_tail:
subs x26, x26, #256
b.ne .Lsg_loop256
```
```asm
; kernel: rev_scan_g_cold_128 items=1 charged=63.0 key=c_scan_g_cold_static
.Lhit128_g:
stp x21, x22, [x27, #16]
stp x24, x26, [x27, #32]
str x23, [x27, #48]
ldp x21, x22, [x27]
mov x24, #0x8000000000000000
mov x26, #0
tst x17, x25
b.eq .Lhg_q0
stp x7, x6, [x27, #64]
stp x5, x4, [x27, #80]
mov x23, #0
b .Lhg_qloop
.Lhg_q0:
tst x16, x25
b.eq .Lhg_mid
stp x3, x2, [x27, #64]
stp x1, x0, [x27, #80]
mov x23, #32
b .Lhg_qloop
.Lhg_mid:
mov x0, x26
mov x26, #0
tst x20, x25
b.eq .Lhg_q2
stp x15, x14, [x27, #64]
stp x13, x12, [x27, #80]
mov x23, #64
b .Lhg_qloop
.Lhg_q2:
tst x19, x25
b.eq .Lhg_done
stp x11, x10, [x27, #64]
stp x9, x8, [x27, #80]
mov x23, #96
.Lhg_qloop:
add x4, x27, #64
.Lhg_wlp:
ldr x5, [x4], #8
and x6, x5, x25
cbz x6, .Lhg_wnxt
.Lhg_blp:
clz x7, x6
lsl x28, x5, x7
lsr x28, x28, #56
sub w28, w28, w22
strb w28, [x21], #1
lsr x28, x24, x7
bic x6, x6, x28
lsr x7, x7, #3
and x28, x23, #63
add x7, x7, x28
lsr x28, x24, x7
orr x26, x26, x28
cbnz x6, .Lhg_blp
.Lhg_wnxt:
add x23, x23, #8
tst x23, #31
b.ne .Lhg_wlp
cmp x23, #64
b.lo .Lhg_q0
b.eq .Lhg_mid
cmp x23, #96
b.eq .Lhg_q2
.Lhg_done:
str x21, [x27]
ldr x23, [x27, #48]
stp x0, x26, [x23], #16
ldp x21, x22, [x27, #16]
ldp x24, x26, [x27, #32]
tbnz x23, #4, .Lsg_mid
b .Lsg_tail
```
```asm
; kernel: rev_scan_f_u8_512 items=512 charged=0.155 key=c_scan_f_raw
.Lsf_loop512:
ldp x0, x1, [x23], #16
orr x0, x0, x1
cbnz x0, .Lhit128_f
.Lsf_ref0:
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
ldp x0, x1, [x23], #16
orr x0, x0, x1
cbnz x0, .Lhit128_f
.Lsf_ref1:
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
ldp x0, x1, [x23], #16
orr x0, x0, x1
cbnz x0, .Lhit128_f
.Lsf_ref2:
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
ldp x0, x1, [x23], #16
orr x0, x0, x1
cbnz x0, .Lhit128_f
.Lsf_ref3:
stp x24, x24, [x21, #16]
stp x24, x24, [x21, #32]
stp x24, x24, [x21, #48]
stp x24, x24, [x21, #64]
stp x24, x24, [x21, #80]
stp x24, x24, [x21, #96]
stp x24, x24, [x21, #112]
stp x24, x24, [x21], #128
stp xzr, xzr, [x22, #16]
stp xzr, xzr, [x22, #32]
stp xzr, xzr, [x22, #48]
stp xzr, xzr, [x22, #64]
stp xzr, xzr, [x22, #80]
stp xzr, xzr, [x22, #96]
stp xzr, xzr, [x22, #112]
stp xzr, xzr, [x22], #128
subs x26, x26, #512
b.ne .Lsf_loop512
```
```asm
; kernel: rev_scan_f_cold_128 items=1 charged=37.0 key=c_scan_f_cold_static
.Lhit128_f:
ldp x0, x1, [x23, #-16]
sub x7, x21, x25
cbz x0, .Lhf_hi
.Lhf_lo_lp:
clz x2, x0
sub x3, x14, x2
ldrb w4, [x19], #1
ldrb w5, [x21, x3]
ldrb w6, [x22, x3]
add w5, w5, w6
add w5, w5, w4
cmp w5, w27
b.lo .Lhf_lo_clr
add x4, x7, x3
str x4, [x20], #8
.Lhf_lo_clr:
lsr x4, x28, x2
bic x0, x0, x4
cbnz x0, .Lhf_lo_lp
.Lhf_hi:
cbz x1, .Lhf_ret
.Lhf_hi_lp:
clz x2, x1
sub x3, x15, x2
ldrb w4, [x19], #1
ldrb w5, [x21, x3]
ldrb w6, [x22, x3]
add w5, w5, w6
add w5, w5, w4
cmp w5, w27
b.lo .Lhf_hi_clr
add x4, x7, x3
str x4, [x20], #8
.Lhf_hi_clr:
lsr x4, x28, x2
bic x1, x1, x4
cbnz x1, .Lhf_hi_lp
.Lhf_ret:
tbnz x23, #5, .Lhf_r23
tbnz x23, #4, .Lsf_ref0
b .Lsf_ref3
.Lhf_r23:
tbnz x23, #4, .Lsf_ref2
b .Lsf_ref1
```
```asm
; kernel: rev_pad_zero_u8_128 items=128 charged=0.080 key=c_pad_zero_raw
.Lpad128:
stp xzr, xzr, [x_b2, #16]
stp xzr, xzr, [x_b2, #32]
stp xzr, xzr, [x_b2, #48]
stp xzr, xzr, [x_b2, #64]
stp xzr, xzr, [x_b2, #80]
stp xzr, xzr, [x_b2, #96]
stp xzr, xzr, [x_b2, #112]
stp xzr, xzr, [x_b2], #128
subs x_cnt, x_cnt, #128
b.ne .Lpad128
```

*(Section 1.7 table-build kernels and Sections 2–5 Montgomery CIOS, NTT/Lehmer, `Z[zeta_12]` Block Wiedemann SpMV, and Descent kernels are in `KERNELS_LA.md`; `audit_cm.py` concatenates `KERNELS.md` + `KERNELS_LA.md` and verifies all 74 ASM blocks.)*
