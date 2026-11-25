# Kernels (Raw-Opcode-Verified Pure-Scalar ARMv8.1-A with Disclosed `FEAT_LSE STADDB`/`STADD`: `Arch3-K2-Lq75-AdaptTg-M15-SplitRot-W1W8Dir32-Bit58-U8Fused-Shoup61-Hex6`)

All hot-loop kernels use pure-scalar general-purpose-register ARMv8.1-A integer instructions (`-march=armv8.1-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `DC ZVA`, zero `fp64`, zero `x18`, zero `x29`, zero `bl`/`blr` — all kernels are inline macros or leaf routines returning via `ret`), plus **two disclosed scalar `FEAT_LSE` atomic integer instructions** from ARMv8.1-A, each priced at `1` scalar instruction = `1 / 20,000` units under `COST_MODEL.md`:
- **`STADDB Ws, [Xn]`** (32-bit byte alias of `LDADDB Ws, WZR, [Xn]`; scalar registers; base-register `[Xn]` addressing only), used for `uint8_t` score updates on fused `S_b1` (`0.25 EiB`) and `S_b2[0..31]` (`8.00 EiB` per slab owner across `M = 15` collision-free chord slab owners) with zero wraparound (deterministically proved `<= 246 <= 255` on `g` at `3.50`b `round` with per-column odd-norm `B_{g,odd,max}(b, c)`-adaptive `T_g_u8(b, c)` and `<= 238 <= 255` on `f` at `2.75`b `round`, verified live in `audit_cm.py`).
- **`STADD Xs, [Xn]`** (64-bit alias of `LDADD Xs, XZR, [Xn]`; `[Xn]` addressing only), used exclusively for in-place `uint64_t pos[1..B]` byte-count accumulation during dense table build across ALL primes `2^{33} <= p < 2^{61.3}` (`8 * N_{a,k} < 2^{61} << 2^{64}` -> overflow mathematically impossible!). All other `FEAT_LSE` opcodes (`LDADD*`, `STADDH`, `SWP*`, `CAS*`) remain forbidden.

## 1. Collection Kernels (`K = 2` Equal-`V_a` Radial Shells `bp = [1, 1]` + Complex-Place `L_{3/4}(b)` Ball + `M = 15` Full-Chord Slab Owners `480` Cols + Split-`p_rot` / 3-Tier Band-2 `W_B=1`/`W_B=8`/`Direct1-W32` + Fused `u8` 128-Cell Scan)

Each shell `k in {0, 1}` owns `V_a,k = 0.9460 * 2^56 < 0.95 * 2^56` (`N_a,k = 4 * V_a,k = 2^57.920`, `stride = 56`, `bp = 1` = all 4 Galois banks inside ONE slab owner). Under the Complex-Place `L_{3/4}(b)` quasi-ball `K_{3/4}(R_{q,k}) = {b in O_H : sum_{m=1}^4 |sigma_m(b)|^{3/4} <= R_{q,k}^{3/4}}` (`R_{q,0} = 1062.99`, `L_max = 2^{-4/3} R_{q,0} = 421.85` proved unconditionally by subadditivity `4 |t|^{3/4} <= 2 R_{q,k}^{3/4}`), `max_m |sigma_m(b)| <= R_{q,0}` and the exact Vandermonde inverse `C_1 = 0.652889` certifies all power-basis limbs satisfy `|b^3_j| < 2^29.55`, `|g_j| < 2^93.26 << 2^127`, and `|f_j| < 2^90.08 << 2^127` (zero 128-bit norm-evaluation overflow, verified live in `unitdomain.py`). To hold an ENTIRE chord of length `<= 422` (plus up to `31` un-trimmed `g`-side overhang columns: `422 + 31 = 453 <= 480`) with **ZERO cross-column ring collisions**, `M = 15` slab owners `s in {0, 1, ..., 14}` each own `32` consecutive chord columns `[32*s .. 32*s + 31]` (`15 * 32 = 480` columns total). Each slab owner's `9.250 EiB < 16.000 EiB` virtual address space partitions cleanly on guard bit `58`: `S_b2[0..31]` (`8.00 EiB`) occupies `bit 58 == 0` at `(col_local << 59) | (0 << 58) | (m << 56) | oid` (`col_local in 0..31`); `S_b1` (`uint8_t`, `0.25 EiB`) occupies `(0 << 59) | (1 << 58) | (m << 56) | oid`; and `M_g[0..31]` (`1.00 EiB`) occupies `bits [63:58] in {0b000011, 0b000101, 0b000111, 0b001001}` (100% disjoint from `S_b2` and `S_b1`). Across all `15` slab owners, total active slab RAM is `138.75 EiB << 2^{120} B`. End-to-end pass ordering: all Band-2 `g`-primes emit the full chord into unique `(slab_id = col >> 5, col_local = col & 31)` slots via `512 KiB`-aligned TX buffers `tx_tbl[-1..15]` (`tx_tbl[-1] = sink_lo`, `tx_tbl[15] = sink_hi` discard out-of-chord overhang without touching any slab); only AFTER all Band-2 `g`-primes finish the chord does each slab owner `s` sweep its `32` columns applying Band-1a `g`-primes into `S_b1`, running `rev_scan_g_u8_fused_128`, setting `M_g[d]`, and clearing `S_b2[d]` (and clearing `<= 31` `g`-side padding columns `[L_chord .. L_chord + 30]` via `rev_pad_zero_u8_128`). Pass 2 (`f`) repeats symmetrically with ZERO padding columns (since `p_f = 60.4 == p_w1`, every `f`-side Band-2 prime uses exact `W_B = 1` chord slicing). Shells execute sequentially, keeping peak global table memory `2^{118.656} B < 2^{120} B` and peak 1-owner build RAM `14.395 EiB < 16.000 EiB`.

### 1.1a `rev_b2_row_k0_step4` (Band-2 4-prime orbit incremental `D_(omega, sigma)` row-step `k_0` init: `33/4 = 8.25` Ir/prime/row + `<= 0.606` outer-chord reset <= `C_B2_ROW_K0 = 8.90` Ir/prime/row)
Because `sigma(omega) = omega` (`7*5 === 4*5 === 13*5 === 5 mod 15`), all `4` primes `s_k = s_0^(7^k)` (`k = 0..3`) in a Band-2 prime orbit share the EXACT same modulus `p` (`2^53.0 <= p < 2^64.0`), the EXACT same `b_0`-coefficient `pwk[0] = 1`, and the EXACT same `b_1`-coefficient `pwk[1] = w = s_0^5 mod p`. Moreover, the `D_(omega, sigma)` fundamental-domain predicate `in_D12(b)` depends ONLY on `(b_2..b_7)`, and on every `(b_1..b_7)` fiber Band-2 sieves the **single contiguous bounding span `[b_{0,min}, b_{0,max}]`** (whose length is unconditionally `<= floor(2^{-4/3} R_{q,k}) + 1 <= 422 < 449` even on the `<= 0.08%` non-convex 2-component fibers such as `b = t + 85*sqrt(5)` with span `399 <= 422`, verified live in `unitdomain.py`!). Thus every fiber executes as ONE uninterrupted contiguous integer chord (`ZERO` multi-component resets in `rev_b2_row_k0_step4`; any interior gap column with `||b||_{3/4} > R_{q,k}` is simply cleared via `rev_pad_zero_u8_128` without running Band-1a or cofactorization, with total gap volume `< 0.05% << 0.5%` covered by `H2_ALLOW = 1.005`, live `3-SE UB <= 0.002` asserted in `unitdomain.py`). Maintaining `k_tail[0..3]` (`32 B` in `[x_st]`) across consecutive `b_1 -> b_1 + 1` rows and converting signed `x_b0` (`|x_b0| <= L_max = 421.85 << 2^53.0 <= p`, precondition `-p < x_b0 < p`) ONCE per 4-prime orbit computes overflow-proof `0 <= x10..x13 < x_p` (`k_0[0..3]`) in `33/4 = 8.25` Ir/prime/row (disjoint GPR map `x_st=x19, x_p=x20, x_d1=x21, x_pms=x22, x_b0=x23`, zero labels, verified bit-exact at `p = 2^64 - 59` and on `Q(zeta_15)` 4-prime orbits in `cm_models.py`). Outer `(b_2..b_7)` chord resets (`<= 25` ops/prime per chord: shared `(b1_lo - 1)*w mod p` + `4` modular adds from `kt_lvl[2..7][0..3]`, `416 B` total per owner) are conservatively bounded at `64 / min(L_b_shells) = 64 / 105.64 = 0.606` Ir/prime/row -> `8.25 + 0.606 = 8.856 <= 8.90 = C_B2_ROW_K0` Ir/prime/row:
```asm
; kernel: rev_b2_row_k0_step4 items=4 charged=8.90 key=c_b2_row_k0
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
Each Band-1a prime `[2^33, 2^53.0)` owns its own dedicated pre-based flat table with `(1 << 58) | (((m - k_c) & 3) << 56) | oid` pre-baked across all 4 Galois conjugates `k_c in {0..3}` (ZERO runtime rotation ops in `rev_b1a_apply_f_u64`, ONE contiguous slice per residue `k_t`, exact on all 4/4 conjugates!). Pre-adding `flat_buf_base` into `off_ptr[k] = flat_buf_base + 8 * prefix_sum[k]` at table-build time (ZERO extra RAM, ZERO extra build ops) eliminates per-bank base adjustments and `lsr #3` shifts (all byte counts `x_cnt = x_end - x_p = 8 * k` remain exact multiples of 8).
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
### 1.2 Shared `uint8_t S_b1` Band-1a (`[2^33, 2^53.0)`): `64x` + Straight-Line `32-16-8-4-2-1` Binary Tail `STADDB` (`2.5625 <= 2.57` Ir/hit, Exact Worst-Case `b1a_long_lk = 31.44 <= 31.5`)

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

### 1.4 `M = 15` Full-Chord Slab-Routed Band-2 (Pre-Rotated `W_B = 1` on `[2^53.0, 2^56.5)`, Shared `W_B = 1` on `[2^56.5, 2^60.4)`, `W_B = 8` on `[2^60.4, 2^61.3)`, Direct 1-Slot `W_SB = 32` on `[2^61.3, p_{g,k})`)

All Band-2 records store `rec = ((k & 31) << 59) | (0 << 58) | (m << 56) | oid` via **`ubfiz x3, x0, #59, #5`** (exact 5-bit residue `r = k & 31`, guard bit `58 == 0`). Each of the `M = 15` slab owners `s in {0..14}` owns `32` contiguous chord columns `[32*s .. 32*s + 31]` (`480` columns total `> L_max + 31 = 452.85`, **ZERO ring collisions**). **Weight-Banded TX-Buffer Aggregation:** Because `uint8_t` sieve weights are quantized at `3.50`b on `g` and `2.75`b on `f`, ALL Band-2 primes span ONLY **`4` constant-weight bands on `g`** (`w_lp in {15, 16, 17, 18}` on `[53.0, 63.6)` — covering the ENTIRE Direct-1 tier `[61.3, 63.6)` with ONE constant `w_lp = 18`!) and ONLY **`4` constant-weight bands on `f`** (`w_lp in {19, 20, 21, 22}` on `[53.0, 60.4)`). Within each weight band, ALL primes share the SAME constant `w_lp` and append directly into the SAME 15 `512 KiB`-aligned TX buffers `tx_tbl[0..14]` (preserving exact weights with ZERO per-record weight field; `tx_tbl[-1] = sink_lo`, `tx_tbl[15] = sink_hi` discard overhang). Whenever any `tx_tbl[s]` reaches `128 KiB = 16,384` entries (`tbnz #17`, leaving `384 KiB` headroom), it executes the dedicated inline drain **`rev_b2_rx_drain_pre64 -> rev_b2_rx_drain_u64 -> rev_b2_rx_drain_bin32`** on `(x_rx = r15, x_fcnt = r14, w_lp = r28)` using ONLY dead scratch **`{x0, x4, x5, x6, x7, x9, x16, x17}`** — programmatically asserted in `cm_models.py` to be **100% disjoint from ALL live Band-2 caller registers `{x1, x2, x3, x8, x10, x11, x12, x13, x19..x28}`** (`98/64 + 109/16384 = 1.5379 <= 1.54`, `<= 1.53125 + 0.015`). At the end of each of the `<= 4` weight bands per side (BEFORE switching `w_lp` and BEFORE scanning the chord), each slab executes ONE final partial drain (`<= 4 * 15 = 60` partial drains/side/chord `<= 60 * 109 = 6,540` ops/side/chord, strictly below the `20,700` ops/side/chord end-of-band drain allowance charged in `audit_cm.py`!).

1. **Tier 2a (`W_B = 1`: Pre-Rotated `[2^53.0, 2^56.5)` + Shared `[2^56.5, 2^60.4)` on BOTH `g` and `f`):** Exact per-residue `off[0..p]` (built in-place via `pos_buf[0] = flat_base` with ZERO copy). On the split pre-rotated sub-tier `[2^53.0, 2^56.5)` (`p_rot = 56.5`), each of the 4 Galois conjugates stores its own pre-rotated table (`m_rot = (m - k_c) & 3` baked at build time, `<= 4.504 EiB < 16 EiB`), so ALL 4/4 conjugates call the existing identity emitter **`rev_b2_emit_id_u16`** (`35/16 = 2.1875` + `1.53125` + `0.015` = **`3.73375 <= 3.735 = c_hit_hi_rot`** Ir/hit, ZERO new ASM!). On `[2^56.5, 2^60.4)`, one shared table per 4-prime orbit serves `1/4` identity + `3/4` rotated conjugates (`(2.1875 + 3*3.1875)/4 + 1.53125 + 0.015 = 4.48375 <= 4.485 = c_hit_hi_u4`). `rev_w1_probe_fast` (`12` ops) sets `x_oc = &off[k_0]` and `x_oe = &off[k_end]` (`k_end <= p`), and `rev_b2_w1_wrap` (`11` ops static; `5` ops dynamic on Range-1 `col0 == 0`) computes the EXACT full-slab count `x21 = (x_oe - x_oc) >> 8` (`floor(L_rem / 32)` -> ZERO over-read beyond `off[k_end]`, ZERO extra slab!). If `col0 & 31 != 0` (Range-2 wraparound mod `p`), `rev_b2_w1_prefix_emit` emits the clipped first partial slab `[x_oc, min(x_oc + 8*(32-(col0&31)), x_oe))` before `.Lw1_aligned`. `rev_b2_w1_slab_hdr` (`11` ops) processes strictly the `x21` full 32-column slabs, and `rev_b2_w1_tail_slab` (`9` ops static; `2` ops exit when `x_oc == x_oe`) emits any clipped final partial slab `[x_oc, x_oe)` (Range-1 chord overhead `12 + 5 + 2 = 19 <= 20.0 = c_tile_hdr_w1`; exact `cand_bd = 1.000`, ZERO overhang, `0` padding cols on `f`-pass since `p_f = 60.4 == p_w1`!). Each slice executes `rev_b2_emit_id_pre16`/`_c123_pre16` (`2` ops) -> `rev_b2_emit_id_u16`/`_c123_u16` (`35/16 = 2.1875` / `51/16 = 3.1875`, including in-loop `tbnz x10, #17` every `16` records = `128 B`!) -> `rev_b2_emit_id_bin8421`/`_c123_bin8421` (`2`-op `.Lb2i_fix16`/`.Lb2c_fix16` + straight-line `8-4-2-1` binary tail, exact worst-case excess `<= +8.8125 <= +8.875` ops across ALL `n >= 0` -> `11 + 8.875 = 19.875 <= 19.9 = c_sb_hdr_u4`).
2. **Tier 2b (`W_B = 8`, `[2^60.4, 2^61.3)` on `g`-side ONLY):** Standard `off[0..B]` (`B = ceil(p/8)`, ZERO sentinel padding needed!). `rev_w8_probe_fast` (`13` ops) + `rev_b2_w8_wrap` (`6` ops dynamic on Range-1 `col0 == 0, delta != 0` -> `19 <= 20.0 = c_tile_hdr_w8`) compute `x_off = &off[q_0]` (`q_0 = k_0 >> 3` in `x21`), `x_oe = &off[q_end]` (`q_end = (k_end + 7) >> 3 <= B`), and the EXACT full 4-bucket slab count `x21 = (x_oe - x_off) >> 5` (`cbz x21, .Lw8_tail_slab` -> ZERO over-read beyond `off[q_end] <= off[B]`, even when `k_0 = p - 1`!). Range-2 wraparound prefixes (`col0 > 0`) and pure slab-aligned chords (`delta == 0`, zero split!) are handled by `rev_b2_w8_prefix_pure` (all suffix-live registers `{r11, r19, r21, r22, r23, r25, r26, r28}` asserted disjoint from BOTH `rev_b2_rx_drain_*` and `rev_b2_emit_*`!). `rev_b2_w8_slab_hdr` (`16` ops) loads `ldp x11, x10, [x_tx]` (**`x11 = tx[s-1]` lower slab, `x10 = tx[s]` upper slab!**), flushes `x10`/`x11` via register-disjoint `rev_b2_rx_drain` (preserving `x2, x3, x8, x10, x11`), splits the single boundary bucket `[off[q], off[q+1])` (`1/4` of entries) via `rev_b2_emit_split_id_pre4/u4/bin21` (`_c123_*`, including in-loop `tbnz x10, #17` + `tbnz x11, #17` every `4` records = `32 B` -> `22/4 = 5.50` / `26/4 = 6.50`!) routing `C = 0` to `x11 = tx[s-1]` and `C = 1` to `x10 = tx[s]` while advancing `x2` to `x3` and preserving `(x3, x8)`, emits the `3` interior buckets `[off[q+1], off[q+4])` (`3/4` of entries, `x2 == x3`) via the **LITERAL** `rev_b2_emit_id_*`/`_c123_*` directly into **`x10 = tx[s]`** (leaving `x11 = tx[s-1]` 100% untouched, asserted live in `cm_models.py`!), and writes back `stp x11, x10, [x_tx]`. `rev_b2_w8_tail_slab` (`15` ops static; `2` ops exit when `x_off == x_oe`) emits any clipped final `< 4`-bucket slab `[x_off, x_oe)` using `ldp x2, x3, [x_off]`, `ldr x8, [x_oe]`, and `ldp/stp x11, x10, [x_tx]` (at `k_0 = p - 1`, `x3 == x8 == off[B]` skips the interior emitter with ZERO out-of-bounds read!). Combined rate: `0.75 * 4.48375 + 0.25 * 7.79625 = 5.311875 <= 5.313 = c_hit_hi_u2`; `16 + 8.875 + 6.0 = 30.875 <= 30.9 = c_sb_hdr_u2`; exact `cand_bd = 1 + 8 / Lb_k`. **Deterministic Arbitrary-Length TX-Buffer Safety Invariant (`W1` + `W8` + `Direct1`):** Every `512 KiB`-aligned buffer `tx[s]` starts at `occ = 0`. Every `W1`/`W8` header tests `tbnz #17` before emitting (entering with `occ <= 131,064 B`). Inside `split_u4`, `tbnz x10/x11, #17` fires every `<= 32 B`, so after `split_bin21` (`<= 24 B`), `occ <= 131,088 B`. Inside `emit_u16`, `tbnz x10, #17` fires every `128 B` (peak `<= 131,088 + 128 = 131,216 B = 0x20090`), and after `bin8421` (`<= 120 B`), `occ <= 131,208 B = 0x20088`. On handoff into `Direct1`, every single write (`+8 B`, peak `<= 131,216 B`) immediately tests `tbnz #17`. Thus for ANY arbitrary slice lengths `0 <= n < 2^{64}` (zero tail assumption!), `occ(tx[s]) <= 131,216 B << 262,144 B` (`256 KiB`) strictly, so bit 17 (`0x20000`) NEVER wraps to `0` and buffer overflow (`512 KiB`) is unreachable (verified adversarially at `n = 20,000` in `w8_verifiers.py`).
3. **Tier 2c (`Direct 1-Slot Sparse u2`, `W_SB = 32`, `[2^61.3, p_{g,k})` on `g`-side ONLY):** `rev_sb_probe_fast` (`13` ops) + `rev_b2_sb_wrap` (`5` ops, including `subs x21, x21, #1; b.lo .Lsd1_tail1` so any `N_bkt == 1` segment — including a range confined to the last bucket before a modulus boundary `p` — branches directly to `.Lsd1_tail1` with ZERO out-of-bounds read!) = `18.0 = c_tile_hdr`. `rev_b2_sieve_direct1_u2` skips empty pairs in `7/2 = 3.5` ops/bucket (`add x_tx, x_tx, #16`) and falls straight through to `.Lsd1_tail1 -> .Lsd1_done: ret`. `rev_b2_sieve_direct1_hit2` (out-of-line past `ret`) uses `cset x6, hs` (bucket 0 -> `tx[s]/tx[s+1]`) and `csinc x6, x_one, x_one, lo` (bucket 1 -> `tx[s+1]/tx[s+2]`) with per-singleton & per-node `tbnz #17` flush guards calling register-disjoint `rev_b2_rx_drain`: deterministic worst-case across ALL 4 pair states & ALL chain lengths `m >= 1` (including `b .Lsd1_next2` on `(0, 1)` = `13` ops and per-element `tbnz #17` on EVERY singleton, loop node, AND terminal overflow record!) is `<= 13.0 * n_ne + 11.0 * m_ovf` (ZERO pair-independence assumption, ZERO unguarded appends!) + `1.54 * lam` receive drain + `<= 18.0` odd-bucket chord tail (`rev_b2_sieve_direct1_tail1`).

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
**Direct 1-Slot Sparse Band-2 (`p >= 2^{61.3}`, `W_SB = 32`):** Every prime inserts `N_a,k = 0.9460 * 2^{58} < 2^{58}` records total, so every overflow node has a unique 58-bit `node_idx in [0, 2^{58})` (`bits [63:58] == 0`, supporting up to `4.000 EiB` of contiguous 16-byte overflow nodes without ever colliding with bit 58!). Overflow heads store `tag = (1 << 58) | node_idx`; with `x_ovf_adj = (x_ovf_base - 2^{62}) mod 2^{64}` precomputed ONCE per prime, `add x4, x_ovf_adj, x1, lsl #4` decodes `x_ovf_base + 16 * node_idx` bit-exact in `1` instruction across all `0 .. 4.000 EiB` (asserted live across bit-58/59/60/61 boundaries in `w8_verifiers.py`). `rev_b2_sb_wrap` computes `x21 = N_bkt - 2` via `subs x21, x21, #1; b.lo .Lsd1_tail1`, jumping directly to `.Lsd1_tail1` whenever `N_bkt == 1` (preventing any out-of-bounds pair load past the last bucket before a modulus boundary `p`), while `N_bkt >= 2` enters `rev_b2_sieve_direct1_u2` (`7/2 = 3.5` ops/empty bucket). On any non-empty pair, `rev_b2_sieve_direct1_hit2` (placed out-of-line past `.Lsd1_done: ret` so `b.hs .Lsd1_u2` falls directly into `.Lsd1_tail1`) emits singletons via `cset x6, hs` (bucket 0 -> `tx[s]/tx[s+1]`) and `csinc x6, x_one, x_one, lo` (bucket 1 -> `tx[s+1]/tx[s+2]`) with per-singleton `tbnz x7, #17` flush (`12` ops for `(1,0)`, `13` ops for `(0,1)` incl. `b .Lsd1_next2`, `21/2 = 10.5` ops/bucket for `(1,1)` -> deterministic `<= 13.0` ops/non-empty bucket with zero pair-independence assumption), follows overflow chains in-register (`ldp x10, x11, [x_tx]` / `[x_tx, #8]`, `9` ops/node worst-case with per-node AND terminal-record `tbnz #17` flush -> `15 + 9*m <= 13.0 + 11.0*m` for all `m >= 1`, ZERO unguarded appends!), and flushes via dedicated register-disjoint `rev_b2_rx_drain_*` (`(x_rx = r15, x_fcnt = r14, w_lp = r28)`, dead scratch `{x0, x4, x5, x6, x7, x9, x16, x17}`, `1.53125 + 109/16384 = 1.5379 <= 1.54` ops/entry); `rev_b2_sieve_direct1_tail1` handles `N_bkt == 1` and odd `N_bkt >= 3` (`+18.0` ops/chord). Total deterministic per-bucket price: `c_sb_d = 3.5 + 13.0*p_ne + 11.0*ovf_m + 1.54*lam`.

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

### 1.5 Fused `uint8_t S_b1 + uint8_t S_b2` 128-Cell SWAR Scan (`68/128 = 0.53125` + `21/128 = 0.16406` + `pad/cold <= c_scan_zero = 0.799`)

On each of the `M = 15` slab owners, `S_b1` (`0.25 EiB`), `S_b2[0..31]` (`8.00 EiB`), and `M_g[0..31]` (`1.00 EiB`) total `9.250 EiB < 16.000 EiB`. Using round-to-nearest weights (`g`: `3.50`b, per-column odd-norm `B_{g,odd,max}(b, c)`-adaptive `T_g_u8(b, c) in [60, 127]` with deterministic no-overflow floor, `bias_g(b, c) = 128 - T_g_u8(b, c)`; `f`: `2.75`b, `T_f_u8 = 49`, `bias_f = 0`), deterministic odd-norm bounds prove `S_b1[i] + S_b2[i] <= 246 <= 255` on `g` (asserted explicitly across quarter-bit `[520..700]` in `cm_models.py`) and `<= 238 <= 255` on `f` (ZERO `uint8_t` wraparound -> ZERO inter-byte carry in 64-bit SWAR `add`). Closed periodic `S_b1` state machine switches `x_b8 = bias_g(b_{r+1}, c) -> 0` on the last active column of Pass 1 (`g`) and reloads `x_f8 = bias_g(b^{(t+1)}_0, c)` (the **next active chord `t + 1`'s** first-column bias for slab `s`, rather than chord `t`'s own first column!) on the final active column of Pass 2 (`f`) via the existing `S_b1` refill store with zero extra instructions (verified live across 3 consecutive chords of lengths `(16, 11, 16)` with distinct per-chord and per-column thresholds `_tg_ch[0..2]` in `cm_models.py`). Exact no-hit path of `rev_scan_g_u8_fused_128`: `32` load/add + `16` clearing stores + `15` ORs + `5` test/branch/bitmap/loop = **`68` opcodes / 128 cells** (`0.53125 <= 0.540` charged; `g`-pass clears `<= 31` padding columns `[L_chord .. L_chord + 30]` via `rev_pad_zero_u8_128`, `f`-pass has ZERO overhang -> total scan+zero+cold `0.7964 <= 0.799 = c_scan_zero`).

```asm
; kernel: rev_scan_zero_wrap items=1 charged=2.0 key=c_scan_zero_wrap
cbz x_cnt, .Lscan_done
.Lscan_done:
ret
```
```asm
; kernel: rev_scan_g_u8_fused_128 items=128 charged=0.540 key=c_scan_g_raw
.Lscan128_g:
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
stp x_b8, x_b8, [x_b1, #16]
stp x_b8, x_b8, [x_b1, #32]
stp x_b8, x_b8, [x_b1, #48]
stp x_b8, x_b8, [x_b1, #64]
stp x_b8, x_b8, [x_b1, #80]
stp x_b8, x_b8, [x_b1, #96]
stp x_b8, x_b8, [x_b1, #112]
stp x_b8, x_b8, [x_b1], #128
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
tst x16, x_8080
b.ne .Lhit128_g
stp xzr, xzr, [x_mg], #16
subs x_cnt, x_cnt, #128
b.ne .Lscan128_g
.Lhit128_g:
```
```asm
; kernel: rev_scan_f_u8_128 items=128 charged=0.165 key=c_scan_f_raw
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
stp x_init_f, x_init_f, [x_b1], #128
stp xzr, xzr, [x_b2, #16]
stp xzr, xzr, [x_b2, #32]
stp xzr, xzr, [x_b2, #48]
stp xzr, xzr, [x_b2, #64]
stp xzr, xzr, [x_b2, #80]
stp xzr, xzr, [x_b2, #96]
stp xzr, xzr, [x_b2, #112]
stp xzr, xzr, [x_b2], #128
subs x_cnt, x_cnt, #128
b.ne .Lscan128_f
.Lhit128_f:
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

*(Section 1.7 table-build kernels and Sections 2–5 Montgomery CIOS, NTT/Lehmer, `Z[zeta_12]` Block Wiedemann SpMV, and Descent kernels are in `KERNELS_LA.md`; `audit_cm.py` concatenates `KERNELS.md` + `KERNELS_LA.md` and verifies all 72 ASM blocks.)*
