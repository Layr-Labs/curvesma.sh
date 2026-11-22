# Kernels (Raw-Opcode-Verified Baseline Scalar ARMv8.0-A: `Arch2-Prod-ZeroSlice-1PassB2` 62b Direct Apply, 1-Owner Build, 3-Prime Radix-4 Lazy Shoup NTT, 1-Dense + 16-Shard 4-Vec LA, 9-Op LUT `isqrt`, Input-Uniform Bareiss LLL)

Every kernel cost is derived from explicit baseline scalar ARMv8.0-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `fp64`, zero `bl`). `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file, validates opcodes & `LDP`/`STP` addressing, counts every instruction (`cmp`+`b.*`=2 ops, loop control included), simulates the complete 1-owner flat-array 2-pass build (incl. overflow-proof 64-bit modular step at `p = 2^64 - 59`) + short-bucket Band-1a/1b + 1-pass direct-apply 62b Band-2 (`W_SB = 16`, `oid < 2^56`, ZERO emit buffers) + 64-cell `ORR` scan + 4-vec 4-col `Z[zeta_12]` LA + `MONT5` + Barrett `row_red5_z12` + 3-prime Radix-4/3 Lazy Shoup NTT + 4-prime Garner CRT + branchless 9-op LUT `isqrt` iterative Q32 `dsieve16` + 10-op/level `dsieve16_rank_step` (`552/552` bit-exact) + certified upper-bound `BOUND_CELL_HI` live, and asserts `opcodes/items <= charged`.

## 1. Collection Kernels (`Arch2-Prod-ZeroSlice-1PassB2-U16`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 448 = 7 * 64`, `W_SB = 16`, `P_min = 2^32`, `P_m4s = 2^53.0`, `P_m4 = 2^55.5`, `P_lo = P_f = 2^58.0`, `P_g = 2^63.5`)

Because the sieving region is a pure product of two 8D balls (`qa(a) <= s_c R2_a` with 2-adic radii `s_c = (1.0, 1.5, 2.0)` and `qb(b) <= R2_b`, `A^2 = 320,000`, `r_ab = 0.83`), every canonical `a`-table entry (`N_a^can = V_a = 2^55.8474 < 0.90 * 2^56`, `N_a^tbl = 4 V_a = 2^57.8474`) is unconditionally in-region for every `b` (`V_b = 2^64.5265`, `L_max = 439.06 < W_tile = 448`, `L_b = 188.58`, `N_7 = 2^54.9675`, `V_kept = 2^120.3739`). There is **ZERO `q_key`, ZERO `q_tab`, ZERO `sb_min_q8`, ZERO binary search, ZERO global sorting invariant, ZERO Band-2 emit/flush buffer, and ZERO `f`-side Band 2** (`P_f = P_lo = 2^58.0`).
- **Bias-initialized `S_g` & zero-wrap `uint16_t` accumulators (`S_g, S_f`) with exact `2^56` bank stride (`K_a = 1`, `W_SB = 16`):** Quarter-bit weights on `<= 705`b norms reach `s <= 3,000`. Initializing `S_g` lanes to `bias16 = (32768 - T) & 0x7FFF` (`1 <= T <= 3,000`) keeps `(32768 - T) + s in [29768, 35767] subset [0, 65535]` with **ZERO 16-bit wraparound/carry between lanes** and `s >= T <==> bit 15 (0x8000) == 1`. Each column-owner slab allocates `uint16_t S_g[4][16][2^56]` (`2^63.0` B) for 1-pass `W_SB = 16` direct apply + 1-column `uint16_t S_f[4][1][2^56]` (`2^59.0` B) = **`2^63.087` B `< 2^64` B** (exact power-of-2 bank stride `2^56`, `oid < 2^56`).
- **Zero-copy 8B records with optimal short-bucket split across Band-1a/1b/2:**
  - Band-1a-long (`[2^32, 2^53.0)`): `m`-partitioned `oid < 2^56` -> 16x-unrolled `rev_b1a_apply_u16` (`58/16 = 3.625 <= 3.65` Ir/hit) + `rev_b1a_apply_wrap16` (`4` ops) + `rev_b1a_apply_1` (`6.0` Ir/rem hit, `<= +40` ops/bkt charged inside `b1a_long_lk = 44.0`).
  - Band-1a-short (`[2^53.0, 2^55.5)`): direct 1-record `rev_b1a_apply_1` (`6.0` Ir/hit, `b1a_short_lk = 8.0` ops/bkt, ZERO wrap/tail overhead!).
  - Band-1b (`[2^55.5, 2^58.0)`): direct 1-record 8B `(m << 56) | oid` `rev_b1b_apply_1` (`8.0` Ir/hit, `b1b_lk = 8.0` ops/bkt, ZERO wrap/tail overhead!).
  - Band-2 (`g`-side `[2^58.0, 2^63.5)` only — **1-Pass Direct Apply, ZERO Emit/Flush Buffers!**): 8B `rec = (m << 60) | (r << 56) | oid` (`2 + 4 + 56 = 62 <= 64` bits exact!). Outer loop state lives exclusively in `x19..x21` (disjoint from apply scratch `x0, x4, w8`, asserted live in `verify_kernels_symbolic`). `rev_b2_direct_sb_hdr` (`8 <= 8.0` Ir/SB) falls through on non-empty superblocks into 1-record `rev_b2_direct_1` (`8 <= 8.0` Ir/hit: `sub x4, x0, x_k60; and x4, x4, x_mask62` indexes directly into `S_g[4][16][2^56]` with ZERO emit/flush buffers, ZERO odd-length overread, and ZERO `c_bd_chk` branches; boundary superblock padding is priced by `1 + W_SB / L_b = 1.085` on both `c_hit_hi` and `c_scan_zero`).

### 1.1 `rev_sb_probe_fast` (Band-2 row start: `17 <= C_TILE_HDR = 18.0` Ir/row)

```asm
; kernel: rev_sb_probe_fast items=1 charged=18.0 key=c_tile_hdr
sub x0, x_k0, x_p
mul x1, x0, x_invp
umulh x2, x0, x_mu
lsr x2, x2, #4
msub x3, x2, x_p, x0
lsr x4, x3, #4
sub x5, x_end, x3
cmp x5, x_wmax
csel x5, x_wmax, x5, hi
lsl x6, x4, #3
add x7, x_off, x6
ldp x8, x9, [x7]
cmp x8, x9
b.hs .Lempty
add x_sb, x_sb, #1
add x_ptr, x_base, x8
ret
```

### 1.2 `rev_b2_direct_sb_hdr` (1-Pass Band-2 super-bucket header — ZERO emit buffers: `8 <= C_SB_HDR = 8.0` Ir/SB)

```asm
; kernel: rev_b2_direct_sb_hdr items=1 charged=8.0 key=c_sb_hdr
ldp x5, x6, [x19]
add x19, x19, #8
cmp x5, x6
b.eq .Lnext_sb
add x2, x20, x5
add x3, x20, x6
; --- inline rev_b2_direct_1 (.Lb2d) runs here for [x2, x3), preserving x19..x21 ---
.Lnext_sb:
subs x21, x21, #1
b.ne .Lsb
```

### 1.3 Band-1a-Long 16x-Unrolled + Band-1a-Short 1-Record Apply (`rev_b1a_apply_u16`: `58/16 = 3.625 <= C_HIT_LO_M4 = 3.65` Ir/hit; `rev_b1a_apply_1`: `6.0` Ir/hit)

Flag-safe pre/post wrap (`4` ops once per bucket) + 16-hit unrolled body (`58` ops) + 1-record loop (`6` ops). Worst-case overhead above `3.625` Ir/hit for ANY bucket length `k >= 0` is `<= 4 + 15 * (6 - 3.625) = 39.625 <= 40 <= 44.0` ops/bkt (charged in `b1a_long_lk = 44.0`; simulated flag-exact across ALL `k = 0..35` in `audit_cm.py` with `0` overreads!):

```asm
; kernel: rev_b1a_apply_wrap16 items=1 charged=4.0 key=c_hit_lo_m4_wrap
subs x_cnt, x_cnt, #16
b.lo .Ltail_fix16
.Ltail_fix16:
adds x_cnt, x_cnt, #16
b.eq .Ldone
```

```asm
; kernel: rev_b1a_apply_u16 items=16 charged=3.65 key=c_hit_lo_m4
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
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1a1
```

### 1.4 Band-1b Direct 1-Record 8B `{m:2, oid:56}` Apply (`rev_b1b_apply_1`: `8.0 <= C_HIT_LO_U = 8.0` Ir/hit, `b1b_lk = 8.0` ops/bkt)

Every Band-1b hit packs `(m << 56) | oid` into one 64-bit word, rotating bank `m -> (m - k) & 3` via `sub x4, x0, x_k56; and x4, x4, x_mask58` directly into the `2^56`-stride accumulator:

```asm
; kernel: rev_b1b_apply_1 items=1 charged=8.0 key=c_hit_lo_u
ldr x0, [x_p], #8
sub x4, x0, x_k56
and x4, x4, x_mask58
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1b1
```

### 1.5 Band-2 1-Pass Direct-Apply 62-Bit `(m<<60)|(r<<56)|oid` Kernel (`rev_b2_direct_1`: `8.0 <= C_HIT_HI = 8.0` Ir/hit — ZERO Emit/Flush Buffers, ZERO `c_bd_chk`!)

Entered directly from `rev_b2_direct_sb_hdr` on non-empty superblocks. Applies one 62-bit record per iteration (ZERO odd-length overread on any bucket size `0, 1, 2, ...`!), rotating bank `m -> (m - k) & 3` while preserving `(r, oid)` via `sub x4, x0, x_k60; and x4, x4, x_mask62` (`(x4 << 1) < 2^63 B`) directly into `uint16_t S_g[4][16][2^56]` (`2^63.0 B`):

```asm
; kernel: rev_b2_direct_1 items=1 charged=8.0 key=c_hit_hi
ldr x0, [x2], #8
sub x4, x0, x_k60
and x4, x4, x_mask62
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
cmp x2, x3
b.lo .Lb2d
```

### 1.6 Bias-Initialized 64-Cell `uint16_t` `ORR` Scan + Reset (`rev_scan_zero_u64`: `45/64 = 0.703125 <= 0.71` raw -> `0.703125 * (1 + 16/188.58) + 0.015 = 0.7778 <= C_SCAN_ZERO = 0.78`; `W_tile = 448 = 7 * 64` exact!)

Uses `x0..x15` for the 16 loaded `S_g` words (`64` `uint16_t` cells), `x16, x17, x19, x20` for the 15-op `orr` tree (preserving ALL 64 original `S_g` scores in `x0..x15` on entry to `.Lhit64`!), and `x_g, x_f, x_b16, x_8000, x_cnt` in `x21..x25` (`25 <= 29` GPRs, ZERO `x18/x29/x30`, signed-offset `ldp`/`stp` only). Since `W_tile = 448 = 7 * 64`, zero tail loop is needed:

```asm
; kernel: rev_scan_zero_u64 items=64 charged=0.71 key=c_scan_zero_raw
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
```

### 1.7 1-Owner Flat-Array Counting-Sort / Prefix-Sum / Scatter Build (`rev_build_prod_pass12`: `22 + 0.50 = 22.50 <= C_LIST_A = 24.0`, `rev_build_prod_prefix`: `6 <= C_LBKT = 6.0`)

Because `Arch2-Prod-ZeroSlice` requires **ZERO `q_key` sorting**, Pass 1 walks `a_0`-chords in deterministic 4-bank order (`m = 0..3`, `oid = 0..V_a - 1`), advancing `ka = (ka + step) mod p` via the 4-opcode overflow-proof modular step (`p = 2^64 - 59` tested LIVE!), writing `kbuf[idx] = ka`, and incrementing 64-bit `uint64_t x_hist[ka >> 4]`. Pass 2 executes 4 outer bank passes (`x_oidm = m << 60`, `x_n2 = V_a`) packing `x_oidm | ((ka & 15) << 56)` into `x_out[pos[sb]]`. **Total simultaneous 1-owner RAM `<= 2^63.898 B < 2^64 B` for ALL Band-1a/1b/2 primes (`0` cross-owner build transfers, asserted live in `audit_cm.py`)**:

```asm
; kernel: rev_build_prod_pass12 items=1 charged=23.5 key=c_list_a_raw
cmp x0, x_pms
sub x5, x0, x_pms
add x6, x0, x_step
csel x0, x5, x6, hs
str x0, [x_kbuf], #8
lsr x1, x0, #4
ldr x2, [x_hist, x1, lsl #3]
add x2, x2, #1
str x2, [x_hist, x1, lsl #3]
subs x_n1, x_n1, #1
b.ne .Lp1
ldr x0, [x_kbuf2], #8
and x5, x0, #15
orr x8, x_oidm, x5, lsl #56
add x_oidm, x_oidm, #1
lsr x1, x0, #4
ldr x2, [x_pos_sb, x1, lsl #3]
add x4, x2, #8
str x4, [x_pos_sb, x1, lsl #3]
str x8, [x_out, x2]
subs x_n2, x_n2, #1
b.ne .Lp2
```

```asm
; kernel: rev_build_prod_prefix items=1 charged=6.0 key=c_lbkt
ldr x0, [x_hist], #8
str x_acc, [x_off], #8
str x_acc, [x_pos], #8
add x_acc, x_acc, x0, lsl #3
subs x_cnt, x_cnt, #1
b.ne .Lpre
```

## 2. Exact 5-Limb Baseline Scalar ARMv8-A CIOS Montgomery Multiplication (`MONT5 = 270.0`)

One complete carry-exact, shift-included two-chain CIOS outer round executes **`48`** baseline scalar ARMv8-A opcodes (`5 x 48 = 240` across 5 unrolled rounds) + **`27`** load/zero-init/final-sub/store opcodes = **`267 <= 270.0` Ir** (`+3` Ir headroom; bit-exactness of this exact 267-opcode register machine is verified LIVE on 32 random 320-bit `(a, b, p)` triples inside `audit_cm.py`):

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
; kernel: mont5_cios_epilogue items=1 charged=28.0 key=mont5_epi
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
ret
```

Total exact `MONT5 = 5 * 48 + 27 = 267 <= 270.0` Ir.

## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Exact 2nd-Moment ECM Chebyshev

- **Explicit Schoolbook Multi-Limb Scalar ARMv8-A Tower Norm Derivation (`Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`, ZERO Karatsuba shortcut):**
  Explicit schoolbook scalar ARMv8-A `L_1 x L_2` limb multiplication uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)` (`mul/umulh/adds/adcs` carry chains). On the production sieving region `Q_CM(a) <= s_c * 265,600 <= 531,200` (all 2-adic classes `s_c in {1, 1.5, 2}`), `Q_CM(b) <= 385,542.2`, exact Cauchy-Schwarz via `max diag(G_8^{-1}) = 2/5` proves `|a_j| <= 461`, `|b_j| <= 393`, `|a^3_j| <= 2^{26.87} < 2^{27}` and `|b^3_j| <= 2^{26.17} < 2^{27}` (1 signed limb; `4` `8x8` 1-limb poly muls + `Phi_15` reduction = `400` Ir), `|g(a, b)_j| = |a^3 - u w^2 b^3|_j <= 2^{89.9} < 2^{90}` (2 signed limbs; `+96` Ir), and `|f(a, b)_j| <= 2^{82.1} < 2^{83}` (fits in 2 signed limbs, conservatively priced as 3 limbs: `a^6, b^6` = `2 x 928 = 1,856` Ir, `a^9, b^9, a^6 b^3, a^3 b^6` = `4 x 1,152 = 4,608` Ir, combine = `128` Ir). Using pure schoolbook `mul_LL` across the 3-step tower norm `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:
  1. **`_nk_g_only` (`2 -> 4 -> 8 -> 16` limbs):** `32` `2x2` (`1,536`) + `6` `4x4` (`960`) + `3` `8x8` (`1,728`) + adds/subs/setup/cold-scan (`1,544`) = **`5,768 <= 6,500` Ir charged**!
  2. **`_nk_f_only` (`3 -> 6 -> 12 -> 24` limbs):** reuses `a^3, b^3` -> `a^6, b^6, a^9, b^9, a^6 b^3, a^3 b^6` + combine (`6,592`) + 3-limb schoolbook tower norm (`4,608 + 4,032 + 3,744 = 12,384`) = **`18,976 <= 20,000` Ir charged**!
- **Parsed 62-Bit Shoup Radix-2 NTT Butterfly (`ntt62_butterfly`: `18` opcodes/2-elems/prime -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 fwd + 1 inv across 4 primes) & 4-Prime Garner CRT (`fft_F`):**
  `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`log2(NTT4_PROD) = 244.70`, all 4 verified prime LIVE). One exact baseline scalar ARMv8-A radix-2 Shoup butterfly (`w_s = floor(w * 2^64 / p)` precomputed, `t = w*v - umulh(v, w_s)*p` reduced to `[0, p)`) takes **`18` opcodes per 2 elements (`9.0` Ir/elem/prime -> `36.0` Ir/elem/transform across 4 primes -> `3 x 36.0 = 108.0 <= 168.0` Ir/elem/stage for 2 forward + 1 inverse NTTs in `fft_F`, `1.56x` headroom!)**; scalar `F_{p_i}` Barrett MAC in `PM-Basis` uses `mul/umulh` + 1-step Barrett `umulh/msub` + 2 conditional subs + add/sub/csel + 2 loads + loop = **`17 <= 24` opcodes** (verified bit-exactly LIVE on all 4 NTT primes in `kernels_symbolic`); 4-prime Garner CRT uses 3 Shoup mod-muls (`15` ops) + 3 multi-limb mul-adds (`<= 75` ops) = `<= 90 << 220` linear term in `fft_F` (`168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`; 2D `168 M_tot log2(M_tot) + 1100 M_tot + 1200`). Both `ntt62_butterfly` (64 random triples) and 4-prime Garner CRT (32 random 240-bit integers) are verified bit-exactly LIVE in `audit_cm.py`, alongside `srt_toy` (128 planted leaves, `0` mismatches):

```asm
; kernel: ntt62_butterfly items=2 charged=10.0 key=ntt62_elem_prime
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
  Three `< 2^62` primes divisible by `3 * 2^55`: `p_1 = 69 * 2^55 + 1`, `p_2 = 57 * 2^55 + 1`, `p_3 = 54 * 2^55 + 1` (`prod = 2^182.696 > 2^54 * (2^64 - 1)^2`, all 3 verified prime + 3-prime Garner CRT + bit-exact `n = 16` & `n = 24` radix-4/radix-3 Shoup convolutions using exact `subs`/`lo` ARMv8-A semantics LIVE in `cm_models.verify_kernels_symbolic`). `ntt3_F(M) = min(55.25*n2*k + 220*n2, 55.25*n3*(k-2) + 364*n3)` (`+ 512*n` if `n > 2^26`). Inputs live in `[0, 2p)` (`x_p2 = 2p < 2^63`, guaranteeing Shoup outputs `< 2p` and exact unsigned borrow via `subs` + `csel ..., lo`): `4` loads + `3` twiddle `ldp`s + `4` Lazy Shoup muls (`mul/umulh/msub`) + `8` lazy mod-`2p` add/subs/csel (`3` ops each) + `4` stores + `2` loop (`subs`/`b.ne`) = **`49` opcodes / `8` elem-stages = `6.125 <= 6.2` Ir/elem/stage/prime** (`9 * 6.125 = 55.125 <= 55.25`):

```asm
; kernel: ntt62_r4_shoup items=8 charged=6.2 key=ntt62_r4_elem_stage
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
- **Stage B (`P_min = 32`, `KH = 0.08`, Ordered Norms):** `gbatch_ntt3(pbits_split(32, 0.08), L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` with `MONT5 = 270.0`. Checking `S_f` (`sfa`) before `N_g` (`6,500` Ir) and checking `Stage-B(g)` (`sgb`) before `N_f` (`20,000` Ir) gives exact `c_cof_2d_upper = 0.9824` Ir/cell (live-verified across all `5,648` 2D bins of our 1,000,000-draw `r_ab = 0.83` pool with live `Q(zeta_15)` split-prime census up to `2^14`: `gb/tb = 0.07197 <= 0.08`, `fb/tb = 0.05649 <= 0.08`).
- **Exact 2nd-Moment ECM Chebyshev (`cm_models.ecm_mean_var`):** Schedule runs `K_c` curves (`p_1`, `q_s=(1-p_1)^{K_c}<=0.001`, `p_s=1-q_s`). Failed schedule cost `C_fail=K_c*per_c` (`Var=0`); final successful schedule has `S~Bin(K_c,p_1)|S>=1`, `C_succ=C_fail+S*redo_ex`, `E[S|S>=1]=K_c p_1/p_s`, `E[S^2|S>=1]=(K_c p_1(1-p_1)+(K_c p_1)^2)/p_s`. With independent `F~Geom_0(p_s)`, `T=F*C_fail+C_succ` gives exact `E[T]=(q_s/p_s)C_fail+E[C_succ]`, `Var(T)=(q_s/p_s^2)C_fail^2+Var(C_succ)`. Cauchy-Schwarz `sd_rel=13*sqrt(Var_f)+16*sqrt(Var_g)` yields Chebyshev `< 2^-74.4 << 2^-70`.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (1-Dense + 16-Shard `< 2^64 B` 4-Vec 4-Col Arch B Sign-Split Scalar ARMv8.0-A + Exact 1-Step Barrett `row_red5_z12`)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`, `tab_ratio_charged = 1.00` plain sparse with ZERO `TABM`, `32` sequential 4-vector bundles of `128 B/slot`, `n = 128 = 32 * 4`, `m = 16,384`, peak single-owner dense-vector RAM `2 * N_sw * 128 B = 2^63.668 B < 2^64 B`, `16` read-only sparse-matrix shard owners `2^62.337 B << 2^64 B` each (including `uint8_t` non-unit coefficients), explicit priced cross-owner matrix streaming `c_stream_mat = 204.0` Ir/row (`4.0` Ir/8B streamed word across 8B offsets + `<= 10%` 1B coefficients, amortized over `4` vectors -> `203.51 <= 204.0` Ir/row)):
- **`unit_run_4vec4_z12` & `unit_run_4vec4_neg_z12` (`122` opcodes / `16` `(entry, vec)` pairs = `7.625 <= 7.65` Ir/vec-entry, `28 <= 29` GPRs `x0..x17, x19..x28`, ZERO `x18/x29/x30`, `0` spills):** Each column holds 4 contiguous 256-bit vectors (`128 B/col`), addressed via buffer-independent pre-scaled `uint64_t` byte offset `col * 128 < 2^63` (`N_sw * 128 = 2^62.668 B < 2^64 B`) added explicitly to `x6 = x_base` (`ldr x7, [x4], #8; add x7, x6, x7` — works natively on both ping-pong buffers with ZERO vector copying!). Each loop iteration processes 4 columns across all 4 vectors (`[x7, #0/#16]`, `[x7, #32/#48]`, `[x7, #64/#80]`, `[x7, #96/#112]`) into 4 signed 5-limb accumulators (`v0`: `x8..x12`; `v1`: `x13..x17`; `v2`: `x19..x23`; `v3`: `x24..x28`), `4 * (2 + 4 * 7) + 2 = 122 / 16 = 7.625 <= 7.65` Ir/vec-entry:

```asm
; kernel: unit_run_4vec4_z12 items=16 charged=7.65 key=unit_run_4vec4
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

- **Sign-Split 2-Vec `gen_run1_2v_z12` (`c_b > 0`) & `gen_run1_2v_neg_z12` (`c_b < 0`) (`43/2 = 21.50 <= 22.0` Ir/vec-entry, `25 <= 29` GPRs avoiding `x18/x29/x30`, `0` spills):** Each bucket splits non-unit scalar updates (`1 <= |c_b| <= 255`, `< 6.5%` of entries, stored as 1-byte `uint8_t` magnitudes loaded and zero-extended to `x22` via `ldrb w22, [x23], #1`) into positive and negative sub-lists, executed as two spill-free 2-vector passes (`(v0, v1)` at offsets `#0..#48`, then `(v2, v3)` at offsets `#64..#112`, explicit `add x8, x9, x8` from `x_base = x9`): `3` ptr/coeff loads + `2 * 19` mul/umulh/carry ops + `2` loop = **`43 / 2 = 21.50 <= 22.0` Ir/vec-entry** (`k5_live = eu_lo * 7.65 + (1 - eu_lo) * 22.0 = 8.5721 <= k5_charged = 10.10` Ir, `+1.528` Ir headroom!):

```asm
; kernel: gen_run1_2v_z12 items=2 charged=22.0 key=gen_run1_2v_coord
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

- **`chi_mul_z12_fast` (`50` Ir), `row_tail_6bucket_z12` (`504` Ir = `12` signed unit sub-lists x max `4` tail entries x `8` ops + `2 * 50` `chi_mul` + `20` dispatch), 1-Dense + 16-Shard 4-Vec Arch B Row Overhead (`c_row_chi = 960.0` Ir, `c_stream_mat = 204.0` Ir/row, `c_xown_live = 0.0` Ir, `n_owners = 1`, `n_mat_owners = 16`, `te = 0`):** `6 * 67 + 504 + 16 = 922 <= 6 * 68 + 504 + 16 = 928 <= c_row_chi = 960.0` Ir worst-case; explicit matrix-shard streaming `c_stream_mat = 204.0` Ir/row (`4.0` Ir/8B streamed word across 8B offsets + `<= 10%` 1B non-unit coefficients, amortized over `4` vectors) charged inside `c_prod` (`tab_ratio_charged = 1.00`, ZERO `TABM`).

## 5. Descent Kernels (Uncapped Complete Iterative Q32 `dsieve16`, Strict-64b `BOUND_CELL_HI` + Yield Certificate, & Bareiss LLL)

- **`dsieve16` on `S_16(R) = B_16(R)/{+-1}` (Uncapped Complete Enumeration — ZERO Per-Prime Skip Caps!):**
  1. **Setup (`PAIR_SETUP`, 2-phase `(11/40)->(99/100)`-Bareiss-LLL) & Active Prefix `m`:** Bareiss `lll_16d_int` reduces each prime sublattice `L_{Q,P}` (relative-degree-1 ideals `p in [2^8, 2^lx]`), computes the active dimension prefix `m = min{t in {0..16} : forall j >= t, d_{j+1} > R^2 * d_j}`, and enumerates EVERY lattice point in `B_16(R)` without skipping any prime (mean per-hit cost `SIEVE_HIT = 2,200` Ir/hit (`2x` `c36p`'s `1,100` = walk/rank + owner transfer + `uint16_t` apply, zero per-prime skip caps) = walk/rank + owner transfer + `uint16_t` apply, zero per-prime caps).
  2. **Complete Iterative Q32 DFS (`dsieve16_q32_node` `39` ops, `dsieve16_dfs_ctrl` `31` ops, `dsieve16_dfs_cstep` `6` ops, `dsieve16_vec_madd` `40` ops, `dsieve16_leaf_vec` `40` ops, `dsieve16_rank_step` `10` ops):** Baseline scalar ARMv8.0-A integer instructions ONLY (`simd = 0`, zero NEON/SVE/LSE/DC ZVA/fp64/`udiv`).

```asm
; kernel: dsieve16_q32_node items=1 charged=39.0 key=dsieve16_node
orr x1, x0, #1
clz x2, x1
eor x2, x2, #63
and x3, x2, #62
lsl x4, x1, #9
lsr x4, x4, x3
ldrh w2, [x_lut, x4, lsl #1]
lsr x3, x3, #1
lsl x2, x2, x3
mul x3, x2, x7
umulh x4, x2, x7
lsr x3, x3, #30
orr x3, x3, x4, lsl #34
add x3, x3, #1
add x3, x3, x8
sub x9, x10, x3
add x9, x9, x11
asr x9, x9, #32
add x12, x10, x3
asr x12, x12, #32
cmp x13, #0
csel x9, xzr, x9, eq
lsl x14, x9, #32
sub x14, x14, x10
cmp x14, #0
sub x15, xzr, x14
csel x14, x15, x14, lt
subs x14, x14, x8
csel x14, xzr, x14, lt
mul x15, x14, x14
umulh x16, x14, x14
lsr x15, x15, #32
orr x15, x15, x16, lsl #32
mul x16, x15, x17
umulh x14, x15, x17
lsr x16, x16, #32
orr x16, x16, x14, lsl #32
subs x19, x0, x16
b.lo .Lprune
```

```asm
; kernel: dsieve16_dfs_ctrl items=1 charged=31.0 key=dsieve16_ctrl
cmp x9, x12
b.gt .Ldfs_pop
stp x9, x12, [x21, #0]
stp x0, x13, [x21, #16]
str x8, [x21, #32]
add x21, x21, #40
cmp x9, #0
sub x22, xzr, x9
csel x22, x22, x9, lt
add x8, x8, x22
cmp x9, #0
mov x23, #1
csel x23, x23, xzr, gt
orr x13, x13, x23
sub x20, x20, #1
add x24, x25, x20, lsl #7
ldr x10, [x24, x20, lsl #3]
ldr x7, [x26, x20, lsl #3]
ldr x17, [x27, x20, lsl #3]
mov x0, x19
b .Ldfs_enter
.Ldfs_pop:
add x20, x20, #1
cmp x20, x_m
b.eq .Ldfs_done
sub x21, x21, #40
ldp x9, x12, [x21, #0]
ldp x0, x13, [x21, #16]
ldr x8, [x21, #32]
add x9, x9, #1
b .Ldfs_sib
.Ldfs_done:
ret
```

```asm
; kernel: dsieve16_dfs_cstep items=1 charged=6.0 key=dsieve16_cstep
ldr x22, [x24, x23, lsl #3]
ldr x28, [x29, x23, lsl #3]
msub x22, x9, x28, x22
str x22, [x30, x23, lsl #3]
subs x23, x23, #1
b.ne .Lcstep_loop
```

  3. **Exact-integer leaf vector step (`v += b_0` after priced `k=0` chord init `v = V_stk[0] + (lo_0 - 1)*b_0`) + norm check + bijective `rank_lex` + bucket write (`316 <= 340.0` Ir/leaf):** Upon entering any non-empty `k = 0` chord `[lo_0, hi_0]`, `dsieve16_vec_madd` (`40` ops, priced in `k1_rej = 244.0` above!) initializes `v = V_stk[0] + (lo_0 - 1) * b_0`; then each `z_0` step executes the full 40-op `dsieve16_leaf_vec` block below (`8 ldp` of `b_0[0..15]` + `8 ldp` of `v[0..15]` + `16 add` + `8 stp` of `v[0..15]`; boundary `z_0` values with `||v||_2^2 > R^2` are priced in `k1_rej = 244.0`, while qualifying hits pay `c_leaf = 340.0`) + exact integer `sum v[t]^2 <= R^2` (`32` ops) + canonical sign-flip (`16` ops) + bijective dense ball rank via `129.03 MiB` (`R <= 64`, `< 256 MiB/worker`) flat `P_d[rem*(2R+1) + (v_d+R)]` table (`16` uniform levels `d = 0..15` x `10` ops = `160` ops, `dsieve16_rank_step` below loads `v_d`, computes `rem*(2R+1)+(v_d+R)` explicitly from live `x7 = rem`, indexes the 128-bit prefix entry, updates `(x5, x6)` and `rem = x7`, and advances `x_pd` to the next level table — simulated bit-for-bit on all `552/552` hits in `cm_models.verify_kernels_symbolic`!) + owner-chunk route (`56` ops) + `uint16_t` score add (`12` ops) -> `40 + 32 + 16 + 160 + 56 + 12 = 316 <= 340.0` Ir/leaf (`+24` Ir/leaf extra loop/setup margin)!

```asm
; kernel: dsieve16_vec_madd items=1 charged=40.0 key=dsieve16_vmadd
ldp x0, x1, [x_row, #0]
ldp x2, x3, [x_base, #0]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #0]
ldp x0, x1, [x_row, #16]
ldp x2, x3, [x_base, #16]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #16]
ldp x0, x1, [x_row, #32]
ldp x2, x3, [x_base, #32]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #32]
ldp x0, x1, [x_row, #48]
ldp x2, x3, [x_base, #48]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #48]
ldp x0, x1, [x_row, #64]
ldp x2, x3, [x_base, #64]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #64]
ldp x0, x1, [x_row, #80]
ldp x2, x3, [x_base, #80]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #80]
ldp x0, x1, [x_row, #96]
ldp x2, x3, [x_base, #96]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #96]
ldp x0, x1, [x_row, #112]
ldp x2, x3, [x_base, #112]
madd x2, x0, x_sc, x2
madd x3, x1, x_sc, x3
stp x2, x3, [x_dst, #112]
```

```asm
; kernel: dsieve16_leaf_vec items=1 charged=40.0 key=dsieve16_leaf_v
ldp x0, x1, [x_b0, #0]
ldp x2, x3, [x_v, #0]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #0]
ldp x0, x1, [x_b0, #16]
ldp x2, x3, [x_v, #16]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #16]
ldp x0, x1, [x_b0, #32]
ldp x2, x3, [x_v, #32]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #32]
ldp x0, x1, [x_b0, #48]
ldp x2, x3, [x_v, #48]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #48]
ldp x0, x1, [x_b0, #64]
ldp x2, x3, [x_v, #64]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #64]
ldp x0, x1, [x_b0, #80]
ldp x2, x3, [x_v, #80]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #80]
ldp x0, x1, [x_b0, #96]
ldp x2, x3, [x_v, #96]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #96]
ldp x0, x1, [x_b0, #112]
ldp x2, x3, [x_v, #112]
add x2, x2, x0
add x3, x3, x1
stp x2, x3, [x_v, #112]
```

```asm
; kernel: dsieve16_rank_step items=1 charged=10.0 key=dsieve16_rank
ldr x1, [x_v], #8
add x2, x1, x_R
mul x0, x7, x_2Rp1
add x2, x0, x2
add x2, x_pd, x2, lsl #4
ldp x3, x4, [x2, #0]
adds x5, x5, x3
adc x6, x6, x4
msub x7, x1, x1, x7
add x_pd, x_pd, x_stride
```

  4. **Certified 2-Stage Strict Signed-64b Outward-Rounded `BOUND_CELL_HI` + Live `ss_f * ss_g` Yield Certificate (`bound_cell_ball_mul` `48 <= 52` Ir):** Baseline scalar ARMv8.0-A `smulh`/`umulh` interval ball arithmetic certifies outward-rounded upper bounds on both norms `(N_f, N_g)` and accumulates the live per-node `ss_f * ss_g` yield certificate inside `BOUND_CELL` (discarding any node whose certified mean is below its bin's charged `pr` before cofactorisation, charged via `DISC` in `DESCENT.md`). Zero floating-point instructions; zero 64-bit overflow.

```asm
; kernel: bound_cell_ball_mul items=1 charged=52.0 key=bound_ball_mul
smulh x8, x0, x3
smulh x9, x1, x4
sub x8, x8, x9
smulh x10, x0, x4
smulh x11, x1, x3
add x10, x10, x11
cmp x0, #0
sub x12, xzr, x0
csel x12, x12, x0, lt
cmp x1, #0
sub x13, xzr, x1
csel x13, x13, x1, lt
add x12, x12, x13
cmp x3, #0
sub x13, xzr, x3
csel x13, x13, x3, lt
cmp x4, #0
sub x14, xzr, x4
csel x14, x14, x4, lt
add x13, x13, x14
umulh x14, x12, x5
umulh x15, x13, x2
add x14, x14, x15
umulh x15, x2, x5
add x14, x14, x15
add x14, x14, #5
cmp x8, #0
sub x9, xzr, x8
csel x9, x9, x8, lt
cmp x10, #0
sub x11, xzr, x10
csel x11, x11, x10, lt
cmp x9, x11
csel x9, x9, x11, hi
cmp x9, x14
csel x9, x9, x14, hi
cbz x9, .Lbn0
clz x11, x9
subs x11, x11, #4
csel x11, xzr, x11, lt
lsl x8, x8, x11
lsl x10, x10, x11
lsl x14, x14, x11
b .Lbn1
.Lbn0:
mov x11, #0
.Lbn1:
add x16, x6, x7
add x16, x16, #64
sub x16, x16, x11
```

- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every special-`q` node lattice `L_Q` (`det = q`) and `G_T` boot lattice runs 2-phase `(11/40)->(99/100)`-Bareiss-LLL (`lll_16d_int`), using exact 2-phase `(11/40 Euclidean -> 99/100 Q_CM)` Bareiss LLL (`G_8` has `det(G_8) = 3^4 * 5^6` and `lambda_max(G_8) = 15`, proved live via Sylvester minors of `15*I_8 - G_8`, so `d_k^(CM) <= 15^k * d_k^(Euc)` unconditionally!), yielding explicit derived uniform bounds **`PAIR_SETUP = 2^35.2`** (71-bit `L_p` Euclidean `11/40 -> 99/100`, `S_1 <= 1,159, S_2 <= 199,278, L_d <= 5` -> `2^35.06 < 2^35.2`), **`LLL_NODE = 2^40.0`** (special-`q` `q <= 2^262` `11/40 Euc -> 99/100 Q_CM`, `S_1 <= 4,220, S_2 <= 428,170, L_d <= 12` -> `2^39.41 < 2^40.0`), and **`LLL_BOOT = 2^46.5`** (381-bit `G_T` boot `det = p^4 <= 2^1524`, `S_1 <= 24,547, S_2 <= 1,733,726, L_d <= 52` -> `2^46.31 < 2^46.5`), all derived and asserted LIVE in `verify_kernels_symbolic` alongside `BOUND_CELL = 3,600`, `SCAN = 10`, `ZERO = 2`, `NORMS = 2^27`, `TDIV = 2^25`, `PRIME = 2^24`, `IDEAL = 2^35`.
