# Kernels (Raw-Opcode-Verified Baseline Scalar ARMv8.0-A: `Arch2-Prod` 8B Zero-Copy, 1-Owner Build, `uint16_t`, `(m-k)&3` Bank Rotation, 1-Step Barrett `row_red5_z12`, Input-Uniform Bareiss LLL)

Every kernel cost is derived from explicit baseline scalar ARMv8.0-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd=0`, zero NEON/SVE, zero `FEAT_LSE`, zero `DC ZVA`, zero `bl`). `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file, validates opcodes & `LDP`/`STP` addressing, counts every instruction (`cmp`+`b.*`=2 ops, loop control included), simulates the complete 1-owner flat-array 2-pass build (incl. overflow-proof 64-bit modular step at `p = 2^64 - 59`) + 1-rec-at-a-time Band-2 multi-SB replay (incl. empty/odd/even buckets) + `(m-k)&3` + SWAR + `MONT5` + Barrett `row_red5_z12` + Shoup NTT + 4-prime Garner CRT + iterative Q32 `dsieve16` + 10-op/level `dsieve16_rank_step` + certified upper-bound `BOUND_CELL_HI` live, and asserts `opcodes/items <= charged`.

## 1. Collection Kernels (`Arch2-Prod-ZeroCopy-U16`, Pure Product `qa <= s_c R2_a` & `qb <= R2_b`, `W_tile = 400`, `W_SB = 32`, `P_m4 = 2^56.5`, `P_lo = P_f = 2^59`, `P_g = 2^64`)

Because the sieving region is a pure product of two 8D balls (`qa(a) <= s_c R2_a` with 2-adic radii `s_c = (1.0, 1.5, 2.0)` and `qb(b) <= R2_b`, `A^2 = 310,000`, `r_ab = 1.00`), every canonical `a`-table entry (`N_a^can = V_a = 2^56.7395 < 2^57`, `N_a^tbl = 4 V_a = 2^58.7395`) is unconditionally in-region for every `b` (`V_b = 2^63.2680`, `L_max = 393.70 < W_tile = 400`, `L_b = 169.10`, `N_7 = 2^53.8663`). There is **ZERO `q_key`, ZERO `q_tab`, ZERO `sb_min_q8`, ZERO binary search, ZERO global sorting invariant, and ZERO `f`-side Band 2** (`P_f = P_lo = 2^59`).
- **Zero-wrap `uint16_t` accumulators (`S_g, S_f`) with exact `2^57` bank stride:** Quarter-bit weights on `<= 700`b norms reach `<= 2,800 << 32,767` (`0` wraparound). Each column owner allocates `4 * 2^57` `uint16_t`s for `S_g` (`2^60` B) + `S_f` (`2^60` B) = `2^61.0` B `<< 2^64` B (exact power-of-2 bank stride `2^57` -> ZERO stride mismatch for `((m - k) & 3) << 57 | oid`).
- **Zero-copy 8B records across Band-1a/1b/2 with exact `(m - k) & 3` rotation:**
  - Band-1a (`[2^24, 2^56.5)`): `m`-partitioned `oid < 2^57` -> `rev_b1a_apply_u4` (`4.00 <= 4.25`).
  - Band-1b (`[2^56.5, 2^59)`): 8B `(m << 57) | oid` (`sub x4,x0,x_k57; and x4,x4,x_mask59`) -> 1-record-at-a-time `rev_b1b_apply_1` (`8.00 <= 8.00`, ZERO remainder path).
  - Band-2 (`g`-side `[2^59, 2^64)` only): 8B `(r << 59) | (m << 57) | oid` (`r = ka & 31`, `5 + 2 + 57 = 64` bits exact!). Outer loop state lives exclusively in callee-saved `x19..x22` (`x19=off_ptr`, `x20=ent_base`, `x21=rem_sb`, `x22=bptr_base`, advancing `add x22, x22, #256` per SB), read-only flush limit in `x18`, bank-rotate constants in `x3, x4` — physically disjoint from emit scratch `x9..x11, x15` (asserted automatically in `verify_kernels_symbolic`). 1-record-at-a-time `rev_b2_prod_emit_1` (`12 <= 12.5` Ir/hit, safe on every bucket length `0, 1, 2, ...`, ZERO odd-length overread) + `rev_b1a_apply_u4` (`4.00 <= 4.25`) = **`16.00 <= C_HIT_HI = 16.75`**.

### 1.1 `rev_sb_probe_fast` (Band-2 row start: `17 <= C_TILE_HDR = 18.0` Ir/row)

```asm
; kernel: rev_sb_probe_fast items=1 charged=18.0 key=c_tile_hdr
sub x0, x_k0, x_p
mul x1, x0, x_invp
umulh x2, x0, x_mu
lsr x2, x2, #4
msub x3, x2, x_p, x0
lsr x4, x3, #5
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

### 1.2 `rev_b2_prod_sb_hdr` (Band-2 super-bucket header with uniform byte offsets: `11 <= C_SB_HDR = 12.0` Ir/SB)

Outer loop state lives exclusively in callee-saved `x19` (`off` ptr), `x20` (`ent` base), `x21` (remaining SB count), `x22` (`bptr` base, explicitly advanced `+256` B = 32 column pointers per SB on `.Lnext_sb`); flush limit lives in dedicated `x18` (never clobbered by `x5`!). Bucket boundaries `off[j], off[j+1]` are uniform 8B-aligned byte offsets:

```asm
; kernel: rev_b2_prod_sb_hdr items=1 charged=12.0 key=c_sb_hdr
ldp x5, x6, [x19]
add x19, x19, #8
cmp x5, x6
b.eq .Lnext_sb
add x0, x20, x5
add x1, x20, x6
mov x2, x22
b .Lemit1
.Lnext_sb:
add x22, x22, #256
subs x21, x21, #1
b.ne .Lsb
```

### 1.3 Band-1a & Band-2-apply zero-copy `uint16_t` 4x-unrolled + 1-record tail direct apply (`rev_b1a_apply_u4`: `16/4 = 4.00 <= C_HIT_LO_M4 = 4.25` Ir/hit; `rev_b1a_apply_tail1`: `6` Ir/rem hit, `<= +8` ops/bkt charged in `b1_lk = 32.0`)

Applies either a Band-1a `m`-partitioned `oid < 2^57` or a pre-rotated Band-2 59-bit word `(((m - k) & 3) << 57) | oid < 2^59` directly into contiguous `uint16_t acc[4 * 2^57]`. Before entering `.Lb1a`, `cmp x_cnt, #4; b.lo .Ltail1` branches directly to the explicit 1-record tail loop `rev_b1a_apply_tail1` for any bucket length `k = 0, 1, 2, 3` or remainder `k mod 4 in {1, 2, 3}` (and after `.Lb1a` exits with `x_cnt in {0..3}`). Total worst-case overhead above `4.00` Ir/hit for ANY bucket length `k >= 0` is `<= 2 + 3 * (6 - 4) = 8` ops per bucket (explicitly charged in `b1_lk = 32.0` for Band-1a; for Band-2 flushes of `>= 256` hits, `<= 8 / 256 = 0.03125` Ir/hit fits inside `4.00 + 0.03125 <= 4.25` Ir/hit; simulated bit-for-bit across all lengths `k = 0..9` in `audit_cm.py`!):

```asm
; kernel: rev_b1a_apply_u4 items=4 charged=4.25 key=c_hit_lo_m4
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
subs x_cnt, x_cnt, #4
b.hs .Lb1a
```

```asm
; kernel: rev_b1a_apply_tail1 items=1 charged=6.0 key=c_hit_lo_m4_tail
ldr x0, [x_p], #8
ldrh w4, [x_acc, x0, lsl #1]
add w4, w4, w_lp
strh w4, [x_acc, x0, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Ltail1
```

### 1.4 Band-1b zero-copy `uint16_t` 1-record-at-a-time 8B packed `{m:2, oid:57}` direct apply (`rev_b1b_apply_1`: `8.00 <= C_HIT_LO_U = 8.00` Ir/hit, ZERO remainder path needed)

Every Band-1b hit packs `(m << 57) | oid` into one 64-bit word (`oid < 2^57`, `m in {0..3}` -> 59 active bits). Because Band-1b covers `[2^56.5, 2^59]` where bucket lengths are frequently `1..3`, `rev_b1b_apply_1` processes **1 record per iteration** (handling EVERY bucket length `k = 0, 1, 2, 3, ...` identically with ZERO unroll-remainder path!), rotating bank `m -> (m - k) & 3` via `sub x4, x0, x_k57; and x4, x4, x_mask59` (`x_k57 = k << 57`, `x_mask59 = (1 << 59) - 1`, verified LIVE on all 16 `(m, k)` pairs and all lengths `k = 0..9`!):

```asm
; kernel: rev_b1b_apply_1 items=1 charged=8.0 key=c_hit_lo_u
ldr x0, [x_p], #8
sub x4, x0, x_k57
and x4, x4, x_mask59
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
subs x_cnt, x_cnt, #1
b.ne .Lb1b
```

### 1.5 Band-2 1-record-at-a-time 8B `(r<<59)|(m<<57)|oid` emit (`rev_b2_prod_emit_1`: `12 <= 12.5` Ir/hit; `+ 4.25 = 16.75` total)

Entered from `rev_b2_prod_sb_hdr` (`b .Lemit1`). Emits **one 8-byte record per iteration** (ZERO odd-length overread bug on any bucket size `0, 1, 2, ...`!), extracts `r = rec >> 59`, rotates bank `m -> (m - k) & 3` via `x3 = k << 57`, `x4 = (1 << 59) - 1`, appends to `bptr[r]`, checks flush limit `x18` (disjoint from `x5` and `x19..x22`!), and loops `cmp x0, x1; b.lo .Lemit1` falling through to `.Lnext_sb`. Boundary super-buckets (`w_sb / L_b = 18.9%` excess) add `3` range-check ops (`sub/cmp/b.hs`, charged `c_bd_chk = 3.0 * (1 + w_sb / L_b)`):

```asm
; kernel: rev_b2_prod_emit_1 items=1 charged=12.5 key=c_b2_emit
ldr x15, [x0], #8
lsr x9, x15, #59
sub x15, x15, x3
and x15, x15, x4
add x10, x2, x9, lsl #3
ldr x11, [x10]
str x15, [x11], #8
str x11, [x10]
cmp x11, x18
b.hs .Lflush
cmp x0, x1
b.lo .Lemit1
```

### 1.6 8-Cell `uint16_t` 64-Bit Borrow-Free SWAR Scan + Zero (`rev_scan_zero_u16`: `14/8 = 1.75` Ir/cell -> `1.75 * 1.0276 + 0.15 = 1.948 <= C_SCAN_ZERO = 2.00`)

Each 64-bit word holds 4 `uint16_t` quarter-bit accumulators (`<= 2,800 << 32,768`). For any threshold `T in [1..32767]`, set `bias16 = (32768 - T) & 0x7fff` broadcast across 4 16-bit lanes (`x_bias16`), `x_7fff = 0x7fff_7fff_7fff_7fff`, `x_8000 = 0x8000_8000_8000_8000`. Then for EVERY `uint16_t` lane `x in [0..65535]`: `x >= T` iff `(((x & 0x7fff) + bias16) | x) & 0x8000 != 0` (verified LIVE across full threshold sweep & boundary edges inside `audit_cm.py`, ZERO borrow between 16-bit lanes!):

```asm
; kernel: rev_scan_zero_u16 items=8 charged=1.80 key=c_scan_zero_raw
ldp x0, x1, [x_g]
stp xzr, xzr, [x_g], #16
and x2, x0, x_7fff
add x2, x2, x_bias16
orr x2, x2, x0
and x3, x1, x_7fff
add x3, x3, x_bias16
orr x3, x3, x1
orr x4, x2, x3
tst x4, x_8000
b.ne .Lhit8
stp xzr, xzr, [x_f], #16
subs x_cnt, x_cnt, #8
b.ne .Lscan8
```

Including tile padding `p2_mult = 1.0276` + setup `0.10` + hit-branch amortization `< 0.05`: `1.75 * 1.0276 + 0.15 = 1.9483 <= C_SCAN_ZERO = 2.00`.

### 1.7 1-Owner Flat-Array Counting-Sort / Prefix-Sum / Scatter Build (`rev_build_prod_pass12`: `22 + 0.50 = 22.50 <= C_LIST_A = 24.0`, `rev_build_prod_prefix`: `6 <= C_LBKT = 6.0`)

Because `Arch2-Prod` requires **ZERO `q_key` sorting**, Pass 1 walks `a_0`-chords in deterministic 4-bank order (`m = 0..3`, `oid = 0..V_a - 1`), advancing `ka = (ka + step) mod p` via the **4-opcode overflow-proof modular step** (`x_pms = p - step` precomputed once per prime; exact for ALL `0 <= ka, step < p < 2^64` with ZERO carry loss even for `p = 2^64 - 59`, tested LIVE!), writing `kbuf[idx] = ka`, and incrementing **64-bit `uint64_t` `x_hist[ka >> 5]`** (`ZERO` 32-bit overflow even at `p = 2^24`!). `rev_build_prod_prefix` writes BOTH immutable byte offsets `x_off[j]` AND mutable scatter cursors `x_pos[j]` in one 6-op pass (plus terminal sentinel `off[n_sb] = 8 * N_a_tbl`). Pass 2 executes **4 outer bank passes (`m = 0, 1, 2, 3`)**, initializing `x_oidm = m << 57` and `x_n2 = V_a` at the start of each bank `m` (so the bank transition after `oid = V_a - 1` in bank `m` explicitly resets `x_oidm = (m + 1) << 57` rather than overflowing into `(m << 57) | V_a`; `4` outer bank-init ops across `N_a^tbl = 2^58.74` entries = `< 10^-16` Ir/entry, verified LIVE across all 4 banks with `V_a = 40 << 2^57` in both `verify_kernels_symbolic` and `check_sigmaB_exactness`!), streaming `x_kbuf` in the identical order (ZERO `x_meta`, ZERO `x_imap`!), packing `((ka & 31) << 59) | x_oidm`, and scattering 8B into `x_out[pos[sb]]`. **Total simultaneous 1-owner RAM (`x_kbuf + x_out + x_off + x_pos + 400 MiB`) `<= 2^63.876 B < 2^64 B` for ALL Band-1a/1b/2 primes (`0` cross-owner build transfers, asserted live in `audit_cm.py`)**:

```asm
; kernel: rev_build_prod_pass12 items=1 charged=23.5 key=c_list_a_raw
cmp x0, x_pms
sub x5, x0, x_pms
add x6, x0, x_step
csel x0, x5, x6, hs
str x0, [x_kbuf], #8
lsr x1, x0, #5
ldr x2, [x_hist, x1, lsl #3]
add x2, x2, #1
str x2, [x_hist, x1, lsl #3]
subs x_n1, x_n1, #1
b.ne .Lp1
ldr x0, [x_kbuf2], #8
and x5, x0, #31
orr x8, x_oidm, x5, lsl #59
add x_oidm, x_oidm, #1
lsr x1, x0, #5
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

Exactness of 4-orbit `sigma`-shared pure-product replay with `D_(omega,sigma)` (`239/239` exact-1), 8B `(r<<59)|(m<<57)|oid` flat-array build, overflow-proof 64-bit modular step, and odd/even/empty bucket replay is verified LIVE in `audit_cm.py` (`sigmaB_exactness`: 48 real BLS12-381 `w' = omega^2` ideals, `0` mismatches).

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
adcs x18, x18, x13
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
adcs x18, x18, x12
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
adcs x17, x18, x13
adcs x18, x19, x14
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
adcs x18, x18, x13
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
mov x18, #0
mov x19, #0
mov x20, #0
subs x10, x15, x22
sbcs x11, x16, x23
sbcs x12, x17, x24
sbcs x13, x18, x25
sbcs x14, x19, x26
sbcs xzr, x20, xzr
csel x15, x10, x15, hs
csel x16, x11, x16, hs
csel x17, x12, x17, hs
csel x18, x13, x18, hs
csel x19, x14, x19, hs
stp x15, x16, [x_out, #0]
stp x17, x18, [x_out, #16]
str x19, [x_out, #32]
ret
```

Total exact `MONT5 = 5 * 48 + 27 = 267 <= 270.0` Ir.

## 3. Explicit Multi-Limb Cyclotomic Tower Norm Kernels, Parsed 4-Prime 62-Bit NTT (`fft_F`) & Garner CRT, Stage B, 2D Cofactorization, & Exact 2nd-Moment ECM Chebyshev

- **Explicit Schoolbook Multi-Limb Scalar ARMv8-A Tower Norm Derivation (`Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`, ZERO Karatsuba shortcut):**
  Explicit schoolbook scalar ARMv8-A `L_1 x L_2` limb multiplication uses `mul_LL(L_1, L_2) = 8 L_1 L_2 + 4 (L_1 + L_2)` (`mul/umulh/adds/adcs` carry chains). On the production sieving region (`|a_j|, |b_j| <= 2^8`), `a^3, b^3 in Z[y]` fit in 1 limb (`<= 27`b; `4` `8x8` 1-limb poly muls + `Phi_15` reduction = `400` Ir), `g(a, b) = a^3 - u w^2 b^3` fits in 2 limbs (`<= 94`b; `+96` Ir), and `f(a, b) = a^9 + (w+1) a^6 b^3 + a^3 b^6 - w b^9` fits in 3 limbs (`<= 128`b; `a^6, b^6` = `2 x 928 = 1,856` Ir, `a^9, b^9, a^6 b^3, a^3 b^6` = `4 x 1,152 = 4,608` Ir, combine = `128` Ir). Using pure schoolbook `mul_LL` across the 3-step tower norm `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:
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
- **Analytical Stage B (`cm_models.stage_b_analytic`):** `gbatch_cost(24197143, L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` with `MONT5 = 270.0`.
- **Live 2D-Histogram Cofactorization (`Q_F_SAFETY = 2.0`):** Under `H1`, `2.0 * rho_exact((bf + alpha_f) / lf)` at left endpoint `bf` gives `>= 2.27x` bin-midpoint `rho_db`, integrated LIVE across 401 bins up to `bg_max = 700`b (`c_cof_2d_upper = 2.8950` ipc).
- **Exact 2nd-Moment ECM Chebyshev (`cm_models.ecm_mean_var`):** Schedule runs `K_c` curves (`p_1`, `q_s=(1-p_1)^{K_c}<=0.001`, `p_s=1-q_s`). Failed schedule cost `C_fail=K_c*per_c` (`Var=0`); final successful schedule has `S~Bin(K_c,p_1)|S>=1`, `C_succ=C_fail+S*redo_ex`, `E[S|S>=1]=K_c p_1/p_s`, `E[S^2|S>=1]=(K_c p_1(1-p_1)+(K_c p_1)^2)/p_s`. With independent `F~Geom_0(p_s)`, `T=F*C_fail+C_succ` gives exact `E[T]=(q_s/p_s)C_fail+E[C_succ]`, `Var(T)=(q_s/p_s^2)C_fail^2+Var(C_succ)`. Cauchy-Schwarz `sd_rel=13*sqrt(Var_f)+16*sqrt(Var_g)` yields Chebyshev `< 2^-74.4 << 2^-70`.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (Straight-Line 4-Limb Scalar ARMv8-A + Exact 1-Step Barrett `row_red5_z12`)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`):
- **`unit_run5_z12` (`43` opcodes / `5` entries = `8.60 <= 8.80` Ir/vec-entry):** Unrolled 5-entry run adding 4-limb entries into 5-limb accumulator (`x10..x14`, carrying cleanly into the 5th limb `x14` via `adc x14, x14, xzr`):

```asm
; kernel: unit_run5_z12 items=5 charged=8.80 key=unit_run5
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v1, [x_idx], #8
ldp x0, x1, [x_v1, #0]
ldp x2, x3, [x_v1, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v2, [x_idx], #8
ldp x0, x1, [x_v2, #0]
ldp x2, x3, [x_v2, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v3, [x_idx], #8
ldp x0, x1, [x_v3, #0]
ldp x2, x3, [x_v3, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v4, [x_idx], #8
ldp x0, x1, [x_v4, #0]
ldp x2, x3, [x_v4, #16]
adds x10, x10, x0
adcs x11, x11, x1
adcs x12, x12, x2
adcs x13, x13, x3
adc x14, x14, xzr
ldr x_v0, [x_idx], #8
add x_ptr, x_ptr, #40
subs x_cnt, x_cnt, #5
b.gt .Lrun5
```

- **`row_red5_z12` (`67 <= 68.0` Ir/bucket-vec — Exact 1-Step Signed 5-Limb Barrett Reduction for ALL `|a_4| < 2^60`):**
  Handles ANY signed 5-limb bucket accumulator `A = (x10..x14)` with `|a_4| < 2^60` (`|A| < 2^316`, covering EVERY `|c| <= 255` row up to `+-2^59 * (R_BLS - 1) > 10^17 * (R_BLS - 1)` with ZERO carry overflow!). Computes branchless `U = |A|`, `top = (u_4 << 2) | (u_3 >> 62)`, `q_hat = umulh(top, mu_r)` where `mu_r = floor(2^318 / R_BLS) = 0x8d54253b7fb78ddf < 2^64` (so `0 <= floor(U / R_BLS) - q_hat <= 2` unconditionally!), `rem = U - q_hat * R_BLS in [0, 3 * R_BLS)`, two branchless conditional `- R_BLS` steps, and branchless sign fix (`if A < 0 and rem != 0: rem = R_BLS - rem`). Verified bit-for-bit LIVE in `audit_cm.py` on `0, +-1, +-3*(R-1), +-200*(R-1), +-51,000*(R-1), +-1,000,000*(R-1), +-2^59*(R-1)`, and 30 random `|A| <= 2^315` inputs:

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
umulh x18, x5, x24
umulh x19, x5, x25
adds x7, x7, x16
adcs x8, x8, x17
adcs x9, x9, x18
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

- **`gen_run1_z12` (`23 <= 29.0` Ir/coord -> `4 * 29.0 = 116.0` Ir/non-unit `Z[zeta_12]` entry):** Because a general non-unit element `c = c_0 + c_1 chi + c_2 chi^2 + c_3 chi^3 in Z[zeta_12]` (`phi(12) = 4`, `|c_j| <= 255`) decomposes into up to 4 signed basis-direction scalar updates `(c_0 * 1, c_1 * chi, c_2 * chi^2, c_3 * chi^3)` into the existing buckets `0..3` (with `c_j > 0` using `adds/adcs` and `c_j < 0` using the identical 23-opcode `subs/sbcs` run on `|c_j|` into the same bucket — ZERO extra buckets beyond the 6 `{chi^0..chi^5}` buckets and zero assumption on `nnz(c) < 4`!), we explicitly charge ALL 4 basis coordinates per non-unit entry: `4 * 23 = 92 <= 4 * 29.0 = 116.0` Ir/non-unit entry (verified live in `audit_cm.py`!). Each coordinate update executes `23 <= 29.0` baseline scalar ARMv8-A opcodes:

```asm
; kernel: gen_run1_z12 items=1 charged=29.0 key=gen_run1_coord
ldp x0, x1, [x_v0, #0]
ldp x2, x3, [x_v0, #16]
mul x4, x0, x_c
mul x5, x1, x_c
mul x6, x2, x_c
mul x7, x3, x_c
adds x10, x10, x4
adcs x11, x11, x5
adcs x12, x12, x6
adcs x13, x13, x7
adc x14, x14, xzr
umulh x4, x0, x_c
umulh x5, x1, x_c
umulh x6, x2, x_c
umulh x7, x3, x_c
adds x11, x11, x4
adcs x12, x12, x5
adcs x13, x13, x6
adcs x14, x14, x7
ldr x_v0, [x_idx], #8
ldr x_c, [x_coef], #8
subs x_cnt, x_cnt, #1
b.gt .Lgen1
```

- **`chi_mul_z12_fast` (`50` Ir), `row_tail_6bucket_z12` (`504` Ir), 1-Owner Arch B Row Overhead (`c_row_chi = 1,000.0` Ir):** `6 * 67 + 504 + 16 = 922 <= 6 * 68 + 504 + 16 = 928 <= c_row_chi = 1,000.0` Ir worst-case (`0` cross-owner sends/combines; `k5_live = eu_lo * 8.80 + (1 - eu_lo) * 116.0 = 16.316 <= k5_charged = 17.50`; `red_row_u_z12 = 18.6` Ir/elem).

## 5. Descent Kernels (Zero-Skip Complete Iterative Q32 `dsieve16`, Strict-64b `BOUND_CELL_HI`, & Bareiss LLL)

- **`dsieve16` on `S_16(R) = B_16(R)/{+-1}` (`R in [60, 128]`, ZERO skipped primes!):**
  1. **Setup (`PAIR_SETUP = 2^35.0` Ir, proved) & ZERO-SKIP `1.00 * WC_LLL99` Expected-Work Budget:** Bareiss `lll_16d_int(9,10)` on `L_{q,Q}` (`det=p<=2^70`, `d_k=p^2<=2^140`, `<=13,815` swaps -> `2^31.20 << 2^35.0` Ir), then Q32 tables `MU, L2_LO, INV_L_HI`. **ZERO primes are skipped**: every prime ideal `Q` with `N(Q) <= 2^{lx}` is sieved completely by iterative Q32 `dsieve16` (100% lossless score accumulation into `S`, ZERO missing score mass before `BOUND_CELL_HI` + `S >= L_hi - T`!). Under `H2` (residue equidistribution of `L_{q,Q}` mod `p`, enclosed inside the unified hard aborting meter `t_run = 2^123.0`), expected internal DFS nodes per prime are bounded by `1.00 * WC_LLL99(R, p)` (`8/8` live proved production `Q(zeta_15)` primes `16b..70b` have `B_int / WC_LLL99 in [0.181, 0.373] < 0.45 << 1.00`, i.e., `2.7x..5.5x` safety headroom, and actual simulated `nodes_q32 / B_int <= 0.0008`!).
  2. **Complete Iterative Q32 DFS (`w_k = 132 + 6*k` Ir at level `k`, weighted max `143.3 <= c_node = 144.0` Ir/node, ZERO floats, opcode-counted live):** Three parsed ARMv8.0-A blocks cover EVERY DFS instruction: `dsieve16_q32_node` (`55` ops: `cbz`+`clz`+5-iter `udiv` NR `isqrt(x0)+1`, proved for all `x0 in [0,2^46]` + `[lo,hi]` + `sq`), `dsieve16_dfs_ctrl` (`31` ops: sibling step, `stp/ldp` stack frame, `asum`/`nz`, child loads), `dsieve16_dfs_cstep` (`6*k+4` ops: incremental `C_stk[k-1][0..k-1]` update via `msub`/`sub`, exact by convexity of `sq`), plus `<=2` convex boundary rejects (`<=40` ops) -> `55+31+40+6+6*k = 132+6*k` Ir/level (`W_int_wc/B_int_wc <= 143.3 <= 144.0` asserted in `verify_kernels_symbolic`; live `R=2,p=31` & `R=64` 84b primes assert `ops_sim <= sum (132+6k) N_k <= 144*B_int`).

```asm
; kernel: dsieve16_q32_node items=1 charged=55.0 key=dsieve16_node
cbz x0, .Lnr_zero
clz x1, x0
mov x2, #65
sub x1, x2, x1
lsr x1, x1, #1
mov x2, #1
lsl x2, x2, x1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
udiv x3, x0, x2
add x2, x2, x3
lsr x2, x2, #1
add x2, x2, #1
b .Lnr_done
.Lnr_zero:
mov x2, #1
.Lnr_done:
mul x3, x2, x7
umulh x4, x2, x7
lsr x3, x3, #16
orr x3, x3, x4, lsl #48
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
umulh x18, x15, x17
lsr x16, x16, #32
orr x16, x16, x18, lsl #32
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
cmp x20, #16
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

  3. **Exact-integer leaf check + bijective `rank_lex` + bucket write (`292 <= 296.0` Ir/leaf):** At `k = 0`, stepping `z_0 -> z_0 + 1` updates `v[0..15] += b_0[0..15]` (`16` adds) + exact integer `sum v[t]^2 <= R^2` (`32` ops) + canonical sign-flip (`16` ops) + bijective dense ball rank via `0.94 GiB` flat `P_d[rem*(2R+1) + (v_d+R)]` table (`16` uniform levels `d = 0..15` x `10` ops = `160` ops, `dsieve16_rank_step` below loads `v_d`, computes `rem*(2R+1)+(v_d+R)` explicitly from live `x7 = rem`, indexes the 128-bit prefix entry, updates `(x5, x6)` and `rem = x7`, and advances `x_pd` to the next level table — simulated bit-for-bit on all `552/552` hits in `audit_cm.py`!) + owner-chunk route (`56` ops) + `uint16_t` score add (`12` ops) -> `16 + 32 + 16 + 160 + 56 + 12 = 292 <= 296.0` Ir/leaf!

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

  4. **Certified 2-stage strict signed-64b 60-bit-mantissa (`smulh`/`umulh`) UPPER-BOUND ball `BOUND_CELL_HI` (`BOUND_CELL_g = 1,800` Ir on all cells + `BOUND_CELL_f = 1,800` Ir on `g`-survivors, ZERO `p_nc` heuristic surcharge, ZERO floats, ZERO 64-bit overflow):** 16-op pre-shift scales every input `|a_j|,|b_j|<2^62` (`+31`b above prod `<2^31`) into 60-bit signed mantissas `<2^60` against Q58 `Z_TAB` (`err <= 0.51 * 2^-58`). Every disk maintains `|re|,|im|,rad < 2^60` with left-normalization (`clz`), so `bound_cell_ball_mul` (`48 <= 52` Ir using baseline scalar ARMv8.0-A `smulh`/`umulh`) preserves `>= 50` bits of relative precision (`rad <= 2^-48 M_poly`) across all 11 Horner steps with ZERO 64-bit overflow. At output, shifting right by 30 gives `r30, i30 < 2^30`, `e30 = ((rad+(1<<30)-1)>>30)+2`, and `m_hi = isqrt(r30^2+i30^2) + 1 + e30`, yielding the **UNCONDITIONAL upper bound** `L_hi = sum_{k=0}^3 2 * log2_ceil(m_hi * 2^{sh+30}) >= log2|N_exact|` (via pure-Taylor `ROUND_CEILING` `_LUT_HI`) for **100% of ALL cells — including near-cancellation cells (`m_lo <= 0`)!** **Consequence 1 (Unconditional Survivor Bound, ZERO `2^-24` Surcharge):** Gating every cell with `S >= L_hi - T` (`T = kA*nA + kB*n_c + 96.0`) **UNCONDITIONALLY** enforces `S >= log2|N_exact| - T` across 100% of all discrete lattice cells, so survivor counts are deterministically `<=` the exact `T`-survivor model (`sv`) with **ZERO extra survivors and ZERO heuristic `p_near_cancel` surcharge!** **Consequence 2 (Proved Two-Sided Tightness on the Entire `QBINS` Domain):** Whenever `m_lo = isqrt(r30^2+i30^2) - e30 >= 64 * e30` across all 4 embeddings (i.e. `|sigma_k(val)| >= 2^-41 M_poly`, covering 100% of the generic-norm `QBINS` domain where norms lie within a few bits of generic and never rely on `-82`-bit near-root cancellation), `m_hi <= (133/128) |z_k|`, proving **`log2|N_exact| <= L_hi <= log2|N_exact| + 0.6212 < log2|N_exact| + 1.0`**! Since `T` budgets `1.0` bit above cofactor/prime-power slack, 100% of well-conditioned smooth relations pass deterministically (`0/151` live cells fail `m_lo >= 64 * e30`, `max_gap_g < 0.25`b, `max_gap_f < 0.25`b `< 0.6212`b!). **Stage 1 (`g_3`, 5 muls/emb):** `16+4*(104+5*48+32+40)=1,680<=1,800` Ir. **Stage 2 (`f_9` on `g`-survivors, 6 muls/emb):** `4*(6*48+96+40)=1,696<=1,800` Ir. Exact integer norm caps `|N_g| <= 2^{ng_max}, |N_f| <= 2^{nf_max}` are enforced inside `NORMS` (`pre_gate`) at zero extra cost (`36/36` production lattices lie `36..51`b below!). SWAR `SCAN+ZERO=6`, `NORMS=2^27`, `TDIV=2^25`, `IDEAL=2^35`, `MONT5=270`, MR `<=2^-80`.

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

- **100% Input-Uniform 2-Stage Pure-Int Bareiss LLL (`11/40 -> 99/100`, ZERO Float LLL):** Every reduction runs `_dw_int(B, 11, 40)` -> `_dw_int(B_1, 99, 100)` with full pre-Lovasz size-reduction (`|2*lam[i][j]| <= d[j+1]`). Swaps shrink `D = prod_{k=1}^{15} d_k` by `< delta` with `max_k d_k <= D_max` non-increasing; `L_lam = L_b + L_d`, `L_num = 2 L_d + L_b`, Stage-1 bound `d_k <= 40^{k(16-k)/2} d_{16}^{k/16}` (`||b_i||^2 <= 5 D_{max,2}`, asserted in `verify_kernels_symbolic`):
  1. **Special-`q` (`q <= 2^265`, `LLL_NODE = 2^40`):** Stage 1 (`S_1 <= 4,381`): `<= 2^32.41` Ir; Stage 2 (`S_2 <= 398,940`): `<= 2^39.00` Ir -> **`2^39.01 << 2^40`** Ir.
  2. **381b `G_T` Boot (`d_{16} <= 2^3384`, `LLL_BOOT = 2^48`):** Stage 1 (`S_1 <= 36,076`): `<= 2^43.16` Ir; Stage 2 (`S_2 <= 1,875,188`): `<= 2^46.30` Ir -> **`2^46.45 << 2^48`** Ir.
