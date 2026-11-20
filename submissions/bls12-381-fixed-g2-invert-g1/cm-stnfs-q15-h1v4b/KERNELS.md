# Kernels (100% Raw-Opcode-Verified Baseline Scalar ARMv8-A: `uint16_t` Accumulators, Zero SIMD, Zero `FEAT_LSE`, Zero External C Measurements)

Every kernel cost constant in this package is derived directly from explicit line-by-line baseline scalar ARMv8-A assembly (`-march=armv8-a -mgeneral-regs-only`, `simd = 0`, zero NEON/SVE registers `v0..v31`/`q0..q31`, zero `FEAT_LSE` atomics, zero `DC ZVA`, zero `bl` libc calls, zero external unbundled C compilation measurements). `audit_cm.py` (`kernels_symbolic`) parses every fenced `asm` block directly from this file, validates every mnemonic and operand against the strict baseline scalar ARMv8-A allowlist, counts EVERY actual instruction line (counting `cmp` and conditional branch `b.*` as 2 separate instructions and including `subs`/`b.ne` loop control), and asserts `counted_opcodes / items <= charged_constant`.

## 1. Collection Kernels (`Arch2-Rev-SigmaB-v3-ZeroCopy-U16`, `W_tile = 552`, `W_SB = 64`, `P_lo = 2^55`, `P_m4 = 2^48`, `P_g = 2^63`, `P_f = 2^60`)

- **Zero-wrap `uint16_t` accumulators (`S_g, S_f`):** Quarter-bit log weights on `<= 700`-bit norms reach at most `4 * 700 = 2,800 << 32,767 < 65,535`, fitting in `uint16_t` with ZERO wraparound (`2 B` for `S_g` + `2 B` for `S_f` = `4 B` per `(a_oid, m)` cell; total column-owner arena `4 * N_a_tbl = 2^61.829 B < 2^64 B`). All hit-apply kernels use native baseline scalar `ldrh / add / strh` with shifted halfword offset `[x_acc, x_idx, lsl #1]`.
- **Zero-copy Band-1a/1b/2:** Band-1a (`[2^24, 2^48)`) and Band-1b (`[2^48, 2^55)`) apply hits directly in-place from the contiguous `[ptr_lo, ptr_hi)` bucket slice into `acc[4][2^58]` (zero intermediate copy buffer). Because `log2(V_a) = 57.8294 < 58`, each Band-1b entry packs `(m << 58) | oid` into ONE 64-bit word (`8 B`, 60 active bits). Band-2 (`[2^55, P_hi)`) emits `(t, oid)` into per-row buckets `bucket[t]` and applies each bucket directly via `rev_b1a_apply_u4` (zero copy).

### 1.1 `rev_sb_probe_fast` (Band-2 row start: `18` Ir/row <= `C_TILE_HDR = 20.0`)

```asm
; kernel: rev_sb_probe_fast items=1 charged=20.0 key=c_tile_hdr
sub x0, x_k0, x_p
mul x1, x0, x_invp
umulh x2, x0, x_mu
lsr x2, x2, #4
msub x3, x2, x_p, x0
add x4, x3, x_sb_m1
lsr x4, x4, #6
sub x5, x_end, x3
cmp x5, x_wmax
csel x5, x_wmax, x5, hi
lsl x6, x4, #3
add x7, x_off, x6
ldp x8, x9, [x7]
cmp x8, x9
b.gt .Lempty
add x_sb, x_sb, #1
add x_ptr, x_base, x8
ret
```

### 1.2 `rev_b2_minkey_gate` skip (`6` Ir/skip <= `C_SB_SKIP = 6.0`) & enter (`18` Ir/enter <= `C_SB_ENTER = 20.0`)

`cmp` and `b.ls` are counted as 2 separate instructions:

```asm
; kernel: rev_b2_minkey_skip items=1 charged=6.0 key=c_sb_skip
ldrb w0, [x1], #1
ldrb w2, [x3], #1
cmp w0, w2
b.ls .Lenter
subs x4, x4, #1
b.ne .Lloop
```

```asm
; kernel: rev_b2_minkey_enter items=1 charged=20.0 key=c_sb_enter
ldr x5, [x_off, x_j, lsl #3]
ldr x6, [x_off, x_j1, lsl #3]
subs x7, x6, x5
b.eq .Lnext_sb
add x_p, x_ent, x5, lsl #4
add x_pe, x_ent, x6, lsl #4
mov w_lim, w2
ldr x_qtab, [x_sp, #0]
ldp x_bptr, x_blim, [x_sp, #8]
ldp x_base, x_len, [x_sp, #24]
ldp x_lp, x_tmp, [x_sp, #40]
mov w_cnt, #0
cmp x_p, x_pe
b.hs .Lnext_sb
subs x4, x4, #1
b.ne .Lloop
add x_j, x_j, #1
b .Ldone
```

### 1.3 `rev_b2_scan_emit_8b` candidate check (`11` Ir/cand <= `C_CAND = 11.5`)

```asm
; kernel: rev_b2_scan_cand items=1 charged=11.5 key=c_cand
ldp w0, w1, [x_p], #16
cmp w0, w_lim
b.hi .Ldone
sub w2, w1, w_base
cmp w2, w_len
b.hs .Lnext
ldr w3, [x_qtab, x2, lsl #2]
cmp w0, w3
b.hi .Lnext
cmp x_p, x_pe
b.lo .Lcand
```

### 1.4 Band-1a zero-copy `uint16_t` 4x-unrolled direct apply (`rev_b1a_apply_u4`: `16/4 = 4.00` Ir/hit <= `C_HIT_LO_M4 = 4.25`)

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
b.ne .Lb1a
```

### 1.5 Band-1b zero-copy `uint16_t` 4x-unrolled 8B packed `{m:2, oid:58}` direct apply (`rev_b1b_apply_u4`: `24/4 = 6.00` Ir/hit <= `C_HIT_LO_U = 6.25`)

Every Band-1b hit packs `(m << 58) | oid` into one 64-bit word (`V_a < 2^58`, `m in {0..3}` -> 60 active bits <= 64). Rotating bank `m -> (m - m_b) & 3` in the contiguous 4-bank `uint16_t` array `acc[4 * 2^58]` is `sub` + `and #mask60`, directly followed by `ldrh/add/strh` with ZERO intermediate copy:

```asm
; kernel: rev_b1b_apply_u4 items=4 charged=6.25 key=c_hit_lo_u
ldp x0, x1, [x_p], #16
ldp x2, x3, [x_p], #16
sub x4, x0, x_mb58
and x4, x4, x_mask60
ldrh w8, [x_acc, x4, lsl #1]
add w8, w8, w_lp
strh w8, [x_acc, x4, lsl #1]
sub x5, x1, x_mb58
and x5, x5, x_mask60
ldrh w9, [x_acc, x5, lsl #1]
add w9, w9, w_lp
strh w9, [x_acc, x5, lsl #1]
sub x6, x2, x_mb58
and x6, x6, x_mask60
ldrh w10, [x_acc, x6, lsl #1]
add w10, w10, w_lp
strh w10, [x_acc, x6, lsl #1]
sub x7, x3, x_mb58
and x7, x7, x_mask60
ldrh w11, [x_acc, x7, lsl #1]
add w11, w11, w_lp
strh w11, [x_acc, x7, lsl #1]
subs x_cnt, x_cnt, #4
b.ne .Lb1b
```

### 1.6 Band-2 zero-copy emit + direct `uint16_t` apply (`8 + 4.00 = 12.00` Ir/hit <= `C_HIT_HI = 12.5`)

Emit appends the 8B `oid` directly to `bucket[t]` (`8` opcodes), and row processing applies `bucket[t]` directly via `rev_b1a_apply_u4` (`4.00` opcodes/hit, ZERO copy):

```asm
; kernel: rev_b2_emit items=1 charged=8.5 key=c_hit_hi_emit
ldur x4, [x_p, #-8]
add x5, x_bptr, x2, lsl #3
ldr x6, [x5]
str x4, [x6], #8
str x6, [x5]
cmp x6, x_blim
b.hs .Lflush
b .Lnext
```

### 1.7 8-Cell `uint16_t` 64-Bit Borrow-Free SWAR Scan + Zero (`rev_scan_zero_u16`: `14/8 = 1.75` Ir/cell -> `1.75 * 1.0276 + 0.15 = 1.948 <= C_SCAN_ZERO = 2.00`)

Each 64-bit word holds 4 `uint16_t` quarter-bit accumulators (`<= 2,800 << 32,768`). For any threshold `T in [1..32767]`, set `bias16 = (32768 - T) & 0x7fff` broadcast across 4 16-bit lanes (`x_bias16`), `x_7fff = 0x7fff_7fff_7fff_7fff`, `x_8000 = 0x8000_8000_8000_8000`. Then for EVERY `uint16_t` lane `x in [0..65535]`: `x >= T` iff `(((x & 0x7fff) + bias16) | x) & 0x8000 != 0` (verified LIVE across full threshold sweep & boundary edges inside `audit_cm.py`, ZERO borrow between 16-bit lanes!). Processing 8 `uint16_t` cells (`2` 64-bit words) per loop iteration with unconditional zeroing of `S_g` and cold-path zeroing of `S_f`:

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

### 1.8 2-Pass `a_0`-Chord Residue Build (`rev_build_pass12`: `28` opcodes + `0.50` chord-start amort. = `28.50 <= C_LIST_A = 30.0`, `4` Ir/SB <= `C_LBKT = 4.0`)

Every Band-2 16-byte entry uses the exact packed little-endian layout `struct.pack("<IIQ", q_key, ka_lo, oid_m)`: bytes `[0..3] = q_key` (`uint32_t`, top byte `q8 = q_key >> 24`), bytes `[4..7] = ka_lo = ka & 0xffffffff` (`uint32_t`), bytes `[8..15] = oid_m = (m << 58) | oid` (`uint64_t`). Because `sigma(1) = 1` fixes `a_0` on all 4 `sigma`-conjugates, Pass 1 walks `a_0`-chords (mean length `L_a = 233.88`), advancing `ka = (ka + r^{-1}) mod p` (`4` opcodes; chord-start dot product `60 / 233.88 = 0.26 <= 0.50` Ir/entry), storing `ka` to `ka_by_pos[pos]`, and incrementing `hist[ka >> 6]` (`12` opcodes). Pass 2 walks `pos = 0..N_a_tbl-1` in global `key_a`-sorted order, loads `(q_key, oid_m)` via legal immediate post-index `ldp x7, x3, [x_meta], #16`, packs `x8 = q_key | (ka << 32)`, writes the 16-byte entry via `add x_tmp, x_out, x2` + legal zero-immediate `stp x8, x3, [x_tmp]` (matching `rev_b2_scan_cand`'s `ldp w0, w1, [x_p], #16` on bytes `[0..7]` and `rev_b2_emit`'s `ldur x4, [x_p, #-8]` on bytes `[8..15]`!), and updates `sb_min_q8[sbi]` (`16` opcodes -> `12 + 16 = 28` opcodes + `0.50` = **`28.50 <= 30.0`** Ir/entry; builder/scanner/emit 16-byte layout agreement verified LIVE in `audit_cm.py`):

```asm
; kernel: rev_build_pass12 items=1 charged=29.5 key=c_list_a_raw
add x0, x0, x_step
sub x5, x0, x_p
cmp x0, x_p
csel x0, x5, x0, hs
ldr x_pos, [x_imap], #8
str x0, [x_kbuf, x_pos, lsl #3]
lsr x1, x0, #6
ldr w2, [x_hist, x1, lsl #2]
add w2, w2, #1
str w2, [x_hist, x1, lsl #2]
subs x_n1, x_n1, #1
b.ne .Lp1
ldr x0, [x_kbuf], #8
ldp x7, x3, [x_meta], #16
orr x8, x7, x0, lsl #32
lsr x1, x0, #6
ldr x2, [x_pos_sb, x1, lsl #3]
add x4, x2, #16
str x4, [x_pos_sb, x1, lsl #3]
add x_tmp, x_out, x2
stp x8, x3, [x_tmp]
lsr w5, w7, #24
ldrb w6, [x_minq, x1]
cmp w5, w6
csel w6, w5, w6, lo
strb w6, [x_minq, x1]
subs x_n2, x_n2, #1
b.ne .Lp2
```

```asm
; kernel: rev_build_sb_min items=1 charged=4.0 key=c_lbkt
ldrb w0, [x_minq], #1
cmp w0, w_cur
csel w_cur, w0, w_cur, lo
strb w_cur, [x_outq], #1
```

Band-1b lookup (`rev_b1_early4_find`, `[2^48, 2^55)`): empty early-exit `keys[0] > neg_qb` = `11` Ir (`ldr, add, ldr, cmp, b.gt` + 6 header/return); `+14 + cnt` 4-strided scan when `cnt > 0`; binary-search fallback `24 + 8 log2(mu + 1)`. Exactness of 4-orbit `sigma`-shared replay with `D_(omega,sigma)` (`239/239` exact-1), per-SB `q_min_j`, and 1-byte `sb_min_q8` gate is verified live in `audit_cm.py` (`sigmaB_exactness`: 48 real BLS12-381 `w' = omega^2` ideals, `0` mismatches).

## 2. Exact 5-Limb Baseline Scalar ARMv8-A CIOS Montgomery Multiplication (`MONT5 = 270.0`)

One complete carry-exact, shift-included two-chain CIOS outer round (`ldr bi` + Pass A low halves `a[j]*bi` + Pass B high halves `a[j]*bi` + `m = t0*n0inv` + Pass C low halves `m*p[j]` shifting `t[1..6] -> t[0..5]` in place + Pass D high halves `m*p[j]`) executes **`48`** baseline scalar ARMv8-A opcodes (`5 x 48 = 240` across 5 unrolled rounds) + **`27`** load/zero-init/final-sub/store opcodes = **`267 <= 270.0` Ir** (`+3` Ir headroom; bit-exactness of this exact 267-opcode register machine is verified LIVE on 32 random 320-bit `(a, b, p)` triples inside `audit_cm.py`):

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

## 3. Exact Cyclotomic Tower Norm Kernels, 4-Prime 62-Bit NTT (`fft_F`), Analytical Stage B, 2D Cofactorization, & Exact 2nd-Moment ECM Chebyshev

- **Tower `Q(zeta_15) -> Q(omega, sqrt(5)) -> Q(omega) -> Q`:** `_norm_z15_2limb` = `3,357` Ir; `_nk_g_only` = `4,204 <= 4,500` Ir; `_nk_f_only` = `6,167 <= 6,520` Ir.
- **Closed-form non-recursive 4-prime 62-bit radix-2 NTT (`cm_models.fft_F`):** `p_i = c_i * 2^53 + 1 < 2^62` for `c_i in {269, 276, 305, 309}` (`log2(NTT4_PROD) = 244.70`, all 4 verified prime LIVE). 1D (`2 N <= 2^53`): `168 M log2(M) + 220 M + 1200 + (504 M if N > 2^26 else 0)`. 2D (`2 N > 2^53`): `168 M_tot log2(M_tot) + 1100 M_tot + 1200`. `srt_toy` verifies exact scaled-remainder-tree remaindering on 128 planted leaves (`0` mismatches).
- **Analytical Stage B (`cm_models.stage_b_analytic`):** `gbatch_cost(24197143, L) + 10 * MONT5 * (L/5)^2 + 2 * nb * (12 L + 16)` with `MONT5 = 270.0`.
- **Live 2D-Histogram Cofactorization (`Q_F_SAFETY = 2.0`):** Under `H1`, evaluating `2.0 * rho_exact((bf + alpha_f) / lf)` at the left endpoint `bf` of each 4-bit bin `[bf, bf + 4]` gives `>= 2.27x` the bin-midpoint de Bruijn smoothness density `rho_db` (matching the 6x Yukon-accepted `Q_F_SAFETY = 2.0`). Integrated LIVE across all 337 bins up to `bg_max = 700`b, giving **`c_cof_2d_upper = 2.8407`** ipc.
- **Exact Conditional-Dependence 2nd-Moment ECM Chebyshev Bound (`cm_models.ecm_mean_var`):** Each schedule runs `K_c = curves` independent curves (`p_1` per-curve success, schedule failure `q_s = (1 - p_1)^{K_c} <= 0.001`, schedule success `p_s = 1 - q_s >= 0.999`). Failed schedules have `0` curve successes and deterministic cost `C_fail = K_c * per_c` (`Var = 0`). The final successful schedule has `S ~ Bin(K_c, p_1) | S >= 1` successes and cost `C_succ = C_fail + S * redo_ex`, with exact truncated-binomial moments `E[S | S >= 1] = K_c p_1 / p_s`, `E[S^2 | S >= 1] = (K_c p_1 (1 - p_1) + (K_c p_1)^2) / p_s`, `Var(C_succ) = (E[S^2 | S >= 1] - E[S | S >= 1]^2) * redo_ex^2`. Since the number of failed schedules `F ~ Geom_0(p_s)` before first success is strictly independent of `C_succ`, total work `T = F * C_fail + C_succ` has exact `E[T] = (q_s / p_s) C_fail + E[C_succ]` and exact `Var(T) = (q_s / p_s^2) C_fail^2 + Var(C_succ)`. Allowing worst-case intra-relation correlation via Cauchy-Schwarz (`sd_rel = 13*sqrt(Var_f) + 16*sqrt(Var_g)`) across independent relations yields exact Chebyshev tail `< 2^-74.4 << 2^-70`.

## 4. `Z[zeta_12]` Block Wiedemann Kernels (Straight-Line 4-Limb Scalar ARMv8-A)

Over `Z/rZ` (`r` is 255-bit -> 4 `uint64_t` limbs, `phi(12) = 4` coordinates in `Z[zeta_12]`):
- **`unit_run5_z12` (`44` opcodes / `5` entries = `8.80` Ir/vec-entry):** Unrolled 5-entry run adding 4-limb entries into 5-limb accumulator (`5 x (2 ldp + 4 adds/adcs + 2 carry/ptr) + 4 loop/control = 44` opcodes -> `8.80` Ir/entry):

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

- **`gen_run1_z12` (`25.0 <= 33.0` Ir/vec-entry):** `2 ldp` + `4 x (mul, umulh, adds, adcs)` (16) + carry/loop (7) = `25.0` Ir (charged `33.0`).
- **`chi_mul_z12_fast` (`50` Ir), `row_red5_z12` (`53` Ir), `row_tail_6bucket_z12` (`504` Ir), 1-Owner Arch B Row Overhead (`785` Ir):** `504 + 265 + 16 = 785` Ir worst-case (`0` cross-owner sends/combines; `k5_live = eu_lo * 8.80 + (1 - eu_lo) * 33.0 = 10.494 <= k5_charged = 11.50`; `red_row_u_z12 = 18.6` Ir/elem).

## 5. Descent Kernels (`16D` Rank-2 `dsieve16`, Explicit 381-Bit & 265-Bit Bareiss LLL Bounds, & Limb-Scaled 100-Round Miller-Rabin)

- **Analytical `dsieve16` counts (`sieve_mult = 1.0`, each op already carrying per-instruction headroom):** `SIEVE_HIT = 1,100`, `SCAN + ZERO = 12`, `BOUND_CELL = 2,400`, `PAIR_SETUP = 1.33 x 2^25.1`, `NORMS = 2^27`, `TDIV = 2^25`, `IDEAL = 2^35`, `ecm_instr` (`MONT5 = 270.0`), and exact limb-scaled 100-round Miller-Rabin **`mr_instr(nb) = 100 * 2 * nb * MONT5 * (ceil(nb/64)/5)^2`** (`<= 4^(-100) = 2^(-200)` error/test, `<= 2^(-80)` over `< 2^120` run calls).
- **Explicit input-uniform 16D Bareiss LLL limb + swap bounds (`lll_16d_int`, verified LIVE in `audit_cm.py`):**
  1. **Special-`q` lattice (`q <= 2^265`, `det = q`, `LLL_NODE = 2^38` Ir):** By Hadamard's inequality on `B_sq`, the initial Gram-determinant potential satisfies `log2(D_init) <= 16 * 530 = 8,480` bits -> input-uniform worst-case swaps `S_wc <= 8,480 / log2(4/3) = 20,432`; uniform de Weger/Bareiss bound `|d_k|, |lam_{i,j}| <= D_max` has `L_wc <= 10` limbs (640 bits), with explicit `(16`-coord + `15`-coeff`)` size-reduction and swap update `<= 16 * 80 * L_wc^2 = 128,000` scalar ARMv8-A opcodes/swap -> input-uniform bound **`<= 20,432 * 128,000 = 2^31.28 << LLL_NODE = 2^38` Ir** (verified live on a genuine 265-bit special-`q` lattice: `L = 9` limbs, `2^27.2 << 2^38`).
  2. **381-bit `G_T` Boot lattice (`det = p^8 * S^16 = 2^3384`, `LLL_BOOT = 2^42` Ir):** Two-stage `D55` 280-digit (`15`-limb) pre-reduction (`S_dec_wc <= 91,440 / log2(1/0.80) = 284,043` swaps x `250,000` Ir/swap = `2^36.05`) followed by exact `delta = 99/100` integer Bareiss LLL (input-uniform `L_wc <= 64` limbs after `D55`, post-`D55` potential `<= 1,563` bits -> `S_dw_wc <= 107,800` swaps x `16 * 80 * 64^2 = 5.24e6` Ir/swap = `2^39.04`) -> input-uniform bound **`<= 2^39.21 << LLL_BOOT = 2^42` Ir** (verified live on all 8/8 genuine production `G_T` challenges: `L = 53` limbs, `~2^30.8 << 2^42`, and every opened run is also bounded pointwise by the unified hard work meter `t_run = 2^123.0`).
