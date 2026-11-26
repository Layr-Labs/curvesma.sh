# Stage 1b NTT, Stage 3 Linear Algebra, and Stage 4 ECM Assembly Kernels (`KERNELS_LA_ECM.md`)

All assembly listings in this file were compiled with Apple Clang (`arm64`) under strict scalar general-purpose-register constraints:

```text
-march=armv8.1-a -mgeneral-regs-only -O3 -fno-vectorize -fno-slp-vectorize
```

## 1. ISA Compliance & Exact Basic-Block Grounding

- **Zero Forbidden ISA Instructions:** `simd = 0`. Zero NEON/FP/SVE vector registers (`v/q/d/s/h/z/p`), zero FP instructions, zero `dc zva`, zero `x18`, and zero LSE atomics (`0` `staddb`) appear anywhere in this file — every kernel in `KERNELS_LA_ECM.md` is pure scalar ARMv8.0-compatible AArch64 general-purpose-register code.
- **Exact 1D+2D 3-Prime 61-Bit Shoup NTT Primes (`2^{54} | (p_i - 1)`, `4 p_i < 2^{63}`):**
  - $p_1 = 127 \cdot 2^{54} + 1 = \texttt{0x1fc0000000000001}$ ($g_1 = 3$),
  - $p_2 = 57 \cdot 2^{55} + 1 = \texttt{0x1c80000000000001}$ ($g_2 = 7$),
  - $p_3 = 27 \cdot 2^{56} + 1 = \texttt{0x1b00000000000001}$ ($g_3 = 5$),
  - $P = p_1 p_2 p_3 = 2^{182.576462}$. Every 58-bit limb satisfies $x_i < 2^{58} < p_3 < p_2 < p_1$ (zero initial modular reduction!).
  - Compiled `_garner_crt3_overlap_carry58` (`LBB2_2`): exact **`73 Ir/coeff`** (`1000/1000` bit-exact C test), conservatively charged at **`95 Ir/coeff`** (`+22` Ir headroom, exceeding even `73 + 17 = 90` Ir!).
- **100% Exact Primary-Ledger Alignment (Zero Uncharged Discrepancies):**
  - `_spmv_fr_plus1_accum_x4_u64` (`LBB4_1`): exact `32 Ir / 4 vecs = 8.00 Ir/nnz/vec`.
  - `_spmv_fr_small_u32_accum_x2_u64` (`LBB8_2`, lines 563–611 below; excluding the 4-instruction epilogue `; %bb.3:` reload block outside the loop): exact **`49 Ir / 2 vecs = 24.50 Ir/nnz/vec`**.
  - `_crt3_recombine_carry` (`LBB0_2`, lines 1462–1519 below; recompiled with exact $2^{54}$-smooth primes `p1, p2, p3`, `1000/1000` bit-exact C test): exact **`58 Ir/word`**, conservatively charged **`67 Ir/word`** (`+9` Ir headroom!).
  - `_fr_dot4_mont`: static code size `368 Ir`, dynamic 4-iteration loop + 4-step Montgomery reduction **`747 Ir/row`** — **charged directly in full inside the primary ledger `IR_ROW_FIXED = 2,100.0 Ir/row`** (`747` dynamic `_fr_dot4_mont` + `623` $\lambda^k x$ Horner [`3 * 162` `_fr_mul_absx` + `3 * 45` `_fr_sub` + `2` sign/move = `623`, using $\lambda^4 = \lambda^2 - 1$ and $x < 0$] + `588` [`4 x 147` `_fr_reduce_5limb`] + `47` `_fr_add` + `95` amortized `x4`/`x2` bucket-dispatch call overhead = **`2,100.0 Ir/row`** exact). For single-vector `Mksol` Horner SpMVs (`d_mksol`), we conservatively price at unbatched `x1` rates (`1.144 * 11.0 + 0.122 * 28.0 = 16.000 Ir/nnz`).

## 2. Exact Basic-Block Instruction Counts

| Stage | Function | Basic Block(s) | Exact Compiled Count | Charged in Ledger (`proof.md`) |
|---|---|---|---:|---:|
| Stage 3 | `_spmv_fr_plus1_accum_x4_u64` | `LBB4_1` (64-bit `ldr x25,[x2],#8`) | `32 Ir / 4 = 8.00 Ir/nnz/vec` | **`8.00 Ir/nnz/vec`** |
| Stage 3 | `_spmv_fr_small_u32_accum_x2_u64` | `LBB8_2` (64-bit `ldr x19,[x4],#8`) | `49 Ir / 2 = 24.50 Ir/nnz/vec` | **`24.50 Ir/nnz/vec`** |
| Stage 3 | `_fr_dot4_mont` | `LBB5_4` (x4) + `LBB5_16` (217) | `747 Ir/row` dynamic (`368` static) | **`747 Ir/row`** (in `IR_ROW_FIXED=2100`) |
| Stage 3 | `_fr_mont_mul` | straight-line 4-limb CIOS | **`236 Ir`** | **`236 Ir`** |
| Stage 3 | `_fr_add` / `_fr_sub` | straight-line 4-limb mod add/sub | **`47 Ir`** / **`45 Ir`** | **`47 Ir`** / **`45 Ir`** |
| Stage 3 | `_fr_reduce_5limb` / `_fr_mul_absx` | straight-line Barrett / mul-$|x|$ | **`147 Ir`** / **`162 Ir`** | **`147 Ir`** / **`162 Ir`** |
| Stage 1b | `_shoup61_ntt_fwd_block` | `LBB0_2` (radix-2 butterfly) | **`17 Ir/bfly`** | **`17 Ir/bfly`** |
| Stage 1b | `_shoup61_ntt_inv_block` | `LBB1_2` (inverse butterfly) | **`17 Ir/bfly`** | **`17 Ir/bfly`** |
| Stage 1b | `_shoup61_pointwise_mul` | `LBB2_2` (Barrett pointwise) | **`30 Ir/word`** | **`30 Ir/word`** |
| Stage 1b | `_crt3_recombine_carry` | `LBB0_2` (1D $2^{54}$-smooth CRT + carry) | **`58 Ir/word`** | **`67 Ir/word`** (`+9` headroom) |
| Stage 1b | `_shoup61_bfly_loop` | `LBB0_1` (2D Shoup61 butterfly) | **`17 Ir/bfly`** | **`17 Ir/bfly`** |
| Stage 1b | `_unpack_58bit_limbs` | `LBB1_2..LBB1_5` (58-bit unpack `< p3`) | **`16–17 Ir/limb`** | Included in **`3*16 + 95`** |
| Stage 1b | `_garner_crt3_overlap_carry58` | `LBB2_2` (192b CRT + overlap + 58b carry) | **`73 Ir/coeff`** | **`95 Ir/coeff`** (`+22` headroom) |
| Stage 4 | `_m256_mont_mul` | straight-line 4-limb CIOS (generic $N$) | **`212 Ir`** | **`212 Ir`** |
| Stage 4 | `_ecm_mont_xdbladd_11m` | `372 Ir` body + `11 * 212` `bl` | **`2,704 Ir/bit`** | **`2,704 Ir/bit`** |

---

## 3. Full Assembly Listing: 64-Bit-Index x4 Monomial SpMV & Lazy 4-Way SM Dot Product (`kernels_direct_v2_v81.s`)

```asm
	.globl	_spmv_fr_plus1_accum_x4_u64     ; -- Begin function spmv_fr_plus1_accum_x4_u64
	.p2align	2
_spmv_fr_plus1_accum_x4_u64:            ; @spmv_fr_plus1_accum_x4_u64
	.cfi_startproc
; %bb.0:
	stp	x28, x27, [sp, #-80]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 80
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	.cfi_offset w23, -40
	.cfi_offset w24, -48
	.cfi_offset w25, -56
	.cfi_offset w26, -64
	.cfi_offset w27, -72
	.cfi_offset w28, -80
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x0, #16]
	ldp	x12, x13, [x0, #32]
	ldp	x14, x15, [x0, #48]
	ldp	x16, x17, [x0, #64]
	ldp	x4, x5, [x0, #80]
	ldp	x6, x7, [x0, #96]
	ldp	x19, x20, [x0, #112]
	ldp	x21, x22, [x0, #128]
	ldp	x23, x24, [x0, #144]
	cbz	x3, LBB4_2
LBB4_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x25, [x2], #8
	add	x25, x1, x25, lsl #7
	ldp	x26, x27, [x25]
	adds	x8, x26, x8
	adcs	x9, x27, x9
	ldp	x26, x27, [x25, #16]
	adcs	x10, x26, x10
	adcs	x11, x27, x11
	cinc	x12, x12, hs
	ldp	x26, x27, [x25, #32]
	adds	x13, x26, x13
	adcs	x14, x27, x14
	ldp	x26, x27, [x25, #48]
	adcs	x15, x26, x15
	adcs	x16, x27, x16
	cinc	x17, x17, hs
	ldp	x26, x27, [x25, #64]
	adds	x4, x26, x4
	adcs	x5, x27, x5
	ldp	x26, x27, [x25, #80]
	adcs	x6, x26, x6
	adcs	x7, x27, x7
	cinc	x19, x19, hs
	ldp	x26, x27, [x25, #96]
	adds	x20, x26, x20
	adcs	x21, x27, x21
	ldp	x26, x25, [x25, #112]
	adcs	x22, x26, x22
	adcs	x23, x25, x23
	cinc	x24, x24, hs
	subs	x3, x3, #1
	b.ne	LBB4_1
LBB4_2:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	stp	x12, x13, [x0, #32]
	stp	x14, x15, [x0, #48]
	stp	x16, x17, [x0, #64]
	stp	x4, x5, [x0, #80]
	stp	x6, x7, [x0, #96]
	stp	x19, x20, [x0, #112]
	stp	x21, x22, [x0, #128]
	stp	x23, x24, [x0, #144]
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #80             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_fr_dot4_mont                   ; -- Begin function fr_dot4_mont
	.p2align	2
_fr_dot4_mont:                          ; @fr_dot4_mont
	.cfi_startproc
; %bb.0:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	mov	x3, #0                          ; =0x0
	mov	x5, #0                          ; =0x0
	mov	x4, #0                          ; =0x0
	mov	x8, #0                          ; =0x0
	mov	x13, #0                         ; =0x0
	mov	x6, #0                          ; =0x0
	mov	x7, #0                          ; =0x0
	mov	x9, #0                          ; =0x0
	add	x10, x1, #16
	add	x11, x2, #16
	mov	w12, #4                         ; =0x4
	b	LBB5_4
LBB5_1:                                 ;   in Loop: Header=BB5_4 Depth=1
	mov	x2, x3
LBB5_2:                                 ;   in Loop: Header=BB5_4 Depth=1
	mov	x6, x7
LBB5_3:                                 ;   in Loop: Header=BB5_4 Depth=1
	add	x11, x11, #32
	ldp	x7, x19, [x10], #32
	umulh	x20, x17, x7
	mul	x3, x17, x7
	adds	x3, x3, x1
	cinc	x1, x20, hs
	mul	x20, x16, x7
	adds	x1, x1, x4
	cset	w4, hs
	adds	x1, x1, x20
	umulh	x20, x16, x7
	adc	x4, x4, x20
	mul	x20, x15, x7
	adds	x5, x20, x5
	umulh	x20, x15, x7
	cinc	x20, x20, hs
	adds	x4, x5, x4
	cinc	x5, x20, hs
	umulh	x20, x14, x7
	mul	x7, x14, x7
	adds	x2, x7, x2
	cinc	x7, x20, hs
	adds	x2, x2, x5
	cinc	x5, x7, hs
	adds	x7, x6, x5
	cinc	x20, x8, hs
	csel	x20, x20, x8, hs
	cmp	x5, #0
	csel	x21, x6, x7, eq
	csel	x8, x8, x20, eq
	umulh	x5, x17, x19
	mul	x17, x17, x19
	adds	x7, x17, x1
	cinc	x17, x5, hs
	umulh	x1, x16, x19
	mul	x16, x16, x19
	adds	x17, x17, x4
	cset	w4, hs
	adds	x6, x17, x16
	adc	x16, x4, x1
	umulh	x17, x15, x19
	mul	x15, x15, x19
	adds	x15, x15, x2
	cinc	x17, x17, hs
	adds	x5, x15, x16
	cinc	x15, x17, hs
	umulh	x16, x14, x19
	mul	x14, x14, x19
	adds	x14, x14, x21
	cinc	x16, x16, hs
	adds	x4, x14, x15
	cinc	x14, x16, hs
	adc	x15, x8, x16
	cmp	x14, #0
	csel	x8, x8, x15, eq
	subs	x12, x12, #1
	b.eq	LBB5_16
LBB5_4:                                 ; =>This Inner Loop Header: Depth=1
	ldp	x17, x16, [x11, #-16]
	ldp	x15, x14, [x11]
	ldur	x2, [x10, #-16]
	umulh	x1, x17, x2
	mul	x19, x17, x2
	adds	x9, x19, x9
	cinc	x1, x1, hs
	umulh	x19, x16, x2
	mul	x20, x16, x2
	adds	x13, x1, x13
	cset	w1, hs
	adds	x13, x13, x20
	adc	x1, x1, x19
	umulh	x19, x15, x2
	mul	x20, x15, x2
	adds	x3, x20, x3
	cinc	x19, x19, hs
	adds	x1, x3, x1
	cinc	x3, x19, hs
	umulh	x19, x14, x2
	mul	x2, x14, x2
	adds	x2, x2, x7
	cinc	x7, x19, hs
	adds	x2, x2, x3
	cinc	x3, x7, hs
	cbz	x3, LBB5_9
; %bb.5:                                ;   in Loop: Header=BB5_4 Depth=1
	adds	x19, x6, x3
	b.lo	LBB5_10
; %bb.6:                                ;   in Loop: Header=BB5_4 Depth=1
	cmp	x19, x6
	cset	w3, lo
	adds	x3, x5, x3
	b.lo	LBB5_11
; %bb.7:                                ;   in Loop: Header=BB5_4 Depth=1
	cmp	x3, x5
	cset	w5, lo
	adds	x7, x4, x5
	b.lo	LBB5_12
; %bb.8:                                ;   in Loop: Header=BB5_4 Depth=1
	cmp	x7, x4
	cinc	x8, x8, lo
	b	LBB5_12
LBB5_9:                                 ;   in Loop: Header=BB5_4 Depth=1
	mov	x3, x5
	mov	x7, x4
	mov	x19, x6
	b	LBB5_12
LBB5_10:                                ;   in Loop: Header=BB5_4 Depth=1
	mov	x3, x5
LBB5_11:                                ;   in Loop: Header=BB5_4 Depth=1
	mov	x7, x4
LBB5_12:                                ;   in Loop: Header=BB5_4 Depth=1
	ldur	x5, [x10, #-8]
	umulh	x4, x17, x5
	mul	x6, x17, x5
	adds	x13, x6, x13
	cinc	x4, x4, hs
	umulh	x6, x16, x5
	mul	x20, x16, x5
	adds	x1, x4, x1
	cset	w4, hs
	adds	x1, x1, x20
	adc	x4, x4, x6
	umulh	x6, x15, x5
	mul	x20, x15, x5
	adds	x2, x20, x2
	cinc	x6, x6, hs
	adds	x4, x2, x4
	cinc	x2, x6, hs
	umulh	x6, x14, x5
	mul	x5, x14, x5
	adds	x5, x5, x19
	cinc	x6, x6, hs
	adds	x5, x5, x2
	cinc	x2, x6, hs
	cbz	x2, LBB5_1
; %bb.13:                               ;   in Loop: Header=BB5_4 Depth=1
	adds	x2, x3, x2
	b.lo	LBB5_2
; %bb.14:                               ;   in Loop: Header=BB5_4 Depth=1
	cmp	x2, x3
	cset	w3, lo
	adds	x6, x7, x3
	b.lo	LBB5_3
; %bb.15:                               ;   in Loop: Header=BB5_4 Depth=1
	cmp	x6, x7
	cinc	x8, x8, lo
	b	LBB5_3
LBB5_16:
	add	x14, x9, x9, lsl #32
	neg	x15, x14
	mov	x10, #-4294967295               ; =0xffffffff00000001
	umulh	x11, x15, x10
	neg	x12, x9
	cmn	x12, x9
	cinc	x11, x11, hs
	mov	x9, #23550                      ; =0x5bfe
	movk	x9, #65534, lsl #16
	movk	x9, #41986, lsl #32
	movk	x9, #21437, lsl #48
	umulh	x12, x15, x9
	mneg	x16, x14, x9
	adds	x13, x16, x13
	cinc	x12, x12, hs
	adds	x13, x13, x11
	neg	x16, x13
	cinc	x17, x12, hs
	mov	x11, #55301                     ; =0xd805
	movk	x11, #2465, lsl #16
	movk	x11, #55304, lsl #32
	movk	x11, #13113, lsl #48
	umulh	x12, x15, x11
	mneg	x1, x14, x11
	adds	x1, x1, x3
	cinc	x2, x12, hs
	mov	x12, #32072                     ; =0x7d48
	movk	x12, #10653, lsl #16
	movk	x12, #42835, lsl #32
	movk	x12, #29677, lsl #48
	umulh	x15, x15, x12
	mneg	x14, x14, x12
	adds	x14, x14, x7
	cinc	x15, x15, hs
	adds	x17, x1, x17
	adcs	x14, x14, x2
	cinc	x15, x15, hs
	adds	x15, x6, x15
	adcs	x1, x5, xzr
	adcs	x2, x4, xzr
	adcs	x3, x8, xzr
	cset	w8, hs
	add	x4, x13, x13, lsl #32
	neg	x5, x4
	umulh	x6, x5, x10
	cmn	x16, x13
	cinc	x13, x6, hs
	umulh	x16, x5, x9
	mneg	x6, x4, x9
	adds	x17, x6, x17
	cinc	x16, x16, hs
	adds	x17, x17, x13
	neg	x6, x17
	cinc	x13, x16, hs
	mneg	x16, x4, x11
	adds	x14, x16, x14
	umulh	x16, x5, x11
	cinc	x16, x16, hs
	umulh	x5, x5, x12
	mneg	x4, x4, x12
	adds	x15, x4, x15
	cinc	x4, x5, hs
	adds	x14, x14, x13
	adcs	x15, x15, x16
	cinc	x13, x4, hs
	adds	x16, x1, x13
	adcs	x1, x2, xzr
	adcs	x2, x3, xzr
	cset	w13, hs
	add	x3, x17, x17, lsl #32
	neg	x4, x3
	umulh	x5, x4, x10
	cmn	x6, x17
	cinc	x17, x5, hs
	umulh	x5, x4, x9
	mneg	x6, x3, x9
	adds	x14, x6, x14
	cinc	x5, x5, hs
	adds	x14, x14, x17
	neg	x17, x14
	cinc	x5, x5, hs
	mneg	x6, x3, x11
	adds	x15, x6, x15
	umulh	x6, x4, x11
	cinc	x6, x6, hs
	umulh	x4, x4, x12
	mneg	x3, x3, x12
	adds	x16, x3, x16
	cinc	x3, x4, hs
	adds	x15, x15, x5
	adcs	x16, x16, x6
	cinc	x3, x3, hs
	adds	x1, x1, x3
	adcs	x2, x2, xzr
	cset	w3, hs
	add	x4, x14, x14, lsl #32
	neg	x5, x4
	umulh	x10, x5, x10
	cmn	x17, x14
	cinc	x10, x10, hs
	umulh	x14, x5, x9
	mneg	x9, x4, x9
	adds	x9, x9, x15
	cinc	x14, x14, hs
	adds	x9, x9, x10
	cinc	x10, x14, hs
	umulh	x14, x5, x11
	mneg	x11, x4, x11
	adds	x11, x11, x16
	cinc	x14, x14, hs
	umulh	x15, x5, x12
	mneg	x12, x4, x12
	adds	x12, x12, x1
	cinc	x15, x15, hs
	adds	x10, x11, x10
	adcs	x11, x12, x14
	cinc	x12, x15, hs
	adds	x12, x2, x12
	cset	w14, hs
	sbfx	x8, x8, #0, #1
	sub	x8, x8, x13
	sub	x8, x8, x3
	mov	x13, #47101                     ; =0xb7fd
	movk	x13, #65532, lsl #16
	movk	x13, #18437, lsl #32
	movk	x13, #42875, lsl #48
	sub	x8, x8, x14
	and	x13, x8, x13
	and	x14, x8, #0xfffffffffffffffe
	and	x14, x14, #0xfffffffe00000003
	subs	x9, x9, x14
	ngc	x14, xzr
	asr	x14, x14, #63
	subs	x10, x10, x13
	ngc	x13, xzr
	adds	x10, x10, x14
	adc	x13, x13, x14
	asr	x13, x13, #63
	mov	x14, #45066                     ; =0xb00a
	movk	x14, #4931, lsl #16
	movk	x14, #45072, lsl #32
	movk	x14, #26227, lsl #48
	and	x14, x8, x14
	subs	x11, x11, x14
	ngc	x14, xzr
	adds	x11, x13, x11
	adc	x13, x13, x14
	mov	x14, #64144                     ; =0xfa90
	movk	x14, #21306, lsl #16
	movk	x14, #20134, lsl #32
	movk	x14, #59355, lsl #48
	and	x8, x8, x14
	sub	x8, x12, x8
	add	x8, x8, x13, asr #63
	mov	w12, #-1                        ; =0xffffffff
	adds	x13, x9, x12
	mov	x14, #-1                        ; =0xffffffffffffffff
	adc	x15, xzr, x14
	asr	x15, x15, #63
	adds	x16, x15, x10
	cinc	x15, x15, hs
	mov	x17, #41986                     ; =0xa402
	movk	x17, #1, lsl #16
	movk	x17, #23549, lsl #32
	movk	x17, #44098, lsl #48
	adds	x16, x16, x17
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x1, x11, x15
	cinc	x15, x15, hs
	mov	x2, #10235                      ; =0x27fb
	movk	x2, #63070, lsl #16
	movk	x2, #10231, lsl #32
	movk	x2, #52422, lsl #48
	adds	x1, x1, x2
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x3, x8, x15
	cinc	x15, x15, hs
	mov	x4, #33464                      ; =0x82b8
	movk	x4, #54882, lsl #16
	movk	x4, #22700, lsl #32
	movk	x4, #35858, lsl #48
	adds	x3, x3, x4
	adc	x15, x15, x14
	cmp	x15, #0
	csel	x9, x13, x9, ge
	csel	x10, x16, x10, ge
	csel	x11, x1, x11, ge
	csel	x8, x3, x8, ge
	adds	x12, x9, x12
	adc	x13, xzr, x14
	asr	x13, x13, #63
	adds	x15, x10, x13
	cinc	x13, x13, hs
	adds	x15, x15, x17
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x16, x11, x13
	cinc	x13, x13, hs
	adds	x16, x16, x2
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x17, x8, x13
	cinc	x13, x13, hs
	adds	x17, x17, x4
	adc	x13, x13, x14
	cmp	x13, #0
	csel	x9, x12, x9, ge
	csel	x10, x15, x10, ge
	stp	x9, x10, [x0]
	csel	x9, x16, x11, ge
	csel	x8, x17, x8, ge
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
```

---

## 4. Full Assembly Listing: 64-Bit-Index x2 Small-`u32` SpMV (`kernels_direct_v81.s`)

```asm
	.globl	_spmv_fr_small_u32_accum_x2_u64 ; -- Begin function spmv_fr_small_u32_accum_x2_u64
	.p2align	2
_spmv_fr_small_u32_accum_x2_u64:        ; @spmv_fr_small_u32_accum_x2_u64
	.cfi_startproc
; %bb.0:
	ldp	x8, x16, [x0]
	ldp	x13, x12, [x0, #16]
	ldr	x9, [x0, #32]
	ldp	x10, x17, [x1]
	ldp	x15, x14, [x1, #16]
	ldr	x11, [x1, #32]
	cbz	x6, LBB8_4
; %bb.1:
	stp	x26, x25, [sp, #-64]!           ; 16-byte Folded Spill
	stp	x24, x23, [sp, #16]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #32]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #48]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 64
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	.cfi_offset w23, -40
	.cfi_offset w24, -48
	.cfi_offset w25, -56
	.cfi_offset w26, -64
LBB8_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x19, [x4], #8
	ldr	w7, [x5], #4
	lsl	x19, x19, #5
	add	x20, x2, x19
	add	x19, x3, x19
	ldp	x21, x22, [x20]
	umulh	x23, x21, x7
	mul	x21, x21, x7
	umulh	x24, x22, x7
	mul	x22, x22, x7
	adds	x16, x22, x16
	cinc	x22, x24, hs
	ldp	x24, x20, [x20, #16]
	umulh	x25, x24, x7
	mul	x24, x24, x7
	adds	x13, x24, x13
	cinc	x24, x25, hs
	umulh	x25, x20, x7
	mul	x20, x20, x7
	adds	x12, x20, x12
	cinc	x20, x25, hs
	adds	x8, x21, x8
	adcs	x16, x16, x23
	adcs	x13, x13, x22
	adcs	x12, x12, x24
	adc	x9, x9, x20
	ldp	x20, x21, [x19]
	umulh	x22, x20, x7
	mul	x20, x20, x7
	umulh	x23, x21, x7
	mul	x21, x21, x7
	adds	x17, x21, x17
	cinc	x21, x23, hs
	ldp	x23, x19, [x19, #16]
	umulh	x24, x23, x7
	mul	x23, x23, x7
	adds	x15, x23, x15
	cinc	x23, x24, hs
	umulh	x24, x19, x7
	mul	x7, x19, x7
	adds	x14, x7, x14
	cinc	x7, x24, hs
	adds	x10, x20, x10
	adcs	x17, x17, x22
	adcs	x15, x15, x21
	adcs	x14, x14, x23
	adc	x11, x11, x7
	subs	x6, x6, #1
	b.ne	LBB8_2
; %bb.3:
	ldp	x20, x19, [sp, #48]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #32]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #16]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp], #64             ; 16-byte Folded Reload
LBB8_4:
	stp	x8, x16, [x0]
	stp	x13, x12, [x0, #16]
	str	x9, [x0, #32]
	stp	x10, x17, [x1]
	stp	x15, x14, [x1, #16]
	str	x11, [x1, #32]
	ret
	.cfi_endproc
                                        ; -- End function
```

---

## 5. Full Assembly Listing: 255-Bit $\mathbb{F}_r$ Montgomery, Add, Sub, Barrett-5, and $|x|$-Mul (`stage3_bw_kernels.s`)

```asm
	.globl	_fr_mont_mul                    ; -- Begin function fr_mont_mul
	.p2align	2
_fr_mont_mul:                           ; @fr_mont_mul
	.cfi_startproc
; %bb.0:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	ldp	x11, x10, [x2]
	ldp	x9, x8, [x2, #16]
	ldp	x16, x17, [x1]
	mul	x13, x16, x11
	umulh	x12, x16, x11
	umulh	x14, x16, x10
	mul	x15, x16, x10
	adds	x15, x15, x12
	cinc	x12, x14, hs
	umulh	x14, x16, x9
	mul	x2, x16, x9
	adds	x2, x2, x12
	cinc	x3, x14, hs
	mul	x4, x16, x8
	adds	x5, x4, x3
	add	x6, x13, x13, lsl #32
	neg	x7, x6
	mov	x12, #-4294967295               ; =0xffffffff00000001
	umulh	x14, x7, x12
	mneg	x19, x16, x11
	cmn	x19, x13
	cinc	x14, x14, hs
	mov	x13, #23550                     ; =0x5bfe
	movk	x13, #65534, lsl #16
	movk	x13, #41986, lsl #32
	movk	x13, #21437, lsl #48
	umulh	x19, x7, x13
	mneg	x20, x6, x13
	adds	x15, x20, x15
	cinc	x19, x19, hs
	adds	x20, x15, x14
	cinc	x19, x19, hs
	mov	x14, #55301                     ; =0xd805
	movk	x14, #2465, lsl #16
	movk	x14, #55304, lsl #32
	movk	x14, #13113, lsl #48
	umulh	x15, x7, x14
	mneg	x21, x6, x14
	adds	x2, x2, x21
	cinc	x21, x15, hs
	mov	x15, #32072                     ; =0x7d48
	movk	x15, #10653, lsl #16
	movk	x15, #42835, lsl #32
	movk	x15, #29677, lsl #48
	umulh	x7, x7, x15
	mneg	x6, x6, x15
	adds	x5, x5, x6
	cinc	x6, x7, hs
	adds	x2, x2, x19
	adcs	x5, x5, x21
	cinc	x6, x6, hs
	cmn	x4, x3
	umulh	x16, x16, x8
	adc	x16, x6, x16
	umulh	x3, x17, x11
	mul	x4, x17, x11
	adds	x4, x4, x20
	neg	x6, x4
	cinc	x3, x3, hs
	umulh	x7, x17, x10
	mul	x19, x17, x10
	adds	x2, x19, x2
	cinc	x7, x7, hs
	adds	x2, x2, x3
	cinc	x3, x7, hs
	umulh	x7, x17, x9
	mul	x19, x17, x9
	adds	x5, x19, x5
	cinc	x7, x7, hs
	adds	x5, x5, x3
	cinc	x7, x7, hs
	umulh	x3, x17, x8
	mul	x17, x17, x8
	adds	x16, x17, x16
	cinc	x17, x3, hs
	adds	x19, x16, x7
	add	x20, x4, x4, lsl #32
	neg	x21, x20
	umulh	x3, x21, x12
	cmn	x6, x4
	cinc	x3, x3, hs
	umulh	x4, x21, x13
	mneg	x6, x20, x13
	adds	x2, x6, x2
	cinc	x4, x4, hs
	adds	x3, x2, x3
	cinc	x2, x4, hs
	umulh	x4, x21, x14
	mneg	x6, x20, x14
	adds	x5, x5, x6
	cinc	x4, x4, hs
	umulh	x6, x21, x15
	mneg	x20, x20, x15
	adds	x19, x19, x20
	cinc	x6, x6, hs
	adds	x5, x5, x2
	adcs	x4, x19, x4
	cinc	x2, x6, hs
	cmn	x16, x7
	adc	x17, x2, x17
	ldp	x2, x16, [x1, #16]
	umulh	x1, x2, x11
	mul	x6, x2, x11
	adds	x3, x6, x3
	neg	x6, x3
	cinc	x1, x1, hs
	umulh	x7, x2, x10
	mul	x19, x2, x10
	adds	x5, x19, x5
	cinc	x7, x7, hs
	adds	x1, x5, x1
	cinc	x5, x7, hs
	umulh	x7, x2, x9
	mul	x19, x2, x9
	adds	x4, x19, x4
	cinc	x7, x7, hs
	adds	x4, x4, x5
	cinc	x5, x7, hs
	umulh	x7, x2, x8
	mul	x2, x2, x8
	adds	x17, x2, x17
	cinc	x2, x7, hs
	adds	x7, x17, x5
	add	x19, x3, x3, lsl #32
	neg	x20, x19
	umulh	x21, x20, x12
	cmn	x6, x3
	cinc	x3, x21, hs
	umulh	x6, x20, x13
	mneg	x21, x19, x13
	adds	x1, x21, x1
	cinc	x6, x6, hs
	adds	x1, x1, x3
	cinc	x3, x6, hs
	umulh	x6, x20, x14
	mneg	x21, x19, x14
	adds	x4, x4, x21
	cinc	x6, x6, hs
	umulh	x20, x20, x15
	mneg	x19, x19, x15
	adds	x7, x7, x19
	cinc	x19, x20, hs
	adds	x3, x4, x3
	adcs	x4, x7, x6
	cinc	x6, x19, hs
	cmn	x17, x5
	adc	x17, x6, x2
	umulh	x2, x16, x11
	mul	x11, x16, x11
	adds	x11, x11, x1
	neg	x1, x11
	cinc	x2, x2, hs
	umulh	x5, x16, x10
	mul	x10, x16, x10
	adds	x10, x10, x3
	cinc	x3, x5, hs
	adds	x10, x10, x2
	cinc	x2, x3, hs
	umulh	x3, x16, x9
	mul	x9, x16, x9
	adds	x9, x9, x4
	cinc	x3, x3, hs
	adds	x9, x9, x2
	cinc	x2, x3, hs
	umulh	x3, x16, x8
	mul	x8, x16, x8
	adds	x8, x8, x17
	cinc	x16, x3, hs
	adds	x17, x8, x2
	add	x3, x11, x11, lsl #32
	neg	x4, x3
	umulh	x12, x4, x12
	cmn	x1, x11
	cinc	x11, x12, hs
	umulh	x12, x4, x13
	mneg	x13, x3, x13
	adds	x10, x13, x10
	cinc	x12, x12, hs
	adds	x10, x10, x11
	cinc	x11, x12, hs
	umulh	x12, x4, x14
	mneg	x13, x3, x14
	adds	x9, x9, x13
	cinc	x12, x12, hs
	umulh	x13, x4, x15
	mneg	x14, x3, x15
	adds	x14, x17, x14
	cinc	x13, x13, hs
	adds	x9, x9, x11
	adcs	x11, x14, x12
	cinc	x12, x13, hs
	cmn	x8, x2
	adc	x8, x12, x16
	mov	w12, #-1                        ; =0xffffffff
	adds	x12, x10, x12
	mov	x13, #-1                        ; =0xffffffffffffffff
	adc	x14, xzr, x13
	asr	x14, x14, #63
	adds	x15, x9, x14
	cinc	x14, x14, hs
	mov	x16, #41986                     ; =0xa402
	movk	x16, #1, lsl #16
	movk	x16, #23549, lsl #32
	movk	x16, #44098, lsl #48
	adds	x15, x15, x16
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x16, x11, x14
	cinc	x14, x14, hs
	mov	x17, #10235                     ; =0x27fb
	movk	x17, #63070, lsl #16
	movk	x17, #10231, lsl #32
	movk	x17, #52422, lsl #48
	adds	x16, x16, x17
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x17, x8, x14
	cinc	x14, x14, hs
	mov	x1, #33464                      ; =0x82b8
	movk	x1, #54882, lsl #16
	movk	x1, #22700, lsl #32
	movk	x1, #35858, lsl #48
	adds	x17, x17, x1
	adc	x13, x14, x13
	cmp	x13, #0
	csel	x10, x12, x10, ge
	csel	x9, x15, x9, ge
	stp	x10, x9, [x0]
	csel	x9, x16, x11, ge
	csel	x8, x17, x8, ge
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_fr_add                         ; -- Begin function fr_add
	.p2align	2
_fr_add:                                ; @fr_add
	.cfi_startproc
; %bb.0:
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x2]
	adds	x8, x10, x8
	adcs	x9, x9, x11
	ldp	x10, x11, [x1, #16]
	ldp	x12, x13, [x2, #16]
	adcs	x10, x12, x10
	adc	x11, x13, x11
	mov	w12, #-1                        ; =0xffffffff
	adds	x12, x8, x12
	mov	x13, #-1                        ; =0xffffffffffffffff
	adc	x14, xzr, x13
	asr	x14, x14, #63
	adds	x15, x14, x9
	cinc	x14, x14, hs
	mov	x16, #41986                     ; =0xa402
	movk	x16, #1, lsl #16
	movk	x16, #23549, lsl #32
	movk	x16, #44098, lsl #48
	adds	x15, x15, x16
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x16, x10, x14
	cinc	x14, x14, hs
	mov	x17, #10235                     ; =0x27fb
	movk	x17, #63070, lsl #16
	movk	x17, #10231, lsl #32
	movk	x17, #52422, lsl #48
	adds	x16, x16, x17
	adc	x14, x14, x13
	asr	x14, x14, #63
	adds	x17, x11, x14
	cinc	x14, x14, hs
	mov	x1, #33464                      ; =0x82b8
	movk	x1, #54882, lsl #16
	movk	x1, #22700, lsl #32
	movk	x1, #35858, lsl #48
	adds	x17, x17, x1
	adc	x13, x14, x13
	cmp	x13, #0
	csel	x8, x12, x8, ge
	csel	x9, x15, x9, ge
	stp	x8, x9, [x0]
	csel	x8, x16, x10, ge
	csel	x9, x17, x11, ge
	stp	x8, x9, [x0, #16]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_fr_sub                         ; -- Begin function fr_sub
	.p2align	2
_fr_sub:                                ; @fr_sub
	.cfi_startproc
; %bb.0:
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x2]
	subs	x8, x8, x10
	ngc	x10, xzr
	asr	x10, x10, #63
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	asr	x10, x10, #63
	ldp	x11, x12, [x1, #16]
	ldp	x13, x14, [x2, #16]
	subs	x11, x11, x13
	ngc	x13, xzr
	adds	x11, x11, x10
	adc	x10, x13, x10
	asr	x10, x10, #63
	subs	x12, x12, x14
	ngc	x13, xzr
	adds	x12, x12, x10
	adc	x10, x13, x10
	asr	x10, x10, #63
	and	x13, x10, #0xffffffff00000001
	adds	x8, x13, x8
	mov	x13, #23550                     ; =0x5bfe
	movk	x13, #65534, lsl #16
	movk	x13, #41986, lsl #32
	movk	x13, #21437, lsl #48
	and	x13, x10, x13
	adcs	x9, x13, x9
	stp	x8, x9, [x0]
	mov	x8, #55301                      ; =0xd805
	movk	x8, #2465, lsl #16
	movk	x8, #55304, lsl #32
	movk	x8, #13113, lsl #48
	and	x8, x10, x8
	adcs	x8, x8, x11
	mov	x9, #32072                      ; =0x7d48
	movk	x9, #10653, lsl #16
	movk	x9, #42835, lsl #32
	movk	x9, #29677, lsl #48
	and	x9, x10, x9
	adc	x9, x9, x12
	stp	x8, x9, [x0, #16]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_fr_reduce_5limb                ; -- Begin function fr_reduce_5limb
	.p2align	2
_fr_reduce_5limb:                       ; @fr_reduce_5limb
	.cfi_startproc
; %bb.0:
	ldp	x9, x8, [x1, #24]
	ldr	x10, [x1, #16]
	extr	x11, x8, x9, #63
	mov	x12, #36319                     ; =0x8ddf
	movk	x12, #32695, lsl #16
	movk	x12, #9531, lsl #32
	movk	x12, #36180, lsl #48
	mul	x13, x11, x12
	umulh	x11, x11, x12
	extr	x11, x11, x13, #63
	sub	x12, x11, x11, lsl #32
	mov	x13, #-4294967295               ; =0xffffffff00000001
	umulh	x13, x11, x13
	mov	x14, #23550                     ; =0x5bfe
	movk	x14, #65534, lsl #16
	movk	x14, #41986, lsl #32
	movk	x14, #21437, lsl #48
	umulh	x15, x11, x14
	mul	x14, x11, x14
	adds	x13, x13, x14
	mov	x14, #55301                     ; =0xd805
	movk	x14, #2465, lsl #16
	movk	x14, #55304, lsl #32
	movk	x14, #13113, lsl #48
	umulh	x16, x11, x14
	mul	x14, x11, x14
	adcs	x14, x14, x15
	mov	x15, #32072                     ; =0x7d48
	movk	x15, #10653, lsl #16
	movk	x15, #42835, lsl #32
	movk	x15, #29677, lsl #48
	umulh	x17, x11, x15
	mul	x11, x11, x15
	adcs	x11, x11, x16
	cinc	x15, x17, hs
	ldp	x16, x17, [x1]
	subs	x12, x16, x12
	ngc	x16, xzr
	asr	x16, x16, #63
	subs	x13, x17, x13
	ngc	x17, xzr
	adds	x13, x13, x16
	adc	x16, x17, x16
	asr	x16, x16, #63
	subs	x10, x10, x14
	ngc	x14, xzr
	adds	x10, x16, x10
	adc	x14, x16, x14
	asr	x14, x14, #63
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x14, x9
	adc	x11, x14, x11
	sub	x8, x15, x8
	add	x8, x8, x11, lsr #63
	and	x11, x8, #0xfffffffffffffffe
	and	x11, x11, #0xfffffffe00000003
	subs	x11, x12, x11
	ngc	x12, xzr
	asr	x12, x12, #63
	mov	x14, #47101                     ; =0xb7fd
	movk	x14, #65532, lsl #16
	movk	x14, #18437, lsl #32
	movk	x14, #42875, lsl #48
	and	x14, x8, x14
	subs	x13, x13, x14
	ngc	x14, xzr
	adds	x13, x12, x13
	adc	x12, x12, x14
	asr	x12, x12, #63
	mov	x14, #45066                     ; =0xb00a
	movk	x14, #4931, lsl #16
	movk	x14, #45072, lsl #32
	movk	x14, #26227, lsl #48
	and	x14, x8, x14
	subs	x10, x10, x14
	ngc	x14, xzr
	adds	x10, x12, x10
	adc	x12, x12, x14
	mov	x14, #64144                     ; =0xfa90
	movk	x14, #21306, lsl #16
	movk	x14, #20134, lsl #32
	movk	x14, #59355, lsl #48
	and	x8, x8, x14
	sub	x8, x9, x8
	add	x8, x8, x12, asr #63
	mov	w9, #-1                         ; =0xffffffff
	adds	x12, x11, x9
	mov	x14, #-1                        ; =0xffffffffffffffff
	adc	x15, xzr, x14
	asr	x15, x15, #63
	adds	x16, x13, x15
	cinc	x15, x15, hs
	mov	x17, #41986                     ; =0xa402
	movk	x17, #1, lsl #16
	movk	x17, #23549, lsl #32
	movk	x17, #44098, lsl #48
	adds	x16, x16, x17
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x1, x10, x15
	cinc	x15, x15, hs
	mov	x2, #10235                      ; =0x27fb
	movk	x2, #63070, lsl #16
	movk	x2, #10231, lsl #32
	movk	x2, #52422, lsl #48
	adds	x1, x1, x2
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x3, x8, x15
	cinc	x15, x15, hs
	mov	x4, #33464                      ; =0x82b8
	movk	x4, #54882, lsl #16
	movk	x4, #22700, lsl #32
	movk	x4, #35858, lsl #48
	adds	x3, x3, x4
	adc	x15, x15, x14
	cmp	x15, #0
	csel	x11, x12, x11, ge
	csel	x12, x16, x13, ge
	csel	x10, x1, x10, ge
	csel	x8, x3, x8, ge
	adds	x9, x11, x9
	adc	x13, xzr, x14
	asr	x13, x13, #63
	adds	x15, x12, x13
	cinc	x13, x13, hs
	adds	x15, x15, x17
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x16, x10, x13
	cinc	x13, x13, hs
	adds	x16, x16, x2
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x17, x8, x13
	cinc	x13, x13, hs
	adds	x17, x17, x4
	adc	x13, x13, x14
	cmp	x13, #0
	csel	x9, x9, x11, ge
	csel	x11, x15, x12, ge
	stp	x9, x11, [x0]
	csel	x9, x16, x10, ge
	csel	x8, x17, x8, ge
	stp	x9, x8, [x0, #16]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_fr_mul_absx                    ; -- Begin function fr_mul_absx
	.p2align	2
_fr_mul_absx:                           ; @fr_mul_absx
	.cfi_startproc
; %bb.0:
	mov	x8, #65536                      ; =0x10000
	movk	x8, #53761, lsl #48
	ldp	x9, x10, [x1]
	mul	x11, x9, x8
	umulh	x9, x9, x8
	umulh	x12, x10, x8
	mul	x10, x10, x8
	adds	x9, x10, x9
	cinc	x10, x12, hs
	ldp	x12, x13, [x1, #16]
	umulh	x14, x12, x8
	mul	x12, x12, x8
	adds	x10, x10, x12
	cinc	x12, x14, hs
	umulh	x14, x13, x8
	mul	x8, x13, x8
	adds	x8, x12, x8
	cinc	x12, x14, hs
	extr	x13, x12, x8, #63
	mov	x14, #36319                     ; =0x8ddf
	movk	x14, #32695, lsl #16
	movk	x14, #9531, lsl #32
	movk	x14, #36180, lsl #48
	mul	x15, x13, x14
	umulh	x13, x13, x14
	extr	x13, x13, x15, #63
	sub	x14, x13, x13, lsl #32
	mov	x15, #-4294967295               ; =0xffffffff00000001
	umulh	x15, x13, x15
	mov	x16, #23550                     ; =0x5bfe
	movk	x16, #65534, lsl #16
	movk	x16, #41986, lsl #32
	movk	x16, #21437, lsl #48
	umulh	x17, x13, x16
	mul	x16, x13, x16
	adds	x15, x15, x16
	mov	x16, #55301                     ; =0xd805
	movk	x16, #2465, lsl #16
	movk	x16, #55304, lsl #32
	movk	x16, #13113, lsl #48
	umulh	x1, x13, x16
	mul	x16, x13, x16
	adcs	x16, x16, x17
	mov	x17, #32072                     ; =0x7d48
	movk	x17, #10653, lsl #16
	movk	x17, #42835, lsl #32
	movk	x17, #29677, lsl #48
	umulh	x2, x13, x17
	mul	x13, x13, x17
	adcs	x13, x13, x1
	cinc	x17, x2, hs
	subs	x11, x11, x14
	ngc	x14, xzr
	asr	x14, x14, #63
	subs	x9, x9, x15
	ngc	x15, xzr
	adds	x9, x9, x14
	adc	x14, x15, x14
	asr	x14, x14, #63
	subs	x10, x10, x16
	ngc	x15, xzr
	adds	x10, x14, x10
	adc	x14, x14, x15
	asr	x14, x14, #63
	subs	x8, x8, x13
	ngc	x13, xzr
	adds	x8, x14, x8
	adc	x13, x14, x13
	sub	x12, x17, x12
	add	x12, x12, x13, lsr #63
	and	x13, x12, #0xfffffffffffffffe
	and	x13, x13, #0xfffffffe00000003
	subs	x11, x11, x13
	ngc	x13, xzr
	asr	x13, x13, #63
	mov	x14, #47101                     ; =0xb7fd
	movk	x14, #65532, lsl #16
	movk	x14, #18437, lsl #32
	movk	x14, #42875, lsl #48
	and	x14, x12, x14
	subs	x9, x9, x14
	ngc	x14, xzr
	adds	x9, x13, x9
	adc	x13, x13, x14
	asr	x13, x13, #63
	mov	x14, #45066                     ; =0xb00a
	movk	x14, #4931, lsl #16
	movk	x14, #45072, lsl #32
	movk	x14, #26227, lsl #48
	and	x14, x12, x14
	subs	x10, x10, x14
	ngc	x14, xzr
	adds	x10, x13, x10
	adc	x13, x13, x14
	mov	x14, #64144                     ; =0xfa90
	movk	x14, #21306, lsl #16
	movk	x14, #20134, lsl #32
	movk	x14, #59355, lsl #48
	and	x12, x12, x14
	sub	x8, x8, x12
	add	x8, x8, x13, asr #63
	mov	w12, #-1                        ; =0xffffffff
	adds	x13, x11, x12
	mov	x14, #-1                        ; =0xffffffffffffffff
	adc	x15, xzr, x14
	asr	x15, x15, #63
	adds	x16, x9, x15
	cinc	x15, x15, hs
	mov	x17, #41986                     ; =0xa402
	movk	x17, #1, lsl #16
	movk	x17, #23549, lsl #32
	movk	x17, #44098, lsl #48
	adds	x16, x16, x17
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x1, x10, x15
	cinc	x15, x15, hs
	mov	x2, #10235                      ; =0x27fb
	movk	x2, #63070, lsl #16
	movk	x2, #10231, lsl #32
	movk	x2, #52422, lsl #48
	adds	x1, x1, x2
	adc	x15, x15, x14
	asr	x15, x15, #63
	adds	x3, x8, x15
	cinc	x15, x15, hs
	mov	x4, #33464                      ; =0x82b8
	movk	x4, #54882, lsl #16
	movk	x4, #22700, lsl #32
	movk	x4, #35858, lsl #48
	adds	x3, x3, x4
	adc	x15, x15, x14
	cmp	x15, #0
	csel	x11, x13, x11, ge
	csel	x9, x16, x9, ge
	csel	x10, x1, x10, ge
	csel	x8, x3, x8, ge
	adds	x12, x11, x12
	adc	x13, xzr, x14
	asr	x13, x13, #63
	adds	x15, x9, x13
	cinc	x13, x13, hs
	adds	x15, x15, x17
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x16, x10, x13
	cinc	x13, x13, hs
	adds	x16, x16, x2
	adc	x13, x13, x14
	asr	x13, x13, #63
	adds	x17, x8, x13
	cinc	x13, x13, hs
	adds	x17, x17, x4
	adc	x13, x13, x14
	cmp	x13, #0
	csel	x11, x12, x11, ge
	csel	x9, x15, x9, ge
	stp	x11, x9, [x0]
	csel	x9, x16, x10, ge
	csel	x8, x17, x8, ge
	stp	x9, x8, [x0, #16]
	ret
	.cfi_endproc
                                        ; -- End function
```

---

## 6. Full Assembly Listing: Stage 1b Shoup61 NTT & 3-Prime CRT Recombination (`stage1b_cofactor_ntt.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_shoup61_ntt_fwd_block          ; -- Begin function shoup61_ntt_fwd_block
	.p2align	2
_shoup61_ntt_fwd_block:                 ; @shoup61_ntt_fwd_block
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB0_3
; %bb.1:
	lsl	x8, x5, #1
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x9, [x0]
	ldr	x10, [x2], #8
	ldr	x11, [x1]
	ldr	x12, [x3], #8
	umulh	x12, x12, x11
	mul	x10, x10, x11
	msub	x10, x12, x5, x10
	cmp	x9, x8
	csel	x11, xzr, x8, lo
	sub	x9, x9, x11
	add	x11, x10, x9
	str	x11, [x0], #8
	add	x9, x9, x8
	sub	x9, x9, x10
	str	x9, [x1], #8
	subs	x4, x4, #1
	b.ne	LBB0_2
LBB0_3:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_shoup61_ntt_inv_block          ; -- Begin function shoup61_ntt_inv_block
	.p2align	2
_shoup61_ntt_inv_block:                 ; @shoup61_ntt_inv_block
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB1_3
; %bb.1:
	lsl	x8, x5, #1
LBB1_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x9, [x0]
	ldr	x10, [x1]
	add	x11, x10, x9
	cmp	x11, x8
	csel	x12, xzr, x8, lo
	sub	x11, x11, x12
	add	x9, x9, x8
	sub	x9, x9, x10
	str	x11, [x0], #8
	ldr	x10, [x2], #8
	ldr	x11, [x3], #8
	umulh	x11, x11, x9
	mul	x9, x10, x9
	msub	x9, x11, x5, x9
	str	x9, [x1], #8
	subs	x4, x4, #1
	b.ne	LBB1_2
LBB1_3:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_shoup61_pointwise_mul          ; -- Begin function shoup61_pointwise_mul
	.p2align	2
_shoup61_pointwise_mul:                 ; @shoup61_pointwise_mul
	.cfi_startproc
; %bb.0:
	cbz	x3, LBB2_3
; %bb.1:
	lsl	x8, x4, #1
LBB2_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x9, [x1], #8
	ldr	x10, [x2], #8
	cmp	x9, x8
	csel	x11, xzr, x8, lo
	sub	x9, x9, x11
	cmp	x9, x4
	csel	x11, xzr, x4, lo
	sub	x9, x9, x11
	cmp	x10, x8
	csel	x11, xzr, x8, lo
	sub	x10, x10, x11
	cmp	x10, x4
	csel	x11, xzr, x4, lo
	sub	x10, x10, x11
	mul	x11, x10, x9
	umulh	x9, x10, x9
	extr	x9, x9, x11, #60
	mul	x10, x9, x5
	umulh	x9, x9, x5
	extr	x9, x9, x10, #62
	msub	x9, x9, x4, x11
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	cmp	x9, x4
	csel	x10, xzr, x4, lo
	sub	x9, x9, x10
	str	x9, [x0], #8
	subs	x3, x3, #1
	b.ne	LBB2_2
LBB2_3:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_crt3_recombine_carry           ; -- Begin function crt3_recombine_carry
	.p2align	2
_crt3_recombine_carry:                  ; @crt3_recombine_carry
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB0_4
; %bb.1:
	stp	x28, x27, [sp, #-80]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 80
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	.cfi_offset w23, -40
	.cfi_offset w24, -48
	.cfi_offset w25, -56
	.cfi_offset w26, -64
	.cfi_offset w27, -72
	.cfi_offset w28, -80
	mov	x9, #0                          ; =0x0
	mov	x14, #0                         ; =0x0
	ldp	x12, x8, [sp, #104]
	mov	x10, #-2053641430080946177      ; =0xe37fffffffffffff
	mov	x11, #2053641430080946176       ; =0x1c80000000000000
	ldp	x21, x13, [sp, #88]
	mov	x15, #1                         ; =0x1
	movk	x15, #7296, lsl #48
	mov	x16, #-1945555039024054273      ; =0xe4ffffffffffffff
	mov	x17, #1945555039024054272       ; =0x1b00000000000000
	mov	x19, #1                         ; =0x1
	movk	x19, #6912, lsl #48
	mov	x20, #1                         ; =0x1
	movk	x20, #8128, lsl #48
	ldr	x22, [sp, #80]
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x23, [x0], #8
	ldr	x24, [x1], #8
	add	x25, x23, x10
	cmp	x23, x11
	csel	x25, x25, x23, hi
	add	x26, x24, x15
	cmp	x24, x25
	csel	x24, x26, x24, lo
	sub	x24, x24, x25
	umulh	x25, x24, x6
	mul	x24, x24, x5
	madd	x24, x25, x10, x24
	add	x25, x24, x10
	cmp	x24, x11
	csel	x24, x25, x24, hi
	ldr	x25, [x2], #8
	umulh	x26, x24, x8
	mul	x27, x24, x12
	madd	x26, x26, x16, x27
	add	x27, x26, x16
	cmp	x26, x17
	csel	x26, x27, x26, hi
	add	x27, x23, x16
	cmp	x23, x17
	csel	x27, x27, x23, hi
	add	x26, x26, x27
	add	x27, x26, x16
	cmp	x26, x17
	csel	x26, x27, x26, hi
	add	x27, x25, x19
	cmp	x25, x26
	csel	x25, x27, x25, lo
	sub	x25, x25, x26
	umulh	x26, x25, x22
	mul	x25, x25, x7
	madd	x25, x26, x16, x25
	add	x26, x25, x16
	umulh	x27, x24, x20
	cmp	x25, x17
	csel	x25, x26, x25, hi
	mul	x24, x24, x20
	umulh	x26, x25, x21
	mul	x28, x25, x21
	adds	x14, x23, x14
	cset	w23, hs
	adds	x14, x14, x24
	adc	x23, x23, x27
	adds	x24, x14, x28
	adc	x14, x23, x26
	umulh	x23, x25, x13
	mul	x25, x25, x13
	adds	x9, x25, x9
	cinc	x23, x23, hs
	adds	x14, x9, x14
	cinc	x9, x23, hs
	str	x24, [x3], #8
	subs	x4, x4, #1
	b.ne	LBB0_2
; %bb.3:
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #80             ; 16-byte Folded Reload
LBB0_4:
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```

---

## 7. Full Assembly Listing: Stage 4 Torsion-16 Edwards ECM (`stage4_ecm_torsion16.s`)

```asm
	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_m256_mont_mul                  ; -- Begin function m256_mont_mul
	.p2align	2
_m256_mont_mul:                         ; @m256_mont_mul
	.cfi_startproc
; %bb.0:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
	ldp	x16, x15, [x2]
	ldp	x14, x12, [x2, #16]
	ldp	x11, x10, [x3]
	ldp	x9, x8, [x3, #16]
	ldr	x13, [x3, #32]
	ldp	x17, x2, [x1]
	mul	x3, x17, x16
	umulh	x4, x17, x16
	umulh	x5, x17, x15
	mul	x6, x17, x15
	adds	x4, x6, x4
	cinc	x5, x5, hs
	umulh	x6, x17, x14
	mul	x7, x17, x14
	adds	x5, x7, x5
	cinc	x6, x6, hs
	mul	x7, x17, x12
	adds	x19, x7, x6
	mul	x20, x13, x3
	mul	x21, x20, x11
	cmn	x21, x3
	umulh	x3, x20, x11
	cinc	x3, x3, hs
	mul	x21, x20, x10
	adds	x4, x21, x4
	umulh	x21, x20, x10
	cinc	x21, x21, hs
	adds	x3, x4, x3
	cinc	x4, x21, hs
	mul	x21, x20, x9
	adds	x5, x5, x21
	umulh	x21, x20, x9
	cinc	x21, x21, hs
	adds	x4, x5, x4
	cinc	x5, x21, hs
	umulh	x21, x20, x8
	mul	x20, x20, x8
	adds	x19, x19, x20
	cinc	x20, x21, hs
	adds	x5, x19, x5
	cinc	x19, x20, hs
	cmn	x7, x6
	umulh	x17, x17, x12
	adc	x17, x19, x17
	umulh	x6, x2, x16
	mul	x7, x2, x16
	adds	x3, x7, x3
	cinc	x6, x6, hs
	umulh	x7, x2, x15
	mul	x19, x2, x15
	adds	x4, x19, x4
	cinc	x7, x7, hs
	adds	x4, x4, x6
	cinc	x6, x7, hs
	umulh	x7, x2, x14
	mul	x19, x2, x14
	adds	x5, x19, x5
	cinc	x7, x7, hs
	adds	x5, x5, x6
	cinc	x6, x7, hs
	umulh	x7, x2, x12
	mul	x2, x2, x12
	adds	x17, x2, x17
	cinc	x7, x7, hs
	adds	x19, x17, x6
	mul	x20, x13, x3
	umulh	x2, x20, x11
	mul	x21, x20, x11
	cmn	x21, x3
	cinc	x2, x2, hs
	umulh	x3, x20, x10
	mul	x21, x20, x10
	adds	x4, x21, x4
	cinc	x3, x3, hs
	adds	x2, x4, x2
	cinc	x3, x3, hs
	umulh	x4, x20, x9
	mul	x21, x20, x9
	adds	x5, x5, x21
	cinc	x21, x4, hs
	adds	x4, x5, x3
	cinc	x3, x21, hs
	umulh	x5, x20, x8
	mul	x20, x20, x8
	adds	x19, x19, x20
	cinc	x20, x5, hs
	adds	x5, x19, x3
	cinc	x3, x20, hs
	cmn	x17, x6
	adc	x17, x3, x7
	ldp	x3, x1, [x1, #16]
	umulh	x6, x3, x16
	mul	x7, x3, x16
	adds	x2, x7, x2
	cinc	x6, x6, hs
	umulh	x7, x3, x15
	mul	x19, x3, x15
	adds	x4, x19, x4
	cinc	x7, x7, hs
	adds	x4, x4, x6
	cinc	x6, x7, hs
	umulh	x7, x3, x14
	mul	x19, x3, x14
	adds	x5, x19, x5
	cinc	x7, x7, hs
	adds	x5, x5, x6
	cinc	x6, x7, hs
	umulh	x7, x3, x12
	mul	x3, x3, x12
	adds	x17, x3, x17
	cinc	x3, x7, hs
	adds	x7, x17, x6
	mul	x19, x13, x2
	umulh	x20, x19, x11
	mul	x21, x19, x11
	cmn	x21, x2
	cinc	x2, x20, hs
	umulh	x20, x19, x10
	mul	x21, x19, x10
	adds	x4, x21, x4
	cinc	x20, x20, hs
	adds	x2, x4, x2
	cinc	x4, x20, hs
	umulh	x20, x19, x9
	mul	x21, x19, x9
	adds	x5, x5, x21
	cinc	x20, x20, hs
	adds	x4, x5, x4
	cinc	x5, x20, hs
	umulh	x20, x19, x8
	mul	x19, x19, x8
	adds	x7, x7, x19
	cinc	x19, x20, hs
	adds	x5, x7, x5
	cinc	x7, x19, hs
	cmn	x17, x6
	adc	x17, x7, x3
	umulh	x3, x1, x16
	mul	x16, x1, x16
	adds	x16, x16, x2
	cinc	x2, x3, hs
	umulh	x3, x1, x15
	mul	x15, x1, x15
	adds	x15, x15, x4
	cinc	x3, x3, hs
	adds	x15, x15, x2
	cinc	x2, x3, hs
	umulh	x3, x1, x14
	mul	x14, x1, x14
	adds	x14, x14, x5
	cinc	x3, x3, hs
	adds	x14, x14, x2
	cinc	x2, x3, hs
	umulh	x3, x1, x12
	mul	x12, x1, x12
	adds	x12, x12, x17
	cinc	x17, x3, hs
	adds	x1, x12, x2
	mul	x13, x13, x16
	umulh	x3, x13, x11
	mul	x4, x13, x11
	cmn	x4, x16
	cinc	x16, x3, hs
	umulh	x3, x13, x10
	mul	x4, x13, x10
	adds	x15, x4, x15
	cinc	x3, x3, hs
	adds	x15, x15, x16
	cinc	x16, x3, hs
	umulh	x3, x13, x9
	mul	x4, x13, x9
	adds	x14, x14, x4
	cinc	x3, x3, hs
	adds	x14, x14, x16
	cinc	x16, x3, hs
	umulh	x3, x13, x8
	mul	x13, x13, x8
	adds	x13, x1, x13
	cinc	x1, x3, hs
	adds	x13, x13, x16
	cinc	x16, x1, hs
	cmn	x12, x2
	adc	x12, x16, x17
	subs	x11, x15, x11
	ngc	x16, xzr
	asr	x16, x16, #63
	subs	x10, x14, x10
	ngc	x17, xzr
	adds	x10, x16, x10
	adc	x16, x16, x17
	asr	x16, x16, #63
	subs	x9, x13, x9
	ngc	x17, xzr
	adds	x9, x16, x9
	adc	x16, x16, x17
	asr	x16, x16, #63
	subs	x8, x12, x8
	ngc	x17, xzr
	adds	x8, x16, x8
	adc	x16, x16, x17
	cmp	x16, #0
	csel	x11, x11, x15, ge
	csel	x10, x10, x14, ge
	stp	x11, x10, [x0]
	csel	x9, x9, x13, ge
	csel	x8, x8, x12, ge
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_ecm_mont_xdbladd_11m           ; -- Begin function ecm_mont_xdbladd_11m
	.p2align	2
_ecm_mont_xdbladd_11m:                  ; @ecm_mont_xdbladd_11m
	.cfi_startproc
; %bb.0:
	stp	x28, x27, [sp, #-96]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #80]             ; 16-byte Folded Spill
	add	x29, sp, #80
	sub	sp, sp, #464
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	mov	x19, x7
	mov	x22, x6
	str	x5, [sp]                        ; 8-byte Folded Spill
	mov	x23, x4
	mov	x25, x3
	mov	x24, x1
	mov	x20, x0
	ldp	x9, x10, [x2]
	ldp	x14, x12, [x2, #32]
	adds	x16, x14, x9
	adcs	x17, x12, x10
	ldp	x13, x8, [x2, #16]
	ldp	x15, x11, [x2, #48]
	adcs	x0, x15, x13
	adcs	x1, x11, x8
	cset	w2, hs
	ldp	x26, x27, [x7]
	subs	x3, x16, x26
	ngc	x4, xzr
	asr	x4, x4, #63
	subs	x5, x17, x27
	ngc	x6, xzr
	adds	x5, x4, x5
	adc	x4, x4, x6
	asr	x4, x4, #63
	ldp	x28, x21, [x7, #16]
	subs	x6, x0, x28
	ngc	x7, xzr
	adds	x6, x4, x6
	adc	x4, x4, x7
	asr	x4, x4, #63
	subs	x7, x1, x21
	ngc	x30, xzr
	adds	x7, x4, x7
	adc	x4, x4, x30
	lsr	x4, x4, #63
	bic	x2, x4, x2
	neg	x4, x2
	and	x16, x16, x4
	sub	x2, x2, #1
	and	x3, x2, x3
	orr	x16, x3, x16
	and	x17, x4, x17
	and	x3, x2, x5
	orr	x17, x3, x17
	stp	x16, x17, [x29, #-120]
	and	x16, x4, x0
	and	x17, x2, x6
	orr	x16, x17, x16
	and	x17, x4, x1
	and	x0, x2, x7
	orr	x17, x0, x17
	stp	x16, x17, [x29, #-104]
	subs	x9, x9, x14
	ngc	x14, xzr
	asr	x14, x14, #63
	subs	x10, x10, x12
	ngc	x12, xzr
	adds	x10, x10, x14
	adc	x12, x12, x14
	asr	x12, x12, #63
	subs	x13, x13, x15
	ngc	x14, xzr
	adds	x13, x12, x13
	adc	x12, x12, x14
	asr	x12, x12, #63
	subs	x8, x8, x11
	ngc	x11, xzr
	adds	x8, x12, x8
	adc	x11, x12, x11
	asr	x11, x11, #63
	and	x12, x26, x11
	adds	x9, x9, x12
	and	x12, x27, x11
	adcs	x10, x10, x12
	stp	x9, x10, [x29, #-152]
	and	x9, x28, x11
	adcs	x9, x13, x9
	and	x10, x21, x11
	adc	x8, x10, x8
	stp	x9, x8, [x29, #-136]
	sub	x0, x29, #248
	sub	x1, x29, #120
	sub	x2, x29, #120
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, sp, #264
	sub	x1, x29, #152
	sub	x2, x29, #152
	mov	x3, x19
	bl	_m256_mont_mul
	ldp	x8, x9, [x29, #-248]
	ldp	x10, x11, [sp, #264]
	subs	x8, x8, x10
	ngc	x10, xzr
	asr	x10, x10, #63
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	asr	x10, x10, #63
	ldp	x11, x12, [x29, #-232]
	ldp	x13, x14, [sp, #280]
	subs	x11, x11, x13
	ngc	x13, xzr
	adds	x11, x10, x11
	adc	x10, x10, x13
	asr	x10, x10, #63
	subs	x12, x12, x14
	ngc	x13, xzr
	adds	x12, x10, x12
	adc	x10, x10, x13
	asr	x10, x10, #63
	and	x13, x26, x10
	adds	x8, x8, x13
	and	x13, x27, x10
	adcs	x9, x9, x13
	stp	x8, x9, [sp, #232]
	and	x8, x28, x10
	adcs	x8, x11, x8
	and	x9, x21, x10
	adc	x9, x9, x12
	stp	x8, x9, [sp, #248]
	ldp	x9, x10, [x25]
	ldp	x11, x12, [x25, #32]
	adds	x13, x11, x9
	adcs	x14, x12, x10
	ldp	x15, x8, [x25, #16]
	ldp	x16, x17, [x25, #48]
	adcs	x0, x16, x15
	adcs	x1, x17, x8
	cset	w2, hs
	subs	x3, x13, x26
	ngc	x4, xzr
	asr	x4, x4, #63
	subs	x5, x14, x27
	ngc	x6, xzr
	adds	x5, x5, x4
	adc	x4, x6, x4
	asr	x4, x4, #63
	subs	x6, x0, x28
	ngc	x7, xzr
	adds	x6, x4, x6
	adc	x4, x4, x7
	asr	x4, x4, #63
	subs	x7, x1, x21
	ngc	x25, xzr
	adds	x7, x4, x7
	adc	x4, x4, x25
	lsr	x4, x4, #63
	bic	x2, x4, x2
	neg	x4, x2
	and	x13, x13, x4
	sub	x2, x2, #1
	and	x3, x2, x3
	orr	x13, x3, x13
	and	x14, x4, x14
	and	x3, x2, x5
	orr	x14, x3, x14
	stp	x13, x14, [x29, #-184]
	and	x13, x4, x0
	and	x14, x2, x6
	orr	x13, x14, x13
	and	x14, x4, x1
	and	x0, x2, x7
	orr	x14, x0, x14
	stp	x13, x14, [x29, #-168]
	subs	x9, x9, x11
	ngc	x11, xzr
	asr	x11, x11, #63
	subs	x10, x10, x12
	ngc	x12, xzr
	adds	x10, x10, x11
	adc	x11, x12, x11
	asr	x11, x11, #63
	subs	x12, x15, x16
	ngc	x13, xzr
	adds	x12, x11, x12
	adc	x11, x11, x13
	asr	x11, x11, #63
	subs	x8, x8, x17
	ngc	x13, xzr
	adds	x8, x11, x8
	adc	x11, x11, x13
	asr	x11, x11, #63
	and	x13, x26, x11
	adds	x9, x9, x13
	and	x13, x27, x11
	adcs	x10, x10, x13
	stp	x9, x10, [x29, #-216]
	and	x9, x28, x11
	adcs	x9, x12, x9
	and	x10, x21, x11
	adc	x8, x10, x8
	stp	x9, x8, [x29, #-200]
	add	x0, sp, #200
	sub	x1, x29, #216
	sub	x2, x29, #120
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, sp, #168
	sub	x1, x29, #184
	sub	x2, x29, #152
	mov	x3, x19
	bl	_m256_mont_mul
	ldp	x9, x10, [sp, #200]
	ldp	x11, x12, [sp, #168]
	adds	x13, x11, x9
	adcs	x14, x12, x10
	ldp	x15, x8, [sp, #216]
	ldp	x16, x17, [sp, #184]
	adcs	x0, x16, x15
	adcs	x1, x17, x8
	cset	w2, hs
	subs	x3, x13, x26
	ngc	x4, xzr
	asr	x4, x4, #63
	subs	x5, x14, x27
	ngc	x6, xzr
	adds	x5, x5, x4
	adc	x4, x6, x4
	asr	x4, x4, #63
	subs	x6, x0, x28
	ngc	x7, xzr
	adds	x6, x4, x6
	adc	x4, x4, x7
	asr	x4, x4, #63
	subs	x7, x1, x21
	ngc	x25, xzr
	adds	x7, x4, x7
	adc	x4, x4, x25
	lsr	x4, x4, #63
	bic	x2, x4, x2
	neg	x4, x2
	and	x13, x13, x4
	sub	x2, x2, #1
	and	x3, x2, x3
	orr	x13, x3, x13
	and	x14, x4, x14
	and	x3, x2, x5
	orr	x14, x3, x14
	stp	x13, x14, [sp, #136]
	and	x13, x4, x0
	and	x14, x2, x6
	orr	x13, x14, x13
	and	x14, x4, x1
	and	x0, x2, x7
	orr	x14, x0, x14
	stp	x13, x14, [sp, #152]
	subs	x9, x9, x11
	ngc	x11, xzr
	asr	x11, x11, #63
	subs	x10, x10, x12
	ngc	x12, xzr
	adds	x10, x10, x11
	adc	x11, x12, x11
	asr	x11, x11, #63
	subs	x12, x15, x16
	ngc	x13, xzr
	adds	x12, x11, x12
	adc	x11, x11, x13
	asr	x11, x11, #63
	subs	x8, x8, x17
	ngc	x13, xzr
	adds	x8, x11, x8
	adc	x11, x11, x13
	asr	x11, x11, #63
	and	x13, x26, x11
	adds	x9, x9, x13
	and	x13, x27, x11
	adcs	x10, x10, x13
	stp	x9, x10, [sp, #104]
	and	x9, x28, x11
	adcs	x9, x12, x9
	and	x10, x21, x11
	adc	x8, x10, x8
	stp	x9, x8, [sp, #120]
	add	x1, sp, #136
	add	x2, sp, #136
	mov	x0, x24
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, sp, #72
	add	x1, sp, #104
	add	x2, sp, #104
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, x24, #32
	add	x1, sp, #72
	mov	x2, x23
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, sp, #40
	add	x2, sp, #264
	mov	x1, x22
	mov	x3, x19
	bl	_m256_mont_mul
	sub	x1, x29, #248
	add	x2, sp, #40
	mov	x0, x20
	mov	x3, x19
	bl	_m256_mont_mul
	add	x0, sp, #8
	add	x2, sp, #232
	ldr	x1, [sp]                        ; 8-byte Folded Reload
	mov	x3, x19
	bl	_m256_mont_mul
	ldp	x8, x9, [sp, #40]
	ldp	x10, x11, [sp, #8]
	adds	x8, x10, x8
	adcs	x9, x11, x9
	ldp	x10, x11, [sp, #56]
	ldp	x12, x13, [sp, #24]
	adcs	x10, x12, x10
	adcs	x11, x13, x11
	cset	w12, hs
	subs	x13, x8, x26
	ngc	x14, xzr
	asr	x14, x14, #63
	subs	x15, x9, x27
	ngc	x16, xzr
	adds	x15, x15, x14
	adc	x14, x16, x14
	asr	x14, x14, #63
	subs	x16, x10, x28
	ngc	x17, xzr
	adds	x16, x14, x16
	adc	x14, x14, x17
	asr	x14, x14, #63
	subs	x17, x11, x21
	ngc	x0, xzr
	adds	x17, x14, x17
	adc	x14, x14, x0
	lsr	x14, x14, #63
	bic	x12, x14, x12
	neg	x14, x12
	and	x8, x8, x14
	sub	x12, x12, #1
	and	x13, x12, x13
	orr	x8, x13, x8
	and	x9, x14, x9
	and	x13, x12, x15
	orr	x9, x13, x9
	stp	x8, x9, [sp, #40]
	and	x8, x14, x10
	and	x9, x12, x16
	orr	x8, x9, x8
	and	x9, x14, x11
	and	x10, x12, x17
	orr	x9, x10, x9
	stp	x8, x9, [sp, #56]
	add	x0, x20, #32
	add	x1, sp, #232
	add	x2, sp, #40
	mov	x3, x19
	bl	_m256_mont_mul
	add	sp, sp, #464
	ldp	x29, x30, [sp, #80]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #96             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```

---

## 8. Full Assembly Listing: 2D Bivariate 3-Prime 61-Bit Shoup NTT & 192-Bit Garner CRT + Overlap + 58-Bit Carry (`kernels_ntt2d_shoup61_v81.s`)

```asm
.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 5	sdk_version 15, 5
	.globl	_shoup61_bfly_loop              ; -- Begin function shoup61_bfly_loop
	.p2align	2
_shoup61_bfly_loop:                     ; @shoup61_bfly_loop
	.cfi_startproc
; %bb.0:
	cbz	x4, LBB0_2
LBB0_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x8, [x0]
	ldr	x9, [x3], #8
	ldr	x10, [x1]
	umulh	x9, x9, x10
	ldr	x11, [x2], #8
	mul	x10, x11, x10
	cmp	x8, x6
	csel	x11, xzr, x6, lo
	msub	x9, x9, x5, x10
	sub	x8, x8, x11
	add	x10, x9, x8
	str	x10, [x0], #8
	add	x8, x8, x6
	sub	x8, x8, x9
	str	x8, [x1], #8
	subs	x4, x4, #1
	b.ne	LBB0_1
LBB0_2:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_unpack_58bit_limbs             ; -- Begin function unpack_58bit_limbs
	.p2align	2
_unpack_58bit_limbs:                    ; @unpack_58bit_limbs
	.cfi_startproc
; %bb.0:
	cbz	x2, LBB1_6
; %bb.1:
	mov	x8, #0                          ; =0x0
	b	LBB1_4
LBB1_2:                                 ;   in Loop: Header=BB1_4 Depth=1
	ldr	x11, [x9, #8]
	neg	w10, w10
	lsl	x10, x11, x10
LBB1_3:                                 ;   in Loop: Header=BB1_4 Depth=1
	ldr	x9, [x9]
	and	x11, x8, #0x3e
	lsr	x9, x9, x11
	orr	x9, x10, x9
	and	x9, x9, #0x3ffffffffffffff
	str	x9, [x1], #8
	add	x8, x8, #58
	subs	x2, x2, #1
	b.eq	LBB1_6
LBB1_4:                                 ; =>This Inner Loop Header: Depth=1
	and	w10, w8, #0x3e
	lsr	x9, x8, #6
	add	x9, x0, x9, lsl #3
	cmp	w10, #7
	b.hs	LBB1_2
; %bb.5:                                ;   in Loop: Header=BB1_4 Depth=1
	mov	x10, #0                         ; =0x0
	b	LBB1_3
LBB1_6:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_garner_crt3_overlap_carry58    ; -- Begin function garner_crt3_overlap_carry58
	.p2align	2
_garner_crt3_overlap_carry58:           ; @garner_crt3_overlap_carry58
	.cfi_startproc
; %bb.0:
	cbz	x7, LBB2_4
; %bb.1:
	stp	x28, x27, [sp, #-96]!           ; 16-byte Folded Spill
	stp	x26, x25, [sp, #16]             ; 16-byte Folded Spill
	stp	x24, x23, [sp, #32]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #48]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #64]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #80]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 96
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	.cfi_offset w25, -72
	.cfi_offset w26, -80
	.cfi_offset w27, -88
	.cfi_offset w28, -96
	mov	x22, #0                         ; =0x0
	mov	x25, #0                         ; =0x0
	ldp	x9, x8, [sp, #144]
	ldp	x11, x10, [sp, #128]
	mov	x12, #-2053641430080946177      ; =0xe37fffffffffffff
	mov	x13, #2053641430080946176       ; =0x1c80000000000000
	ldp	x15, x14, [sp, #112]
	mov	x17, #-1945555039024054273      ; =0xe4ffffffffffffff
	mov	x19, #1945555039024054272       ; =0x1b00000000000000
	mov	x21, #1                         ; =0x1
	movk	x21, #8128, lsl #48
	ldp	x24, x23, [sp, #96]
LBB2_2:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x26, [x0], #8
	ldr	x27, [x1], #8
	add	x28, x26, x12
	cmp	x26, x13
	csel	x28, x28, x26, hi
	mov	x16, #1                         ; =0x1
	movk	x16, #7296, lsl #48
	add	x30, x27, x16
	cmp	x27, x28
	csel	x27, x30, x27, lo
	sub	x27, x27, x28
	umulh	x28, x27, x23
	mul	x27, x27, x24
	madd	x27, x28, x12, x27
	add	x28, x27, x12
	ldr	x30, [x2], #8
	cmp	x27, x13
	csel	x27, x28, x27, hi
	umulh	x28, x27, x8
	mul	x16, x27, x9
	madd	x16, x28, x17, x16
	add	x28, x16, x17
	cmp	x16, x19
	csel	x16, x28, x16, hi
	add	x28, x26, x17
	cmp	x26, x19
	csel	x28, x28, x26, hi
	add	x16, x16, x28
	add	x28, x16, x17
	cmp	x16, x19
	csel	x16, x28, x16, hi
	mov	x20, #1                         ; =0x1
	movk	x20, #6912, lsl #48
	add	x28, x30, x20
	cmp	x30, x16
	csel	x28, x28, x30, lo
	sub	x16, x28, x16
	umulh	x28, x16, x14
	mul	x16, x16, x15
	madd	x16, x28, x17, x16
	add	x28, x16, x17
	cmp	x16, x19
	csel	x16, x28, x16, hi
	umulh	x28, x27, x21
	mul	x27, x27, x21
	adds	x26, x27, x26
	cinc	x27, x28, hs
	umulh	x28, x16, x11
	mul	x30, x16, x11
	adds	x26, x26, x30
	umulh	x30, x16, x10
	ldr	x20, [x3], #8
	adc	x27, x27, x28
	mul	x16, x16, x10
	adds	x20, x20, x25
	cset	w25, hs
	ldr	x28, [x4], #8
	adds	x22, x28, x22
	cset	w28, hs
	adds	x16, x27, x16
	ldr	x27, [x5], #8
	adc	x27, x27, x30
	adds	x16, x22, x16
	cinc	x22, x28, hs
	adds	x20, x20, x26
	and	x26, x20, #0x3ffffffffffffff
	str	x26, [x6], #8
	adcs	x16, x16, x25
	adc	x22, x27, x22
	extr	x25, x16, x20, #58
	extr	x22, x22, x16, #58
	subs	x7, x7, #1
	b.ne	LBB2_2
; %bb.3:
	ldp	x29, x30, [sp, #80]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #48]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #32]             ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #16]             ; 16-byte Folded Reload
	ldp	x28, x27, [sp], #96             ; 16-byte Folded Reload
LBB2_4:
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
```
