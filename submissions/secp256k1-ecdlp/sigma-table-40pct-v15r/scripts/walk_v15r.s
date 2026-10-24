	.build_version macos, 27, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_v15o_align                     ; -- Begin function v15o_align
	.p2align	2
_v15o_align:                            ; @v15o_align
	.cfi_startproc
	sub	sp, sp, #128
	stp	x22, x21, [sp, #80]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #96]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #112]            ; 16-byte Folded Spill
	add	x29, sp, #112
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	mov	x20, x1
	mov	x19, x0
Lloh0:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh1:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh2:
	ldr	x8, [x8]
	stur	x8, [x29, #-40]
Lloh3:
	adrp	x21, _BETA@PAGE
Lloh4:
	add	x21, x21, _BETA@PAGEOFF
	add	x0, sp, #40
	mov	x1, x19
	mov	x2, x21
	bl	_fmul
	add	x0, sp, #8
	add	x1, sp, #40
	mov	x2, x21
	bl	_fmul
	ldr	x10, [x19, #24]
	ldr	x12, [sp, #64]
	cmp	x12, x10
	b.hs	LBB0_3
	ldr	x13, [sp, #56]
LBB0_2:
	ldp	x8, x9, [sp, #40]
	mov	w0, #1                          ; =0x1
	mov	x10, x12
	mov	x11, x13
	b	LBB0_5
LBB0_3:
	ldp	x8, x9, [x19]
	ldr	x11, [x19, #16]
	b.ls	LBB0_13
LBB0_4:
	mov	x0, #0                          ; =0x0
LBB0_5:
	ldr	x12, [sp, #32]
	cmp	x12, x10
	b.hs	LBB0_8
	ldr	x13, [sp, #24]
LBB0_7:
	ldp	x8, x9, [sp, #8]
	mov	w0, #2                          ; =0x2
	mov	x10, x12
	mov	x11, x13
	b	LBB0_9
LBB0_8:
	b.ls	LBB0_18
LBB0_9:
	stp	x8, x9, [x19]
	stp	x11, x10, [x19, #16]
	mov	w10, #15381                     ; =0x3c15
	mul	w9, w9, w10
	eor	w8, w8, w9
	ubfx	x8, x8, #3, #11
	str	x8, [x19, #64]
	ldr	x8, [x19, #32]
	and	w9, w8, #0x1
	str	w9, [x20]
	tbz	w8, #0, LBB0_11
	mov	x9, #-977                       ; =0xfffffffffffffc2f
	movk	x9, #65534, lsl #32
	subs	x8, x9, x8
	ngc	x9, xzr
	ldp	x10, x11, [x19, #40]
	mvn	x10, x10
	adds	x10, x9, x10
	mvn	x11, x11
	adc	x9, x11, x9
	ldr	x11, [x19, #56]
	mvn	x11, x11
	stp	x8, x10, [x19, #32]
	stp	x9, x11, [x19, #48]
LBB0_11:
	ldur	x8, [x29, #-40]
Lloh5:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh6:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh7:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB0_23
	ldp	x29, x30, [sp, #112]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #96]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #128
	ret
LBB0_13:
	ldr	x13, [sp, #56]
	cmp	x13, x11
	b.lo	LBB0_2
	b.hi	LBB0_4
	ldr	x14, [sp, #48]
	cmp	x14, x9
	b.lo	LBB0_2
	mov	x0, #0                          ; =0x0
	cmp	x14, x9
	b.hi	LBB0_5
	ldr	x14, [sp, #40]
	cmp	x14, x8
	b.lo	LBB0_2
	b	LBB0_5
LBB0_18:
	ldr	x13, [sp, #24]
	cmp	x13, x11
	b.lo	LBB0_7
	b.hi	LBB0_9
	ldr	x14, [sp, #16]
	cmp	x14, x9
	b.lo	LBB0_7
	b.hi	LBB0_9
	ldr	x14, [sp, #8]
	cmp	x14, x8
	b.lo	LBB0_7
	b	LBB0_9
LBB0_23:
	bl	___stack_chk_fail
	.loh AdrpAdd	Lloh3, Lloh4
	.loh AdrpLdrGotLdr	Lloh0, Lloh1, Lloh2
	.loh AdrpLdrGotLdr	Lloh5, Lloh6, Lloh7
	.cfi_endproc
	.p2align	2                               ; -- Begin function fmul
_fmul:                                  ; @fmul
	.cfi_startproc
	mov	x8, #-1                         ; =0xffffffffffffffff
	mov	x9, #977                        ; =0x3d1
	movk	x9, #1, lsl #32
	ldp	x12, x11, [x2]
	ldp	x13, x16, [x1]
	mul	x10, x12, x13
	umulh	x14, x12, x13
	umulh	x15, x11, x13
	mul	x17, x11, x13
	adds	x17, x14, x17
	cinc	x3, x15, hs
	ldp	x15, x14, [x2, #16]
	umulh	x2, x15, x13
	mul	x4, x15, x13
	adds	x3, x3, x4
	cinc	x2, x2, hs
	umulh	x4, x14, x13
	mul	x13, x14, x13
	adds	x2, x2, x13
	cinc	x13, x4, hs
	cmp	x13, #0
	csel	x4, xzr, x13, eq
	umulh	x5, x12, x16
	mul	x13, x12, x16
	adds	x13, x13, x17
	cinc	x17, x5, hs
	umulh	x5, x11, x16
	mul	x6, x11, x16
	adds	x17, x17, x3
	cset	w3, hs
	adds	x17, x17, x6
	adc	x3, x3, x5
	umulh	x5, x15, x16
	mul	x6, x15, x16
	adds	x2, x3, x2
	cset	w3, hs
	adds	x2, x2, x6
	adc	x3, x3, x5
	umulh	x5, x14, x16
	mul	x16, x14, x16
	adds	x3, x3, x4
	cset	w4, hs
	adds	x16, x3, x16
	adcs	x3, x4, x5
	csel	x3, xzr, x3, eq
	ldp	x4, x1, [x1, #16]
	umulh	x5, x12, x4
	mul	x6, x12, x4
	adds	x17, x6, x17
	cinc	x5, x5, hs
	umulh	x6, x11, x4
	mul	x7, x11, x4
	adds	x2, x5, x2
	cset	w5, hs
	adds	x2, x2, x7
	adc	x5, x5, x6
	umulh	x6, x15, x4
	mul	x7, x15, x4
	adds	x16, x5, x16
	cset	w5, hs
	adds	x16, x16, x7
	adc	x5, x5, x6
	umulh	x6, x14, x4
	mul	x4, x14, x4
	adds	x3, x5, x3
	cset	w5, hs
	adds	x3, x3, x4
	adcs	x4, x5, x6
	csel	x4, xzr, x4, eq
	umulh	x5, x12, x1
	mul	x12, x12, x1
	adds	x12, x12, x2
	cinc	x2, x5, hs
	umulh	x5, x11, x1
	mul	x11, x11, x1
	adds	x16, x2, x16
	cset	w2, hs
	adds	x11, x16, x11
	adc	x16, x2, x5
	umulh	x2, x15, x1
	mul	x15, x15, x1
	adds	x16, x16, x3
	cset	w3, hs
	cmn	x16, x15
	adc	x2, x3, x2
	umulh	x3, x14, x1
	mul	x14, x14, x1
	adds	x1, x2, x4
	cset	w2, hs
	adds	x14, x1, x14
	adc	x1, x2, x3
	umulh	x2, x1, x9
	mul	x1, x1, x9
	adds	x12, x1, x12
	adcs	x11, x11, x2
	adc	x15, x16, x15
	umulh	x16, x14, x9
	mul	x14, x14, x9
	adds	x14, x14, x17
	adcs	x12, x12, x16
	cinc	x11, x11, hs
	umulh	x16, x15, x9
	mul	x15, x15, x9
	adds	x13, x15, x13
	adcs	x14, x14, x16
	cinc	x12, x12, hs
	umulh	x15, x11, x9
	mul	x11, x11, x9
	adds	x10, x11, x10
	adcs	x11, x13, x15
	cinc	x13, x14, hs
	adds	x14, x10, x9
	adc	x15, xzr, x8
	asr	x16, x15, #63
	adds	x15, x15, x11
	cinc	x16, x16, hs
	adds	x15, x15, #1
	adc	x16, x16, x8
	asr	x17, x16, #63
	adds	x16, x16, x13
	cinc	x17, x17, hs
	adds	x16, x16, #1
	adc	x17, x17, x8
	asr	x1, x17, #63
	adds	x17, x17, x12
	cinc	x1, x1, hs
	adds	x17, x17, #1
	adc	x1, x1, x8
	adds	x9, x14, x9
	cinc	x2, x8, hs
	asr	x3, x2, #63
	adds	x2, x2, x15
	cinc	x3, x3, hs
	adds	x2, x2, #1
	adc	x3, x3, x8
	asr	x4, x3, #63
	adds	x3, x3, x16
	cinc	x4, x4, hs
	adds	x3, x3, #1
	adc	x4, x4, x8
	asr	x5, x4, #63
	adds	x4, x4, x17
	cinc	x5, x5, hs
	adds	x4, x4, #1
	adc	x8, x5, x8
	tst	w1, #0x1
	csel	x8, x1, x8, ne
	csel	x4, x17, x4, ne
	csel	x3, x16, x3, ne
	csel	x2, x15, x2, ne
	csel	x9, x14, x9, ne
	tst	w1, #0x1
	csel	x12, x12, x17, ne
	csel	x13, x13, x16, ne
	csel	x11, x11, x15, ne
	csel	x10, x10, x14, ne
	tst	w8, #0x1
	csel	x8, x12, x4, ne
	csel	x11, x11, x2, ne
	csel	x9, x10, x9, ne
	stp	x9, x11, [x0]
	csel	x9, x13, x3, ne
	stp	x9, x8, [x0, #16]
	ret
	.cfi_endproc
	.globl	_v15o_label                     ; -- Begin function v15o_label
	.p2align	2
_v15o_label:                            ; @v15o_label
	.cfi_startproc
	sub	sp, sp, #80
	stp	x20, x19, [sp, #48]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #64]             ; 16-byte Folded Spill
	add	x29, sp, #64
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	mov	x19, x1
Lloh8:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh9:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh10:
	ldr	x8, [x8]
	stur	x8, [x29, #-24]
	cmp	x2, #1
	b.eq	LBB2_3
	cbnz	x2, LBB2_4
Lloh11:
	adrp	x9, _R_MONT_ONE@PAGE
Lloh12:
	add	x9, x9, _R_MONT_ONE@PAGEOFF
	b	LBB2_5
LBB2_3:
Lloh13:
	adrp	x9, _LAMBDA_R@PAGE
Lloh14:
	add	x9, x9, _LAMBDA_R@PAGEOFF
	b	LBB2_5
LBB2_4:
Lloh15:
	adrp	x9, _LAMBDA2_R@PAGE
Lloh16:
	add	x9, x9, _LAMBDA2_R@PAGEOFF
LBB2_5:
	ldp	x13, x16, [x9]
	stp	x13, x16, [sp, #8]
	add	x8, x9, #24
	ldr	x2, [x9, #16]!
	ldr	x17, [x8]
	stp	x2, x17, [sp, #24]
Lloh17:
	adrp	x8, _R_CONST@PAGE
Lloh18:
	add	x8, x8, _R_CONST@PAGEOFF
	ldp	x12, x11, [x8]
	ldp	x10, x9, [x8, #16]
	mov	x14, #0                         ; =0x0
	mov	x15, #0                         ; =0x0
	cbz	w3, LBB2_7
	subs	x13, x12, x13
	ngc	x3, xzr
	subs	x16, x11, x16
	ngc	x6, xzr
	adds	x16, x16, x3
	adc	x3, x6, x3
	subs	x2, x10, x2
	ngc	x6, xzr
	adds	x2, x2, x3
	asr	x3, x3, #63
	adc	x3, x6, x3
	stp	x13, x16, [sp, #8]
	sub	x13, x9, x17
	add	x13, x13, x3
	stp	x2, x13, [sp, #24]
LBB2_7:
	ldp	x13, x16, [x0]
	ldp	x17, x2, [x4]
	adds	x13, x17, x13
	adcs	x16, x16, xzr
	cset	w17, hs
	ldp	x3, x6, [x0, #16]
	ldp	x7, x4, [x4, #16]
	adds	x3, x7, x3
	cset	w7, hs
	adds	x4, x4, x6
	cset	w6, hs
	adds	x16, x16, x2
	adcs	x17, x3, x17
	adcs	x2, x4, x7
	cinc	x3, x6, hs
	subs	x12, x13, x12
	ngc	x4, xzr
	subs	x11, x16, x11
	ngc	x15, x15
	adds	x11, x11, x4
	adc	x15, x15, x4
	asr	x4, x15, #63
	subs	x10, x17, x10
	ngc	x14, x14
	adds	x10, x10, x15
	adc	x14, x14, x4
	asr	x15, x14, #63
	subs	x9, x2, x9
	ngc	x4, xzr
	adds	x9, x9, x14
	adc	x14, x4, x15
	cbnz	x3, LBB2_10
	and	x14, x14, #0x1
	tbz	w14, #0, LBB2_10
	mov	x12, x13
	mov	x11, x16
	mov	x10, x17
	mov	x9, x2
LBB2_10:
	stp	x12, x11, [x0]
	stp	x10, x9, [x0, #16]
	ldp	x9, x10, [x19]
	ldp	x11, x12, [x5]
	adds	x9, x11, x9
	adcs	x10, x10, xzr
	cset	w13, hs
	ldp	x11, x14, [x19, #16]
	ldp	x15, x16, [x5, #16]
	adds	x15, x15, x11
	cset	w17, hs
	adds	x14, x16, x14
	cset	w16, hs
	adds	x11, x10, x12
	adcs	x12, x15, x13
	adcs	x14, x14, x17
	cinc	x17, x16, hs
	ldp	x10, x13, [x8]
	subs	x10, x9, x10
	ngc	x15, xzr
	subs	x13, x11, x13
	ngc	x16, xzr
	adds	x13, x13, x15
	adc	x15, x16, x15
	asr	x16, x15, #63
	ldp	x2, x3, [x8, #16]
	subs	x8, x12, x2
	ngc	x2, xzr
	adds	x8, x8, x15
	adc	x15, x2, x16
	asr	x16, x15, #63
	subs	x2, x14, x3
	ngc	x3, xzr
	adds	x15, x2, x15
	adc	x16, x3, x16
	cbnz	x17, LBB2_13
	and	x16, x16, #0x1
	tbz	w16, #0, LBB2_13
	mov	x10, x9
	mov	x13, x11
	mov	x8, x12
	mov	x15, x14
LBB2_13:
	stp	x10, x13, [x19]
	stp	x8, x15, [x19, #16]
	add	x2, sp, #8
	mov	x1, x0
	bl	_rmont
	add	x2, sp, #8
	mov	x0, x19
	mov	x1, x19
	bl	_rmont
	ldur	x8, [x29, #-24]
Lloh19:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh20:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh21:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB2_15
	ldp	x29, x30, [sp, #64]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #80
	ret
LBB2_15:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh8, Lloh9, Lloh10
	.loh AdrpAdd	Lloh11, Lloh12
	.loh AdrpAdd	Lloh13, Lloh14
	.loh AdrpAdd	Lloh15, Lloh16
	.loh AdrpAdd	Lloh17, Lloh18
	.loh AdrpLdrGotLdr	Lloh19, Lloh20, Lloh21
	.cfi_endproc
	.p2align	2                               ; -- Begin function rmont
_rmont:                                 ; @rmont
	.cfi_startproc
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
Lloh22:
	adrp	x8, _R_N0@PAGE
Lloh23:
	ldr	x8, [x8, _R_N0@PAGEOFF]
Lloh24:
	adrp	x12, _R_CONST@PAGE
Lloh25:
	add	x12, x12, _R_CONST@PAGEOFF
	ldp	x10, x9, [x1]
	ldp	x11, x17, [x2]
	mul	x14, x11, x10
	umulh	x13, x11, x10
	umulh	x15, x11, x9
	mul	x16, x11, x9
	adds	x4, x13, x16
	cinc	x13, x15, hs
	ldp	x16, x15, [x1, #16]
	umulh	x1, x11, x16
	mul	x3, x11, x16
	adds	x5, x13, x3
	cinc	x13, x1, hs
	umulh	x1, x11, x15
	mul	x11, x11, x15
	adds	x6, x13, x11
	cinc	x3, x1, hs
	mul	x1, x8, x14
	ldp	x13, x11, [x12]
	umulh	x7, x1, x13
	mul	x19, x1, x13
	cmn	x19, x14
	cinc	x14, x7, hs
	umulh	x7, x11, x1
	mul	x19, x11, x1
	adds	x14, x4, x14
	cset	w4, hs
	adds	x19, x14, x19
	adc	x4, x4, x7
	ldp	x14, x12, [x12, #16]
	umulh	x7, x14, x1
	mul	x20, x14, x1
	adds	x4, x5, x4
	cset	w5, hs
	adds	x4, x4, x20
	adc	x5, x5, x7
	umulh	x7, x12, x1
	mul	x1, x12, x1
	adds	x5, x6, x5
	cset	w6, hs
	adds	x1, x5, x1
	adc	x7, x6, x7
	adds	x5, x3, x7
	umulh	x6, x17, x10
	mul	x20, x17, x10
	adds	x19, x19, x20
	cinc	x6, x6, hs
	umulh	x20, x17, x9
	mul	x21, x17, x9
	adds	x4, x4, x6
	cset	w6, hs
	adds	x4, x4, x21
	adc	x6, x6, x20
	umulh	x20, x17, x16
	mul	x21, x17, x16
	adds	x1, x1, x6
	cset	w6, hs
	adds	x21, x1, x21
	adc	x1, x6, x20
	umulh	x6, x17, x15
	mul	x17, x17, x15
	adds	x1, x1, x5
	cset	w5, hs
	adds	x20, x1, x17
	adc	x17, x5, x6
	mul	x5, x8, x19
	umulh	x1, x5, x13
	mul	x6, x5, x13
	cmn	x6, x19
	cinc	x1, x1, hs
	umulh	x6, x11, x5
	mul	x19, x11, x5
	adds	x1, x4, x1
	cset	w4, hs
	adds	x1, x1, x19
	adc	x4, x4, x6
	umulh	x6, x14, x5
	mul	x19, x14, x5
	adds	x4, x21, x4
	cset	w21, hs
	adds	x4, x4, x19
	adc	x6, x21, x6
	umulh	x19, x12, x5
	mul	x5, x12, x5
	adds	x6, x20, x6
	cset	w20, hs
	adds	x5, x6, x5
	adc	x6, x20, x19
	cmn	x3, x7
	ldp	x3, x2, [x2, #16]
	umulh	x7, x3, x10
	mul	x19, x3, x10
	adcs	x17, x17, x6
	mrs	x6, NZCV
	adds	x1, x1, x19
	umulh	x19, x3, x9
	mul	x20, x3, x9
	cinc	x7, x7, hs
	adds	x4, x4, x7
	cset	w7, hs
	adds	x4, x4, x20
	umulh	x20, x3, x16
	mul	x21, x3, x16
	adc	x7, x7, x19
	adds	x5, x5, x7
	cset	w7, hs
	adds	x5, x5, x21
	umulh	x19, x3, x15
	mul	x3, x3, x15
	adc	x7, x7, x20
	adds	x17, x7, x17
	cset	w7, hs
	adds	x17, x17, x3
	adc	x3, x7, x19
	msr	NZCV, x6
	adcs	x3, x3, xzr
	cset	w6, hs
	mul	x7, x8, x1
	umulh	x19, x7, x13
	mul	x20, x7, x13
	cmn	x20, x1
	cinc	x1, x19, hs
	umulh	x19, x11, x7
	mul	x20, x11, x7
	adds	x1, x4, x1
	cset	w4, hs
	adds	x1, x1, x20
	adc	x4, x4, x19
	umulh	x19, x14, x7
	mul	x20, x14, x7
	adds	x4, x5, x4
	cset	w5, hs
	adds	x4, x4, x20
	adc	x5, x5, x19
	umulh	x19, x12, x7
	mul	x7, x12, x7
	adds	x17, x17, x5
	cset	w5, hs
	adds	x17, x17, x7
	adc	x5, x5, x19
	adds	x7, x3, x5
	umulh	x19, x2, x10
	mul	x10, x2, x10
	adds	x10, x1, x10
	cinc	x1, x19, hs
	umulh	x19, x2, x9
	mul	x9, x2, x9
	adds	x1, x4, x1
	cset	w4, hs
	adds	x9, x1, x9
	adc	x1, x4, x19
	umulh	x4, x2, x16
	mul	x16, x2, x16
	adds	x17, x17, x1
	cset	w1, hs
	adds	x16, x17, x16
	adc	x17, x1, x4
	umulh	x1, x2, x15
	mul	x15, x2, x15
	adds	x17, x17, x7
	cset	w2, hs
	adds	x15, x17, x15
	adc	x17, x2, x1
	cmn	x3, x5
	adcs	x17, x17, x6
	cset	w1, hs
	mul	x2, x8, x10
	umulh	x8, x2, x13
	mul	x3, x2, x13
	cmn	x3, x10
	cinc	x8, x8, hs
	umulh	x10, x11, x2
	mul	x3, x11, x2
	adds	x8, x9, x8
	cset	w9, hs
	adds	x8, x8, x3
	adc	x9, x9, x10
	umulh	x10, x14, x2
	mul	x3, x14, x2
	adds	x9, x16, x9
	cset	w16, hs
	adds	x9, x9, x3
	adc	x10, x16, x10
	umulh	x16, x12, x2
	mul	x2, x12, x2
	adds	x10, x15, x10
	cset	w15, hs
	adds	x10, x10, x2
	adc	x16, x15, x16
	adds	x2, x17, x16
	subs	x13, x8, x13
	ngc	x15, xzr
	and	x15, x15, #0x1
	adds	x11, x15, x11
	cset	w3, hs
	subs	x15, x9, x11
	ngc	x11, x3
	and	x11, x11, #0x1
	adds	x11, x11, x14
	cset	w3, hs
	subs	x14, x10, x11
	ngc	x11, x3
	and	x11, x11, #0x1
	adds	x11, x11, x12
	cset	w3, hs
	subs	x12, x2, x11
	ngc	x2, x3
	adds	x11, x17, x16
	asr	x16, x2, #63
	adcs	xzr, x2, x1
	cinc	w16, w16, hs
	tbnz	w16, #0, LBB3_2
	mov	x11, x12
	mov	x10, x14
	mov	x9, x15
	mov	x8, x13
LBB3_2:
	stp	x8, x9, [x0]
	stp	x10, x11, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.loh AdrpAdd	Lloh24, Lloh25
	.loh AdrpLdr	Lloh22, Lloh23
	.cfi_endproc
	.globl	_v15o_detect                    ; -- Begin function v15o_detect
	.p2align	2
_v15o_detect:                           ; @v15o_detect
	.cfi_startproc
	ldp	x8, x11, [x0]
	ldp	x12, x13, [x1]
	ldp	x10, x9, [x0, #16]
	ldp	x14, x15, [x1, #16]
	cmp	x8, x12
	ccmp	x11, x13, #0, eq
	ccmp	x10, x14, #0, eq
	ccmp	x9, x15, #0, eq
	b.eq	LBB4_5
	ldr	x12, [x2]
	cmp	x8, x12
	b.ne	LBB4_6
	ldr	x12, [x2, #8]
	cmp	x11, x12
	b.ne	LBB4_6
	ldr	x11, [x2, #16]
	cmp	x10, x11
	b.ne	LBB4_6
	ldr	x10, [x2, #24]
	cmp	x9, x10
	b.ne	LBB4_6
LBB4_5:
	ldr	x9, [x3]
	add	x9, x9, #1
	str	x9, [x3]
	ldr	x9, [x0, #24]
LBB4_6:
	ldr	x10, [x4, #24]
	cmp	x9, x10
	b.hs	LBB4_8
LBB4_7:
	ldr	x9, [x0]
	str	x9, [x4]
	ldr	x9, [x0, #8]
	str	x9, [x4, #8]
	ldr	x9, [x0, #16]
	str	x9, [x4, #16]
	ldr	x9, [x0, #24]
	str	x9, [x4, #24]
	str	xzr, [x5]
	b	LBB4_9
LBB4_8:
	b.ls	LBB4_10
LBB4_9:
	ldr	x9, [x5]
	add	x10, x9, #1
	cmp	x10, #67
	mov	w10, #67                        ; =0x43
	csinc	x9, x10, x9, hs
	str	x9, [x5]
	mvn	w8, w8
	ldr	x9, [x1]
	str	x9, [x2]
	ldr	x9, [x1, #8]
	str	x9, [x2, #8]
	ldr	x9, [x1, #16]
	str	x9, [x2, #16]
	ldr	x9, [x1, #24]
	str	x9, [x2, #24]
	ldr	x9, [x0]
	str	x9, [x1]
	ldr	x9, [x0, #8]
	str	x9, [x1, #8]
	ldr	x9, [x0, #16]
	str	x9, [x1, #16]
	ldr	x9, [x0, #24]
	str	x9, [x1, #24]
	and	x0, x8, #0x1
	ret
LBB4_10:
	ldr	x9, [x0, #16]
	ldr	x10, [x4, #16]
	cmp	x9, x10
	b.lo	LBB4_7
	b.hi	LBB4_9
	ldr	x9, [x0, #8]
	ldr	x10, [x4, #8]
	cmp	x9, x10
	b.lo	LBB4_7
	b.hi	LBB4_9
	ldr	x9, [x0]
	ldr	x10, [x4]
	cmp	x9, x10
	b.lo	LBB4_7
	b	LBB4_9
	.cfi_endproc
	.globl	_v15o_report                    ; -- Begin function v15o_report
	.p2align	2
_v15o_report:                           ; @v15o_report
	.cfi_startproc
	ldp	x9, x10, [x2]
	mov	x8, #31765                      ; =0x7c15
	movk	x8, #32586, lsl #16
	movk	x8, #31161, lsl #32
	movk	x8, #40503, lsl #48
	mul	x8, x10, x8
	eor	x8, x8, x9
	and	x11, x8, x1
	mov	w8, #96                         ; =0x60
	madd	x8, x11, x8, x0
	ldp	x14, x15, [x8]
	ldp	x13, x12, [x8, #16]
	orr	x16, x15, x14
	orr	x17, x13, x12
	orr	x16, x16, x17
	cbz	x16, LBB5_17
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_4
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_4
	ldr	x13, [x2, #24]
	cmp	x12, x13
	b.eq	LBB5_16
LBB5_4:
	add	x8, x11, #1
	and	x8, x8, x1
	mov	w12, #96                        ; =0x60
	madd	x8, x8, x12, x0
	ldp	x14, x15, [x8]
	ldp	x13, x12, [x8, #16]
	orr	x16, x15, x14
	orr	x17, x13, x12
	orr	x16, x16, x17
	cbz	x16, LBB5_17
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_8
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_8
	ldr	x13, [x2, #24]
	cmp	x12, x13
	b.eq	LBB5_16
LBB5_8:
	add	x8, x11, #2
	and	x8, x8, x1
	mov	w12, #96                        ; =0x60
	madd	x8, x8, x12, x0
	ldp	x14, x15, [x8]
	ldp	x13, x12, [x8, #16]
	orr	x16, x15, x14
	orr	x17, x13, x12
	orr	x16, x16, x17
	cbz	x16, LBB5_17
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_12
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_12
	ldr	x13, [x2, #24]
	cmp	x12, x13
	b.eq	LBB5_16
LBB5_12:
	add	x8, x11, #3
	and	x8, x8, x1
	mov	w11, #96                        ; =0x60
	madd	x8, x8, x11, x0
	ldp	x13, x14, [x8]
	ldp	x12, x11, [x8, #16]
	orr	x15, x14, x13
	orr	x16, x12, x11
	orr	x15, x15, x16
	cbz	x15, LBB5_17
	cmp	x13, x9
	ccmp	x14, x10, #0, eq
	b.ne	LBB5_17
	ldr	x10, [x2, #16]
	cmp	x12, x10
	b.ne	LBB5_17
	ldr	x10, [x2, #24]
	cmp	x11, x10
	b.ne	LBB5_17
LBB5_16:
	ldr	x9, [x8, #32]
	str	x9, [x5]
	ldr	x9, [x8, #40]
	str	x9, [x5, #8]
	ldr	x9, [x8, #48]
	str	x9, [x5, #16]
	ldr	x9, [x8, #56]
	str	x9, [x5, #24]
	ldr	x9, [x8, #64]
	str	x9, [x6]
	ldr	x9, [x8, #72]
	str	x9, [x6, #8]
	ldr	x9, [x8, #80]
	str	x9, [x6, #16]
	ldr	x8, [x8, #88]
	str	x8, [x6, #24]
	mov	w0, #1                          ; =0x1
	ret
LBB5_17:
	mov	x0, #0                          ; =0x0
	str	x9, [x8]
	ldr	x9, [x2, #8]
	str	x9, [x8, #8]
	ldr	x9, [x2, #16]
	str	x9, [x8, #16]
	ldr	x9, [x2, #24]
	str	x9, [x8, #24]
	ldr	x9, [x3]
	str	x9, [x8, #32]
	ldr	x9, [x3, #8]
	str	x9, [x8, #40]
	ldr	x9, [x3, #16]
	str	x9, [x8, #48]
	ldr	x9, [x3, #24]
	str	x9, [x8, #56]
	ldr	x9, [x4]
	str	x9, [x8, #64]
	ldr	x9, [x4, #8]
	str	x9, [x8, #72]
	ldr	x9, [x4, #16]
	str	x9, [x8, #80]
	ldr	x9, [x4, #24]
	str	x9, [x8, #88]
	ret
	.cfi_endproc
	.globl	_v15r_init_for_Q                ; -- Begin function v15r_init_for_Q
	.p2align	2
_v15r_init_for_Q:                       ; @v15r_init_for_Q
	.cfi_startproc
	sub	sp, sp, #144
	stp	x24, x23, [sp, #80]             ; 16-byte Folded Spill
	stp	x22, x21, [sp, #96]             ; 16-byte Folded Spill
	stp	x20, x19, [sp, #112]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #128]            ; 16-byte Folded Spill
	add	x29, sp, #128
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	.cfi_offset w23, -56
	.cfi_offset w24, -64
	mov	x20, x0
	mov	x19, #0                         ; =0x0
Lloh26:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh27:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh28:
	ldr	x8, [x8]
	stur	x8, [x29, #-56]
	ldr	x8, [x1]
Lloh29:
	adrp	x22, _SMALL@GOTPAGE
Lloh30:
	ldr	x22, [x22, _SMALL@GOTPAGEOFF]
	str	x8, [x22, #64]
	ldr	x8, [x2]
	mov	x21, x22
	str	x8, [x21, #96]!
	ldr	x8, [x1, #8]
	stur	x8, [x21, #-24]
	ldr	x8, [x2, #8]
	str	x8, [x21, #8]
	ldr	x8, [x1, #16]
	stur	x8, [x21, #-16]
	ldr	x8, [x2, #16]
	str	x8, [x21, #16]
	ldr	x8, [x1, #24]
	stur	x8, [x21, #-8]
	ldr	x8, [x2, #24]
	str	x8, [x21, #24]
	mov	w23, #14                        ; =0xe
	b	LBB6_3
LBB6_1:                                 ;   in Loop: Header=BB6_3 Depth=1
	add	x19, x19, #1
LBB6_2:                                 ;   in Loop: Header=BB6_3 Depth=1
	add	x21, x21, #64
	subs	x23, x23, #1
	b.eq	LBB6_5
LBB6_3:                                 ; =>This Inner Loop Header: Depth=1
	sub	x2, x21, #32
	add	x0, sp, #40
	add	x1, sp, #8
	add	x4, x22, #64
	add	x5, x22, #96
	mov	x3, x21
	bl	_v15q_add
	cbnz	w0, LBB6_1
	ldr	x8, [sp, #40]
	str	x8, [x21, #32]
	ldr	x8, [sp, #8]
	ldur	q0, [sp, #48]
	stur	q0, [x21, #40]
	ldur	q0, [sp, #16]
	stur	q0, [x21, #72]
	ldr	x9, [sp, #64]
	stp	x9, x8, [x21, #56]
	ldr	x8, [sp, #32]
	str	x8, [x21, #88]
	b	LBB6_2
LBB6_5:
	mov	x21, #0                         ; =0x0
	add	x23, x20, #32
Lloh31:
	adrp	x20, _UBASE@GOTPAGE
Lloh32:
	ldr	x20, [x20, _UBASE@GOTPAGEOFF]
	mov	w24, #34953                     ; =0x8889
	b	LBB6_7
LBB6_6:                                 ;   in Loop: Header=BB6_7 Depth=1
	sub	x8, x23, #32
	ldur	q0, [sp, #40]
	ldur	q1, [sp, #56]
	stp	q0, q1, [x8]
	ldur	q0, [sp, #8]
	ldur	q1, [sp, #24]
	stp	q0, q1, [x23], #64
	add	x21, x21, #1
	add	x20, x20, #64
	cmp	x21, #2048
	b.eq	LBB6_9
LBB6_7:                                 ; =>This Inner Loop Header: Depth=1
	mul	w8, w21, w24
	lsr	w8, w8, #19
	sub	w8, w8, w8, lsl #4
	add	w8, w21, w8
	ubfiz	x8, x8, #6, #10
	add	x8, x22, x8
	add	x0, sp, #40
	add	x1, sp, #8
	add	x3, x20, #32
	add	x4, x8, #64
	add	x5, x8, #96
	mov	x2, x20
	bl	_v15q_add
	cbz	w0, LBB6_6
	add	x19, x19, #1
	str	xzr, [sp, #8]
	str	xzr, [sp, #40]
	b	LBB6_6
LBB6_9:
	ldur	x8, [x29, #-56]
Lloh33:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh34:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh35:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB6_11
	mov	x0, x19
	ldp	x29, x30, [sp, #128]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #112]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #96]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #144
	ret
LBB6_11:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh29, Lloh30
	.loh AdrpLdrGotLdr	Lloh26, Lloh27, Lloh28
	.loh AdrpLdrGot	Lloh31, Lloh32
	.loh AdrpLdrGotLdr	Lloh33, Lloh34, Lloh35
	.cfi_endproc
	.globl	_v15q_add                       ; -- Begin function v15q_add
	.p2align	2
_v15q_add:                              ; @v15q_add
	.cfi_startproc
	sub	sp, sp, #320
	stp	x28, x27, [sp, #224]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #240]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #256]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #272]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #288]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #304]            ; 16-byte Folded Spill
	add	x29, sp, #304
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
	mov	x23, x4
	mov	x20, x3
	mov	x21, x2
	mov	x19, x1
	mov	x22, x0
Lloh36:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh37:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh38:
	ldr	x8, [x8]
	stur	x8, [x29, #-96]
	ldp	x26, x27, [x2]
	ldp	x11, x10, [x4]
	ldp	x25, x24, [x2, #16]
	ldp	x9, x8, [x4, #16]
	cmp	x26, x11
	ccmp	x27, x10, #0, eq
	ccmp	x25, x9, #0, eq
	ccmp	x24, x8, #0, eq
	b.eq	LBB7_2
	subs	x12, x11, x26
	ngc	x11, xzr
	subs	x10, x10, x27
	ngc	x13, xzr
	adds	x10, x10, x11
	adc	x11, x13, x11
	asr	x13, x11, #63
	subs	x9, x9, x25
	ngc	x14, xzr
	adds	x9, x9, x11
	adc	x11, x14, x13
	asr	x13, x11, #63
	subs	x8, x8, x24
	ngc	x14, xzr
	adds	x8, x8, x11
	adc	x11, x14, x13
	mov	x13, #-977                      ; =0xfffffffffffffc2f
	movk	x13, #65534, lsl #32
	adds	x14, x12, x13
	adcs	x15, x10, xzr
	cset	w16, hs
	subs	x15, x15, #1
	adcs	x16, x9, x16
	cset	w17, hs
	subs	x16, x16, #1
	adc	x17, x8, x17
	sub	x17, x17, #1
	tst	w11, #0x1
	csel	x11, x8, x17, eq
	csel	x9, x9, x16, eq
	csel	x10, x10, x15, eq
	csel	x8, x12, x14, eq
	stp	x8, x10, [sp, #80]
	stp	x9, x11, [sp, #96]
	ldp	x12, x14, [x5]
	ldp	x15, x16, [x20]
	subs	x12, x12, x15
	ngc	x15, xzr
	subs	x14, x14, x16
	ngc	x16, xzr
	adds	x14, x14, x15
	adc	x15, x16, x15
	asr	x16, x15, #63
	ldp	x17, x0, [x5, #16]
	ldp	x1, x2, [x20, #16]
	subs	x17, x17, x1
	ngc	x1, xzr
	adds	x15, x17, x15
	adc	x16, x1, x16
	asr	x17, x16, #63
	subs	x0, x0, x2
	ngc	x1, xzr
	adds	x16, x0, x16
	adc	x17, x1, x17
	adds	x13, x12, x13
	adcs	x0, x14, xzr
	cset	w1, hs
	subs	x0, x0, #1
	adcs	x1, x15, x1
	cset	w2, hs
	subs	x1, x1, #1
	adc	x2, x16, x2
	sub	x2, x2, #1
	tst	w17, #0x1
	csel	x16, x16, x2, eq
	csel	x14, x14, x0, eq
	csel	x12, x12, x13, eq
	stp	x12, x14, [sp, #48]
	csel	x12, x15, x1, eq
	stp	x12, x16, [sp, #64]
	b	LBB7_7
LBB7_2:
	ldr	x8, [x20]
	ldr	x9, [x5]
	mov	w0, #-1                         ; =0xffffffff
	cmp	x8, x9
	b.ne	LBB7_22
	ldr	x8, [x20, #8]
	ldr	x9, [x5, #8]
	cmp	x8, x9
	b.ne	LBB7_22
	ldr	x8, [x20, #16]
	ldr	x9, [x5, #16]
	cmp	x8, x9
	b.ne	LBB7_22
	ldr	x8, [x20, #24]
	ldr	x9, [x5, #24]
	cmp	x8, x9
	b.ne	LBB7_22
	add	x0, sp, #16
	mov	x1, x21
	mov	x2, x21
	bl	_fmul
Lloh39:
	adrp	x2, _v15q_add.THREE@PAGE
Lloh40:
	add	x2, x2, _v15q_add.THREE@PAGEOFF
	add	x0, sp, #48
	add	x1, sp, #16
	bl	_fmul
Lloh41:
	adrp	x2, _v15q_add.TWO@PAGE
Lloh42:
	add	x2, x2, _v15q_add.TWO@PAGEOFF
	add	x0, sp, #80
	mov	x1, x20
	bl	_fmul
	ldp	x8, x10, [sp, #80]
	ldp	x9, x11, [sp, #96]
LBB7_7:
	orr	x8, x10, x8
	orr	x9, x9, x11
	orr	x8, x8, x9
	cbz	x8, LBB7_21
	str	x19, [sp, #8]                   ; 8-byte Folded Spill
Lloh43:
	adrp	x8, l___const.rinv.one@PAGE
Lloh44:
	add	x8, x8, l___const.rinv.one@PAGEOFF
	ldp	q0, q1, [x8]
	stp	q0, q1, [x29, #-128]
	mov	w28, #255                       ; =0xff
Lloh45:
	adrp	x19, _P_EXP@PAGE
Lloh46:
	add	x19, x19, _P_EXP@PAGEOFF
	b	LBB7_10
LBB7_9:                                 ;   in Loop: Header=BB7_10 Depth=1
	sub	w28, w28, #1
	cmn	w28, #1
	b.eq	LBB7_12
LBB7_10:                                ; =>This Inner Loop Header: Depth=1
	add	x0, sp, #144
	sub	x1, x29, #128
	sub	x2, x29, #128
	bl	_fmul
	ldp	x8, x9, [sp, #144]
	stp	x8, x9, [x29, #-128]
	ldp	x10, x11, [sp, #160]
	stp	x10, x11, [x29, #-112]
	lsr	w12, w28, #6
	ldr	x12, [x19, w12, uxtw #3]
	lsr	x12, x12, x28
	tbz	w12, #0, LBB7_9
	add	x0, sp, #112
	sub	x1, x29, #128
	add	x2, sp, #80
	bl	_fmul
	ldp	x8, x9, [sp, #112]
	stp	x8, x9, [x29, #-128]
	ldp	x10, x11, [sp, #128]
	stp	x10, x11, [x29, #-112]
	b	LBB7_9
LBB7_12:
	stp	x8, x9, [sp, #144]
	stp	x10, x11, [sp, #160]
	sub	x0, x29, #128
	add	x1, sp, #48
	add	x2, sp, #144
	bl	_fmul
	add	x0, sp, #16
	sub	x1, x29, #128
	sub	x2, x29, #128
	bl	_fmul
	ldp	x8, x9, [sp, #16]
	subs	x8, x8, x26
	ngc	x10, xzr
	subs	x9, x9, x27
	ngc	x12, xzr
	adds	x11, x9, x10
	adc	x9, x12, x10
	ldp	x10, x12, [sp, #32]
	subs	x10, x10, x25
	ngc	x13, xzr
	adds	x10, x10, x9
	asr	x9, x9, #63
	adc	x9, x13, x9
	asr	x13, x9, #63
	subs	x12, x12, x24
	ngc	x14, xzr
	adds	x9, x12, x9
	adc	x12, x14, x13
	tbz	w12, #0, LBB7_14
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x11, x11, xzr
	cset	w12, hs
	subs	x11, x11, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x9, x9, x12
	sub	x9, x9, #1
LBB7_14:
	ldr	x19, [sp, #8]                   ; 8-byte Folded Reload
	stp	x8, x11, [x22]
	stp	x10, x9, [x22, #16]
	ldp	x12, x13, [x23]
	subs	x8, x8, x12
	ngc	x12, xzr
	subs	x11, x11, x13
	ngc	x13, xzr
	adds	x11, x11, x12
	adc	x12, x13, x12
	asr	x13, x12, #63
	ldp	x14, x15, [x23, #16]
	subs	x10, x10, x14
	ngc	x14, xzr
	adds	x10, x10, x12
	adc	x12, x14, x13
	asr	x13, x12, #63
	subs	x9, x9, x15
	ngc	x14, xzr
	adds	x9, x9, x12
	adc	x12, x14, x13
	tbz	w12, #0, LBB7_16
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x11, x11, xzr
	cset	w12, hs
	subs	x11, x11, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x9, x9, x12
	sub	x9, x9, #1
LBB7_16:
	stp	x8, x11, [x22]
	stp	x10, x9, [x22, #16]
	ldp	x12, x13, [x21]
	subs	x8, x12, x8
	ngc	x12, xzr
	subs	x11, x13, x11
	ngc	x13, xzr
	adds	x11, x11, x12
	adc	x12, x13, x12
	asr	x13, x12, #63
	ldp	x14, x15, [x21, #16]
	subs	x10, x14, x10
	ngc	x14, xzr
	adds	x10, x10, x12
	adc	x12, x14, x13
	asr	x13, x12, #63
	subs	x9, x15, x9
	ngc	x14, xzr
	adds	x9, x9, x12
	adc	x12, x14, x13
	tbz	w12, #0, LBB7_18
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x11, x11, xzr
	cset	w12, hs
	subs	x11, x11, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x9, x9, x12
	sub	x9, x9, #1
LBB7_18:
	stp	x8, x11, [sp, #112]
	stp	x10, x9, [sp, #128]
	sub	x1, x29, #128
	add	x2, sp, #112
	mov	x0, x19
	bl	_fmul
	ldp	x8, x9, [x19]
	ldp	x10, x11, [x20]
	subs	x8, x8, x10
	ngc	x10, xzr
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x11, x11, x10
	ldp	x10, x12, [x19, #16]
	ldp	x13, x14, [x20, #16]
	subs	x10, x10, x13
	ngc	x13, xzr
	adds	x10, x10, x11
	asr	x11, x11, #63
	adc	x11, x13, x11
	asr	x13, x11, #63
	subs	x12, x12, x14
	ngc	x14, xzr
	adds	x11, x12, x11
	adc	x12, x14, x13
	tbz	w12, #0, LBB7_20
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x9, x9, xzr
	cset	w12, hs
	subs	x9, x9, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x11, x11, x12
	sub	x11, x11, #1
LBB7_20:
	mov	w0, #0                          ; =0x0
	stp	x8, x9, [x19]
	stp	x10, x11, [x19, #16]
	b	LBB7_22
LBB7_21:
	mov	w0, #-1                         ; =0xffffffff
LBB7_22:
	ldur	x8, [x29, #-96]
Lloh47:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh48:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh49:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB7_24
	ldp	x29, x30, [sp, #304]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #288]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #272]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #256]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #240]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #224]            ; 16-byte Folded Reload
	add	sp, sp, #320
	ret
LBB7_24:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh36, Lloh37, Lloh38
	.loh AdrpAdd	Lloh41, Lloh42
	.loh AdrpAdd	Lloh39, Lloh40
	.loh AdrpAdd	Lloh45, Lloh46
	.loh AdrpAdd	Lloh43, Lloh44
	.loh AdrpLdrGotLdr	Lloh47, Lloh48, Lloh49
	.cfi_endproc
	.globl	_v15r_action                    ; -- Begin function v15r_action
	.p2align	2
_v15r_action:                           ; @v15r_action
	.cfi_startproc
	ldr	x8, [x0]
	cmp	x8, #15
	b.ls	LBB8_2
	mov	w0, #1                          ; =0x1
	ret
LBB8_2:
	ldr	x8, [x1]
	cmp	x8, #66
	cset	w0, hi
	ret
	.cfi_endproc
	.globl	_v15r_reset_hist                ; -- Begin function v15r_reset_hist
	.p2align	2
_v15r_reset_hist:                       ; @v15r_reset_hist
	.cfi_startproc
	str	xzr, [x0]
	str	xzr, [x1]
	mov	x8, #-1                         ; =0xffffffffffffffff
	str	x8, [x3]
	str	xzr, [x0, #8]
	str	xzr, [x1, #8]
	str	x8, [x3, #8]
	str	xzr, [x0, #16]
	str	xzr, [x1, #16]
	str	x8, [x3, #16]
	str	xzr, [x0, #24]
	str	xzr, [x1, #24]
	str	x8, [x3, #24]
	str	xzr, [x2]
	str	xzr, [x4]
	ret
	.cfi_endproc
	.globl	_v15r_escape                    ; -- Begin function v15r_escape
	.p2align	2
_v15r_escape:                           ; @v15r_escape
	.cfi_startproc
	sub	sp, sp, #464
	stp	x28, x27, [sp, #368]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #384]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #400]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #416]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #432]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #448]            ; 16-byte Folded Spill
	add	x29, sp, #448
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
	mov	x28, x7
	mov	x20, x6
	mov	x26, x5
	mov	x27, x4
	stp	x3, x2, [sp, #48]               ; 16-byte Folded Spill
	stp	x0, x1, [sp, #64]               ; 16-byte Folded Spill
	ldr	x8, [x29, #64]
	str	x8, [sp, #40]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #56]
	str	x8, [sp, #32]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #48]
	str	x8, [sp, #24]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #40]
	str	x8, [sp, #16]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #32]
	str	x8, [sp, #8]                    ; 8-byte Folded Spill
	ldp	x24, x23, [x29, #16]
Lloh50:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh51:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh52:
	ldr	x8, [x8]
	stur	x8, [x29, #-104]
	mov	w25, #1                         ; =0x1
	mov	w8, #15                         ; =0xf
LBB10_1:                                ; =>This Loop Header: Depth=1
	mov	x19, x8
	b	LBB10_3
LBB10_2:                                ;   in Loop: Header=BB10_3 Depth=2
	sub	x19, x19, #1
LBB10_3:                                ;   Parent Loop BB10_1 Depth=1
	tbz	w25, #0, LBB10_6
	lsr	x8, x24, x19
	tbnz	w8, #0, LBB10_9
	cbnz	x19, LBB10_2
	b	LBB10_13
LBB10_6:                                ;   in Loop: Header=BB10_3 Depth=2
	sub	x0, x29, #208
	add	x1, sp, #208
	sub	x2, x29, #144
	sub	x3, x29, #176
	sub	x4, x29, #144
	sub	x5, x29, #176
	bl	_v15q_add
	cbnz	w0, LBB10_30
	ldp	q0, q1, [x29, #-208]
	stp	q0, q1, [x29, #-144]
	ldp	q0, q1, [sp, #208]
	stp	q0, q1, [x29, #-176]
	lsr	x8, x24, x19
	tbnz	w8, #0, LBB10_10
	cbnz	x19, LBB10_2
	b	LBB10_12
LBB10_9:                                ;   in Loop: Header=BB10_1 Depth=1
	mov	x22, x27
	mov	x21, x26
	b	LBB10_11
LBB10_10:                               ;   in Loop: Header=BB10_1 Depth=1
	sub	x22, x29, #208
	add	x21, sp, #208
	sub	x0, x29, #208
	add	x1, sp, #208
	sub	x2, x29, #144
	sub	x3, x29, #176
	mov	x4, x27
	mov	x5, x26
	bl	_v15q_add
	cbnz	w0, LBB10_30
LBB10_11:                               ;   in Loop: Header=BB10_1 Depth=1
	mov	w25, #0                         ; =0x0
	ldp	q0, q1, [x22]
	stp	q0, q1, [x29, #-144]
	ldp	q0, q1, [x21]
	stp	q0, q1, [x29, #-176]
	sub	x8, x19, #1
	cbnz	x19, LBB10_1
LBB10_12:
	ldur	x8, [x29, #-144]
	ldur	x9, [x29, #-120]
	str	x8, [sp, #176]
	ldur	q0, [x29, #-176]
	str	q0, [sp, #144]
	sub	x8, x29, #144
	ldur	q0, [x8, #8]
	stur	q0, [sp, #184]
	ldp	x10, x8, [x29, #-160]
	str	x10, [sp, #160]
	str	x9, [sp, #200]
	b	LBB10_14
LBB10_13:
	mov	x8, #0                          ; =0x0
	stp	xzr, xzr, [sp, #176]
	stp	xzr, xzr, [sp, #152]
	str	xzr, [sp, #144]
	stp	xzr, xzr, [sp, #192]
LBB10_14:
	str	x8, [sp, #168]
	mov	w21, #1                         ; =0x1
	mov	w8, #15                         ; =0xf
LBB10_15:                               ; =>This Loop Header: Depth=1
	mov	x19, x8
	b	LBB10_17
LBB10_16:                               ;   in Loop: Header=BB10_17 Depth=2
	sub	x19, x19, #1
LBB10_17:                               ;   Parent Loop BB10_15 Depth=1
	tbz	w21, #0, LBB10_20
	lsr	x8, x23, x19
	tbnz	w8, #0, LBB10_23
	cbnz	x19, LBB10_16
	b	LBB10_27
LBB10_20:                               ;   in Loop: Header=BB10_17 Depth=2
	sub	x0, x29, #208
	add	x1, sp, #208
	sub	x2, x29, #144
	sub	x3, x29, #176
	sub	x4, x29, #144
	sub	x5, x29, #176
	bl	_v15q_add
	cbnz	w0, LBB10_30
	ldp	q0, q1, [x29, #-208]
	stp	q0, q1, [x29, #-144]
	ldp	q0, q1, [sp, #208]
	stp	q0, q1, [x29, #-176]
	lsr	x8, x23, x19
	tbnz	w8, #0, LBB10_24
	cbnz	x19, LBB10_16
	b	LBB10_26
LBB10_23:                               ;   in Loop: Header=BB10_15 Depth=1
	mov	x25, x20
	mov	x22, x28
	b	LBB10_25
LBB10_24:                               ;   in Loop: Header=BB10_15 Depth=1
	sub	x25, x29, #208
	add	x22, sp, #208
	sub	x0, x29, #208
	add	x1, sp, #208
	sub	x2, x29, #144
	sub	x3, x29, #176
	mov	x4, x20
	mov	x5, x28
	bl	_v15q_add
	cbnz	w0, LBB10_30
LBB10_25:                               ;   in Loop: Header=BB10_15 Depth=1
	mov	w21, #0                         ; =0x0
	ldp	q0, q1, [x25]
	stp	q0, q1, [x29, #-144]
	ldp	q0, q1, [x22]
	stp	q0, q1, [x29, #-176]
	sub	x8, x19, #1
	cbnz	x19, LBB10_15
LBB10_26:
	ldur	x8, [x29, #-144]
	ldur	x9, [x29, #-120]
	str	x8, [sp, #112]
	ldur	q0, [x29, #-176]
	str	q0, [sp, #80]
	sub	x8, x29, #144
	ldur	q0, [x8, #8]
	stur	q0, [sp, #120]
	ldp	x10, x8, [x29, #-160]
	str	x10, [sp, #96]
	str	x9, [sp, #136]
	b	LBB10_28
LBB10_27:
	mov	x8, #0                          ; =0x0
	stp	xzr, xzr, [sp, #112]
	stp	xzr, xzr, [sp, #88]
	str	xzr, [sp, #80]
	stp	xzr, xzr, [sp, #128]
LBB10_28:
	str	x8, [sp, #104]
	sub	x0, x29, #144
	sub	x1, x29, #176
	add	x4, sp, #176
	add	x5, sp, #144
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	mov	x2, x20
	mov	x3, x19
	bl	_v15q_add
	cbnz	w0, LBB10_30
	sub	x0, x29, #208
	add	x1, sp, #208
	sub	x2, x29, #144
	sub	x3, x29, #176
	add	x4, sp, #112
	add	x5, sp, #80
	bl	_v15q_add
	cbz	w0, LBB10_33
LBB10_30:
	mov	w0, #-1                         ; =0xffffffff
LBB10_31:
	ldur	x8, [x29, #-104]
Lloh53:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh54:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh55:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB10_40
	ldp	x29, x30, [sp, #448]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #432]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #416]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #400]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #384]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #368]            ; 16-byte Folded Reload
	add	sp, sp, #464
	ret
LBB10_33:
	ldp	q0, q1, [x29, #-208]
	stp	q0, q1, [x29, #-144]
	ldp	q0, q1, [sp, #208]
	stp	q0, q1, [x29, #-176]
	ldp	x8, x9, [x29, #-144]
	str	x8, [x20]
	ldp	x8, x10, [x29, #-176]
	str	x8, [x19]
	str	x9, [x20, #8]
	str	x10, [x19, #8]
	ldp	x8, x9, [x29, #-128]
	str	x8, [x20, #16]
	ldp	x8, x10, [x29, #-160]
	str	x8, [x19, #16]
	str	x9, [x20, #24]
	str	x10, [x19, #24]
	ldr	x11, [sp, #56]                  ; 8-byte Folded Reload
	ldp	x8, x10, [x11]
	adds	x9, x24, x8
	adcs	x10, x10, xzr
	ldp	x8, x11, [x11, #16]
	adcs	x12, x8, xzr
	adcs	x13, x11, xzr
	cset	w0, hs
Lloh56:
	adrp	x8, _R_CONST@PAGE
Lloh57:
	add	x8, x8, _R_CONST@PAGEOFF
	ldp	x11, x14, [x8]
	subs	x11, x9, x11
	ngc	x15, xzr
	subs	x14, x10, x14
	ngc	x16, xzr
	adds	x14, x14, x15
	adc	x15, x16, x15
	asr	x16, x15, #63
	ldp	x17, x1, [x8, #16]
	subs	x17, x12, x17
	ngc	x2, xzr
	adds	x15, x17, x15
	adc	x16, x2, x16
	asr	x17, x16, #63
	subs	x1, x13, x1
	ngc	x2, xzr
	adds	x16, x1, x16
	adc	x17, x2, x17
	tbnz	w0, #0, LBB10_36
	and	x17, x17, #0x1
	tbz	w17, #0, LBB10_36
	mov	x11, x9
	mov	x14, x10
	mov	x15, x12
	mov	x16, x13
LBB10_36:
	ldp	x5, x9, [sp, #48]               ; 16-byte Folded Reload
	stp	x11, x14, [x9]
	stp	x15, x16, [x9, #16]
	ldp	x9, x10, [x5]
	adds	x9, x23, x9
	adcs	x11, x10, xzr
	ldp	x10, x12, [x5, #16]
	adcs	x13, x10, xzr
	adcs	x14, x12, xzr
	cset	w17, hs
	ldp	x10, x12, [x8]
	subs	x10, x9, x10
	ngc	x15, xzr
	subs	x12, x11, x12
	ngc	x16, xzr
	adds	x12, x12, x15
	adc	x15, x16, x15
	asr	x16, x15, #63
	ldp	x0, x1, [x8, #16]
	subs	x8, x13, x0
	ngc	x0, xzr
	adds	x8, x8, x15
	adc	x15, x0, x16
	asr	x16, x15, #63
	subs	x0, x14, x1
	ngc	x1, xzr
	adds	x15, x0, x15
	adc	x16, x1, x16
	tbnz	w17, #0, LBB10_39
	and	x16, x16, #0x1
	tbz	w16, #0, LBB10_39
	mov	x10, x9
	mov	x12, x11
	mov	x8, x13
	mov	x15, x14
LBB10_39:
	ldp	x3, x4, [sp, #32]               ; 16-byte Folded Reload
	ldp	x1, x2, [sp, #16]               ; 16-byte Folded Reload
	ldr	x0, [sp, #8]                    ; 8-byte Folded Reload
	stp	x10, x12, [x5]
	stp	x8, x15, [x5, #16]
	bl	_v15r_reset_hist
	mov	w0, #0                          ; =0x0
	b	LBB10_31
LBB10_40:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh50, Lloh51, Lloh52
	.loh AdrpLdrGotLdr	Lloh53, Lloh54, Lloh55
	.loh AdrpAdd	Lloh56, Lloh57
	.cfi_endproc
	.globl	_v15q_step                      ; -- Begin function v15q_step
	.p2align	2
_v15q_step:                             ; @v15q_step
	.cfi_startproc
	sub	sp, sp, #384
	stp	x28, x27, [sp, #288]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #304]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #320]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #336]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #352]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #368]            ; 16-byte Folded Spill
	add	x29, sp, #368
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
	mov	x19, x1
	mov	x20, x0
Lloh58:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh59:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh60:
	ldr	x8, [x8]
	stur	x8, [x29, #-96]
	mov	x8, #-6148914691236517206       ; =0xaaaaaaaaaaaaaaaa
	movk	x8, #43691
	umulh	x8, x2, x8
	lsr	x8, x8, #12
	mov	w9, #6144                       ; =0x1800
Lloh61:
	adrp	x10, _TBL@GOTPAGE
Lloh62:
	ldr	x10, [x10, _TBL@GOTPAGEOFF]
	msub	x8, x8, x9, x2
	add	x8, x10, x8, lsl #6
	ldp	x7, x21, [x8]
	ldp	x26, x25, [x0]
	subs	x12, x7, x26
	ngc	x9, xzr
	subs	x10, x21, x25
	ngc	x11, xzr
	adds	x13, x10, x9
	adc	x9, x11, x9
	asr	x10, x9, #63
	ldp	x23, x24, [x8, #16]
	ldp	x28, x27, [x0, #16]
	subs	x11, x23, x28
	ngc	x14, xzr
	adds	x11, x11, x9
	adc	x9, x14, x10
	asr	x10, x9, #63
	subs	x14, x24, x27
	ngc	x15, xzr
	adds	x9, x14, x9
	adc	x10, x15, x10
	mov	x14, #-977                      ; =0xfffffffffffffc2f
	movk	x14, #65534, lsl #32
	adds	x15, x12, x14
	adcs	x16, x13, xzr
	cset	w17, hs
	subs	x16, x16, #1
	adcs	x17, x11, x17
	cset	w0, hs
	subs	x17, x17, #1
	adc	x0, x9, x0
	sub	x0, x0, #1
	tst	w10, #0x1
	csel	x9, x9, x0, eq
	csel	x10, x11, x17, eq
	csel	x11, x13, x16, eq
	csel	x12, x12, x15, eq
	stp	x12, x11, [sp, #144]
	stp	x10, x9, [sp, #160]
	ldp	x13, x15, [x8, #32]
	ldp	x5, x3, [x1]
	subs	x13, x13, x5
	ngc	x16, xzr
	subs	x15, x15, x3
	ngc	x17, xzr
	adds	x15, x15, x16
	adc	x16, x17, x16
	asr	x17, x16, #63
	ldp	x0, x8, [x8, #48]
	ldp	x6, x4, [x1, #16]
	subs	x0, x0, x6
	ngc	x1, xzr
	adds	x16, x0, x16
	adc	x17, x1, x17
	asr	x0, x17, #63
	subs	x8, x8, x4
	ngc	x1, xzr
	adds	x8, x8, x17
	adc	x17, x1, x0
	adds	x14, x13, x14
	adcs	x0, x15, xzr
	cset	w1, hs
	subs	x0, x0, #1
	adcs	x1, x16, x1
	cset	w2, hs
	subs	x1, x1, #1
	adc	x2, x8, x2
	sub	x2, x2, #1
	tst	w17, #0x1
	csel	x8, x8, x2, eq
	csel	x16, x16, x1, eq
	csel	x15, x15, x0, eq
	csel	x13, x13, x14, eq
	stp	x13, x15, [sp, #112]
	stp	x16, x8, [sp, #128]
	orr	x8, x10, x9
	orr	x9, x11, x12
	orr	x8, x8, x9
	cbz	x8, LBB11_14
	stp	x21, x7, [sp]                   ; 16-byte Folded Spill
	stp	x6, x5, [sp, #16]               ; 16-byte Folded Spill
	stp	x4, x3, [sp, #32]               ; 16-byte Folded Spill
Lloh63:
	adrp	x8, l___const.rinv.one@PAGE
Lloh64:
	add	x8, x8, l___const.rinv.one@PAGEOFF
	ldp	q0, q1, [x8]
	stp	q0, q1, [x29, #-128]
	mov	w21, #255                       ; =0xff
Lloh65:
	adrp	x22, _P_EXP@PAGE
Lloh66:
	add	x22, x22, _P_EXP@PAGEOFF
	b	LBB11_3
LBB11_2:                                ;   in Loop: Header=BB11_3 Depth=1
	sub	w21, w21, #1
	cmn	w21, #1
	b.eq	LBB11_5
LBB11_3:                                ; =>This Inner Loop Header: Depth=1
	sub	x0, x29, #160
	sub	x1, x29, #128
	sub	x2, x29, #128
	bl	_fmul
	ldp	x8, x9, [x29, #-160]
	stp	x8, x9, [x29, #-128]
	ldp	x10, x11, [x29, #-144]
	stp	x10, x11, [x29, #-112]
	lsr	w12, w21, #6
	ldr	x12, [x22, w12, uxtw #3]
	lsr	x12, x12, x21
	tbz	w12, #0, LBB11_2
	add	x0, sp, #176
	sub	x1, x29, #128
	add	x2, sp, #144
	bl	_fmul
	ldp	x8, x9, [sp, #176]
	stp	x8, x9, [x29, #-128]
	ldp	x10, x11, [sp, #192]
	stp	x10, x11, [x29, #-112]
	b	LBB11_2
LBB11_5:
	stp	x8, x9, [x29, #-160]
	stp	x10, x11, [x29, #-144]
	sub	x0, x29, #128
	add	x1, sp, #112
	sub	x2, x29, #160
	bl	_fmul
	add	x0, sp, #176
	sub	x1, x29, #128
	sub	x2, x29, #128
	bl	_fmul
	ldp	x8, x9, [sp, #176]
	subs	x8, x8, x26
	ngc	x10, xzr
	subs	x9, x9, x25
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	asr	x12, x10, #63
	ldp	x11, x13, [sp, #192]
	subs	x11, x11, x28
	ngc	x14, xzr
	adds	x11, x11, x10
	adc	x10, x14, x12
	asr	x12, x10, #63
	subs	x13, x13, x27
	ngc	x14, xzr
	adds	x10, x13, x10
	adc	x12, x14, x12
	tbz	w12, #0, LBB11_7
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x9, x9, xzr
	cset	w12, hs
	subs	x9, x9, #1
	adcs	x11, x11, x12
	cset	w12, hs
	subs	x11, x11, #1
	adc	x10, x10, x12
	sub	x10, x10, #1
LBB11_7:
	ldp	x12, x13, [sp]                  ; 16-byte Folded Reload
	subs	x21, x8, x13
	ngc	x8, xzr
	subs	x9, x9, x12
	ngc	x12, xzr
	adds	x22, x9, x8
	adc	x8, x12, x8
	asr	x9, x8, #63
	subs	x11, x11, x23
	ngc	x12, xzr
	adds	x23, x11, x8
	adc	x8, x12, x9
	asr	x9, x8, #63
	subs	x10, x10, x24
	ngc	x11, xzr
	adds	x24, x10, x8
	adc	x8, x11, x9
	tbz	w8, #0, LBB11_9
	mov	x8, #-977                       ; =0xfffffffffffffc2f
	movk	x8, #65534, lsl #32
	adds	x21, x21, x8
	adcs	x8, x22, xzr
	cset	w9, hs
	subs	x22, x8, #1
	adcs	x8, x23, x9
	cset	w9, hs
	subs	x23, x8, #1
	adc	x8, x24, x9
	sub	x24, x8, #1
LBB11_9:
	subs	x8, x26, x21
	ngc	x10, xzr
	subs	x9, x25, x22
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	asr	x11, x10, #63
	subs	x12, x28, x23
	ngc	x13, xzr
	adds	x10, x12, x10
	adc	x11, x13, x11
	asr	x12, x11, #63
	subs	x13, x27, x24
	ngc	x14, xzr
	adds	x11, x13, x11
	adc	x12, x14, x12
	tbz	w12, #0, LBB11_11
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x9, x9, xzr
	cset	w12, hs
	subs	x9, x9, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x11, x11, x12
	sub	x11, x11, #1
LBB11_11:
	stp	x8, x9, [sp, #48]
	stp	x10, x11, [sp, #64]
	add	x0, sp, #80
	sub	x1, x29, #128
	add	x2, sp, #48
	bl	_fmul
	ldp	x8, x9, [sp, #80]
	ldp	x13, x10, [sp, #16]             ; 16-byte Folded Reload
	subs	x8, x8, x10
	ngc	x10, xzr
	ldp	x14, x11, [sp, #32]             ; 16-byte Folded Reload
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x11, x11, x10
	ldp	x10, x12, [sp, #96]
	subs	x10, x10, x13
	ngc	x13, xzr
	adds	x10, x10, x11
	asr	x11, x11, #63
	adc	x11, x13, x11
	asr	x13, x11, #63
	subs	x12, x12, x14
	ngc	x14, xzr
	adds	x11, x12, x11
	adc	x12, x14, x13
	tbz	w12, #0, LBB11_13
	mov	x12, #-977                      ; =0xfffffffffffffc2f
	movk	x12, #65534, lsl #32
	adds	x8, x8, x12
	adcs	x9, x9, xzr
	cset	w12, hs
	subs	x9, x9, #1
	adcs	x10, x10, x12
	cset	w12, hs
	subs	x10, x10, #1
	adc	x11, x11, x12
	sub	x11, x11, #1
LBB11_13:
	mov	w0, #0                          ; =0x0
	stp	x21, x22, [x20]
	stp	x23, x24, [x20, #16]
	stp	x8, x9, [x19]
	stp	x10, x11, [x19, #16]
	b	LBB11_15
LBB11_14:
	mov	w0, #-1                         ; =0xffffffff
LBB11_15:
	ldur	x8, [x29, #-96]
Lloh67:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh68:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh69:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB11_17
	ldp	x29, x30, [sp, #368]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #352]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #336]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #320]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #304]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #288]            ; 16-byte Folded Reload
	add	sp, sp, #384
	ret
LBB11_17:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh61, Lloh62
	.loh AdrpLdrGotLdr	Lloh58, Lloh59, Lloh60
	.loh AdrpAdd	Lloh65, Lloh66
	.loh AdrpAdd	Lloh63, Lloh64
	.loh AdrpLdrGotLdr	Lloh67, Lloh68, Lloh69
	.cfi_endproc
	.globl	_v15q_solve                     ; -- Begin function v15q_solve
	.p2align	2
_v15q_solve:                            ; @v15q_solve
	.cfi_startproc
	sub	sp, sp, #224
	stp	x22, x21, [sp, #176]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #192]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #208]            ; 16-byte Folded Spill
	add	x29, sp, #208
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	mov	x19, x4
Lloh70:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh71:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh72:
	ldr	x8, [x8]
	stur	x8, [x29, #-40]
	ldp	x8, x9, [x0]
	ldp	x10, x11, [x2]
	subs	x8, x8, x10
	ngc	x10, xzr
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	stp	x8, x9, [sp, #32]
	asr	x12, x10, #63
	ldp	x11, x13, [x0, #16]
	ldp	x14, x15, [x2, #16]
	subs	x11, x11, x14
	ngc	x14, xzr
	adds	x11, x11, x10
	adc	x10, x14, x12
	asr	x12, x10, #63
	subs	x13, x13, x15
	ngc	x14, xzr
	adds	x10, x13, x10
	adc	x12, x14, x12
	stp	x11, x10, [sp, #48]
	tbz	w12, #0, LBB12_2
Lloh73:
	adrp	x12, _R_CONST@PAGE
Lloh74:
	add	x12, x12, _R_CONST@PAGEOFF
	ldp	x13, x14, [x12]
	adds	x8, x8, x13
	adcs	x9, x9, xzr
	cset	w13, hs
	ldp	x15, x12, [x12, #16]
	adds	x11, x11, x15
	cset	w15, hs
	adds	x9, x9, x14
	stp	x8, x9, [sp, #32]
	adcs	x8, x11, x13
	adc	x9, x10, x15
	add	x9, x12, x9
	stp	x8, x9, [sp, #48]
LBB12_2:
	ldp	x8, x9, [x1]
	ldp	x10, x11, [x3]
	subs	x8, x8, x10
	ngc	x10, xzr
	subs	x9, x9, x11
	ngc	x11, xzr
	adds	x9, x9, x10
	adc	x10, x11, x10
	stp	x8, x9, [sp]
	asr	x12, x10, #63
	ldp	x11, x13, [x1, #16]
	ldp	x14, x15, [x3, #16]
	subs	x11, x11, x14
	ngc	x14, xzr
	adds	x11, x11, x10
	adc	x10, x14, x12
	asr	x12, x10, #63
	subs	x13, x13, x15
	ngc	x14, xzr
	adds	x10, x13, x10
	adc	x12, x14, x12
	stp	x11, x10, [sp, #16]
	tbz	w12, #0, LBB12_4
Lloh75:
	adrp	x12, _R_CONST@PAGE
Lloh76:
	add	x12, x12, _R_CONST@PAGEOFF
	ldp	x13, x14, [x12]
	adds	x8, x8, x13
	adcs	x9, x9, xzr
	cset	w13, hs
	ldp	x15, x12, [x12, #16]
	adds	x11, x11, x15
	cset	w15, hs
	adds	x9, x9, x14
	stp	x8, x9, [sp]
	adcs	x11, x11, x13
	adc	x10, x10, x15
	add	x10, x12, x10
	stp	x11, x10, [sp, #16]
LBB12_4:
	orr	x9, x11, x9
	orr	x8, x9, x8
	orr	x8, x8, x10
	cbz	x8, LBB12_11
Lloh77:
	adrp	x2, _R2_R@PAGE
Lloh78:
	add	x2, x2, _R2_R@PAGEOFF
	add	x0, sp, #64
	mov	x1, sp
	bl	_rmont
Lloh79:
	adrp	x8, _R_MONT_ONE@PAGE
Lloh80:
	add	x8, x8, _R_MONT_ONE@PAGEOFF
	ldp	q0, q1, [x8]
	stp	q0, q1, [x29, #-80]
	mov	w20, #255                       ; =0xff
Lloh81:
	adrp	x21, _R_EXP@PAGE
Lloh82:
	add	x21, x21, _R_EXP@PAGEOFF
	b	LBB12_7
LBB12_6:                                ;   in Loop: Header=BB12_7 Depth=1
	sub	w20, w20, #1
	cmn	w20, #1
	b.eq	LBB12_9
LBB12_7:                                ; =>This Inner Loop Header: Depth=1
	add	x0, sp, #96
	sub	x1, x29, #80
	sub	x2, x29, #80
	bl	_rmont
	ldp	q0, q1, [sp, #96]
	stp	q0, q1, [x29, #-80]
	lsr	w8, w20, #6
	ldr	x8, [x21, w8, uxtw #3]
	lsr	x8, x8, x20
	tbz	w8, #0, LBB12_6
	sub	x0, x29, #80
	add	x1, sp, #96
	add	x2, sp, #64
	bl	_rmont
	b	LBB12_6
LBB12_9:
Lloh83:
	adrp	x2, l___const.rinv.one@PAGE
Lloh84:
	add	x2, x2, l___const.rinv.one@PAGEOFF
	add	x0, sp, #96
	sub	x1, x29, #80
	bl	_rmont
	sub	x0, x29, #80
	add	x1, sp, #32
	add	x2, sp, #96
	bl	_rmont
Lloh85:
	adrp	x2, _R2_R@PAGE
Lloh86:
	add	x2, x2, _R2_R@PAGEOFF
	add	x0, sp, #64
	sub	x1, x29, #80
	bl	_rmont
	ldp	x11, x10, [sp, #64]
	orr	x12, x10, x11
	ldp	x9, x8, [sp, #80]
	orr	x13, x9, x8
	orr	x12, x12, x13
	cbz	x12, LBB12_11
	mov	w0, #0                          ; =0x0
Lloh87:
	adrp	x12, _R_CONST@PAGE
Lloh88:
	add	x12, x12, _R_CONST@PAGEOFF
	ldp	x13, x14, [x12]
	subs	x11, x13, x11
	ngc	x13, xzr
	subs	x10, x14, x10
	ngc	x14, xzr
	adds	x10, x10, x13
	adc	x13, x14, x13
	ldp	x14, x12, [x12, #16]
	subs	x9, x14, x9
	ngc	x14, xzr
	adds	x9, x9, x13
	asr	x13, x13, #63
	adc	x13, x14, x13
	stp	x11, x10, [x19]
	sub	x8, x12, x8
	add	x8, x8, x13
	stp	x9, x8, [x19, #16]
	b	LBB12_12
LBB12_11:
	mov	w0, #-1                         ; =0xffffffff
LBB12_12:
	ldur	x8, [x29, #-40]
Lloh89:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh90:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh91:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB12_14
	ldp	x29, x30, [sp, #208]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #192]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #176]            ; 16-byte Folded Reload
	add	sp, sp, #224
	ret
LBB12_14:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh70, Lloh71, Lloh72
	.loh AdrpAdd	Lloh73, Lloh74
	.loh AdrpAdd	Lloh75, Lloh76
	.loh AdrpAdd	Lloh81, Lloh82
	.loh AdrpAdd	Lloh79, Lloh80
	.loh AdrpAdd	Lloh77, Lloh78
	.loh AdrpAdd	Lloh85, Lloh86
	.loh AdrpAdd	Lloh83, Lloh84
	.loh AdrpAdd	Lloh87, Lloh88
	.loh AdrpLdrGotLdr	Lloh89, Lloh90, Lloh91
	.cfi_endproc
	.globl	_v15r_dbg                       ; -- Begin function v15r_dbg
	.p2align	2
_v15r_dbg:                              ; @v15r_dbg
	.cfi_startproc
	sub	sp, sp, #80
	stp	x20, x19, [sp, #48]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #64]             ; 16-byte Folded Spill
	add	x29, sp, #64
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	mov	x19, x0
Lloh92:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh93:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh94:
	ldr	x8, [x8]
	stur	x8, [x29, #-24]
	add	x0, sp, #8
	bl	_rmont
Lloh95:
	adrp	x2, _R2_R@PAGE
Lloh96:
	add	x2, x2, _R2_R@PAGEOFF
	add	x1, sp, #8
	mov	x0, x19
	bl	_rmont
	ldur	x8, [x29, #-24]
Lloh97:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh98:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh99:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB13_2
	ldp	x29, x30, [sp, #64]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #80
	ret
LBB13_2:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh97, Lloh98, Lloh99
	.loh AdrpAdd	Lloh95, Lloh96
	.loh AdrpLdrGotLdr	Lloh92, Lloh93, Lloh94
	.cfi_endproc
	.globl	_v15r_dbg_inv                   ; -- Begin function v15r_dbg_inv
	.p2align	2
_v15r_dbg_inv:                          ; @v15r_dbg_inv
	.cfi_startproc
	sub	sp, sp, #160
	stp	x22, x21, [sp, #112]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #128]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #144]            ; 16-byte Folded Spill
	add	x29, sp, #144
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	.cfi_offset w21, -40
	.cfi_offset w22, -48
	mov	x19, x0
Lloh100:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh101:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh102:
	ldr	x8, [x8]
	stur	x8, [x29, #-40]
Lloh103:
	adrp	x2, _R2_R@PAGE
Lloh104:
	add	x2, x2, _R2_R@PAGEOFF
	mov	x0, sp
	bl	_rmont
Lloh105:
	adrp	x8, _R_MONT_ONE@PAGE
Lloh106:
	add	x8, x8, _R_MONT_ONE@PAGEOFF
	ldp	q0, q1, [x8]
	stp	q0, q1, [sp, #64]
	mov	w20, #255                       ; =0xff
Lloh107:
	adrp	x21, _R_EXP@PAGE
Lloh108:
	add	x21, x21, _R_EXP@PAGEOFF
	b	LBB14_2
LBB14_1:                                ;   in Loop: Header=BB14_2 Depth=1
	sub	w20, w20, #1
	cmn	w20, #1
	b.eq	LBB14_4
LBB14_2:                                ; =>This Inner Loop Header: Depth=1
	add	x0, sp, #32
	add	x1, sp, #64
	add	x2, sp, #64
	bl	_rmont
	ldp	q0, q1, [sp, #32]
	stp	q0, q1, [sp, #64]
	lsr	w8, w20, #6
	ldr	x8, [x21, w8, uxtw #3]
	lsr	x8, x8, x20
	tbz	w8, #0, LBB14_1
	add	x0, sp, #64
	add	x1, sp, #32
	mov	x2, sp
	bl	_rmont
	b	LBB14_1
LBB14_4:
Lloh109:
	adrp	x2, l___const.rinv.one@PAGE
Lloh110:
	add	x2, x2, l___const.rinv.one@PAGEOFF
	add	x1, sp, #64
	mov	x0, x19
	bl	_rmont
	ldur	x8, [x29, #-40]
Lloh111:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh112:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh113:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB14_6
	ldp	x29, x30, [sp, #144]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #128]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #112]            ; 16-byte Folded Reload
	add	sp, sp, #160
	ret
LBB14_6:
	bl	___stack_chk_fail
	.loh AdrpAdd	Lloh107, Lloh108
	.loh AdrpAdd	Lloh105, Lloh106
	.loh AdrpAdd	Lloh103, Lloh104
	.loh AdrpLdrGotLdr	Lloh100, Lloh101, Lloh102
	.loh AdrpLdrGotLdr	Lloh111, Lloh112, Lloh113
	.loh AdrpAdd	Lloh109, Lloh110
	.cfi_endproc
	.globl	_v15q_mul                       ; -- Begin function v15q_mul
	.p2align	2
_v15q_mul:                              ; @v15q_mul
	.cfi_startproc
	sub	sp, sp, #240
	stp	x28, x27, [sp, #144]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #160]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #176]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #192]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #208]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #224]            ; 16-byte Folded Spill
	add	x29, sp, #224
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
	mov	x21, x4
	mov	x22, x3
	mov	x23, x2
	mov	x19, x1
	mov	x20, x0
	mov	w25, #0                         ; =0x0
Lloh114:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh115:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh116:
	ldr	x8, [x8]
	stur	x8, [x29, #-88]
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp, #96]
	stp	q0, q0, [sp, #64]
	mov	w24, #255                       ; =0xff
	b	LBB15_4
LBB15_1:                                ;   in Loop: Header=BB15_4 Depth=1
	lsr	w8, w24, #6
	ldr	x8, [x23, w8, uxtw #3]
	lsr	x8, x8, x24
	mov	x27, x22
	mov	x26, x21
	tbz	w8, #0, LBB15_3
LBB15_2:                                ;   in Loop: Header=BB15_4 Depth=1
	ldp	q0, q1, [x27]
	stp	q0, q1, [sp, #96]
	ldp	q0, q1, [x26]
	stp	q0, q1, [sp, #64]
	mov	w25, #1                         ; =0x1
LBB15_3:                                ;   in Loop: Header=BB15_4 Depth=1
	sub	w24, w24, #1
	cmn	w24, #1
	b.eq	LBB15_8
LBB15_4:                                ; =>This Inner Loop Header: Depth=1
	cbz	w25, LBB15_1
	add	x0, sp, #32
	mov	x1, sp
	add	x2, sp, #96
	add	x3, sp, #64
	add	x4, sp, #96
	add	x5, sp, #64
	bl	_v15q_add
	cbnz	w0, LBB15_10
	ldp	q0, q1, [sp, #32]
	stp	q0, q1, [sp, #96]
	ldp	q0, q1, [sp]
	stp	q0, q1, [sp, #64]
	lsr	w8, w24, #6
	ldr	x8, [x23, w8, uxtw #3]
	lsr	x8, x8, x24
	tbz	w8, #0, LBB15_3
	add	x27, sp, #32
	mov	x26, sp
	add	x0, sp, #32
	mov	x1, sp
	add	x2, sp, #96
	add	x3, sp, #64
	mov	x4, x22
	mov	x5, x21
	bl	_v15q_add
	cbz	w0, LBB15_2
	b	LBB15_10
LBB15_8:
	cbz	w25, LBB15_10
	mov	w0, #0                          ; =0x0
	ldr	x8, [sp, #96]
	str	x8, [x20]
	ldur	q0, [sp, #104]
	stur	q0, [x20, #8]
	ldr	x8, [sp, #120]
	str	x8, [x20, #24]
	ldr	x8, [sp, #64]
	str	x8, [x19]
	ldur	q0, [sp, #72]
	stur	q0, [x19, #8]
	ldr	x8, [sp, #88]
	str	x8, [x19, #24]
	b	LBB15_11
LBB15_10:
	mov	w0, #-1                         ; =0xffffffff
LBB15_11:
	ldur	x8, [x29, #-88]
Lloh117:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh118:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh119:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB15_13
	ldp	x29, x30, [sp, #224]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #208]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #192]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #176]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #160]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #144]            ; 16-byte Folded Reload
	add	sp, sp, #240
	ret
LBB15_13:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh114, Lloh115, Lloh116
	.loh AdrpLdrGotLdr	Lloh117, Lloh118, Lloh119
	.cfi_endproc
	.section	__DATA,__data
	.globl	_NODE_MASK                      ; @NODE_MASK
	.p2align	3, 0x0
_NODE_MASK:
	.quad	1023                            ; 0x3ff
	.globl	_R_CONST                        ; @R_CONST
	.p2align	3, 0x0
_R_CONST:
	.quad	-4624529908474429119            ; 0xbfd25e8cd0364141
	.quad	-4994812053365940165            ; 0xbaaedce6af48a03b
	.quad	-2                              ; 0xfffffffffffffffe
	.quad	-1                              ; 0xffffffffffffffff
	.globl	_R_N0                           ; @R_N0
.zerofill __DATA,__common,_R_N0,8,3
	.globl	_LAMBDA_R                       ; @LAMBDA_R
	.p2align	3, 0x0
_LAMBDA_R:
	.quad	-1117478481910666082            ; 0xf07deb3dc9926c9e
	.quad	3212165691671483468             ; 0x2c93e7ad83c6944c
	.quad	8334304762764295569             ; 0x73a9660652697d91
	.quad	5992109773982062137             ; 0x532840178558d639
	.globl	_LAMBDA2_R                      ; @LAMBDA2_R
	.p2align	3, 0x0
_LAMBDA2_R:
	.quad	-8131581335038192156            ; 0x8f26d1dbd6da15e4
	.quad	5244954275306187818             ; 0x48c9d21fdacaac2a
	.quad	-8334304762764295572            ; 0x8c5699f9ad96826c
	.quad	-5992109773982062138            ; 0xacd7bfe87aa729c6
	.globl	_R_MONT_ONE                     ; @R_MONT_ONE
	.p2align	3, 0x0
_R_MONT_ONE:
	.quad	4624529908474429119             ; 0x402da1732fc9bebf
	.quad	4994812053365940164             ; 0x4551231950b75fc4
	.quad	1                               ; 0x1
	.quad	0                               ; 0x0
	.globl	_R2_R                           ; @R2_R
	.p2align	3, 0x0
_R2_R:
	.quad	-8544188223573208768            ; 0x896cf21467d7d140
	.quad	8364476168144746616             ; 0x741496c20e7cf878
	.quad	-1830724362361305146            ; 0xe697f5e45bcd07c6
	.quad	-7104678183822779451            ; 0x9d671cd581c69bc5
	.section	__TEXT,__const
	.globl	_BETA                           ; @BETA
	.p2align	3, 0x0
_BETA:
	.quad	-4523465429756870162            ; 0xc1396c28719501ee
	.quad	-7138124642204153451            ; 0x9cf0497512f58995
	.quad	7954561588662645993             ; 0x6e64479eac3434e9
	.quad	8856726876819556112             ; 0x7ae96a2b657c0710
	.globl	_P_EXP                          ; @P_EXP
	.p2align	3, 0x0
_P_EXP:
	.quad	-4294968275                     ; 0xfffffffefffffc2d
	.quad	-1                              ; 0xffffffffffffffff
	.quad	-1                              ; 0xffffffffffffffff
	.quad	-1                              ; 0xffffffffffffffff
	.globl	_R_EXP                          ; @R_EXP
	.p2align	3, 0x0
_R_EXP:
	.quad	-4624529908474429121            ; 0xbfd25e8cd036413f
	.quad	-4994812053365940165            ; 0xbaaedce6af48a03b
	.quad	-2                              ; 0xfffffffffffffffe
	.quad	-1                              ; 0xffffffffffffffff
	.comm	_SMALL,1024,3                   ; @SMALL
	.comm	_UBASE,131072,3                 ; @UBASE
	.comm	_TBL,393216,3                   ; @TBL
	.p2align	3, 0x0                          ; @v15q_add.TWO
_v15q_add.TWO:
	.quad	2                               ; 0x2
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
	.p2align	3, 0x0                          ; @v15q_add.THREE
_v15q_add.THREE:
	.quad	3                               ; 0x3
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
	.comm	_NODE_HEAD,98304,3              ; @NODE_HEAD
	.p2align	3, 0x0                          ; @__const.rinv.one
l___const.rinv.one:
	.quad	1                               ; 0x1
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
	.quad	0                               ; 0x0
.subsections_via_symbols
