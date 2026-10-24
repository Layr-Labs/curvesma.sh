	.build_version macos, 27, 0	sdk_version 27, 0
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_v15o_align                     ; -- Begin function v15o_align
	.p2align	2
_v15o_align:                            ; @v15o_align
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #112
	stp	x20, x19, [sp, #80]             ; 16-byte Folded Spill
	stp	x29, x30, [sp, #96]             ; 16-byte Folded Spill
	add	x29, sp, #96
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
	.cfi_offset w19, -24
	.cfi_offset w20, -32
	mov	x20, x1
	mov	x19, x0
Lloh0:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh1:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh2:
	ldr	x8, [x8]
	stur	x8, [x29, #-24]
	add	x0, sp, #40
	mov	x1, x19
	bl	_fmul
	add	x0, sp, #8
	add	x1, sp, #40
	bl	_fmul
	ldr	x10, [x19, #24]
	ldr	x12, [sp, #64]
	cmp	x12, x10
	b.hs	LBB0_3
; %bb.1:
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
; %bb.6:
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
; %bb.10:
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
	ldur	x8, [x29, #-24]
Lloh3:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh4:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh5:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB0_23
; %bb.12:
	ldp	x29, x30, [sp, #96]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #112
	ret
LBB0_13:
	ldr	x13, [sp, #56]
	cmp	x13, x11
	b.lo	LBB0_2
; %bb.14:
	b.hi	LBB0_4
; %bb.15:
	ldr	x14, [sp, #48]
	cmp	x14, x9
	b.lo	LBB0_2
; %bb.16:
	mov	x0, #0                          ; =0x0
	cmp	x14, x9
	b.hi	LBB0_5
; %bb.17:
	ldr	x14, [sp, #40]
	cmp	x14, x8
	b.lo	LBB0_2
	b	LBB0_5
LBB0_18:
	ldr	x13, [sp, #24]
	cmp	x13, x11
	b.lo	LBB0_7
; %bb.19:
	b.hi	LBB0_9
; %bb.20:
	ldr	x14, [sp, #16]
	cmp	x14, x9
	b.lo	LBB0_7
; %bb.21:
	b.hi	LBB0_9
; %bb.22:
	ldr	x14, [sp, #8]
	cmp	x14, x8
	b.lo	LBB0_7
	b	LBB0_9
LBB0_23:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh0, Lloh1, Lloh2
	.loh AdrpLdrGotLdr	Lloh3, Lloh4, Lloh5
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function fmul
_fmul:                                  ; @fmul
	.cfi_startproc
; %bb.0:
	mov	x8, #-1                         ; =0xffffffffffffffff
	mov	x9, #977                        ; =0x3d1
	movk	x9, #1, lsl #32
Lloh6:
	adrp	x13, _BETA@GOTPAGE
Lloh7:
	ldr	x13, [x13, _BETA@GOTPAGEOFF]
	ldp	x12, x11, [x13]
	ldp	x16, x17, [x1]
	mul	x10, x12, x16
	umulh	x14, x12, x16
	umulh	x15, x11, x16
	mul	x2, x11, x16
	adds	x2, x14, x2
	cinc	x3, x15, hs
	ldp	x15, x14, [x13, #16]
	umulh	x13, x15, x16
	mul	x4, x15, x16
	adds	x3, x3, x4
	cinc	x13, x13, hs
	umulh	x4, x14, x16
	mul	x16, x14, x16
	adds	x16, x13, x16
	cinc	x13, x4, hs
	cmp	x13, #0
	csel	x4, xzr, x13, eq
	umulh	x5, x12, x17
	mul	x13, x12, x17
	adds	x13, x13, x2
	cinc	x2, x5, hs
	umulh	x5, x11, x17
	mul	x6, x11, x17
	adds	x2, x2, x3
	cset	w3, hs
	adds	x2, x2, x6
	adc	x3, x3, x5
	umulh	x5, x15, x17
	mul	x6, x15, x17
	adds	x16, x3, x16
	cset	w3, hs
	adds	x16, x16, x6
	adc	x3, x3, x5
	umulh	x5, x14, x17
	mul	x17, x14, x17
	adds	x3, x3, x4
	cset	w4, hs
	adds	x17, x3, x17
	adcs	x3, x4, x5
	csel	x3, xzr, x3, eq
	ldp	x4, x1, [x1, #16]
	umulh	x5, x12, x4
	mul	x6, x12, x4
	adds	x2, x6, x2
	cinc	x5, x5, hs
	umulh	x6, x11, x4
	mul	x7, x11, x4
	adds	x16, x5, x16
	cset	w5, hs
	adds	x16, x16, x7
	adc	x5, x5, x6
	umulh	x6, x15, x4
	mul	x7, x15, x4
	adds	x17, x5, x17
	cset	w5, hs
	adds	x17, x17, x7
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
	adds	x12, x12, x16
	cinc	x16, x5, hs
	umulh	x5, x11, x1
	mul	x11, x11, x1
	adds	x16, x16, x17
	cset	w17, hs
	adds	x11, x16, x11
	adc	x16, x17, x5
	umulh	x17, x15, x1
	mul	x15, x15, x1
	adds	x16, x16, x3
	cset	w3, hs
	cmn	x16, x15
	adc	x17, x3, x17
	umulh	x3, x14, x1
	mul	x14, x14, x1
	adds	x17, x17, x4
	cset	w1, hs
	adds	x14, x17, x14
	adc	x17, x1, x3
	umulh	x1, x17, x9
	mul	x17, x17, x9
	adds	x12, x17, x12
	adcs	x11, x11, x1
	adc	x15, x16, x15
	umulh	x16, x14, x9
	mul	x14, x14, x9
	adds	x14, x14, x2
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
	.loh AdrpLdrGot	Lloh6, Lloh7
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_label                     ; -- Begin function v15o_label
	.p2align	2
_v15o_label:                            ; @v15o_label
	.cfi_startproc
; %bb.0:
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
; %bb.1:
	cbnz	x2, LBB2_4
; %bb.2:
Lloh11:
	adrp	x9, _R_MONT_ONE@GOTPAGE
Lloh12:
	ldr	x9, [x9, _R_MONT_ONE@GOTPAGEOFF]
	b	LBB2_5
LBB2_3:
Lloh13:
	adrp	x9, _LAMBDA_R@GOTPAGE
Lloh14:
	ldr	x9, [x9, _LAMBDA_R@GOTPAGEOFF]
	b	LBB2_5
LBB2_4:
Lloh15:
	adrp	x9, _LAMBDA2_R@GOTPAGE
Lloh16:
	ldr	x9, [x9, _LAMBDA2_R@GOTPAGEOFF]
LBB2_5:
	ldp	x13, x16, [x9]
	stp	x13, x16, [sp, #8]
	add	x8, x9, #24
	ldr	x2, [x9, #16]!
	ldr	x17, [x8]
	stp	x2, x17, [sp, #24]
Lloh17:
	adrp	x8, _R_CONST@GOTPAGE
Lloh18:
	ldr	x8, [x8, _R_CONST@GOTPAGEOFF]
	ldp	x12, x11, [x8]
	ldp	x10, x9, [x8, #16]
	mov	x14, #0                         ; =0x0
	mov	x15, #0                         ; =0x0
	cbz	w3, LBB2_7
; %bb.6:
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
; %bb.8:
	and	x14, x14, #0x1
	tbz	w14, #0, LBB2_10
; %bb.9:
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
; %bb.11:
	and	x16, x16, #0x1
	tbz	w16, #0, LBB2_13
; %bb.12:
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
; %bb.14:
	ldp	x29, x30, [sp, #64]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #48]             ; 16-byte Folded Reload
	add	sp, sp, #80
	ret
LBB2_15:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh8, Lloh9, Lloh10
	.loh AdrpLdrGot	Lloh11, Lloh12
	.loh AdrpLdrGot	Lloh13, Lloh14
	.loh AdrpLdrGot	Lloh15, Lloh16
	.loh AdrpLdrGot	Lloh17, Lloh18
	.loh AdrpLdrGotLdr	Lloh19, Lloh20, Lloh21
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function rmont
_rmont:                                 ; @rmont
	.cfi_startproc
; %bb.0:
	stp	x22, x21, [sp, #-32]!           ; 16-byte Folded Spill
	stp	x20, x19, [sp, #16]             ; 16-byte Folded Spill
	.cfi_def_cfa_offset 32
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	.cfi_offset w21, -24
	.cfi_offset w22, -32
Lloh22:
	adrp	x8, _R_N0@GOTPAGE
Lloh23:
	ldr	x8, [x8, _R_N0@GOTPAGEOFF]
Lloh24:
	ldr	x8, [x8]
Lloh25:
	adrp	x12, _R_CONST@GOTPAGE
Lloh26:
	ldr	x12, [x12, _R_CONST@GOTPAGEOFF]
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
; %bb.1:
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
	.loh AdrpLdrGot	Lloh25, Lloh26
	.loh AdrpLdrGotLdr	Lloh22, Lloh23, Lloh24
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_detect                    ; -- Begin function v15o_detect
	.p2align	2
_v15o_detect:                           ; @v15o_detect
	.cfi_startproc
; %bb.0:
	ldp	x8, x11, [x0]
	ldp	x12, x13, [x1]
	ldp	x10, x9, [x0, #16]
	ldp	x14, x15, [x1, #16]
	cmp	x8, x12
	ccmp	x11, x13, #0, eq
	ccmp	x10, x14, #0, eq
	ccmp	x9, x15, #0, eq
	b.eq	LBB4_5
; %bb.1:
	ldr	x12, [x2]
	cmp	x8, x12
	b.ne	LBB4_6
; %bb.2:
	ldr	x12, [x2, #8]
	cmp	x11, x12
	b.ne	LBB4_6
; %bb.3:
	ldr	x11, [x2, #16]
	cmp	x10, x11
	b.ne	LBB4_6
; %bb.4:
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
; %bb.11:
	b.hi	LBB4_9
; %bb.12:
	ldr	x9, [x0, #8]
	ldr	x10, [x4, #8]
	cmp	x9, x10
	b.lo	LBB4_7
; %bb.13:
	b.hi	LBB4_9
; %bb.14:
	ldr	x9, [x0]
	ldr	x10, [x4]
	cmp	x9, x10
	b.lo	LBB4_7
	b	LBB4_9
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_report                    ; -- Begin function v15o_report
	.p2align	2
_v15o_report:                           ; @v15o_report
	.cfi_startproc
; %bb.0:
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
; %bb.1:
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_4
; %bb.2:
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_4
; %bb.3:
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
; %bb.5:
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_8
; %bb.6:
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_8
; %bb.7:
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
; %bb.9:
	cmp	x14, x9
	ccmp	x15, x10, #0, eq
	b.ne	LBB5_12
; %bb.10:
	ldr	x14, [x2, #16]
	cmp	x13, x14
	b.ne	LBB5_12
; %bb.11:
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
; %bb.13:
	cmp	x13, x9
	ccmp	x14, x10, #0, eq
	b.ne	LBB5_17
; %bb.14:
	ldr	x10, [x2, #16]
	cmp	x12, x10
	b.ne	LBB5_17
; %bb.15:
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
                                        ; -- End function
	.globl	_v15r_init_for_Q                ; -- Begin function v15r_init_for_Q
	.p2align	2
_v15r_init_for_Q:                       ; @v15r_init_for_Q
	.cfi_startproc
; %bb.0:
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
Lloh27:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh28:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh29:
	ldr	x8, [x8]
	stur	x8, [x29, #-56]
	ldr	x8, [x1]
Lloh30:
	adrp	x22, _SMALL@GOTPAGE
Lloh31:
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
; %bb.4:                                ;   in Loop: Header=BB6_3 Depth=1
	ldp	x8, x9, [sp, #40]
	stp	x8, x9, [x21, #32]
	ldp	x8, x9, [sp, #8]
	stp	x8, x9, [x21, #64]
	ldp	x8, x9, [sp, #56]
	stp	x8, x9, [x21, #48]
	ldp	x8, x9, [sp, #24]
	stp	x8, x9, [x21, #80]
	b	LBB6_2
LBB6_5:
	mov	x21, #0                         ; =0x0
	add	x23, x20, #32
Lloh32:
	adrp	x20, _UBASE@GOTPAGE
Lloh33:
	ldr	x20, [x20, _UBASE@GOTPAGEOFF]
	mov	w24, #34953                     ; =0x8889
	b	LBB6_7
LBB6_6:                                 ;   in Loop: Header=BB6_7 Depth=1
	sub	x8, x23, #32
	ldp	x9, x10, [sp, #40]
	stp	x9, x10, [x8]
	ldp	x9, x10, [sp, #8]
	stp	x9, x10, [x23]
	ldp	x9, x10, [sp, #24]
	stp	x9, x10, [x23, #16]
	ldp	x9, x10, [sp, #56]
	add	x21, x21, #1
	add	x20, x20, #64
	add	x23, x23, #64
	stp	x9, x10, [x8, #16]
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
; %bb.8:                                ;   in Loop: Header=BB6_7 Depth=1
	add	x19, x19, #1
	str	xzr, [sp, #8]
	str	xzr, [sp, #40]
	b	LBB6_6
LBB6_9:
	ldur	x8, [x29, #-56]
Lloh34:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh35:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh36:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB6_11
; %bb.10:
	mov	x0, x19
	ldp	x29, x30, [sp, #128]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #112]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #96]             ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #144
	ret
LBB6_11:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh30, Lloh31
	.loh AdrpLdrGotLdr	Lloh27, Lloh28, Lloh29
	.loh AdrpLdrGot	Lloh32, Lloh33
	.loh AdrpLdrGotLdr	Lloh34, Lloh35, Lloh36
	.cfi_endproc
                                        ; -- End function
	.globl	_v15r_action                    ; -- Begin function v15r_action
	.p2align	2
_v15r_action:                           ; @v15r_action
	.cfi_startproc
; %bb.0:
	ldr	x8, [x0]
	cmp	x8, #15
	b.ls	LBB7_2
; %bb.1:
	mov	w0, #1                          ; =0x1
	ret
LBB7_2:
	ldr	x8, [x1]
	cmp	x8, #66
	cset	w0, hi
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_v15r_reset_hist                ; -- Begin function v15r_reset_hist
	.p2align	2
_v15r_reset_hist:                       ; @v15r_reset_hist
	.cfi_startproc
; %bb.0:
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
                                        ; -- End function
	.globl	_v15r_escape                    ; -- Begin function v15r_escape
	.p2align	2
_v15r_escape:                           ; @v15r_escape
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #448
	stp	x28, x27, [sp, #352]            ; 16-byte Folded Spill
	stp	x26, x25, [sp, #368]            ; 16-byte Folded Spill
	stp	x24, x23, [sp, #384]            ; 16-byte Folded Spill
	stp	x22, x21, [sp, #400]            ; 16-byte Folded Spill
	stp	x20, x19, [sp, #416]            ; 16-byte Folded Spill
	stp	x29, x30, [sp, #432]            ; 16-byte Folded Spill
	add	x29, sp, #432
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
	mov	x19, x6
	mov	x26, x5
	mov	x27, x4
	str	x3, [sp, #40]                   ; 8-byte Folded Spill
	stp	x2, x0, [sp, #56]               ; 16-byte Folded Spill
	str	x1, [sp, #72]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #64]
	str	x8, [sp, #48]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #56]
	str	x8, [sp, #32]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #48]
	str	x8, [sp, #24]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #40]
	str	x8, [sp, #16]                   ; 8-byte Folded Spill
	ldr	x8, [x29, #32]
	str	x8, [sp, #8]                    ; 8-byte Folded Spill
	ldp	x24, x23, [x29, #16]
Lloh37:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh38:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh39:
	ldr	x8, [x8]
	stur	x8, [x29, #-96]
	mov	w25, #1                         ; =0x1
	mov	w8, #15                         ; =0xf
LBB9_1:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB9_3 Depth 2
	mov	x20, x8
	b	LBB9_3
LBB9_2:                                 ;   in Loop: Header=BB9_3 Depth=2
	sub	x20, x20, #1
LBB9_3:                                 ;   Parent Loop BB9_1 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	tbz	w25, #0, LBB9_6
; %bb.4:                                ;   in Loop: Header=BB9_3 Depth=2
	lsr	x8, x24, x20
	tbnz	w8, #0, LBB9_9
; %bb.5:                                ;   in Loop: Header=BB9_3 Depth=2
	cbnz	x20, LBB9_2
	b	LBB9_13
LBB9_6:                                 ;   in Loop: Header=BB9_3 Depth=2
	sub	x0, x29, #192
	add	x1, sp, #208
	sub	x2, x29, #128
	sub	x3, x29, #160
	sub	x4, x29, #128
	sub	x5, x29, #160
	bl	_v15q_add
	cbnz	w0, LBB9_30
; %bb.7:                                ;   in Loop: Header=BB9_3 Depth=2
	ldp	x8, x9, [x29, #-192]
	stp	x8, x9, [x29, #-128]
	ldp	x8, x9, [x29, #-176]
	stp	x8, x9, [x29, #-112]
	ldp	x8, x9, [sp, #208]
	stp	x8, x9, [x29, #-160]
	ldp	x8, x9, [sp, #224]
	stp	x8, x9, [x29, #-144]
	lsr	x8, x24, x20
	tbnz	w8, #0, LBB9_10
; %bb.8:                                ;   in Loop: Header=BB9_3 Depth=2
	cbnz	x20, LBB9_2
	b	LBB9_12
LBB9_9:                                 ;   in Loop: Header=BB9_1 Depth=1
	mov	x22, x27
	mov	x21, x26
	b	LBB9_11
LBB9_10:                                ;   in Loop: Header=BB9_1 Depth=1
	sub	x22, x29, #192
	add	x21, sp, #208
	sub	x0, x29, #192
	add	x1, sp, #208
	sub	x2, x29, #128
	sub	x3, x29, #160
	mov	x4, x27
	mov	x5, x26
	bl	_v15q_add
	cbnz	w0, LBB9_30
LBB9_11:                                ;   in Loop: Header=BB9_1 Depth=1
	mov	w25, #0                         ; =0x0
	ldp	x8, x9, [x22]
	stp	x8, x9, [x29, #-128]
	ldp	x8, x9, [x22, #16]
	stp	x8, x9, [x29, #-112]
	ldp	x8, x9, [x21]
	stp	x8, x9, [x29, #-160]
	ldp	x9, x10, [x21, #16]
	sub	x8, x20, #1
	stp	x9, x10, [x29, #-144]
	cbnz	x20, LBB9_1
LBB9_12:
	ldp	x8, x9, [x29, #-128]
	stp	x8, x9, [sp, #176]
	ldp	x8, x9, [x29, #-160]
	stp	x8, x9, [sp, #144]
	ldp	x9, x8, [x29, #-144]
	str	x9, [sp, #160]
	ldp	x9, x10, [x29, #-112]
	stp	x9, x10, [sp, #192]
	b	LBB9_14
LBB9_13:
	mov	x8, #0                          ; =0x0
	stp	xzr, xzr, [sp, #176]
	stp	xzr, xzr, [sp, #144]
	str	xzr, [sp, #160]
	stp	xzr, xzr, [sp, #192]
LBB9_14:
	str	x8, [sp, #168]
	mov	w21, #1                         ; =0x1
	mov	w8, #15                         ; =0xf
LBB9_15:                                ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB9_17 Depth 2
	mov	x20, x8
	b	LBB9_17
LBB9_16:                                ;   in Loop: Header=BB9_17 Depth=2
	sub	x20, x20, #1
LBB9_17:                                ;   Parent Loop BB9_15 Depth=1
                                        ; =>  This Inner Loop Header: Depth=2
	tbz	w21, #0, LBB9_20
; %bb.18:                               ;   in Loop: Header=BB9_17 Depth=2
	lsr	x8, x23, x20
	tbnz	w8, #0, LBB9_23
; %bb.19:                               ;   in Loop: Header=BB9_17 Depth=2
	cbnz	x20, LBB9_16
	b	LBB9_27
LBB9_20:                                ;   in Loop: Header=BB9_17 Depth=2
	sub	x0, x29, #192
	add	x1, sp, #208
	sub	x2, x29, #128
	sub	x3, x29, #160
	sub	x4, x29, #128
	sub	x5, x29, #160
	bl	_v15q_add
	cbnz	w0, LBB9_30
; %bb.21:                               ;   in Loop: Header=BB9_17 Depth=2
	ldp	x8, x9, [x29, #-192]
	stp	x8, x9, [x29, #-128]
	ldp	x8, x9, [x29, #-176]
	stp	x8, x9, [x29, #-112]
	ldp	x8, x9, [sp, #208]
	stp	x8, x9, [x29, #-160]
	ldp	x8, x9, [sp, #224]
	stp	x8, x9, [x29, #-144]
	lsr	x8, x23, x20
	tbnz	w8, #0, LBB9_24
; %bb.22:                               ;   in Loop: Header=BB9_17 Depth=2
	cbnz	x20, LBB9_16
	b	LBB9_26
LBB9_23:                                ;   in Loop: Header=BB9_15 Depth=1
	mov	x25, x19
	mov	x22, x28
	b	LBB9_25
LBB9_24:                                ;   in Loop: Header=BB9_15 Depth=1
	sub	x25, x29, #192
	add	x22, sp, #208
	sub	x0, x29, #192
	add	x1, sp, #208
	sub	x2, x29, #128
	sub	x3, x29, #160
	mov	x4, x19
	mov	x5, x28
	bl	_v15q_add
	cbnz	w0, LBB9_30
LBB9_25:                                ;   in Loop: Header=BB9_15 Depth=1
	mov	w21, #0                         ; =0x0
	ldp	x8, x9, [x25]
	stp	x8, x9, [x29, #-128]
	ldp	x8, x9, [x25, #16]
	stp	x8, x9, [x29, #-112]
	ldp	x8, x9, [x22]
	stp	x8, x9, [x29, #-160]
	ldp	x9, x10, [x22, #16]
	sub	x8, x20, #1
	stp	x9, x10, [x29, #-144]
	cbnz	x20, LBB9_15
LBB9_26:
	ldp	x8, x9, [x29, #-128]
	stp	x8, x9, [sp, #112]
	ldp	x8, x9, [x29, #-160]
	stp	x8, x9, [sp, #80]
	ldp	x9, x8, [x29, #-144]
	str	x9, [sp, #96]
	ldp	x9, x10, [x29, #-112]
	stp	x9, x10, [sp, #128]
	b	LBB9_28
LBB9_27:
	mov	x8, #0                          ; =0x0
	stp	xzr, xzr, [sp, #112]
	stp	xzr, xzr, [sp, #80]
	str	xzr, [sp, #96]
	stp	xzr, xzr, [sp, #128]
LBB9_28:
	str	x8, [sp, #104]
	sub	x0, x29, #128
	sub	x1, x29, #160
	add	x4, sp, #176
	add	x5, sp, #144
	ldp	x20, x19, [sp, #64]             ; 16-byte Folded Reload
	mov	x2, x20
	mov	x3, x19
	bl	_v15q_add
	cbnz	w0, LBB9_30
; %bb.29:
	sub	x0, x29, #192
	add	x1, sp, #208
	sub	x2, x29, #128
	sub	x3, x29, #160
	add	x4, sp, #112
	add	x5, sp, #80
	bl	_v15q_add
	cbz	w0, LBB9_33
LBB9_30:
	mov	w0, #-1                         ; =0xffffffff
LBB9_31:
	ldur	x8, [x29, #-96]
Lloh40:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh41:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh42:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB9_40
; %bb.32:
	ldp	x29, x30, [sp, #432]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #416]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #400]            ; 16-byte Folded Reload
	ldp	x24, x23, [sp, #384]            ; 16-byte Folded Reload
	ldp	x26, x25, [sp, #368]            ; 16-byte Folded Reload
	ldp	x28, x27, [sp, #352]            ; 16-byte Folded Reload
	add	sp, sp, #448
	ret
LBB9_33:
	ldp	x8, x9, [x29, #-192]
	stp	x8, x9, [x29, #-128]
	ldp	x9, x10, [x29, #-176]
	stp	x9, x10, [x29, #-112]
	ldp	x9, x10, [sp, #208]
	stp	x9, x10, [x29, #-160]
	ldp	x10, x11, [sp, #224]
	stp	x10, x11, [x29, #-144]
	str	x8, [x20]
	str	x9, [x19]
	ldp	x8, x9, [x29, #-120]
	str	x8, [x20, #8]
	ldur	x8, [x29, #-152]
	str	x8, [x19, #8]
	str	x9, [x20, #16]
	str	x10, [x19, #16]
	ldur	x8, [x29, #-104]
	str	x8, [x20, #24]
	mov	x8, x11
	str	x8, [x19, #24]
	ldr	x11, [sp, #56]                  ; 8-byte Folded Reload
	ldp	x8, x10, [x11]
	adds	x9, x24, x8
	adcs	x10, x10, xzr
	ldp	x8, x11, [x11, #16]
	adcs	x12, x8, xzr
	adcs	x13, x11, xzr
	cset	w0, hs
Lloh43:
	adrp	x8, _R_CONST@GOTPAGE
Lloh44:
	ldr	x8, [x8, _R_CONST@GOTPAGEOFF]
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
	tbnz	w0, #0, LBB9_36
; %bb.34:
	and	x17, x17, #0x1
	tbz	w17, #0, LBB9_36
; %bb.35:
	mov	x11, x9
	mov	x14, x10
	mov	x15, x12
	mov	x16, x13
LBB9_36:
	ldr	x5, [sp, #40]                   ; 8-byte Folded Reload
	ldr	x9, [sp, #56]                   ; 8-byte Folded Reload
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
	tbnz	w17, #0, LBB9_39
; %bb.37:
	and	x16, x16, #0x1
	tbz	w16, #0, LBB9_39
; %bb.38:
	mov	x10, x9
	mov	x12, x11
	mov	x8, x13
	mov	x15, x14
LBB9_39:
	ldr	x4, [sp, #48]                   ; 8-byte Folded Reload
	ldp	x2, x3, [sp, #24]               ; 16-byte Folded Reload
	ldp	x0, x1, [sp, #8]                ; 16-byte Folded Reload
	stp	x10, x12, [x5]
	stp	x8, x15, [x5, #16]
	bl	_v15r_reset_hist
	mov	w0, #0                          ; =0x0
	b	LBB9_31
LBB9_40:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh37, Lloh38, Lloh39
	.loh AdrpLdrGotLdr	Lloh40, Lloh41, Lloh42
	.loh AdrpLdrGot	Lloh43, Lloh44
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
