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
	add	x8, x8, x8, lsl #1
	add	x8, x8, x0
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
	sub	sp, sp, #96
	stp	x29, x30, [sp, #80]             ; 16-byte Folded Spill
	add	x29, sp, #80
	.cfi_def_cfa w29, 16
	.cfi_offset w30, -8
	.cfi_offset w29, -16
Lloh6:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh7:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh8:
	ldr	x8, [x8]
	stur	x8, [x29, #-8]
	mov	x8, sp
	movi.2d	v0, #0000000000000000
	stp	q0, q0, [sp]
Lloh9:
	adrp	x12, _BETA@GOTPAGE
Lloh10:
	ldr	x12, [x12, _BETA@GOTPAGEOFF]
	stp	q0, q0, [sp, #32]
	ldp	x9, x10, [x12]
	ldp	x11, x12, [x12, #16]
	ldr	x13, [x1]
	umulh	x14, x9, x13
	mul	x15, x9, x13
	ldp	x16, x17, [sp]
	adds	x15, x15, x16
	cinc	x14, x14, hs
	mul	x16, x10, x13
	adds	x14, x14, x17
	cset	w17, hs
	adds	x14, x14, x16
	umulh	x16, x10, x13
	adc	x16, x17, x16
	stp	x15, x14, [sp]
	mul	x14, x11, x13
	ldp	x15, x17, [sp, #16]
	adds	x15, x16, x15
	cset	w16, hs
	adds	x14, x15, x14
	umulh	x15, x11, x13
	adc	x15, x16, x15
	umulh	x16, x12, x13
	mul	x13, x12, x13
	adds	x15, x15, x17
	cset	w17, hs
	adds	x15, x15, x13
	adc	x13, x17, x16
	stp	x14, x15, [sp, #16]
	cbz	x13, LBB1_5
; %bb.1:
	ldr	x14, [sp, #32]
	adds	x13, x14, x13
	str	x13, [sp, #32]
	b.lo	LBB1_5
; %bb.2:
	ldr	x13, [sp, #40]
	adds	x13, x13, #1
	str	x13, [sp, #40]
	b.lo	LBB1_5
; %bb.3:
	ldr	x13, [sp, #48]
	adds	x13, x13, #1
	str	x13, [sp, #48]
	b.lo	LBB1_5
; %bb.4:
	ldr	x13, [sp, #56]
	add	x13, x13, #1
	str	x13, [sp, #56]
LBB1_5:
	mov	x13, #0                         ; =0x0
	add	x14, x8, #32
	add	x15, x1, #8
	b	LBB1_7
LBB1_6:                                 ;   in Loop: Header=BB1_7 Depth=1
	add	x13, x13, #1
	add	x14, x14, #8
	cmp	x13, #3
	b.eq	LBB1_13
LBB1_7:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x16, [x15, x13, lsl #3]
	umulh	x17, x9, x16
	mul	x1, x9, x16
	ldp	x2, x3, [x14, #-24]
	adds	x1, x1, x2
	cinc	x17, x17, hs
	mul	x2, x10, x16
	adds	x17, x17, x3
	cset	w3, hs
	adds	x17, x17, x2
	umulh	x2, x10, x16
	adc	x2, x3, x2
	stp	x1, x17, [x14, #-24]
	mul	x17, x11, x16
	ldp	x1, x3, [x14, #-8]
	adds	x1, x2, x1
	cset	w2, hs
	adds	x17, x1, x17
	umulh	x1, x11, x16
	adc	x1, x2, x1
	umulh	x2, x12, x16
	mul	x16, x12, x16
	adds	x1, x1, x3
	cset	w3, hs
	adds	x16, x1, x16
	adc	x1, x3, x2
	stp	x17, x16, [x14, #-8]
	cbz	x1, LBB1_6
; %bb.8:                                ;   in Loop: Header=BB1_7 Depth=1
	add	x16, x13, #5
	ldr	x17, [x14, #8]
	adds	x1, x17, x1
	cset	w2, hs
	str	x1, [x14, #8]
	cmp	x16, #6
	b.hi	LBB1_6
; %bb.9:                                ;   in Loop: Header=BB1_7 Depth=1
	cbz	w2, LBB1_6
; %bb.10:                               ;   in Loop: Header=BB1_7 Depth=1
	cmp	x1, x17
	cset	w1, lo
	ldr	x17, [x14, #16]
	adds	x1, x17, x1
	cset	w2, hs
	str	x1, [x14, #16]
	cmp	x16, #5
	b.hi	LBB1_6
; %bb.11:                               ;   in Loop: Header=BB1_7 Depth=1
	cbz	w2, LBB1_6
; %bb.12:                               ;   in Loop: Header=BB1_7 Depth=1
	cmp	x1, x17
	ldr	x16, [x14, #24]
	cinc	x16, x16, lo
	str	x16, [x14, #24]
	b	LBB1_6
LBB1_13:
	mov	x9, #-1                         ; =0xffffffffffffffff
	mov	x10, #977                       ; =0x3d1
	movk	x10, #1, lsl #32
	ldr	x11, [sp, #56]
	umulh	x12, x11, x10
	mul	x11, x11, x10
	ldp	x13, x14, [sp, #24]
	adds	x11, x11, x13
	adcs	x12, x14, x12
	ldp	x13, x14, [sp, #40]
	cinc	x13, x13, hs
	umulh	x15, x14, x10
	mul	x14, x14, x10
	ldp	x17, x16, [sp, #8]
	adds	x14, x14, x16
	adcs	x11, x11, x15
	cinc	x12, x12, hs
	umulh	x15, x13, x10
	mul	x13, x13, x10
	adds	x13, x13, x17
	adcs	x14, x14, x15
	cinc	x16, x11, hs
	umulh	x11, x12, x10
	mul	x12, x12, x10
	ldr	x8, [x8]
	adds	x17, x12, x8
	adcs	x1, x13, x11
	cinc	x14, x14, hs
	adds	x2, x17, x10
	adc	x8, xzr, x9
	asr	x11, x8, #63
	adds	x8, x8, x1
	cinc	x11, x11, hs
	adds	x3, x8, #1
	adc	x8, x11, x9
	asr	x11, x8, #63
	adds	x8, x8, x14
	cinc	x11, x11, hs
	adds	x4, x8, #1
	adc	x8, x11, x9
	asr	x11, x8, #63
	adds	x8, x8, x16
	cinc	x11, x11, hs
	adds	x8, x8, #1
	adc	x5, x11, x9
	adds	x13, x2, x10
	cinc	x10, x9, hs
	asr	x11, x10, #63
	adds	x10, x10, x3
	cinc	x11, x11, hs
	adds	x12, x10, #1
	adc	x10, x11, x9
	asr	x11, x10, #63
	adds	x10, x10, x4
	cinc	x11, x11, hs
	adds	x6, x10, #1
	adc	x10, x11, x9
	asr	x11, x10, #63
	adds	x10, x10, x8
	cinc	x11, x11, hs
	adds	x10, x10, #1
	adc	x9, x11, x9
	tst	w5, #0x1
	csel	x15, x9, x5, eq
	csel	x10, x10, x8, eq
	csel	x11, x6, x4, eq
	csel	x12, x12, x3, eq
	csel	x13, x13, x2, eq
	tst	w5, #0x1
	csel	x8, x8, x16, eq
	csel	x9, x4, x14, eq
	csel	x14, x3, x1, eq
	csel	x16, x2, x17, eq
	tbnz	w15, #0, LBB1_15
; %bb.14:
	mov	x8, x10
	mov	x9, x11
	mov	x14, x12
	mov	x16, x13
LBB1_15:
	stp	x16, x14, [x0]
	stp	x9, x8, [x0, #16]
	ldur	x8, [x29, #-8]
Lloh11:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh12:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh13:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB1_17
; %bb.16:
	ldp	x29, x30, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #96
	ret
LBB1_17:
	bl	___stack_chk_fail
	.loh AdrpLdrGot	Lloh9, Lloh10
	.loh AdrpLdrGotLdr	Lloh6, Lloh7, Lloh8
	.loh AdrpLdrGotLdr	Lloh11, Lloh12, Lloh13
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
Lloh14:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh15:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh16:
	ldr	x8, [x8]
	stur	x8, [x29, #-24]
	cmp	x2, #1
	b.eq	LBB2_3
; %bb.1:
	cbnz	x2, LBB2_4
; %bb.2:
Lloh17:
	adrp	x9, _R_MONT_ONE@GOTPAGE
Lloh18:
	ldr	x9, [x9, _R_MONT_ONE@GOTPAGEOFF]
	b	LBB2_5
LBB2_3:
Lloh19:
	adrp	x9, _LAMBDA2_R@GOTPAGE
Lloh20:
	ldr	x9, [x9, _LAMBDA2_R@GOTPAGEOFF]
	b	LBB2_5
LBB2_4:
Lloh21:
	adrp	x9, _LAMBDA_R@GOTPAGE
Lloh22:
	ldr	x9, [x9, _LAMBDA_R@GOTPAGEOFF]
LBB2_5:
	ldp	x13, x16, [x9]
	stp	x13, x16, [sp, #8]
	add	x8, x9, #24
	ldr	x2, [x9, #16]!
	ldr	x17, [x8]
	stp	x2, x17, [sp, #24]
Lloh23:
	adrp	x8, _R_CONST@GOTPAGE
Lloh24:
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
Lloh25:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh26:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh27:
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
	.loh AdrpLdrGotLdr	Lloh14, Lloh15, Lloh16
	.loh AdrpLdrGot	Lloh17, Lloh18
	.loh AdrpLdrGot	Lloh19, Lloh20
	.loh AdrpLdrGot	Lloh21, Lloh22
	.loh AdrpLdrGot	Lloh23, Lloh24
	.loh AdrpLdrGotLdr	Lloh25, Lloh26, Lloh27
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
	mov	x16, #0                         ; =0x0
	mov	x10, #0                         ; =0x0
	mov	x12, #0                         ; =0x0
	mov	x15, #0                         ; =0x0
	mov	x14, #0                         ; =0x0
	mov	x6, #0                          ; =0x0
Lloh28:
	adrp	x8, _R_N0@GOTPAGE
Lloh29:
	ldr	x8, [x8, _R_N0@GOTPAGEOFF]
Lloh30:
	ldr	x17, [x8]
Lloh31:
	adrp	x8, _R_CONST@GOTPAGE
Lloh32:
	ldr	x8, [x8, _R_CONST@GOTPAGEOFF]
	ldp	x3, x4, [x1]
	ldp	x5, x1, [x1, #16]
	ldp	x13, x11, [x8]
	ldp	x9, x8, [x8, #16]
LBB3_1:                                 ; =>This Inner Loop Header: Depth=1
	ldr	x7, [x2, x16]
	umulh	x19, x7, x3
	mul	x20, x7, x3
	adds	x10, x20, x10
	cinc	x19, x19, hs
	umulh	x20, x7, x4
	mul	x21, x7, x4
	adds	x12, x19, x12
	cset	w19, hs
	adds	x12, x12, x21
	adc	x19, x19, x20
	umulh	x20, x7, x5
	mul	x21, x7, x5
	adds	x15, x19, x15
	cset	w19, hs
	adds	x15, x15, x21
	adc	x19, x19, x20
	umulh	x20, x7, x1
	mul	x7, x7, x1
	adds	x14, x19, x14
	cset	w19, hs
	adds	x14, x14, x7
	adc	x7, x19, x20
	adds	x6, x6, x7
	cset	w7, hs
	mul	x19, x17, x10
	umulh	x20, x19, x13
	mul	x21, x19, x13
	cmn	x21, x10
	cinc	x10, x20, hs
	umulh	x20, x11, x19
	mul	x21, x11, x19
	adds	x10, x12, x10
	cset	w12, hs
	adds	x10, x10, x21
	adc	x12, x12, x20
	umulh	x20, x9, x19
	mul	x21, x9, x19
	adds	x12, x15, x12
	cset	w15, hs
	adds	x12, x12, x21
	adc	x15, x15, x20
	umulh	x20, x8, x19
	mul	x19, x8, x19
	adds	x14, x14, x15
	cset	w21, hs
	adds	x15, x14, x19
	adc	x14, x21, x20
	adds	x14, x6, x14
	cinc	x6, x7, hs
	add	x16, x16, #8
	cmp	x16, #32
	b.ne	LBB3_1
; %bb.2:
	subs	x16, x10, x13
	ngc	x17, xzr
	and	x17, x17, #0x1
	adds	x17, x17, x11
	cset	w1, hs
	subs	x17, x12, x17
	ngc	x1, x1
	and	x1, x1, #0x1
	adds	x1, x1, x9
	cset	w2, hs
	subs	x1, x15, x1
	ngc	x2, x2
	and	x2, x2, #0x1
	adds	x2, x2, x8
	cset	w3, hs
	subs	x2, x14, x2
	ngc	x3, x3
	subs	x13, x16, x13
	ngc	x4, xzr
	and	x4, x4, #0x1
	adds	x11, x4, x11
	cset	w4, hs
	subs	x11, x17, x11
	ngc	x4, x4
	and	x4, x4, #0x1
	adds	x9, x4, x9
	cset	w4, hs
	subs	x9, x1, x9
	ngc	x4, x4
	and	x4, x4, #0x1
	adds	x8, x4, x8
	cset	w4, hs
	subs	x8, x2, x8
	ngc	x4, x4
	tst	w3, #0x1
	csel	x4, x3, x4, ne
	csel	x8, x2, x8, ne
	csel	x9, x1, x9, ne
	csel	x11, x17, x11, ne
	csel	x13, x16, x13, ne
	tst	w3, #0x1
	csel	x14, x14, x2, ne
	csel	x15, x15, x1, ne
	csel	x12, x12, x17, ne
	csel	x10, x10, x16, ne
	tst	w4, #0x1
	csel	x8, x14, x8, ne
	csel	x11, x12, x11, ne
	csel	x10, x10, x13, ne
	stp	x10, x11, [x0]
	csel	x9, x15, x9, ne
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.loh AdrpLdrGot	Lloh31, Lloh32
	.loh AdrpLdrGotLdr	Lloh28, Lloh29, Lloh30
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_detect                    ; -- Begin function v15o_detect
	.p2align	2
_v15o_detect:                           ; @v15o_detect
	.cfi_startproc
; %bb.0:
	ldp	q0, q1, [x0]
	ldr	x8, [x0]
	ldp	q2, q3, [x1]
	cmeq.2d	v1, v1, v3
	cmeq.2d	v0, v0, v2
	uzp1.4s	v0, v0, v1
	mvn.16b	v0, v0
	umaxv.4s	s0, v0
	fmov	w9, s0
	tbz	w9, #0, LBB4_2
; %bb.1:
	ldp	q0, q1, [x0]
	ldp	q2, q3, [x2]
	cmeq.2d	v1, v1, v3
	cmeq.2d	v0, v0, v2
	uzp1.4s	v0, v0, v1
	mvn.16b	v0, v0
	umaxv.4s	s0, v0
	fmov	w9, s0
	tbnz	w9, #0, LBB4_3
LBB4_2:
	ldr	x9, [x3]
	add	x9, x9, #1
	str	x9, [x3]
LBB4_3:
	ldr	x9, [x0, #24]
	ldr	x10, [x4, #24]
	cmp	x9, x10
	b.hs	LBB4_5
LBB4_4:
	ldr	x9, [x0]
	str	x9, [x4]
	ldr	x9, [x0, #8]
	str	x9, [x4, #8]
	ldr	x9, [x0, #16]
	str	x9, [x4, #16]
	ldr	x9, [x0, #24]
	str	x9, [x4, #24]
	str	xzr, [x5]
	b	LBB4_6
LBB4_5:
	b.ls	LBB4_7
LBB4_6:
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
LBB4_7:
	ldr	x9, [x0, #16]
	ldr	x10, [x4, #16]
	cmp	x9, x10
	b.lo	LBB4_4
; %bb.8:
	b.hi	LBB4_6
; %bb.9:
	ldr	x9, [x0, #8]
	ldr	x10, [x4, #8]
	cmp	x9, x10
	b.lo	LBB4_4
; %bb.10:
	b.hi	LBB4_6
; %bb.11:
	ldr	x9, [x0]
	ldr	x10, [x4]
	cmp	x9, x10
	b.lo	LBB4_4
	b	LBB4_6
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
.subsections_via_symbols
