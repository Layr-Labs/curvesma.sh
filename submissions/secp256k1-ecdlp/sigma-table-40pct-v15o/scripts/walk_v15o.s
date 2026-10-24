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
	b.hs	LBB0_2
LBB0_1:
	ldp	x8, x9, [sp, #40]
	mov	w0, #1                          ; =0x1
	mov	x10, x12
	ldr	x11, [sp, #56]
	b	LBB0_4
LBB0_2:
	ldp	x8, x9, [x19]
	ldr	x11, [x19, #16]
	b.ls	LBB0_11
LBB0_3:
	mov	x0, #0                          ; =0x0
LBB0_4:
	ldr	x12, [sp, #32]
	cmp	x12, x10
	b.hs	LBB0_6
LBB0_5:
	ldp	x8, x9, [sp, #8]
	mov	w0, #2                          ; =0x2
	mov	x10, x12
	ldr	x11, [sp, #24]
	b	LBB0_7
LBB0_6:
	b.ls	LBB0_16
LBB0_7:
	stp	x8, x9, [x19]
	stp	x11, x10, [x19, #16]
	and	x8, x8, #0x7ff
	add	x8, x8, x8, lsl #1
	add	x8, x8, x0
	str	x8, [x19, #64]
	ldr	x8, [x19, #32]
	and	w9, w8, #0x1
	str	w9, [x20]
	tbz	w8, #0, LBB0_9
; %bb.8:
	mov	x9, #-977                       ; =0xfffffffffffffc2f
	movk	x9, #65534, lsl #32
	subs	x8, x9, x8
	ngc	x9, xzr
	ldp	x10, x11, [x19, #40]
	mvn	x10, x10
	adds	x10, x9, x10
	cinc	x9, x9, hs
	mvn	x11, x11
	asr	x12, x9, #63
	adds	x9, x9, x11
	ldr	x11, [x19, #56]
	mvn	x11, x11
	stp	x8, x10, [x19, #32]
	adc	x8, x11, x12
	stp	x9, x8, [x19, #48]
LBB0_9:
	ldur	x8, [x29, #-24]
Lloh3:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh4:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh5:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB0_21
; %bb.10:
	ldp	x29, x30, [sp, #96]             ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #80]             ; 16-byte Folded Reload
	add	sp, sp, #112
	ret
LBB0_11:
	ldr	x13, [sp, #56]
	cmp	x13, x11
	b.lo	LBB0_1
; %bb.12:
	b.hi	LBB0_3
; %bb.13:
	ldr	x13, [sp, #48]
	cmp	x13, x9
	b.lo	LBB0_1
; %bb.14:
	mov	x0, #0                          ; =0x0
	cmp	x13, x9
	b.hi	LBB0_4
; %bb.15:
	ldr	x13, [sp, #40]
	cmp	x13, x8
	b.lo	LBB0_1
	b	LBB0_4
LBB0_16:
	ldr	x13, [sp, #24]
	cmp	x13, x11
	b.lo	LBB0_5
; %bb.17:
	b.hi	LBB0_7
; %bb.18:
	ldr	x13, [sp, #16]
	cmp	x13, x9
	b.lo	LBB0_5
; %bb.19:
	b.hi	LBB0_7
; %bb.20:
	ldr	x13, [sp, #8]
	cmp	x13, x8
	b.lo	LBB0_5
	b	LBB0_7
LBB0_21:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh0, Lloh1, Lloh2
	.loh AdrpLdrGotLdr	Lloh3, Lloh4, Lloh5
	.cfi_endproc
                                        ; -- End function
	.p2align	2                               ; -- Begin function fmul
_fmul:                                  ; @fmul
	.cfi_startproc
; %bb.0:
Lloh6:
	adrp	x9, _BETA@GOTPAGE
Lloh7:
	ldr	x9, [x9, _BETA@GOTPAGEOFF]
	ldp	x13, x10, [x9]
	ldp	x17, x14, [x1]
	mul	x11, x13, x17
	umulh	x8, x13, x17
	umulh	x12, x10, x17
	mul	x15, x10, x17
	adds	x8, x8, x15
	cinc	x15, x12, hs
	umulh	x16, x13, x14
	mul	x12, x13, x14
	adds	x12, x8, x12
	cinc	x16, x16, hs
	ldp	x8, x9, [x9, #16]
	umulh	x2, x8, x17
	mul	x3, x8, x17
	adds	x15, x15, x3
	cinc	x2, x2, hs
	adds	x15, x16, x15
	cset	w16, hs
	umulh	x3, x10, x14
	mul	x4, x10, x14
	adds	x4, x15, x4
	adc	x3, x16, x3
	ldp	x16, x15, [x1, #16]
	umulh	x5, x13, x16
	mul	x1, x13, x16
	adds	x1, x4, x1
	cinc	x4, x5, hs
	umulh	x5, x9, x17
	mul	x17, x9, x17
	adds	x17, x2, x17
	cinc	x2, x5, hs
	adds	x17, x3, x17
	cset	w3, hs
	umulh	x5, x8, x14
	mul	x6, x8, x14
	adds	x17, x17, x6
	adc	x3, x3, x5
	adds	x17, x4, x17
	cset	w4, hs
	umulh	x5, x10, x16
	mul	x6, x10, x16
	adds	x17, x17, x6
	adc	x4, x4, x5
	umulh	x5, x13, x15
	mul	x13, x13, x15
	adds	x13, x17, x13
	cinc	x17, x5, hs
	adds	x2, x3, x2
	cset	w3, hs
	umulh	x5, x9, x14
	mul	x14, x9, x14
	adds	x14, x2, x14
	adc	x2, x3, x5
	adds	x14, x4, x14
	cset	w3, hs
	umulh	x4, x8, x16
	mul	x5, x8, x16
	adds	x14, x14, x5
	adc	x3, x3, x4
	adds	x14, x17, x14
	cset	w17, hs
	umulh	x4, x10, x15
	mul	x10, x10, x15
	adds	x14, x14, x10
	adc	x17, x17, x4
	mov	x10, #977                       ; =0x3d1
	movk	x10, #1, lsl #32
	umulh	x4, x14, x10
	mul	x14, x14, x10
	adds	x11, x11, x14
	adcs	x12, x12, x4
	cset	w14, hs
	adds	x2, x3, x2
	cset	w3, hs
	umulh	x4, x9, x16
	mul	x16, x9, x16
	adds	x16, x2, x16
	adc	x2, x3, x4
	adds	x16, x17, x16
	cset	w17, hs
	umulh	x3, x8, x15
	mul	x8, x8, x15
	adds	x8, x16, x8
	adc	x16, x17, x3
	umulh	x17, x8, x10
	mul	x8, x8, x10
	adds	x8, x12, x8
	adc	x12, x14, x17
	adds	x12, x12, x1
	cset	w14, hs
	adds	x16, x16, x2
	cset	w17, hs
	umulh	x1, x9, x15
	mul	x9, x9, x15
	adds	x9, x16, x9
	adc	x15, x17, x1
	umulh	x16, x9, x10
	mul	x9, x9, x10
	adds	x12, x12, x9
	adc	x9, x14, x16
	adds	x13, x9, x13
	cset	w9, hs
	umulh	x14, x15, x10
	mul	x15, x15, x10
	adds	x16, x13, x15
	adc	x9, x9, x14
	mul	x9, x9, x10
	umulh	x14, x9, x10
	mul	x9, x9, x10
	adds	x9, x11, x9
	cinc	x11, x14, hs
	adcs	xzr, x8, x14
	adcs	xzr, x12, xzr
	adcs	xzr, x16, xzr
	cset	w14, hs
	mul	x14, x14, x10
	adds	x9, x14, x9
	adcs	x11, x11, x8
	adcs	x12, x12, xzr
	adc	x8, x13, x15
	adds	x10, x9, x10
	mov	x15, #-1                        ; =0xffffffffffffffff
	adc	x13, xzr, x15
	adds	x14, x11, x13
	asr	x13, x13, #63
	cinc	x16, x13, hs
	adds	x13, x14, #1
	adc	x14, x16, x15
	adds	x16, x12, x14
	asr	x14, x14, #63
	cinc	x17, x14, hs
	adds	x14, x16, #1
	adc	x15, x17, x15
	adds	x16, x8, x15
	asr	x15, x15, #63
	cinc	w15, w15, hs
	cmn	x16, #1
	mov	w16, #1                         ; =0x1
	adc	w15, w15, w16
	tbnz	w15, #0, LBB1_2
; %bb.1:
	mov	x8, #0                          ; =0x0
	mov	x9, x10
	mov	x11, x13
	mov	x12, x14
LBB1_2:
	stp	x9, x11, [x0]
	stp	x12, x8, [x0, #16]
	ret
	.loh AdrpLdrGot	Lloh6, Lloh7
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_label                     ; -- Begin function v15o_label
	.p2align	2
_v15o_label:                            ; @v15o_label
	.cfi_startproc
; %bb.0:
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
	mov	x21, x3
	mov	x19, x1
	mov	x20, x0
Lloh8:
	adrp	x8, ___stack_chk_guard@GOTPAGE
Lloh9:
	ldr	x8, [x8, ___stack_chk_guard@GOTPAGEOFF]
Lloh10:
	ldr	x8, [x8]
	stur	x8, [x29, #-40]
	cmp	x2, #1
	b.eq	LBB2_3
; %bb.1:
	cbnz	x2, LBB2_4
; %bb.2:
	mov	w8, #1                          ; =0x1
	stp	x8, xzr, [sp, #64]
	stp	xzr, xzr, [sp, #80]
	b	LBB2_5
LBB2_3:
Lloh11:
	adrp	x8, _LAMBDA_R@GOTPAGE
Lloh12:
	ldr	x8, [x8, _LAMBDA_R@GOTPAGEOFF]
	ldp	q0, q1, [x8]
	stp	q0, q1, [sp, #64]
	b	LBB2_5
LBB2_4:
Lloh13:
	adrp	x8, _LAMBDA_R@GOTPAGE
Lloh14:
	ldr	x8, [x8, _LAMBDA_R@GOTPAGEOFF]
	ldp	x9, x10, [x8]
	stp	x9, x10, [sp, #32]
	ldp	x11, x8, [x8, #16]
	stp	x11, x8, [sp, #48]
	stp	x9, x10, [sp]
	stp	x11, x8, [sp, #16]
	add	x0, sp, #64
	add	x1, sp, #32
	mov	x2, sp
	bl	_rmont
LBB2_5:
	add	x2, sp, #64
	mov	x0, x20
	mov	x1, x20
	bl	_rmont
	add	x2, sp, #64
	mov	x0, x19
	mov	x1, x19
	bl	_rmont
	cbz	w21, LBB2_11
; %bb.6:
Lloh15:
	adrp	x16, _R_CONST@GOTPAGE
Lloh16:
	ldr	x16, [x16, _R_CONST@GOTPAGEOFF]
	ldp	x9, x8, [x16]
	ldp	x10, x11, [x20]
	subs	x0, x9, x10
	ngc	x10, xzr
	and	x10, x10, #0x1
	adds	x10, x10, x11
	cset	w11, hs
	subs	x1, x8, x10
	ngc	x10, x11
	and	x10, x10, #0x1
	ldp	x11, x12, [x20, #16]
	adds	x13, x10, x11
	cset	w14, hs
	ldp	x11, x10, [x16, #16]
	subs	x2, x11, x13
	ngc	x13, x14
	and	x13, x13, #0x1
	adds	x12, x13, x12
	cset	w13, hs
	subs	x3, x10, x12
	mov	x12, #0                         ; =0x0
	ngc	x14, x13
	tbnz	w14, #0, LBB2_8
; %bb.7:
	mov	x13, #0                         ; =0x0
	mov	x15, #0                         ; =0x0
	mov	x17, #0                         ; =0x0
	b	LBB2_9
LBB2_8:
	mov	x13, #0                         ; =0x0
	mov	x15, #0                         ; =0x0
	mov	x17, #0                         ; =0x0
	stp	x0, x1, [x20]
	stp	x2, x3, [x20, #16]
	ldp	x9, x8, [x16]
	ldp	x11, x10, [x16, #16]
LBB2_9:
	and	x14, x14, #0x1
	ldp	x16, x0, [x19]
	adds	x14, x14, x16
	cset	w16, hs
	subs	x9, x9, x14
	sbc	x14, x17, x16
	and	x14, x14, #0x1
	adds	x14, x14, x0
	cset	w16, hs
	subs	x8, x8, x14
	sbc	x14, x15, x16
	and	x14, x14, #0x1
	ldp	x15, x16, [x19, #16]
	adds	x14, x14, x15
	cset	w15, hs
	subs	x11, x11, x14
	sbc	x13, x13, x15
	and	x13, x13, #0x1
	adds	x13, x13, x16
	cset	w14, hs
	subs	x10, x10, x13
	sbc	x12, x12, x14
	tbz	w12, #0, LBB2_11
; %bb.10:
	stp	x9, x8, [x19]
	stp	x11, x10, [x19, #16]
LBB2_11:
	ldur	x8, [x29, #-40]
Lloh17:
	adrp	x9, ___stack_chk_guard@GOTPAGE
Lloh18:
	ldr	x9, [x9, ___stack_chk_guard@GOTPAGEOFF]
Lloh19:
	ldr	x9, [x9]
	cmp	x9, x8
	b.ne	LBB2_13
; %bb.12:
	ldp	x29, x30, [sp, #144]            ; 16-byte Folded Reload
	ldp	x20, x19, [sp, #128]            ; 16-byte Folded Reload
	ldp	x22, x21, [sp, #112]            ; 16-byte Folded Reload
	add	sp, sp, #160
	ret
LBB2_13:
	bl	___stack_chk_fail
	.loh AdrpLdrGotLdr	Lloh8, Lloh9, Lloh10
	.loh AdrpLdrGot	Lloh11, Lloh12
	.loh AdrpLdrGot	Lloh13, Lloh14
	.loh AdrpLdrGot	Lloh15, Lloh16
	.loh AdrpLdrGotLdr	Lloh17, Lloh18, Lloh19
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
Lloh20:
	adrp	x11, _R_CONST@GOTPAGE
Lloh21:
	ldr	x11, [x11, _R_CONST@GOTPAGEOFF]
	ldp	x9, x8, [x11, #16]
	ldp	x14, x12, [x2, #16]
	ldp	x16, x15, [x2]
	ldp	x11, x10, [x11]
Lloh22:
	adrp	x13, _R_N0@GOTPAGE
Lloh23:
	ldr	x13, [x13, _R_N0@GOTPAGEOFF]
Lloh24:
	ldr	x13, [x13]
	ldp	x17, x2, [x1]
	mul	x3, x16, x17
	umulh	x4, x16, x17
	umulh	x5, x15, x17
	mul	x6, x15, x17
	adds	x4, x4, x6
	cinc	x5, x5, hs
	umulh	x6, x14, x17
	mul	x7, x14, x17
	adds	x5, x5, x7
	cinc	x6, x6, hs
	mul	x7, x12, x17
	adds	x19, x6, x7
	mul	x20, x13, x3
	mul	x21, x20, x11
	cmn	x21, x3
	umulh	x3, x20, x11
	cinc	x3, x3, hs
	mul	x21, x10, x20
	adds	x3, x4, x3
	cset	w4, hs
	adds	x3, x3, x21
	umulh	x21, x10, x20
	adc	x4, x4, x21
	mul	x21, x9, x20
	adds	x4, x5, x4
	cset	w5, hs
	adds	x4, x4, x21
	umulh	x21, x9, x20
	adc	x5, x5, x21
	umulh	x21, x8, x20
	mul	x20, x8, x20
	adds	x5, x19, x5
	cset	w19, hs
	adds	x5, x5, x20
	adc	x19, x19, x21
	cmn	x6, x7
	umulh	x17, x12, x17
	adc	x17, x19, x17
	umulh	x6, x16, x2
	mul	x7, x16, x2
	adds	x3, x7, x3
	cinc	x6, x6, hs
	umulh	x7, x15, x2
	mul	x19, x15, x2
	adds	x4, x6, x4
	cset	w6, hs
	adds	x4, x4, x19
	adc	x6, x6, x7
	umulh	x7, x14, x2
	mul	x19, x14, x2
	adds	x5, x6, x5
	cset	w6, hs
	adds	x5, x5, x19
	adc	x6, x6, x7
	umulh	x7, x12, x2
	mul	x2, x12, x2
	adds	x17, x6, x17
	cset	w6, hs
	adds	x2, x17, x2
	adc	x17, x6, x7
	mul	x6, x13, x3
	umulh	x7, x6, x11
	mul	x19, x6, x11
	cmn	x19, x3
	cinc	x3, x7, hs
	umulh	x7, x10, x6
	mul	x19, x10, x6
	adds	x3, x4, x3
	cset	w4, hs
	adds	x3, x3, x19
	adc	x4, x4, x7
	umulh	x7, x9, x6
	mul	x19, x9, x6
	adds	x4, x5, x4
	cset	w20, hs
	adds	x5, x4, x19
	adc	x4, x20, x7
	umulh	x7, x8, x6
	mul	x6, x8, x6
	adds	x2, x2, x4
	cset	w4, hs
	adds	x6, x2, x6
	adc	x2, x4, x7
	ldp	x4, x1, [x1, #16]
	umulh	x7, x16, x4
	mul	x19, x16, x4
	adds	x3, x19, x3
	cinc	x7, x7, hs
	umulh	x19, x15, x4
	mul	x20, x15, x4
	adds	x5, x7, x5
	cset	w7, hs
	adds	x5, x5, x20
	adc	x7, x7, x19
	umulh	x19, x14, x4
	mul	x20, x14, x4
	adds	x6, x7, x6
	cset	w7, hs
	adds	x6, x6, x20
	adc	x7, x7, x19
	umulh	x19, x12, x4
	mul	x4, x12, x4
	add	x17, x17, x2
	adds	x17, x7, x17
	cset	w2, hs
	adds	x17, x17, x4
	adc	x2, x2, x19
	mul	x4, x13, x3
	umulh	x7, x4, x11
	mul	x19, x4, x11
	cmn	x19, x3
	cinc	x3, x7, hs
	umulh	x7, x10, x4
	mul	x19, x10, x4
	adds	x3, x5, x3
	cset	w5, hs
	adds	x3, x3, x19
	adc	x5, x5, x7
	umulh	x7, x9, x4
	mul	x19, x9, x4
	adds	x5, x6, x5
	cset	w6, hs
	adds	x5, x5, x19
	adc	x6, x6, x7
	umulh	x7, x8, x4
	mul	x4, x8, x4
	adds	x17, x17, x6
	cset	w6, hs
	adds	x17, x17, x4
	adc	x4, x6, x7
	umulh	x6, x16, x1
	mul	x16, x16, x1
	adds	x16, x16, x3
	cinc	x3, x6, hs
	umulh	x6, x15, x1
	mul	x15, x15, x1
	adds	x3, x3, x5
	cset	w5, hs
	adds	x15, x3, x15
	adc	x3, x5, x6
	umulh	x5, x14, x1
	mul	x14, x14, x1
	adds	x17, x3, x17
	cset	w3, hs
	adds	x14, x17, x14
	adc	x17, x3, x5
	umulh	x3, x12, x1
	mul	x12, x12, x1
	add	x1, x2, x4
	adds	x17, x17, x1
	cset	w1, hs
	adds	x12, x17, x12
	adc	x17, x1, x3
	mul	x13, x13, x16
	umulh	x1, x13, x11
	mul	x2, x13, x11
	cmn	x2, x16
	cinc	x16, x1, hs
	umulh	x1, x10, x13
	mul	x2, x10, x13
	adds	x15, x15, x16
	cset	w16, hs
	adds	x15, x15, x2
	adc	x16, x16, x1
	umulh	x1, x9, x13
	mul	x2, x9, x13
	adds	x14, x14, x16
	cset	w16, hs
	adds	x14, x14, x2
	adc	x16, x16, x1
	umulh	x1, x8, x13
	mul	x13, x8, x13
	adds	x12, x12, x16
	cset	w16, hs
	adds	x12, x12, x13
	adc	x13, x16, x1
	subs	x11, x15, x11
	ngc	x16, xzr
	and	x16, x16, #0x1
	adds	x10, x16, x10
	cset	w16, hs
	subs	x10, x14, x10
	ngc	x16, x16
	and	x16, x16, #0x1
	adds	x9, x16, x9
	cset	w16, hs
	subs	x9, x12, x9
	ngc	x16, x16
	and	x16, x16, #0x1
	adds	x8, x16, x8
	cset	w16, hs
	add	x13, x17, x13
	subs	x8, x13, x8
	ngc	x16, x16
	tst	w16, #0x1
	csel	x11, x15, x11, ne
	csel	x10, x14, x10, ne
	csel	x9, x12, x9, ne
	stp	x11, x10, [x0]
	csel	x8, x13, x8, ne
	stp	x9, x8, [x0, #16]
	ldp	x20, x19, [sp, #16]             ; 16-byte Folded Reload
	ldp	x22, x21, [sp], #32             ; 16-byte Folded Reload
	ret
	.loh AdrpLdrGotLdr	Lloh22, Lloh23, Lloh24
	.loh AdrpLdrGot	Lloh20, Lloh21
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_detect                    ; -- Begin function v15o_detect
	.p2align	2
_v15o_detect:                           ; @v15o_detect
	.cfi_startproc
; %bb.0:
	ldr	x8, [x0]
	ldr	x10, [x1]
	mov	x9, x8
	cmp	x8, x10
	b.ne	LBB4_5
; %bb.1:
	ldr	x10, [x0, #8]
	ldr	x11, [x1, #8]
	mov	x9, x8
	cmp	x10, x11
	b.ne	LBB4_5
; %bb.2:
	ldr	x10, [x0, #16]
	ldr	x11, [x1, #16]
	mov	x9, x8
	cmp	x10, x11
	b.ne	LBB4_5
; %bb.3:
	ldr	x10, [x0, #24]
	ldr	x11, [x1, #24]
	mov	x9, x8
	cmp	x10, x11
	b.ne	LBB4_5
; %bb.4:
	ldr	x9, [x2]
	add	x9, x9, #1
	str	x9, [x2]
	ldr	x9, [x0]
LBB4_5:
	mvn	w8, w8
	str	x9, [x1]
	ldr	x9, [x0, #8]
	str	x9, [x1, #8]
	ldr	x9, [x0, #16]
	str	x9, [x1, #16]
	ldr	x9, [x0, #24]
	str	x9, [x1, #24]
	and	x0, x8, #0x1
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_v15o_report                    ; -- Begin function v15o_report
	.p2align	2
_v15o_report:                           ; @v15o_report
	.cfi_startproc
; %bb.0:
	stp	x20, x19, [sp, #-16]!           ; 16-byte Folded Spill
	.cfi_def_cfa_offset 16
	.cfi_offset w19, -8
	.cfi_offset w20, -16
	mov	x8, x0
	mov	x11, #0                         ; =0x0
	ldp	x10, x12, [x2]
	mov	x9, #31765                      ; =0x7c15
	movk	x9, #32586, lsl #16
	movk	x9, #31161, lsl #32
	movk	x9, #40503, lsl #48
	mul	x9, x12, x9
	eor	x9, x9, x10
	and	x13, x9, x1
	mov	w14, #96                        ; =0x60
	mov	w0, #2                          ; =0x2
	b	LBB5_2
LBB5_1:                                 ;   in Loop: Header=BB5_2 Depth=1
	add	x11, x11, #1
	cmp	x11, x1
	b.hi	LBB5_8
LBB5_2:                                 ; =>This Inner Loop Header: Depth=1
	add	x9, x13, x11
	and	x9, x9, x1
	madd	x9, x9, x14, x8
	ldp	x17, x7, [x9]
	ldp	x16, x15, [x9, #16]
	orr	x19, x7, x17
	orr	x20, x16, x15
	orr	x19, x19, x20
	cbz	x19, LBB5_7
; %bb.3:                                ;   in Loop: Header=BB5_2 Depth=1
	cmp	x17, x10
	ccmp	x7, x12, #0, eq
	b.ne	LBB5_1
; %bb.4:                                ;   in Loop: Header=BB5_2 Depth=1
	ldr	x17, [x2, #16]
	cmp	x16, x17
	b.ne	LBB5_1
; %bb.5:                                ;   in Loop: Header=BB5_2 Depth=1
	ldr	x16, [x2, #24]
	cmp	x15, x16
	b.ne	LBB5_1
; %bb.6:
	ldr	x8, [x9, #32]
	str	x8, [x5]
	ldr	x8, [x9, #40]
	str	x8, [x5, #8]
	ldr	x8, [x9, #48]
	str	x8, [x5, #16]
	ldr	x8, [x9, #56]
	str	x8, [x5, #24]
	ldr	x8, [x9, #64]
	str	x8, [x6]
	ldr	x8, [x9, #72]
	str	x8, [x6, #8]
	ldr	x8, [x9, #80]
	str	x8, [x6, #16]
	ldr	x8, [x9, #88]
	str	x8, [x6, #24]
	mov	w0, #1                          ; =0x1
	ldp	x20, x19, [sp], #16             ; 16-byte Folded Reload
	ret
LBB5_7:
	mov	x0, #0                          ; =0x0
	str	x10, [x9]
	ldr	x8, [x2, #8]
	str	x8, [x9, #8]
	ldr	x8, [x2, #16]
	str	x8, [x9, #16]
	ldr	x8, [x2, #24]
	str	x8, [x9, #24]
	ldr	x8, [x3]
	str	x8, [x9, #32]
	ldr	x8, [x3, #8]
	str	x8, [x9, #40]
	ldr	x8, [x3, #16]
	str	x8, [x9, #48]
	ldr	x8, [x3, #24]
	str	x8, [x9, #56]
	ldr	x8, [x4]
	str	x8, [x9, #64]
	ldr	x8, [x4, #8]
	str	x8, [x9, #72]
	ldr	x8, [x4, #16]
	str	x8, [x9, #80]
	ldr	x8, [x4, #24]
	str	x8, [x9, #88]
LBB5_8:
	ldp	x20, x19, [sp], #16             ; 16-byte Folded Reload
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
